#!/usr/bin/env python3
from __future__ import print_function
import pandas as pd
import numpy as np
from datetime import datetime
import sys
import os
import re
import pickle
import argparse
from collections import defaultdict
import importlib.metadata

sys.path.insert(1, os.path.dirname(__file__))
from core import *
from post import *
# <--- 修改点: 移除了不再需要的 'trans' 和 'susie' 模块导入
import genotypeio, cis
from genotypeio import validate_sample_consistency

def main():
    # <--- 修改点: 更新描述，移除已删除的模式
    parser = argparse.ArgumentParser(description='tensorQTL: GPU-based QTL mapper for cis-mapping modes.')

    parser.add_argument('genotype_path', help='TAD IS scores file (tab-delimited)')
    parser.add_argument('phenotypes', help="Phenotypes in BED format (.bed, .bed.gz, .bed.parquet).")
    parser.add_argument('prefix', help='Prefix for output file names')
    # <--- 修改点: 更新 --mode 参数的选项，仅保留 cis, cis_nominal, cis_independent
    parser.add_argument('--mode', type=str, default='cis', choices=['cis', 'cis_nominal', 'cis_independent'],
                        help='Mapping mode. Default: cis')
    parser.add_argument('--covariates', default=None, help='Covariates file, tab-delimited (covariates x samples)')
    parser.add_argument('--paired_covariate', default=None, help='Single phenotype-specific covariate, tab-delimited (phenotypes x samples)')
    parser.add_argument('--permutations', type=int, default=10000, help='Number of permutations for cis mode. Default: 10000')
    parser.add_argument('--interaction', default=None, type=str, help='Tab-delimited file mapping sample ID to interaction value(s) for cis_nominal mode.')
    parser.add_argument('--cis_output', default=None, type=str, help="Output from 'cis' mode with q-values. Required for independent cis-QTL mapping.")
    parser.add_argument('--phenotype_groups', default=None, type=str, help='Phenotype groups. Header-less TSV with two columns: phenotype_id, group_id')
    parser.add_argument('--window', default=1000000, type=np.int32, help='Cis-window size, in bases. Default: 1000000.')
    parser.add_argument('--logp', action='store_true', help='Compute nominal p-values as -log10(P) for added precision (requires R)')
    
    # <--- 修改点: 移除仅与 trans 和 susie 相关的命令行参数
    # parser.add_argument('--pval_threshold', ...)
    # parser.add_argument('--return_dense', ...)
    # parser.add_argument('--return_r2', ...)
    # parser.add_argument('--batch_size', ...)
    # parser.add_argument('--susie_loci', ...)
    # parser.add_argument('--max_effects', ...)

    parser.add_argument('--best_only', action='store_true', help='Only write lead association for each phenotype (interaction mode in cis_nominal only)')
    parser.add_argument('--chunk_size', default=None, help="For cis-QTL mapping with large data files, load features in chunks.")
    parser.add_argument('--disable_beta_approx', action='store_true', help='Disable Beta-distribution approximation of empirical p-values (not recommended).')
    parser.add_argument('--warn_monomorphic', action='store_true', help='Warn if features with low variance are found.')
    parser.add_argument('--variance_threshold', default=1e-8, type=np.float64, help='Variance threshold for filtering features (e.g., TADs) with low variance. Default: 1e-8')
    parser.add_argument('--fdr', default=0.05, type=np.float64, help='FDR for cis-QTLs')
    parser.add_argument('--qvalue_lambda', default=None, type=np.float64, help='lambda parameter for pi0est in qvalue.')
    parser.add_argument('--seed', default=None, type=int, help='Seed for permutations.')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    # <--- 修改点: 简化输入检查，因为 trans 模式已被移除
    if args.mode == 'cis_independent' and (args.cis_output is None or not os.path.exists(args.cis_output)):
        raise ValueError("Output from 'cis' mode must be provided for 'cis_independent' mode.")
    if args.interaction is not None and args.mode not in ['cis_nominal']:
        raise ValueError("Interactions are only supported in 'cis_nominal' mode.")

    logger = SimpleLogger(os.path.join(args.output_dir, f'{args.prefix}.tensorQTL.{args.mode}.log'))
    logger.write(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Running TensorQTL v{importlib.metadata.version("tensorqtl")}: {args.mode.split("_")[0]}-QTL mapping')
    if torch.cuda.is_available():
        logger.write(f'  * using GPU ({torch.cuda.get_device_name(torch.cuda.current_device())})')
    else:
        logger.write('  * WARNING: using CPU!')
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if args.seed is not None:
        logger.write(f'  * using seed {args.seed}')

    logger.write(f'  * reading phenotypes ({args.phenotypes})')
    # <--- 修改点: 简化表型文件读取逻辑，因为所有剩余模式都需要BED格式
    assert args.phenotypes.lower().endswith(('.bed', '.bed.gz', '.bed.parquet')), "For cis modes, phenotypes must be in BED format."
    phenotype_df, phenotype_pos_df = read_phenotype_bed(args.phenotypes)
    if 'pos' in phenotype_pos_df.columns:
        logger.write(f"  * cis-window detected as position ± {args.window:,}")
    else:
        logger.write(f"  * cis-window detected as [start - {args.window:,}, end + {args.window:,}]")

    if args.covariates is not None:
        logger.write(f'  * reading covariates ({args.covariates})')
        covariates_df = pd.read_csv(args.covariates, sep='\t', index_col=0).T
        assert phenotype_df.columns.equals(covariates_df.index)
    else:
        covariates_df = None

    if args.paired_covariate is not None:
        assert covariates_df is not None, "Covariates matrix must be provided when using paired covariate"
        paired_covariate_df = pd.read_csv(args.paired_covariate, sep='\t', index_col=0)
        assert paired_covariate_df.index.isin(phenotype_df.index).all()
        assert paired_covariate_df.columns.equals(phenotype_df.columns)
    else:
        paired_covariate_df = None

    if args.interaction is not None:
        logger.write(f'  * reading interaction term(s) ({args.interaction})')
        with open(args.interaction) as f:
            f.readline()
            s = f.readline().strip()
        interaction_df = pd.read_csv(args.interaction, sep='\t', index_col=0, header=None if len(s.split('\t')) == 2 else 'infer')
        assert covariates_df.index.isin(interaction_df.index).all()
        interaction_df = interaction_df.loc[covariates_df.index].astype(np.float32)
    else:
        interaction_df = None

    if args.phenotype_groups is not None:
        group_s = pd.read_csv(args.phenotype_groups, sep='\t', index_col=0, header=None).squeeze('columns').rename(None)
        group_dict = group_s.to_dict()
        previous_group = ''
        parsed_groups = 0
        for i in phenotype_df.index:
            if group_dict.get(i) != previous_group:
                parsed_groups += 1
                previous_group = group_dict.get(i)
        if not parsed_groups == len(group_s.unique()):
            raise ValueError('Groups defined in input do not match phenotype file (check sort order).')
    else:
        group_s = None

    if args.chunk_size is None:
        logger.write(f'  * loading TAD IS scores ({args.genotype_path})')
        tad_df, variant_df = genotypeio.load_tad_data(args.genotype_path, select_samples=phenotype_df.columns)
        tad_df, covariates_df = validate_sample_consistency(phenotype_df, tad_df, covariates_df)
        genotype_df = tad_df
        if variant_df is None:
            raise ValueError("Data must include variant positions for cis-mapping modes.")
    else:
        raise NotImplementedError("Chunking logic is not yet adapted for general TAD data files and has been disabled in this version.")

    if args.mode == 'cis':
        res_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df,
                             group_s=group_s, paired_covariate_df=paired_covariate_df, nperm=args.permutations,
                             window=args.window, beta_approx=not args.disable_beta_approx,
                             warn_monomorphic=args.warn_monomorphic, logger=logger, seed=args.seed, verbose=True,
                             variance_threshold=args.variance_threshold)
        logger.write('  * writing output')
        if has_rpy2:
            calculate_qvalues(res_df, fdr=args.fdr, qvalue_lambda=args.qvalue_lambda, logger=logger)
        out_file = os.path.join(args.output_dir, f'{args.prefix}.cis_qtl.txt.gz')
        res_df.to_csv(out_file, sep='\t', float_format='%.6g')

    elif args.mode == 'cis_nominal':
        cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, args.prefix, covariates_df=covariates_df,
                        paired_covariate_df=paired_covariate_df, interaction_df=interaction_df,
                        group_s=group_s, window=args.window, run_eigenmt=True,
                        output_dir=args.output_dir, write_top=True, write_stats=not args.best_only, logger=logger, verbose=True)
        
        if args.cis_output is not None:
            cis_df = pd.read_csv(args.cis_output, sep='\t', index_col=0)
            nominal_prefix = os.path.join(args.output_dir, f'{args.prefix}.cis_qtl_pairs')
            signif_df = get_significant_pairs(cis_df, nominal_prefix, group_s=group_s, fdr=args.fdr)
            signif_df.to_parquet(os.path.join(args.output_dir, f'{args.prefix}.cis_qtl.signif_pairs.parquet'))

    elif args.mode == 'cis_independent':
        summary_df = pd.read_csv(args.cis_output, sep='\t', index_col=0)
        
        res_df = cis.map_independent(genotype_df, variant_df, summary_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df,
                                     group_s=group_s, fdr=args.fdr, nperm=args.permutations, window=args.window,
                                     logger=logger, seed=args.seed, verbose=True, variance_threshold=args.variance_threshold)

        logger.write('  * writing output')
        out_file = os.path.join(args.output_dir, f'{args.prefix}.cis_independent_qtl.txt.gz')
        res_df.to_csv(out_file, sep='\t', index=False, float_format='%.6g')

    # <--- 修改点: 移除了 cis_susie, trans_susie, 和 trans 模式的 elif 代码块

    logger.write(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] Finished mapping')


if __name__ == '__main__':
    main()