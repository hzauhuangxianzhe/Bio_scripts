import torch
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import os
import time
from collections import OrderedDict

sys.path.insert(1, os.path.dirname(__file__))
import genotypeio, eigenmt
from core import *


# in cis.py
def calculate_cis_nominal(genotypes_t, phenotype_t, residualizer=None, return_af=True):
    """
    Calculate nominal associations

    genotypes_t: genotypes x samples
    phenotype_t: single phenotype
   
    residualizer: Residualizer object (see core.py)
    """
    p = phenotype_t.reshape(1,-1)
    r_nominal_t, genotype_var_t, phenotype_var_t = calculate_corr(genotypes_t, p, residualizer=residualizer, return_var=True)
    std_ratio_t = torch.sqrt(phenotype_var_t.reshape(1,-1) / genotype_var_t.reshape(-1,1))
    r_nominal_t = r_nominal_t.squeeze()
    r2_nominal_t = r_nominal_t.double().pow(2)

    if residualizer is not None:
        dof = residualizer.dof
    else:
        dof = p.shape[1] - 2
    slope_t = r_nominal_t * std_ratio_t.squeeze()
    # MODIFIED: Clamp r2 to prevent division by zero if r2 is exactly 1
    tstat_t = r_nominal_t * torch.sqrt(dof / (1 - r2_nominal_t.clamp(max=0.999999)))
    slope_se_t = (slope_t.double() / tstat_t).float()

###### 移除了对等位基因频率的计算  #####
    return tstat_t, slope_t, slope_se_t, r2_nominal_t.float() # <--- 在这里增加r2的返回
#####################################

def calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t,
                               residualizer=None, random_tiebreak=False):
    """
    计算名义关联和置换后的最大关联值。

    返回:
        r_nominal_t[ix] (torch.Tensor): 最佳关联TAD的名义相关系数。
        std_ratio_t[ix] (torch.Tensor): 最佳关联TAD的表型与基因型标准差之比。
        ix (torch.Tensor): 最佳关联TAD在输入中的索引。
        r2_perm_t (torch.Tensor): 每次置换中的最大R^2值向量。
    """
    # 对表型进行置换
    permutations_t = phenotype_t[permutation_ix_t]

    # 计算原始（名义）关联，并获取方差用于计算效应大小
    r_nominal_t, genotype_var_t, phenotype_var_t = calculate_corr(
        genotypes_t, phenotype_t.reshape(1, -1),
        residualizer=residualizer, return_var=True
    )
    std_ratio_t = torch.sqrt(phenotype_var_t.reshape(1, -1) / genotype_var_t.reshape(-1, 1))
    r_nominal_t = r_nominal_t.squeeze(dim=-1)
    std_ratio_t = std_ratio_t.squeeze(dim=-1)

    # 计算基因型与所有置换表型之间的关联的平方(R^2)
    corr_t = calculate_corr(genotypes_t, permutations_t, residualizer=residualizer).pow(2)  # [基因型数量 x 置换次数]
    
    # 过滤掉计算中可能出现的NaN值
    corr_t = corr_t[~torch.isnan(corr_t).any(1), :]
    if corr_t.shape[0] == 0:
        raise ValueError('所有相关性计算结果均为NaN，请检查您的表型数据。')
    
    # 对每一次置换，找到所有基因型中的最大R^2值，构建经验分布
    r2_perm_t, _ = corr_t.max(0)

    # 在名义关联中找到R^2最大的基因型
    r2_nominal_t = r_nominal_t.pow(2)
    r2_nominal_t[torch.isnan(r2_nominal_t)] = -1  # 处理NaN值的技巧，使其不会被argmax选中
    
    if not random_tiebreak:
        # 选择R^2最大的索引
        ix = r2_nominal_t.argmax()
    else:
        # 如果存在多个相同的最大R^2值，随机选择一个
        max_val = r2_nominal_t.max()
        ix_all = torch.nonzero(r2_nominal_t == max_val, as_tuple=True)[0]
        ix = ix_all[torch.randint(0, len(ix_all), [1])[0]]
        
    # 返回四个核心统计值
    return r_nominal_t[ix], std_ratio_t[ix], ix, r2_perm_t


####### 修改后的calculate_association函数（移除了等位基因频率的相关内容） #####
def calculate_association(genotype_df, phenotype_s, covariates_df=None,
                          interaction_s=None,
                          logp=False, window=1000000, verbose=True):
    """
    Standalone helper function for computing the association between
    a set of genotypes and a single phenotype.
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    assert genotype_df.columns.equals(phenotype_s.index)

    # copy to GPU
    phenotype_t = torch.tensor(phenotype_s.values, dtype=torch.float).to(device)
    genotypes_t = torch.tensor(genotype_df.values, dtype=torch.float).to(device)
    impute_mean(genotypes_t)

    dof = phenotype_s.shape[0] - 2

    if covariates_df is not None:
        assert phenotype_s.index.equals(covariates_df.index)
        residualizer = Residualizer(torch.tensor(covariates_df.values, dtype=torch.float32).to(device))
        dof -= covariates_df.shape[1]
    else:
        residualizer = None

    if interaction_s is None:
        # **修改点**: calculate_cis_nominal 不再返回 af, ma_samples, ma_count
        res = calculate_cis_nominal(genotypes_t, phenotype_t, residualizer)
        # **修改点**: 相应地修改解包
        tstat, slope, slope_se, r2 = [i.cpu().numpy() for i in res] ## FIXED ##
        # **修改点**: 从输出的DataFrame中移除af, ma_samples, ma_count
        df = pd.DataFrame({
            'pval_nominal': get_t_pval(tstat, dof, log=logp),
            'slope': slope, 
            'slope_se': slope_se,
            'tstat': tstat,
            'r2': r2, ## ADDED ##
        }, index=genotype_df.index)
    else:
        interaction_t = torch.tensor(interaction_s.values.reshape(1,-1), dtype=torch.float32).to(device)
        
        # **修改点**: 完全移除 maf_threshold_interaction 和 filter_maf_interaction 的逻辑
        # 因为不再进行过滤，所以所有的TAD都会被保留
        mask_t = torch.ones(genotypes_t.shape[0], dtype=torch.bool, device=device)
        
        # 假设 calculate_interaction_nominal 也已被修改，不再返回MAF相关统计量
        res = calculate_interaction_nominal(genotypes_t, phenotype_t.unsqueeze(0), interaction_t, residualizer,
                                            return_sparse=False)
        # **修改点**: 相应地修改解包
        tstat, b, b_se = [i.cpu().numpy() for i in res]
        mask = mask_t.cpu().numpy()
        dof -= 2

        # **修改点**: 从输出的DataFrame中移除af, ma_samples, ma_count
        df = pd.DataFrame({
            'pval_g': get_t_pval(tstat[:,0], dof, log=logp), 'b_g': b[:,0], 'b_g_se': b_se[:,0],
            'pval_i': get_t_pval(tstat[:,1], dof, log=logp), 'b_i': b[:,1], 'b_i_se': b_se[:,1],
            'pval_gi': get_t_pval(tstat[:,2], dof, log=logp), 'b_gi': b[:,2], 'b_gi_se': b_se[:,2],
        }, index=genotype_df.index[mask])

    if df.index.str.startswith('chr').all():  # assume chr_pos_ref_alt_build format
        df['position'] = df.index.map(lambda x: int(x.split('_')[1]))
    return df
##############################################################################



# 注意: 此函数假设 calculate_cis_nominal 和 calculate_interaction_nominal 函数
# 已经被修改，不再返回 af, ma_samples, ma_count。否则会导致解包错误。
def map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix,
                covariates_df=None, paired_covariate_df=None, interaction_df=None,
                group_s=None, window=1000000, run_eigenmt=False, logp=False,
                output_dir='.', write_top=True, write_stats=True, logger=None, verbose=True):
    """
    cis-QTL mapping: nominal associations for all variant-phenotype pairs

    Association results for each chromosome are written to parquet files
    in the format <output_dir>/<prefix>.cis_qtl_pairs.<chr>.parquet

    If interaction_df is provided, the top association per phenotype is
    written to <output_dir>/<prefix>.cis_qtl_top_assoc.txt.gz unless
    write_top is set to False, in which case it is returned as a DataFrame
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if logger is None:
        logger = SimpleLogger()
    if group_s is not None:
        group_dict = group_s.to_dict()

    logger.write('cis-QTL mapping: nominal associations for all variant-phenotype pairs')
    logger.write(f'  * {phenotype_df.shape[1]} samples')
    logger.write(f'  * {phenotype_df.shape[0]} phenotypes')
    if covariates_df is not None:
        assert np.all(phenotype_df.columns == covariates_df.index)
        logger.write(f'  * {covariates_df.shape[1]} covariates')
        residualizer = Residualizer(torch.tensor(covariates_df.values, dtype=torch.float32).to(device))
        dof = phenotype_df.shape[1] - 2 - covariates_df.shape[1]
    else:
        residualizer = None
        dof = phenotype_df.shape[1] - 2
    if paired_covariate_df is not None:
        assert covariates_df is not None
        assert paired_covariate_df.index.isin(phenotype_df.index).all(), f"Paired covariate phenotypes must be present in phenotype matrix."
        assert paired_covariate_df.columns.equals(phenotype_df.columns), f"Paired covariate samples must match samples in phenotype matrix."
        paired_covariate_df = paired_covariate_df.T  # samples x phenotypes
        logger.write(f'  * including phenotype-specific covariate')
    logger.write(f'  * {variant_df.shape[0]} variants')
    if interaction_df is not None:
        assert interaction_df.index.equals(phenotype_df.columns)
        logger.write(f"  * including {interaction_df.shape[1]} interaction term(s)")
        ni = interaction_df.shape[1]
        dof -= 2 * ni
        interaction_t = torch.tensor(interaction_df.values, dtype=torch.float32).to(device)
        
        interaction_mask_t = None
        
        col_order = ['phenotype_id', 'variant_id', 'start_distance']
        if 'pos' not in phenotype_pos_df:
            col_order += ['end_distance']
        col_order += ['pval_g', 'b_g', 'b_g_se']
        if ni == 1:
            col_order += ['pval_i', 'b_i', 'b_i_se', 'pval_gi', 'b_gi', 'b_gi_se']
        else:
            col_order += [k.replace('i', f"i{i+1}") for i in range(0,ni) for k in ['pval_i', 'b_i', 'b_i_se', 'pval_gi', 'b_gi', 'b_gi_se']]

        var_dict = []
        for i,v in enumerate(interaction_df.columns, 1):
            for c in ['pval_i', 'b_i', 'b_i_se']:
                var_dict.append((c.replace('_i', f'_i{i}'), c.replace('_i', f'_{v}')))
            for c in ['pval_gi', 'b_gi', 'b_gi_se']:
                var_dict.append((c.replace('_gi', f'_gi{i}'), c.replace('_gi', f'_g-{v}')))
        var_dict = dict(var_dict)

    logger.write(f'  * cis-window: ±{window:,}')

    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)
    
    igc = genotypeio.InputGeneratorCis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, group_s=group_s, window=window)
    best_assoc = []
    start_time = time.time()
    k = 0
    logger.write('  * Computing associations')
    for chrom in igc.chrs:
        logger.write(f'    Mapping chromosome {chrom}')
        n = 0
        if group_s is None:
            for i in igc.phenotype_pos_df[igc.phenotype_pos_df['chr'] == chrom].index:
                j = igc.cis_ranges[i]
                n += j[1] - j[0] + 1
        else:
            for i in igc.group_s[igc.phenotype_pos_df['chr'] == chrom].drop_duplicates().index:
                j = igc.cis_ranges[i]
                n += j[1] - j[0] + 1

        chr_res = OrderedDict()
        chr_res['phenotype_id'] = []
        chr_res['variant_id'] = []
        chr_res['start_distance'] = np.empty(n, dtype=np.int32)
        if 'pos' not in phenotype_pos_df:
            chr_res['end_distance'] = np.empty(n, dtype=np.int32)
            
        if interaction_df is None:
            chr_res['pval_nominal'] = np.empty(n, dtype=np.float64)
            chr_res['slope'] =        np.empty(n, dtype=np.float32)
            chr_res['slope_se'] =     np.empty(n, dtype=np.float32)
            chr_res['r2'] =           np.empty(n, dtype=np.float32)
        else:
            chr_res['pval_g'] =  np.empty(n, dtype=np.float64)
            chr_res['b_g'] =     np.empty(n, dtype=np.float32)
            chr_res['b_g_se'] =  np.empty(n, dtype=np.float32)
            chr_res['pval_i'] =  np.empty([n, ni], dtype=np.float64)
            chr_res['b_i'] =     np.empty([n, ni], dtype=np.float32)
            chr_res['b_i_se'] =  np.empty([n, ni], dtype=np.float32)
            chr_res['pval_gi'] = np.empty([n, ni], dtype=np.float64)
            chr_res['b_gi'] =    np.empty([n, ni], dtype=np.float32)
            chr_res['b_gi_se'] = np.empty([n, ni], dtype=np.float32)

        start = 0
        if group_s is None:
            for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(igc.generate_data(chrom=chrom, verbose=verbose), k+1):
                phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
                genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
                genotypes_t = genotypes_t[:,genotype_ix_t]
                impute_mean(genotypes_t)

                variant_ids = variant_df.index[genotype_range[0]:genotype_range[-1]+1]
                start_distance = np.int32(variant_df['pos'].values[genotype_range[0]:genotype_range[-1]+1] - igc.phenotype_start[phenotype_id])
                if 'pos' not in phenotype_pos_df:
                    end_distance = np.int32(variant_df['pos'].values[genotype_range[0]:genotype_range[-1]+1] - igc.phenotype_end[phenotype_id])

                var_t = genotypes_t.var(dim=1)
                valid_mask = var_t > 1e-8
                if not torch.any(valid_mask):
                    logger.write(f"    * WARNING: skipping phenotype {phenotype_id} (all TADs have zero variance in window)")
                    continue
                
                genotypes_t = genotypes_t[valid_mask]
                mask = valid_mask.cpu().numpy().astype(bool)
                variant_ids = variant_ids[mask]
                start_distance = start_distance[mask]
                if 'pos' not in phenotype_pos_df:
                    end_distance = end_distance[mask]

                if paired_covariate_df is None or phenotype_id not in paired_covariate_df:
                    iresidualizer = residualizer
                else:
                    iresidualizer = Residualizer(torch.tensor(np.c_[covariates_df, paired_covariate_df[phenotype_id]],
                                                              dtype=torch.float32).to(device))

                if interaction_df is None:
                    res = calculate_cis_nominal(genotypes_t, phenotype_t, residualizer=iresidualizer)
                    tstat, slope, slope_se, r2 = [i.cpu().numpy() for i in res]
                    n = len(variant_ids)
                else:
                    if genotypes_t.shape[0] > 0:
                        res = calculate_interaction_nominal(genotypes_t, phenotype_t.unsqueeze(0), interaction_t,
                                                             residualizer=iresidualizer, return_sparse=False,
                                                             variant_ids=variant_ids)
                        tstat, b, b_se = [i.cpu().numpy() for i in res]
                        n = len(variant_ids)

                        ix = np.nanargmax(np.abs(tstat[:,1+ni:]).max(1))
                        order = [0] + [i if j % 2 == 0 else i+ni for i in range(1,ni+1) for j in range(2)]
                        top_s = [phenotype_id, variant_ids[ix], start_distance[ix]]
                        if 'pos' not in phenotype_pos_df:
                            top_s += [end_distance[ix]]
                        for i in order:
                            top_s += [tstat[ix,i], b[ix,i], b_se[ix,i]]
                        top_s = pd.Series(top_s, index=col_order)
                        if run_eigenmt:
                            top_s['tests_emt'] = eigenmt.compute_tests(genotypes_t, var_thresh=0.99, variant_window=200)
                        best_assoc.append(top_s)
                    else:
                        n=0
                
                if n > 0:
                    chr_res['phenotype_id'].extend([phenotype_id]*n)
                    chr_res['variant_id'].extend(variant_ids)
                    chr_res['start_distance'][start:start+n] = start_distance
                    if 'pos' not in phenotype_pos_df:
                        chr_res['end_distance'][start:start+n] = end_distance
                    if interaction_df is None:
                        chr_res['pval_nominal'][start:start+n] = tstat
                        chr_res['slope'][start:start+n] = slope
                        chr_res['slope_se'][start:start+n] = slope_se
                        chr_res['r2'][start:start+n] = r2
                    else:
                        chr_res['pval_g'][start:start+n]  = tstat[:,0]
                        chr_res['b_g'][start:start+n]     = b[:,0]
                        chr_res['b_g_se'][start:start+n]  = b_se[:,0]
                        chr_res['pval_i'][start:start+n]  = tstat[:,1:1+ni]
                        chr_res['b_i'][start:start+n]     = b[:,1:1+ni]
                        chr_res['b_i_se'][start:start+n]  = b_se[:,1:1+ni]
                        chr_res['pval_gi'][start:start+n] = tstat[:,1+ni:]
                        chr_res['b_gi'][start:start+n]    = b[:,1+ni:]
                        chr_res['b_gi_se'][start:start+n] = b_se[:,1+ni:]

                start += n
        else: # groups
            for k, (phenotypes, genotypes, genotype_range, phenotype_ids, group_id) in enumerate(igc.generate_data(chrom=chrom, verbose=verbose), k+1):
                genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
                genotypes_t = genotypes_t[:,genotype_ix_t]
                impute_mean(genotypes_t)

                variant_ids = variant_df.index[genotype_range[0]:genotype_range[-1]+1]
                start_distance = np.int32(variant_df['pos'].values[genotype_range[0]:genotype_range[-1]+1] - igc.phenotype_start[phenotype_ids[0]])
                if 'pos' not in phenotype_pos_df:
                    end_distance = np.int32(variant_df['pos'].values[genotype_range[0]:genotype_range[-1]+1] - igc.phenotype_end[phenotype_ids[0]])

                var_t = genotypes_t.var(dim=1)
                valid_mask = var_t > 1e-8
                if not torch.any(valid_mask):
                    logger.write(f"    * WARNING: skipping group {group_id} (all TADs have zero variance in window)")
                    continue

                genotypes_t = genotypes_t[valid_mask]
                mask = valid_mask.cpu().numpy().astype(bool)
                variant_ids = variant_ids[mask]
                start_distance = start_distance[mask]
                if 'pos' not in phenotype_pos_df:
                    end_distance = end_distance[mask]
                
                n = len(variant_ids)

                if genotypes_t.shape[0] > 0:
                    phenotype_id = phenotype_ids[0]
                    phenotype_t = torch.tensor(phenotypes[0], dtype=torch.float).to(device)
                    
                    if paired_covariate_df is not None and phenotype_id in paired_covariate_df:
                         iresidualizer = Residualizer(torch.tensor(np.c_[covariates_df, paired_covariate_df[phenotype_id]],
                                                          dtype=torch.float32).to(device))
                    else:
                        iresidualizer = residualizer
                    
                    if interaction_df is None:
                        res = calculate_cis_nominal(genotypes_t, phenotype_t, residualizer=iresidualizer)
                        tstat, slope, slope_se, r2 = [i.cpu().numpy() for i in res]
                    else:
                        res = calculate_interaction_nominal(genotypes_t, phenotype_t.unsqueeze(0), interaction_t,
                                                             residualizer=iresidualizer, return_sparse=False,
                                                             variant_ids=variant_ids)
                        tstat, b, b_se = [i.cpu().numpy() for i in res]
                    px = [phenotype_id]*n

                    for p_idx, (phenotype, phenotype_id) in enumerate(zip(phenotypes, phenotype_ids)):
                        if p_idx == 0: continue
                        phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
                        
                        if paired_covariate_df is not None and phenotype_id in paired_covariate_df:
                             iresidualizer = Residualizer(torch.tensor(np.c_[covariates_df, paired_covariate_df[phenotype_id]],
                                                              dtype=torch.float32).to(device))
                        else:
                            iresidualizer = residualizer

                        if interaction_df is None:
                            res = calculate_cis_nominal(genotypes_t, phenotype_t, residualizer=iresidualizer)
                            tstat0, slope0, slope_se0, r2_0 = [i.cpu().numpy() for i in res]
                        else:
                            res = calculate_interaction_nominal(genotypes_t, phenotype_t.unsqueeze(0), interaction_t,
                                                                 residualizer=iresidualizer, return_sparse=False,
                                                                 variant_ids=variant_ids)
                            tstat0, b0, b_se0 = [i.cpu().numpy() for i in res]

                        if interaction_df is None:
                            ix = np.where(np.abs(tstat0) > np.abs(tstat))[0]
                        else:
                            ix = np.where(np.abs(tstat0[:,2]) > np.abs(tstat[:,2]))[0]

                        for j in ix:
                            px[j] = phenotype_id
                        if interaction_df is None:
                            tstat[ix] = tstat0[ix]
                            slope[ix] = slope0[ix]
                            slope_se[ix] = slope_se0[ix]
                            r2[ix] = r2_0[ix]
                        else:
                            tstat[ix] = tstat0[ix]
                            b[ix] = b0[ix]
                            b_se[ix] = b_se0[ix]

                    chr_res['phenotype_id'].extend(px)
                    chr_res['variant_id'].extend(variant_ids)
                    chr_res['start_distance'][start:start+n] = start_distance
                    if 'pos' not in phenotype_pos_df:
                        chr_res['end_distance'][start:start+n] = end_distance
                    
                    if interaction_df is None:
                        chr_res['pval_nominal'][start:start+n] = tstat
                        chr_res['slope'][start:start+n] = slope
                        chr_res['slope_se'][start:start+n] = slope_se
                        chr_res['r2'][start:start+n] = r2
                    else:
                        chr_res['pval_g'][start:start+n]  = tstat[:,0]
                        chr_res['b_g'][start:start+n]     = b[:,0]
                        chr_res['b_g_se'][start:start+n]  = b_se[:,0]
                        chr_res['pval_i'][start:start+n]  = tstat[:,1:1+ni]
                        chr_res['b_i'][start:start+n]     = b[:,1:1+ni]
                        chr_res['b_i_se'][start:start+n]  = b_se[:,1:1+ni]
                        chr_res['pval_gi'][start:start+n] = tstat[:,1+ni:]
                        chr_res['b_gi'][start:start+n]    = b[:,1+ni:]
                        chr_res['b_gi_se'][start:start+n] = b_se[:,1+ni:]

                start += n

        logger.write(f'    time elapsed: {(time.time()-start_time)/60:.2f} min')

        if start < len(chr_res['start_distance']):
            for x in chr_res:
                if isinstance(chr_res[x], np.ndarray):
                    chr_res[x] = chr_res[x][:start]
                else: # is list
                    chr_res[x] = chr_res[x][:start]

        if write_stats and start > 0:
            if interaction_df is not None:
                cols = ['pval_i', 'b_i', 'b_i_se', 'pval_gi', 'b_gi', 'b_gi_se']
                if ni == 1:
                    for c in cols:
                        chr_res[c] = chr_res[c][:,0]
                else:
                    for i in range(0, ni):
                        for c in cols:
                            chr_res[c.replace('i', f"i{i+1}")] = None
                    for c in cols:
                        for i in range(0, ni):
                            chr_res[c.replace('i', f"i{i+1}")] = chr_res[c][:,i]
                        del chr_res[c]
            
            chr_res_df = pd.DataFrame(chr_res)

            if paired_covariate_df is not None:
                idof = dof - chr_res_df['phenotype_id'].isin(paired_covariate_df.columns).astype(int).values
            else:
                idof = dof

            if interaction_df is None:
                m = chr_res_df['pval_nominal'].notnull()
                m = m[m].index
                if paired_covariate_df is not None:
                    idof = idof[m]
                chr_res_df.loc[m, 'pval_nominal'] = get_t_pval(chr_res_df.loc[m, 'pval_nominal'], idof, log=logp)
            else:
                if ni == 1:
                    m = chr_res_df['pval_gi'].notnull()
                    m = m[m].index
                    if paired_covariate_df is not None:
                        idof = idof[m]
                    chr_res_df.loc[m, 'pval_g'] =  get_t_pval(chr_res_df.loc[m, 'pval_g'], idof, log=logp)
                    chr_res_df.loc[m, 'pval_i'] =  get_t_pval(chr_res_df.loc[m, 'pval_i'], idof, log=logp)
                    chr_res_df.loc[m, 'pval_gi'] = get_t_pval(chr_res_df.loc[m, 'pval_gi'], idof, log=logp)
                else:
                    m = chr_res_df[f'pval_gi{ni}'].notnull()
                    m = m[m].index
                    if paired_covariate_df is not None:
                        idof = idof[m]
                    chr_res_df.loc[m, 'pval_g'] = get_t_pval(chr_res_df.loc[m, 'pval_g'], idof, log=logp)
                    for i in range(1, ni+1):
                        chr_res_df.loc[m, f'pval_i{i}'] =  get_t_pval(chr_res_df.loc[m, f'pval_i{i}'], idof, log=logp)
                        chr_res_df.loc[m, f'pval_gi{i}'] = get_t_pval(chr_res_df.loc[m, f'pval_gi{i}'], idof, log=logp)
                    chr_res_df.rename(columns=var_dict, inplace=True)
                
            print('    * writing output')
            chr_res_df.to_parquet(os.path.join(output_dir, f'{prefix}.cis_qtl_pairs.{chrom}.parquet'))

    if interaction_df is not None and len(best_assoc) > 0:
        best_assoc = pd.concat(best_assoc, axis=1, sort=False).T.set_index('phenotype_id').infer_objects()
        pval_cols = [c for c in best_assoc.columns if c.startswith('pval_')]
        for x in pval_cols:
            best_assoc[x] = best_assoc[x].astype(np.float64)
        m = best_assoc['pval_g'].notnull()
        m = m[m].index
        if paired_covariate_df is not None:
            idof = dof - best_assoc.index.isin(paired_covariate_df.columns).astype(int).values[m]
        else:
            idof = dof
        best_assoc.loc[m, 'pval_g'] = get_t_pval(best_assoc.loc[m, 'pval_g'], idof, log=logp)
        if ni == 1:
            best_assoc.loc[m, 'pval_i'] =  get_t_pval(best_assoc.loc[m, 'pval_i'], idof, log=logp)
            best_assoc.loc[m, 'pval_gi'] = get_t_pval(best_assoc.loc[m, 'pval_gi'], idof, log=logp)
        else:
            for i in range(1, ni+1):
                best_assoc.loc[m, f'pval_i{i}'] =  get_t_pval(best_assoc.loc[m, f'pval_i{i}'], idof, log=logp)
                best_assoc.loc[m, f'pval_gi{i}'] = get_t_pval(best_assoc.loc[m, f'pval_gi{i}'], idof, log=logp)
        if run_eigenmt and ni == 1:
            if group_s is None:
                best_assoc['pval_emt'] = np.minimum(best_assoc['tests_emt']*best_assoc['pval_gi'], 1)
            else:
                best_assoc['pval_emt'] = np.minimum(best_assoc['num_phenotypes']*best_assoc['tests_emt']*best_assoc['pval_gi'], 1)
            best_assoc['pval_adj_bh'] = eigenmt.padjust_bh(best_assoc['pval_emt'])
        if ni > 1:
            best_assoc.rename(columns=var_dict, inplace=True)
        if write_top:
            best_assoc.to_csv(os.path.join(output_dir, f'{prefix}.cis_qtl_top_assoc.txt.gz'),
                               sep='\t', float_format='%.6g')
        else:
            return best_assoc
            
    logger.write('done.')

# in cis.py
def prepare_cis_output(r_nominal, r2_perm, std_ratio, num_var, dof, variant_id,
                       start_distance, end_distance, phenotype_id, nperm=10000, logp=False):
    """Return nominal p-value, etc. as pd.Series"""
    r2_nominal = r_nominal * r_nominal
    r2_nominal = np.clip(r2_nominal, 0, 0.9999999)
    
    pval_perm = (np.sum(r2_perm >= r2_nominal)+1) / (nperm+1)

    slope = r_nominal * std_ratio
    tstat2 = dof * r2_nominal / (1 - r2_nominal)
    slope_se = np.abs(slope) / np.sqrt(tstat2) if tstat2 > 0 else np.inf

    res_s = pd.Series(OrderedDict([
        ('num_var', num_var),
        ('beta_shape1', np.nan),
        ('beta_shape2', np.nan),
        ('true_df', np.nan),
        ('pval_true_df', np.nan),
        ('variant_id', variant_id),
        ('start_distance', start_distance),
        ('end_distance', end_distance),
        ('pval_nominal', pval_from_corr(r2_nominal, dof, logp=logp)),
        ('slope', slope),
        ('slope_se', slope_se),
        ('r2', r2_nominal),  # <--- 在这里添加r2
        ('pval_perm', pval_perm),
        ('pval_beta', np.nan),
    ]), name=phenotype_id)
    return res_s

# <--- [DOF修正] 此函数的参数名是 'dof'，它会接收来自 map_cis 的 'idof' 变量
def _process_group_permutations(buf, variant_df, start_pos, end_pos, dof, group_id, nperm=10000, beta_approx=True, logp=False):
    """
    Merge results for grouped phenotypes
    buf: [r_nominal, std_ratio, var_ix, r2_perm, num_var, phenotype_id]
    """
    max_ix = np.argmax(np.abs([b[0] for b in buf]))
    r_nominal, std_ratio, var_ix = buf[max_ix][:3]
    
    num_var, phenotype_id = buf[max_ix][-2], buf[max_ix][-1]
    
    r2_perm = np.max([b[3] for b in buf], 0)
    variant_id = variant_df.index[var_ix]
    start_distance = variant_df['pos'].values[var_ix] - start_pos
    end_distance = variant_df['pos'].values[var_ix] - end_pos
    
    # <--- [DOF修正] 此处调用 prepare_cis_output 时，传递的 dof 是上游传入的 idof
    res_s = prepare_cis_output(r_nominal, r2_perm, std_ratio, num_var, dof, variant_id, 
                               start_distance, end_distance, phenotype_id, nperm=nperm, logp=logp)
    
    if beta_approx:
        # <--- [DOF修正] 恢复原始代码中的 dof*0.25 因子，以保持统计模型的一致性
        res_s[['pval_beta', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df']] = calculate_beta_approx_pval(r2_perm, r_nominal*r_nominal, dof*0.25)
    res_s['group_id'] = group_id
    res_s['group_size'] = len(buf)
    return res_s
def map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=None,
            group_s=None, paired_covariate_df=None, beta_approx=True, nperm=10000,
            window=1000000, random_tiebreak=False, logger=None, seed=None, logp=False,
            verbose=True, warn_monomorphic=True, variance_threshold=1e-8):
    """Run cis-QTL mapping"""

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if logger is None:
        logger = SimpleLogger()

    logger.write('cis-QTL mapping: empirical p-values for phenotypes')
    logger.write(f'  * {phenotype_df.shape[1]} samples')
    logger.write(f'  * {phenotype_df.shape[0]} phenotypes')
    if group_s is not None:
        logger.write(f'  * {len(group_s.unique())} phenotype groups')
        group_dict = group_s.to_dict()
    if covariates_df is not None:
        assert covariates_df.index.equals(phenotype_df.columns), 'Sample names in phenotype matrix columns and covariate matrix rows do not match!'
        assert ~(covariates_df.isnull().any().any()), f'Missing or null values in covariates matrix, in columns {",".join(covariates_df.columns[covariates_df.isnull().any(axis=0)].astype(str))}'
        logger.write(f'  * {covariates_df.shape[1]} covariates')
        residualizer = Residualizer(torch.tensor(covariates_df.values, dtype=torch.float32).to(device))
        # <--- [DOF修正] 此处 'dof' 是基础自由度
        dof = phenotype_df.shape[1] - 2 - covariates_df.shape[1]
    else:
        residualizer = None
        dof = phenotype_df.shape[1] - 2
    if paired_covariate_df is not None:
        assert covariates_df is not None
        assert paired_covariate_df.index.isin(phenotype_df.index).all(), f"Paired covariate phenotypes must be present in phenotype matrix."
        assert paired_covariate_df.columns.equals(phenotype_df.columns), f"Paired covariate samples must match samples in phenotype matrix."
        paired_covariate_df = paired_covariate_df.T
        logger.write(f'  * including phenotype-specific covariate')
    logger.write(f'  * {genotype_df.shape[0]} variants')
    if random_tiebreak:
        logger.write(f'  * randomly selecting top variant in case of ties')
    logger.write(f'  * cis-window: ±{window:,}')

    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)

    n_samples = phenotype_df.shape[1]
    ix = np.arange(n_samples)
    if seed is not None:
        logger.write(f'  * using seed {seed}')
        np.random.seed(seed)
    permutation_ix_t = torch.LongTensor(np.array([np.random.permutation(ix) for i in range(nperm)])).to(device)

    res_df = []
    igc = genotypeio.InputGeneratorCis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, group_s=group_s, window=window)
    if igc.n_phenotypes == 0:
        raise ValueError('No valid phenotypes found.')
    start_time = time.time()
    logger.write('  * computing permutations')
    if group_s is None:
        for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(igc.generate_data(verbose=verbose), 1):
            genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
            genotypes_t = genotypes_t[:,genotype_ix_t]
            impute_mean(genotypes_t)

            var_t = genotypes_t.var(dim=1)
            mono_t = var_t < variance_threshold
            if mono_t.any():
                genotypes_t = genotypes_t[~mono_t]
                genotype_range = genotype_range[~mono_t.cpu()]
                if warn_monomorphic:
                    logger.write(f'    * WARNING: excluding {mono_t.sum()} low-variance TADs for phenotype {phenotype_id}')

            if genotypes_t.shape[0] == 0:
                logger.write(f'WARNING: skipping {phenotype_id} (no valid variants after filtering)')
                continue

            phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
            # <--- [DOF修正] 'idof' 在此为每个表型单独计算
            if paired_covariate_df is None or phenotype_id not in paired_covariate_df:
                iresidualizer = residualizer
                idof = dof
            else:
                iresidualizer = Residualizer(torch.tensor(np.c_[covariates_df, paired_covariate_df[phenotype_id]],
                                                          dtype=torch.float32).to(device))
                idof = dof - 1
            res = calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t,
                                             residualizer=iresidualizer, random_tiebreak=random_tiebreak)
            r_nominal, std_ratio, var_ix, r2_perm = [i.cpu().numpy() for i in res]
            var_ix = genotype_range[var_ix]
            variant_id = variant_df.index[var_ix]
            start_distance = variant_df['pos'].values[var_ix] - igc.phenotype_start[phenotype_id]
            end_distance = variant_df['pos'].values[var_ix] - igc.phenotype_end[phenotype_id]
            
            res_s = prepare_cis_output(r_nominal, r2_perm, std_ratio, genotypes_t.shape[0], idof, variant_id,
                                       start_distance, end_distance, phenotype_id, nperm=nperm, logp=logp)
            if beta_approx:
                res_s[['pval_beta', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df']] = calculate_beta_approx_pval(r2_perm, r_nominal*r_nominal, idof)
            res_df.append(res_s)
    else:  # grouped mode
        for k, (phenotypes, genotypes, genotype_range, phenotype_ids, group_id) in enumerate(igc.generate_data(verbose=verbose), 1):
            genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
            genotypes_t = genotypes_t[:,genotype_ix_t]
            impute_mean(genotypes_t)

            var_t = genotypes_t.var(dim=1)
            mono_t = var_t < variance_threshold
            if mono_t.any():
                genotypes_t = genotypes_t[~mono_t]
                genotype_range = genotype_range[~mono_t.cpu()]
                if warn_monomorphic:
                    logger.write(f'    * WARNING: excluding {mono_t.sum()} low-variance TADs for group {group_id}')

            if genotypes_t.shape[0] == 0:
                logger.write(f'WARNING: skipping {group_id} (no valid variants after filtering)')
                continue

            buf = []
            # <--- [DOF修正] 'idof' 在循环中会被覆盖，这是遵循原始代码的逻辑
            idof = dof 
            for phenotype, phenotype_id in zip(phenotypes, phenotype_ids):
                phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
                if paired_covariate_df is None or phenotype_id not in paired_covariate_df:
                    iresidualizer = residualizer
                    idof = dof
                else:
                    iresidualizer = Residualizer(torch.tensor(np.c_[covariates_df, paired_covariate_df[phenotype_id]],
                                                              dtype=torch.float32).to(device))
                    idof = dof - 1
                res = calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t,
                                                 residualizer=iresidualizer, random_tiebreak=random_tiebreak)
                res_list = [i.cpu().numpy() for i in res]
                res_list[2] = genotype_range[res_list[2]]
                buf.append(res_list[:4] + [genotypes_t.shape[0], phenotype_id])

            # <--- [DOF修正] 此处传递的 'idof' 是组内最后一个表型的自由度，与原始代码行为一致
            res_s = _process_group_permutations(buf, variant_df, igc.phenotype_start[phenotype_ids[0]],
                                                igc.phenotype_end[phenotype_ids[0]], idof,
                                                group_id, nperm=nperm, beta_approx=beta_approx, logp=logp)
            res_df.append(res_s)

    res_df = pd.concat(res_df, axis=1, sort=False).T
    res_df.index.name = 'phenotype_id'
    logger.write(f'  Time elapsed: {(time.time()-start_time)/60:.2f} min')
    logger.write('done.')
    return res_df.astype(output_dtype_dict).infer_objects()

# 注意: 此函数及其正确运行，强依赖于其调用的下游函数（prepare_cis_output, _process_group_permutations, 
# 和 calculate_cis_permutations）都已被修改，以移除所有与 MAF/af/ma 相关的逻辑。
def map_independent(genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df=None,
                    group_s=None, fdr=0.05, fdr_col='qval', nperm=10000, window=1000000,
                    # <--- [修正-Grok4反馈1.1] 将缺失值默认改为 np.nan，并添加 variance_threshold 参数
                    missing=np.nan, random_tiebreak=False, logger=None, seed=None, logp=False, verbose=True,
                    variance_threshold=1e-8):
    """
    Run independent cis-QTL mapping (forward-backward regression)
    cis_df: output from map_cis, annotated with q-values (calculate_qvalues)
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if logger is None:
        logger = SimpleLogger()

    logger.write('cis-QTL mapping: conditionally independent variants')
    logger.write(f'  * {phenotype_df.shape[1]} samples')

    signif_df = cis_df[cis_df[fdr_col] <= fdr].copy()
    if len(signif_df) == 0:
        raise ValueError(f"No significant phenotypes at FDR ≤ {fdr}.")
    
    cols = ['num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df', 'variant_id',
            'start_distance', 'end_distance',
            'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta']
    if group_s is not None:
        cols += ['group_id', 'group_size']
    signif_df = signif_df[cols]
    
    signif_threshold = signif_df['pval_beta'].max()
    if group_s is None:
        ix = phenotype_df.index[phenotype_df.index.isin(signif_df.index)]
        logger.write(f'  * {signif_df.shape[0]}/{cis_df.shape[0]} significant phenotypes')
    else:
        ix = group_s[phenotype_df.index].loc[group_s[phenotype_df.index].isin(signif_df['group_id'])].index
        logger.write(f'  * {signif_df.shape[0]}/{cis_df.shape[0]} significant groups')
        logger.write(f'    {len(ix)}/{phenotype_df.shape[0]} phenotypes')
        group_dict = group_s.to_dict()

    assert phenotype_df.index.equals(phenotype_pos_df.index)

    if covariates_df is not None:
        assert covariates_df.index.equals(phenotype_df.columns), 'Sample names in phenotype matrix columns and covariate matrix rows do not match!'
        assert ~(covariates_df.isnull().any().any()), f'Missing or null values in covariates matrix.'
        logger.write(f'  * {covariates_df.shape[1]} covariates')

    logger.write(f'  * {genotype_df.shape[0]} variants')
    # <--- [修正-Grok4反馈1.2] 添加日志记录，告知用户正在使用的方差阈值
    logger.write(f'  * using variance_threshold={variance_threshold}')
    if random_tiebreak:
        logger.write(f'  * randomly selecting top variant in case of ties')
    logger.write(f'  * cis-window: ±{window:,}')
    phenotype_df = phenotype_df.loc[ix]
    phenotype_pos_df = phenotype_pos_df.loc[ix]

    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)
    ix_dict = {i:k for k,i in enumerate(genotype_df.index)}

    n_samples = phenotype_df.shape[1]
    ix = np.arange(n_samples)
    if seed is not None:
        logger.write(f'  * using seed {seed}')
        np.random.seed(seed)
    permutation_ix_t = torch.LongTensor(np.array([np.random.permutation(ix) for i in range(nperm)])).to(device)

    res_df = []
    igc = genotypeio.InputGeneratorCis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, group_s=group_s, window=window)
    if igc.n_phenotypes == 0:
        raise ValueError('No valid phenotypes found.')
    logger.write('  * computing independent QTLs')
    start_time = time.time()
    if group_s is None:
        for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(igc.generate_data(verbose=verbose), 1):
            phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
            genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
            genotypes_t = genotypes_t[:,genotype_ix_t]
            impute_mean(genotypes_t)

            # <--- [修正-Grok4反馈3.3] 添加方差过滤
            var_t = genotypes_t.var(dim=1)
            valid_mask = var_t > variance_threshold
            if not torch.any(valid_mask):
                logger.write(f"    * WARNING: skipping phenotype {phenotype_id} (all TADs in window have low variance)")
                continue
            if (~valid_mask).any():
                 logger.write(f'    * filtered {(~valid_mask).sum()} low-variance TADs for phenotype {phenotype_id}')
            genotypes_t = genotypes_t[valid_mask]
            genotype_range = genotype_range[valid_mask.cpu().numpy()]
            
            # 1) forward pass
            forward_df = [signif_df.loc[phenotype_id]]
            covariates = covariates_df.values.copy() if covariates_df is not None else np.empty((phenotype_df.shape[1], 0))
            dosage_dict = {}
            while True:
                variant_id = forward_df[-1]['variant_id']
                ig = genotype_df.values[ix_dict[variant_id], genotype_ix].copy()
                
                # <--- [修正-Grok4反馈1.1, 3.4] 修改缺失值判断和填补逻辑
                if np.isnan(missing):
                    m = np.isnan(ig)
                    if np.any(m):
                        ig[m] = np.nanmean(ig)
                else:
                    m = ig == missing
                    if np.any(m):
                        ig[m] = np.mean(ig[~m])
                
                dosage_dict[variant_id] = ig
                covariates = np.hstack([covariates, ig.reshape(-1,1)]).astype(np.float32)
                dof = phenotype_df.shape[1] - 2 - covariates.shape[1]
                residualizer = Residualizer(torch.tensor(covariates, dtype=torch.float32).to(device))

                res = calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t,
                                                 residualizer=residualizer, random_tiebreak=random_tiebreak)
                
                # <--- [修正-Grok4反馈3.2] 明确解包4个返回值，假设上游函数已修改
                r_nominal, std_ratio, var_ix, r2_perm = [i.cpu().numpy() for i in res[:4]]
                x = calculate_beta_approx_pval(r2_perm, r_nominal*r_nominal, dof)
                
                if x[0] <= signif_threshold:
                    var_ix = genotype_range[var_ix]
                    variant_id = variant_df.index[var_ix]
                    start_distance = variant_df['pos'].values[var_ix] - igc.phenotype_start[phenotype_id]
                    end_distance = variant_df['pos'].values[var_ix] - igc.phenotype_end[phenotype_id]
                    
                    res_s = prepare_cis_output(r_nominal, r2_perm, std_ratio, genotypes.shape[0], dof, variant_id,
                                               start_distance, end_distance, phenotype_id, nperm=nperm, logp=logp)
                    res_s[['pval_beta', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df']] = x
                    forward_df.append(res_s)
                else:
                    break
            forward_df = pd.concat(forward_df, axis=1, sort=False).T
            dosage_df = pd.DataFrame(dosage_dict)

            # 2) backward pass
            if forward_df.shape[0] > 1:
                back_df = []
                variant_set = set()
                for k_b, i in enumerate(forward_df['variant_id'], 1):
                    covariates_b = dosage_df[np.setdiff1d(forward_df['variant_id'], i)].values
                    if covariates_df is not None:
                         covariates_b = np.hstack([covariates_df.values, covariates_b])
                    
                    dof = phenotype_df.shape[1] - 2 - covariates_b.shape[1]
                    residualizer = Residualizer(torch.tensor(covariates_b, dtype=torch.float32).to(device))

                    res = calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t,
                                                     residualizer=residualizer, random_tiebreak=random_tiebreak)
                    r_nominal, std_ratio, var_ix, r2_perm = [i.cpu().numpy() for i in res[:4]]
                    var_ix = genotype_range[var_ix]
                    variant_id = variant_df.index[var_ix]
                    x = calculate_beta_approx_pval(r2_perm, r_nominal*r_nominal, dof)
                    if x[0] <= signif_threshold and variant_id not in variant_set:
                        start_distance = variant_df['pos'].values[var_ix] - igc.phenotype_start[phenotype_id]
                        end_distance = variant_df['pos'].values[var_ix] - igc.phenotype_end[phenotype_id]
                        res_s = prepare_cis_output(r_nominal, r2_perm, std_ratio, genotypes.shape[0], dof, variant_id,
                                                   start_distance, end_distance, phenotype_id, nperm=nperm, logp=logp)
                        res_s[['pval_beta', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df']] = x
                        res_s['rank'] = k_b
                        back_df.append(res_s)
                        variant_set.add(variant_id)
                if len(back_df)>0:
                    res_df.append(pd.concat(back_df, axis=1, sort=False).T)
            else:
                forward_df['rank'] = 1
                res_df.append(forward_df)

    else:  # grouped phenotypes
        for k, (phenotypes, genotypes, genotype_range, phenotype_ids, group_id) in enumerate(igc.generate_data(verbose=verbose), 1):
            genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
            genotypes_t = genotypes_t[:,genotype_ix_t]
            impute_mean(genotypes_t)

            # <--- 修改点: (分组模式) 添加方差过滤
            var_t = genotypes_t.var(dim=1)
            valid_mask = var_t > variance_threshold
            if not torch.any(valid_mask):
                logger.write(f"    * WARNING: skipping group {group_id} (all TADs in window have low variance)")
                continue
            if (~valid_mask).any():
                logger.write(f'    * filtered {(~valid_mask).sum()} low-variance TADs for group {group_id}')
            genotypes_t = genotypes_t[valid_mask]
            genotype_range = genotype_range[valid_mask.cpu().numpy()]
            
            # 1) forward pass
            forward_df = [signif_df[signif_df['group_id'] == group_id].iloc[0]]
            covariates = covariates_df.values.copy() if covariates_df is not None else np.empty((phenotype_df.shape[1], 0))
            dosage_dict = {}
            while True:
                variant_id = forward_df[-1]['variant_id']
                ig = genotype_df.values[ix_dict[variant_id], genotype_ix].copy()
                
                # <--- 修改点: (分组模式) 修改缺失值判断和填补逻辑
                if np.isnan(missing):
                    m = np.isnan(ig)
                    if np.any(m):
                        ig[m] = np.nanmean(ig)
                else:
                    m = ig == missing
                    if np.any(m):
                        ig[m] = np.mean(ig[~m])

                dosage_dict[variant_id] = ig
                covariates = np.hstack([covariates, ig.reshape(-1,1)]).astype(np.float32)
                dof = phenotype_df.shape[1] - 2 - covariates.shape[1]
                residualizer = Residualizer(torch.tensor(covariates, dtype=torch.float32).to(device))

                buf = []
                # <--- [修正-Grok4反馈4.4] 此处dof的计算在分组模式下是一种简化，假设组内无paired covariate差异
                for phenotype, phenotype_id in zip(phenotypes, phenotype_ids):
                    phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
                    res = calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t,
                                                     residualizer=residualizer, random_tiebreak=random_tiebreak)
                    res_list = [i.cpu().numpy() for i in res]
                    res_list[2] = genotype_range[res_list[2]]
                    
                    # <--- 修改点: (分组模式) 修正buf构建，不再包含'g'
                    buf.append(res_list[:4] + [genotypes.shape[0], phenotype_id])
                
                res_s = _process_group_permutations(buf, variant_df, igc.phenotype_start[phenotype_ids[0]],
                                                    igc.phenotype_end[phenotype_ids[0]], dof, group_id, nperm=nperm, logp=logp)

                if res_s['pval_beta'] <= signif_threshold:
                    forward_df.append(res_s)
                else:
                    break
            forward_df = pd.concat(forward_df, axis=1, sort=False).T
            dosage_df = pd.DataFrame(dosage_dict)

            # 2) backward pass
            if forward_df.shape[0] > 1:
                back_df = []
                variant_set = set()
                for k_b, variant_id in enumerate(forward_df['variant_id'], 1):
                    covariates_b = dosage_df[np.setdiff1d(forward_df['variant_id'], variant_id)].values
                    if covariates_df is not None:
                        covariates_b = np.hstack([covariates_df.values, covariates_b])
                    
                    dof = phenotype_df.shape[1] - 2 - covariates_b.shape[1]
                    residualizer = Residualizer(torch.tensor(covariates_b, dtype=torch.float32).to(device))

                    buf = []
                    for phenotype, phenotype_id in zip(phenotypes, phenotype_ids):
                        phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
                        res = calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t,
                                                         residualizer=residualizer, random_tiebreak=random_tiebreak)
                        res_list = [i.cpu().numpy() for i in res]
                        res_list[2] = genotype_range[res_list[2]]
                        buf.append(res_list[:4] + [genotypes.shape[0], phenotype_id])
                    
                    res_s = _process_group_permutations(buf, variant_df, igc.phenotype_start[phenotype_ids[0]],
                                                        igc.phenotype_end[phenotype_ids[0]], dof, group_id, nperm=nperm, logp=logp)

                    if res_s['pval_beta'] <= signif_threshold and variant_id not in variant_set:
                        res_s['rank'] = k_b
                        back_df.append(res_s)
                        variant_set.add(variant_id)
                if len(back_df)>0:
                    res_df.append(pd.concat(back_df, axis=1, sort=False).T)
            else:
                forward_df['rank'] = 1
                res_df.append(forward_df)

    # <--- [修正-Grok4反馈5.1] 增加对空结果集的健壮处理
    res_df = pd.concat(res_df, axis=0, sort=False) if res_df else pd.DataFrame()
    if not res_df.empty:
        res_df.index.name = 'phenotype_id'
        logger.write(f'  Time elapsed: {(time.time()-start_time)/60:.2f} min')
    logger.write('done.')
    return res_df.reset_index().astype(output_dtype_dict) if not res_df.empty else res_df