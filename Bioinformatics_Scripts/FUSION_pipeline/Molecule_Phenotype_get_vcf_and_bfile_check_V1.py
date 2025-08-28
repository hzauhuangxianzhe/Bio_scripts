#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CIWAS主脚本测试模式输出结果核查工具 (最终增强版)
=====================================================

本脚本用于验证主脚本在测试模式（--test_mode）下生成结果的正确性。
它执行四项核心检查：
1.  **VCF范围核查**: 验证提取的SNP是否都落在正确的基因组窗口内。
2.  **SNP数量交叉验证**: 比较VCF内含的SNP数与PLINK日志报告的SNP数是否一致。
3.  **FAM表型值核查**: 验证生成的.fam文件中的表型值是否与原始输入表型文件精确匹配。
4.  **SNP完整性核查 (新增)**: 通过直接查询原始VCF，验证提取出的SNP是否有遗漏或冗余。

它还会为每个表型的表型值生成一张对比散点图，用于可视化验证。

用法:
    python verify_test_mode.py --mode gene \\
        --vcf_file path/to/original_large.vcf.gz \\
        --pos_file path/to/gene_positions_1based.tsv \\
        --pheno_file path/to/gene_phenotypes.tsv \\
        --output_dir ciwas_gene_output \\
        --window_size 250 \\
        --gene_id_col GeneID \\
        --gene_chr_col Chromosome \\
        --gene_tss_col TSS_1based
"""

import os
import sys
import argparse
import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt
import re
import glob
from tqdm import tqdm

def parse_arguments():
    """解析与主脚本一致的命令行参数，用于定位文件"""
    parser = argparse.ArgumentParser(
        description='CIWAS主脚本测试模式输出结果核查工具 (最终增强版)',
        formatter_class=argparse.RawTextHelpFormatter
    )
    # 新增 --vcf_file 参数，用于访问原始VCF
    parser.add_argument('--vcf_file', type=str, required=True, help='原始的、大的VCF文件路径 (用于验证SNP完整性)')
    parser.add_argument('--mode', type=str, required=True, choices=['gene', 'tad', 'loop'], help='选择要核查的分子表型模式')
    parser.add_argument('--pos_file', type=str, required=True, help='原始的位置信息TSV文件路径')
    parser.add_argument('--pheno_file', type=str, required=True, help='原始的表型值TSV文件路径')
    parser.add_argument('--output_dir', type=str, required=True, help='主脚本的输出目录')
    parser.add_argument('--window_size', type=int, required=True, help='主脚本使用的窗口大小 (kb)')
    
    # 为了简化和稳健，将模式相关参数设为可选，然后在代码中检查
    # 基因模式
    parser.add_argument('--gene_id_col', type=str)
    parser.add_argument('--gene_chr_col', type=str)
    parser.add_argument('--gene_tss_col', type=str)
    # TAD模式
    parser.add_argument('--tad_id_col', type=str)
    parser.add_argument('--tad_chr_col', type=str)
    parser.add_argument('--tad_start_col', type=str)
    parser.add_argument('--tad_end_col', type=str)
    # Loop模式
    parser.add_argument('--loop_id_col', type=str)
    parser.add_argument('--loop_chr_col', type=str)
    parser.add_argument('--loop_a1_start_col', type=str)
    parser.add_argument('--loop_a1_end_col', type=str)
    parser.add_argument('--loop_a2_start_col', type=str)
    parser.add_argument('--loop_a2_end_col', type=str)

    args = parser.parse_args()

    # 手动检查模式所需参数是否存在
    if args.mode == 'gene' and not all([args.gene_id_col, args.gene_chr_col, args.gene_tss_col]):
        parser.error("基因模式下，必须提供 --gene_id_col, --gene_chr_col, 和 --gene_tss_col")
    if args.mode == 'tad' and not all([args.tad_id_col, args.tad_chr_col, args.tad_start_col, args.tad_end_col]):
        parser.error("TAD模式下，必须提供 --tad_id_col, --tad_chr_col, --tad_start_col, 和 --tad_end_col")
    if args.mode == 'loop' and not all([args.loop_id_col, args.loop_chr_col, args.loop_a1_start_col, args.loop_a1_end_col, args.loop_a2_start_col, args.loop_a2_end_col]):
        parser.error("Loop模式下，必须提供所有 --loop_*_col 参数")

    return args

def calculate_expected_regions(pheno_row, mode, window_kb):
    """根据主脚本的逻辑，重新计算一个分子表型的理论基因组区域 (0-based)"""
    window_bp = window_kb * 1000
    regions = []
    centers = []
    if mode in ['gene', 'tad']:
        centers.append(pheno_row['center'])
    elif mode == 'loop':
        centers.append(pheno_row['center1'])
        centers.append(pheno_row['center2'])
    
    for center in centers:
        start_0based = max(0, int(center) - window_bp - 1)
        end_0based = int(center) + window_bp
        regions.append({'chrom': str(pheno_row['chrom']), 'start': start_0based, 'end': end_0based})
        
    return regions

def get_variants_from_vcf(vcf_path, regions=None):
    """从VCF文件提取变异，可选地只提取指定区域"""
    variants = set()
    with pysam.VariantFile(vcf_path) as vcf:
        if regions:
            # 从指定区域提取
            for region in regions:
                try:
                    for record in vcf.fetch(region['chrom'], region['start'], region['end']):
                        variants.add((record.chrom, record.pos, record.ref, record.alts))
                except ValueError:
                    # 忽略不存在的染色体
                    pass
        else:
            # 提取整个文件
            for record in vcf.fetch():
                variants.add((record.chrom, record.pos, record.ref, record.alts))
    return variants

def verify_vcf_range(vcf_path, expected_regions):
    """核查VCF文件中的所有SNP是否都落在预期的基因组区域内"""
    min_pos, max_pos, vcf_chrom = float('inf'), float('-inf'), None
    with pysam.VariantFile(vcf_path) as vcf:
        records = list(vcf.fetch())
        if not records: return True, "VCF为空 (区域内无SNP)", "N/A"
        
        for record in records:
            pos_1based = record.pos + 1
            min_pos = min(min_pos, pos_1based)
            max_pos = max(max_pos, pos_1based)
            if vcf_chrom is None: vcf_chrom = record.chrom
    
    is_in_range = False
    for region in expected_regions:
        expected_start_1based = region['start'] + 1
        expected_end_1based = region['end']
        if vcf_chrom == region['chrom'] and min_pos >= expected_start_1based and max_pos <= expected_end_1based:
            is_in_range = True
            break
    
    actual_range_str = f"[{min_pos}-{max_pos}] on {vcf_chrom}"
    expected_range_str = " or ".join([f"[{r['start']+1}-{r['end']}] on {r['chrom']}" for r in expected_regions])
    return is_in_range, actual_range_str, expected_range_str

def verify_snp_count(vcf_path, plink_log_path):
    """交叉验证VCF中的SNP数和PLINK日志中报告的数量"""
    vcf_count = len(get_variants_from_vcf(vcf_path))
    plink_count = "未找到"
    try:
        with open(plink_log_path, 'r') as f:
            log_content = f.read()
            match = re.search(r'(\d+)\s+variants loaded from', log_content) or re.search(r'--vcf: (\d+) variants remaining', log_content)
            if match: plink_count = int(match.group(1))
    except Exception: plink_count = "无法解析"
    return vcf_count == plink_count, vcf_count, plink_count

def verify_fam_phenotypes(fam_path, original_phenos):
    """精确核对.fam文件中的表型值与原始表型值"""
    try:
        fam_df = pd.read_csv(fam_path, sep='\s+', header=None, names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'])
        output_phenos = fam_df['PHENO'].values
        original_phenos_filled = original_phenos.fillna(-9).values
        if len(output_phenos) != len(original_phenos_filled): return False, "样本数量不匹配"
        return np.all(np.isclose(output_phenos, original_phenos_filled)), "通过"
    except Exception as e: return False, f"读取或比较时出错: {e}"

def verify_snp_completeness(extracted_vcf_path, original_vcf_path, expected_regions):
    """【新增】验证提取的SNP是否有遗漏或冗余"""
    try:
        # 从提取的小VCF获取SNP集合
        extracted_snps = get_variants_from_vcf(extracted_vcf_path)
        # 从原始大VCF的对应区域获取基准SNP集合
        ground_truth_snps = get_variants_from_vcf(original_vcf_path, regions=expected_regions)
        
        # 计算差集
        missing_snps = ground_truth_snps - extracted_snps
        extra_snps = extracted_snps - ground_truth_snps
        
        is_ok = not missing_snps and not extra_snps
        return is_ok, len(missing_snps), len(extra_snps)
    except Exception as e:
        return False, -1, -1 # -1 表示检查过程中出错

def plot_pheno_comparison(original_phenos, fam_path, plot_dir, pheno_id):
    """为表型比较结果生成散点图"""
    try:
        fam_df = pd.read_csv(fam_path, sep='\s+', header=None, names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'])
        output_phenos = fam_df['PHENO']
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(original_phenos, output_phenos, alpha=0.6, edgecolors='k', s=50)
        lims = [np.nanmin([ax.get_xlim(), ax.get_ylim()]), np.nanmax([ax.get_xlim(), ax.get_ylim()])]
        ax.plot(lims, lims, 'r--', alpha=0.75, zorder=0, label='y=x (Perfect Match)')
        ax.set_aspect('equal', adjustable='box'); ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_title(f'Phenotype Comparison for {pheno_id}'); ax.set_xlabel('Original Phenotype Value'); ax.set_ylabel('Phenotype Value in .fam File'); ax.grid(True, linestyle='--', alpha=0.5); ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f"{pheno_id}_pheno_comparison.png"))
        plt.close(fig)
    except Exception as e: print(f"    - [可视化] 警告：为ID '{pheno_id}' 生成图像失败: {e}")

def main():
    """主函数，执行所有核查流程"""
    args = parse_arguments()
    
    print("=" * 80 + "\nCIWAS主脚本测试模式输出结果核查工具 (最终增强版)\n" + "=" * 80)

    # --- 1. 加载原始输入文件作为“事实标准” ---
    print(f"正在加载事实标准文件...")
    pos_df = pd.read_csv(args.pos_file, sep='\t')
    if args.mode == 'gene':
        col_map = {args.gene_id_col: 'id', args.gene_chr_col: 'chrom', args.gene_tss_col: 'center'}
    elif args.mode == 'tad':
        col_map = {args.tad_id_col: 'id', args.tad_chr_col: 'chrom', args.tad_start_col: 'start', args.tad_end_col: 'end'}
        pos_df = pos_df.rename(columns=col_map)
        pos_df['center'] = ((pd.to_numeric(pos_df['start']) + pd.to_numeric(pos_df['end'])) / 2).astype(int)
    elif args.mode == 'loop':
        col_map = {args.loop_id_col: 'id', args.loop_chr_col: 'chrom', args.loop_a1_start_col: 'a1_start', args.loop_a1_end_col: 'a1_end', args.loop_a2_start_col: 'a2_start', args.loop_a2_end_col: 'a2_end'}
        pos_df = pos_df.rename(columns=col_map)
        pos_df['center1'] = ((pd.to_numeric(pos_df['a1_start']) + pd.to_numeric(pos_df['a1_end'])) / 2).astype(int)
        pos_df['center2'] = ((pd.to_numeric(pos_df['a2_start']) + pd.to_numeric(pos_df['a2_end'])) / 2).astype(int)
    if args.mode != 'tad' and args.mode != 'loop': # tad and loop already renamed
        pos_df = pos_df.rename(columns=col_map)
    pos_df = pos_df.set_index('id')

    pheno_df = pd.read_csv(args.pheno_file, sep='\t')
    pheno_df = pheno_df.set_index(pheno_df.columns[0])
    
    # --- 2. 查找所有待核查的输出文件 ---
    vcf_dir = os.path.join(args.output_dir, 'vcf')
    plink_dir = os.path.join(args.output_dir, 'plink')
    plot_dir = os.path.join(args.output_dir, 'verification_plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    vcf_files = sorted(glob.glob(os.path.join(vcf_dir, '*.vcf')))
    if not vcf_files:
        print(f"\n错误：在目录 '{vcf_dir}' 中未找到任何 .vcf 文件。"); sys.exit(1)
        
    print(f"\n在输出目录中找到 {len(vcf_files)} 个VCF文件待核查。开始逐一验证...\n")

    all_passed = True
    for vcf_path in tqdm(vcf_files, desc="核查进度"):
        pheno_id = os.path.basename(vcf_path).replace('.vcf', '')
        print(f"\n--- 正在核查ID: {pheno_id} ---")
        
        fam_path = os.path.join(plink_dir, f"{pheno_id}.fam")
        plink_log_path = os.path.join(plink_dir, f"{pheno_id}.log")
        
        if not all(os.path.exists(p) for p in [fam_path, plink_log_path]):
            print(f"  - [错误] 缺少对应的输出文件，跳过对'{pheno_id}'的核查。"); all_passed = False; continue
        try:
            pos_row = pos_df.loc[pheno_id]
            original_phenos = pheno_df.loc[pheno_id]
        except KeyError:
            print(f"  - [错误] 无法在原始输入文件中找到ID '{pheno_id}'，跳过核查。"); all_passed = False; continue

        # --- 执行核查 ---
        expected_regions = calculate_expected_regions(pos_row, args.mode, args.window_size)
        
        # 1. VCF范围核查
        range_ok, actual_r, expected_r = verify_vcf_range(vcf_path, expected_regions)
        print(f"  1. VCF范围核查: {'通过' if range_ok else '失败！'}")
        if not range_ok: print(f"     - 理论范围: {expected_r}\n     - 实际范围: {actual_r}"); all_passed = False

        # 2. SNP数量交叉验证
        count_ok, vcf_c, plink_c = verify_snp_count(vcf_path, plink_log_path)
        print(f"  2. SNP数量核查: {'通过' if count_ok else '失败！'}")
        print(f"     - VCF内数量: {vcf_c}\n     - PLINK日志数量: {plink_c}")
        if not count_ok: all_passed = False

        # 3. FAM表型值核查
        pheno_ok, reason = verify_fam_phenotypes(fam_path, original_phenos)
        print(f"  3. FAM表型值核查: {'通过' if pheno_ok else '失败！'}")
        if not pheno_ok: print(f"     - 失败原因: {reason}"); all_passed = False

        # 4. 【新增】SNP完整性核查
        comp_ok, missing, extra = verify_snp_completeness(vcf_path, args.vcf_file, expected_regions)
        print(f"  4. SNP完整性核查: {'通过' if comp_ok else '失败！'}")
        print(f"     - 遗漏的SNP数: {missing}\n     - 冗余的SNP数: {extra}")
        if not comp_ok: all_passed = False

        # 5. 可视化
        plot_pheno_comparison(original_phenos, fam_path, plot_dir, pheno_id)

    # --- 4. 输出最终总结 ---
    print("\n" + "=" * 80)
    if all_passed: print("核查结论: 所有检查项均已通过！脚本输出结果正确、完整。")
    else: print("核查结论: 发现问题！部分检查项失败，请检查上面的日志。")
    print(f"可视化散点图已保存在: {plot_dir}")
    print("=" * 80)

if __name__ == "__main__":
    main()