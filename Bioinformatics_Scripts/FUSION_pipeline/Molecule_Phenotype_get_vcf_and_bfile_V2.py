#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CIWAS多分子表型bfile文件组提取脚本 (Definitive Production Version)
======================================================================

本脚本为复杂性状全信息关联分析（CIWAS）流程设计，旨在自动化地为多种分子表型
（基因、TAD边界、染色质环）提取对应的基因组区域内的SNP信息，并生成PLINK bfile文件组。

*** 数据坐标系统核心假设 ***
本脚本假设所有用户输入的位置坐标文件 (--pos_file) 中的坐标均为【1-based】，
且区间是闭合的（inclusive, 即包含起始和结束位置）。例如，一个从第100位碱基到
第200位碱基的区域，应表示为 start=100, end=200。TSS等单点坐标也应为1-based。
脚本内部会自动处理向pysam等工具所需的0-based坐标系统的转换。

此版本经过多轮迭代与严格审查，修复了包括坐标系统转换、样本ID校验、数据类型处理
等多个关键问题，确保了数据处理的准确性和流程的稳健性。

核心流程:
1.  **数据加载与严格校验**:
    - 校验VCF与表型文件的样本ID是否完全一致，并按VCF顺序重排表型数据。
    - 检查输入文件中的ID是否有重复。
    - 强制将坐标列转换为数值类型，剔除无效行。
    - 使用内连接（inner join）安全地合并位置与表型数据，并向用户报告因不匹配而丢弃的ID数量。

2.  **区域定义与精确转换**:
    - 基于统一的1-based输入假设，精确计算每个分子表型的分析区域。
    - 将计算出的1-based区域精确转换为pysam等工具所需的0-based、半开半闭区间。

3.  **SNP提取与精确去重**:
    - 从总VCF文件中提取SNP，去重逻辑能正确处理双等位基因的SNP和短Indel。

4.  **格式转换**: 并行将VCF转换为PLINK bfile格式。

5.  **表型信息更新**:
    - 保留表型值为高精度的浮点数，完美支持连续型变量。

用法示例:
    # 基因模式 (输入文件中的坐标为1-based)
    python this_script.py --mode gene \\
        --vcf_file path/to/cotton.vcf.gz \\
        --pos_file path/to/gene_positions_1based.tsv \\
        --pheno_file path/to/gene_phenotypes.tsv \\
        --output_dir ciwas_gene_output \\
        --window_size 250 \\
        --threads 8 \\
        --gene_id_col GeneID \\
        --gene_chr_col Chromosome \\
        --gene_tss_col TSS_1based
"""

import os
import sys
import argparse
import pandas as pd
import subprocess
from tqdm import tqdm
import logging
import pysam
from concurrent.futures import ProcessPoolExecutor
import glob

# --- 日志配置 ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('molecular_phenotype_processing.log', mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='为CIWAS多分子表型提取PLINK bfile文件组 (Definitive Production Version)',
        epilog="""
重要: 本脚本假设所有输入坐标 (--pos_file) 均为 1-based 系统。
例如: TSS位点, 区间的 start 和 end。
""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('--vcf_file', type=str, required=True, help='总的VCF文件路径 (必须是bgzip压缩并已建立tabix索引的.vcf.gz文件)')
    parser.add_argument('--pos_file', type=str, required=True, help='包含分子表型位置信息的TSV文件路径 (坐标须为1-based)')
    parser.add_argument('--pheno_file', type=str, required=True, help='包含分子表型ID和对应样本表型值的TSV文件路径')
    parser.add_argument('--output_dir', type=str, default='ciwas_output', help='所有输出文件的根目录 (默认: ciwas_output)')
    parser.add_argument('--window_size', type=int, default=250, help='从核心位点向单侧扩展的窗口大小 (单位: kb, 默认: 250)')
    parser.add_argument('--threads', type=int, default=4, help='并行处理时使用的CPU核心数 (默认: 4)')
    parser.add_argument('--test_mode', action='store_true', default=False, help='测试模式：开启后，仅处理每个模式下的前20个分子表型')

    # 将模式选择和所有相关参数移至主解析器
    parser.add_argument('--mode', type=str, required=True, choices=['gene', 'tad', 'loop'], help='选择要处理的分子表型模式')

    # 基因模式参数
    parser.add_argument('--gene_id_col', type=str, help='位置文件中代表基因ID的列名')
    parser.add_argument('--gene_chr_col', type=str, help='位置文件中代表染色体的列名')
    parser.add_argument('--gene_tss_col', type=str, help='位置文件中代表基因TSS位点的列名 (1-based)')

    # TAD边界模式参数
    parser.add_argument('--tad_id_col', type=str, help='位置文件中代表TAD边界ID的列名')
    parser.add_argument('--tad_chr_col', type=str, help='位置文件中代表染色体的列名')
    parser.add_argument('--tad_start_col', type=str, help='位置文件中TAD边界起始位置的列名 (1-based)')
    parser.add_argument('--tad_end_col', type=str, help='位置文件中TAD边界结束位置的列名 (1-based)')

    # 染色质环(Loop)模式参数
    parser.add_argument('--loop_id_col', type=str, help='位置文件中代表Loop ID的列名')
    parser.add_argument('--loop_chr_col', type=str, help='位置文件中代表染色体的列名')
    parser.add_argument('--loop_a1_start_col', type=str, help='位置文件中Anchor1起始位置的列名 (1-based)')
    parser.add_argument('--loop_a1_end_col', type=str, help='位置文件中Anchor1结束位置的列名 (1-based)')
    parser.add_argument('--loop_a2_start_col', type=str, help='位置文件中Anchor2起始位置的列名 (1-based)')
    parser.add_argument('--loop_a2_end_col', type=str, help='位置文件中Anchor2结束位置的列名 (1-based)')

    return parser.parse_args()


def check_dependencies():
    """检查PLINK是否已安装"""
    try:
        subprocess.run(['plink', '--help'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logger.info("依赖检查通过: PLINK 可用。")
    except FileNotFoundError:
        logger.error("依赖检查失败: PLINK 未安装或不在系统路径中。")
        sys.exit(1)


def check_vcf_file(vcf_file):
    """检查VCF文件格式和索引"""
    if not vcf_file.endswith('.gz'):
        logger.error(f"VCF文件格式错误: 必须是bgzip压缩格式 (.gz结尾): {vcf_file}")
        return False
    
    index_file = vcf_file + '.tbi'
    if not os.path.exists(index_file):
        logger.warning(f"VCF索引文件 (.tbi) 不存在，正在尝试创建...")
        try:
            pysam.tabix_index(vcf_file, preset="vcf", force=True)
            logger.info(f"成功创建VCF索引: {index_file}")
        except Exception as e:
            logger.error(f"创建VCF索引失败: {e}")
            return False
    return True


def load_data(args):
    """加载、校验和合并所有输入数据"""
    logger.info("开始加载、校验和合并数据...")
    try:
        with pysam.VariantFile(args.vcf_file) as vcf:
            vcf_samples = list(vcf.header.samples)
        if not vcf_samples:
            logger.error("VCF文件中没有找到样本信息。")
            sys.exit(1)
        logger.info(f"VCF文件中找到 {len(vcf_samples)} 个样本作为黄金标准。")

        pheno_df = pd.read_csv(args.pheno_file, sep='\t')
        pheno_id_col = pheno_df.columns[0]
        pheno_samples = list(pheno_df.columns[1:])

        if set(vcf_samples) != set(pheno_samples):
            logger.error("样本ID不匹配！VCF文件中的样本ID与表型文件中的样本ID集合不一致。")
            sys.exit(1)
        
        pheno_df = pheno_df[[pheno_id_col] + vcf_samples]
        logger.info("表型文件样本ID校验通过，并已按VCF样本顺序重排。")

        if pheno_df[pheno_id_col].duplicated().any():
            logger.error(f"表型文件 '{args.pheno_file}' 的ID列 '{pheno_id_col}' 中存在重复值。")
            sys.exit(1)
        pheno_df = pheno_df.set_index(pheno_id_col)
        
        pos_df = pd.read_csv(args.pos_file, sep='\t')
        
        col_map, coord_cols = {}, []
        if args.mode == 'gene':
            # 这里的 col_map 被修改，不再包含 'id' 作为键，而是直接使用原始列名
            pos_id_col = args.gene_id_col
            col_map = {'chrom': args.gene_chr_col, 'center': args.gene_tss_col}
            coord_cols = ['center']
        elif args.mode == 'tad':
            pos_id_col = args.tad_id_col
            col_map = {'chrom': args.tad_chr_col, 'start': args.tad_start_col, 'end': args.tad_end_col}
            coord_cols = ['start', 'end']
        elif args.mode == 'loop':
            pos_id_col = args.loop_id_col
            col_map = {
                'chrom': args.loop_chr_col,
                'a1_start': args.loop_a1_start_col, 'a1_end': args.loop_a1_end_col,
                'a2_start': args.loop_a2_start_col, 'a2_end': args.loop_a2_end_col
            }
            coord_cols = ['a1_start', 'a1_end', 'a2_start', 'a2_end']
        
        # 这一步是关键修正，检查 pos_id_col 是否存在
        if pos_id_col not in pos_df.columns:
            logger.error(f"位置文件中缺少ID列: '{pos_id_col}' (来自模式 '{args.mode}')")
            sys.exit(1)

        # 检查坐标列
        for col in col_map.values():
            if col not in pos_df.columns:
                logger.error(f"位置文件中缺少列: '{col}' (来自模式 '{args.mode}')")
                sys.exit(1)
        
        # 使用 pos_id_col 和 col_map 来重命名列，并确保 pos_id_col 在 col_map 里
        # 修正：将 pos_id_col 映射到 'id'
        col_map_with_id = {pos_id_col: 'id'}
        col_map_with_id.update({v: k for k, v in col_map.items()})
        pos_df = pos_df.rename(columns=col_map_with_id)
        
        # 修正：现在所有操作都使用重命名后的列名
        if pos_df['id'].duplicated().any():
            logger.error(f"位置文件 '{args.pos_file}' 的ID列 '{pos_id_col}' 中存在重复值。")
            sys.exit(1)

        for col in coord_cols:
            if col not in pos_df.columns:
                # 再次检查以防万一
                logger.error(f"内部错误：重命名后缺少坐标列 '{col}'。请检查脚本或输入。")
                sys.exit(1)
            pos_df[col] = pd.to_numeric(pos_df[col], errors='coerce')
        
        invalid_rows = pos_df[coord_cols].isnull().any(axis=1)
        if invalid_rows.any():
            logger.warning(f"位置文件中有 {invalid_rows.sum()} 行包含无效的非数值坐标，这些行将被忽略。")
            pos_df = pos_df.dropna(subset=coord_cols)

        if args.mode == 'tad':
            pos_df['center'] = ((pos_df['start'] + pos_df['end']) / 2).astype(int)
        elif args.mode == 'loop':
            pos_df['center1'] = ((pos_df['a1_start'] + pos_df['a1_end']) / 2).astype(int)
            pos_df['center2'] = ((pos_df['a2_start'] + pos_df['a2_end']) / 2).astype(int)

        logger.info(f"位置文件包含 {len(pos_df)} 个有效ID。表型文件包含 {len(pheno_df)} 个有效ID。")
        # 修正：现在可以安全地使用 'id' 列进行合并了
        merged_df = pd.merge(pos_df, pheno_df, left_on='id', right_index=True, how='inner')
        logger.info(f"数据合并后，共有 {len(merged_df)} 个ID将在后续流程中处理。")
        
        if len(pos_df) > len(merged_df):
            logger.warning(f"{len(pos_df) - len(merged_df)} 个来自位置文件的ID因在表型文件中没有匹配项而被丢弃。")

        if merged_df.empty:
            logger.error("数据合并后为空！没有在位置文件和表型文件中找到任何共同的ID。")
            sys.exit(1)
            
        return merged_df

    except Exception as e:
        logger.error(f"加载和处理数据时发生未知错误: {e}")
        sys.exit(1)

def extract_vcf_regions(input_vcf, output_vcf, regions):
    """
    从VCF提取一个或多个区域的SNP，并精确去重。
    此去重逻辑对于双等位基因SNP和短Indel是稳健的。
    """
    try:
        os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
        unique_variants = {}
        with pysam.VariantFile(input_vcf) as vcf_in:
            for region in regions:
                try:
                    for record in vcf_in.fetch(region['chrom'], region['start'], region['end']):
                        variant_key = (record.chrom, record.pos, record.ref, record.alts)
                        unique_variants[variant_key] = record
                except ValueError:
                    logger.warning(f"在VCF中未找到或无法访问区域 {region['chrom']}:{region['start']}-{region['end']}，已跳过。")

            if not unique_variants:
                logger.warning(f"在指定区域中未提取到任何变异位点。")
                return False

            sorted_variants = sorted(unique_variants.values(), key=lambda r: (r.chrom, r.pos))
            with pysam.VariantFile(output_vcf, 'w', header=vcf_in.header) as vcf_out:
                for record in sorted_variants:
                    vcf_out.write(record)
        return True

    except Exception as e:
        logger.error(f"提取VCF区域时发生严重错误: {e}")
        return False


def process_phenotype(task_info):
    """处理单个分子表型，进行坐标转换和VCF提取"""
    # 任务创建逻辑稳健：'id'和'row'从DataFrame中显式传递
    pheno_id, pheno_row, args = task_info
    chrom = pheno_row['chrom']
    window_bp = args.window_size * 1000

    regions = []
    centers = []
    if args.mode in ['gene', 'tad']:
        centers.append(pheno_row['center'])
    elif args.mode == 'loop':
        centers.append(pheno_row['center1'])
        centers.append(pheno_row['center2'])
        
    for center in centers:
        # 核心修正：将1-based的中心点和窗口，精确转换为pysam所需的0-based半开区间
        start_0based = max(0, int(center) - window_bp - 1)
        end_0based = int(center) + window_bp
        regions.append({'chrom': str(chrom), 'start': start_0based, 'end': end_0based})

    vcf_output_dir = os.path.join(args.output_dir, 'vcf')
    output_vcf = os.path.join(vcf_output_dir, f"{pheno_id}.vcf")
    success = extract_vcf_regions(args.vcf_file, output_vcf, regions)

    if success:
        return {'id': pheno_id, 'success': True, 'vcf_path': output_vcf}
    else:
        if os.path.exists(output_vcf): 
            os.remove(output_vcf)
        return {'id': pheno_id, 'success': False, 'vcf_path': None}


def run_parallel_extraction(pheno_df, args):
    """并行执行VCF区域提取"""
    tasks = [(row['id'], row, args) for _, row in pheno_df.iterrows()]
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        future_to_id = {executor.submit(process_phenotype, task): task[0] for task in tasks}
        for future in tqdm(future_to_id, desc="提取VCF区域", total=len(tasks)):
            try:
                results.append(future.result())
            except Exception as e:
                pheno_id = future_to_id[future]
                logger.error(f"处理ID {pheno_id} 的子进程发生异常: {e}")
                results.append({'id': pheno_id, 'success': False, 'vcf_path': None})
    return results


def vcf_to_plink(task_info):
    """将单个VCF转换为PLINK bfile的工作函数"""
    vcf_file, output_prefix = task_info
    try:
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
        cmd = ['plink', '--vcf', vcf_file, '--make-bed', '--out', output_prefix, '--allow-extra-chr', '--allow-no-sex']
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            logger.error(f"PLINK转换失败 {os.path.basename(vcf_file)}: {result.stderr}")
            return False
        return True
    except Exception as e:
        logger.error(f"运行PLINK命令时发生错误 {os.path.basename(vcf_file)}: {e}")
        return False


def run_parallel_plink_conversion(tasks, threads):
    """并行执行PLINK转换"""
    success_count = 0
    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_task = {executor.submit(vcf_to_plink, task): task for task in tasks}
        for future in tqdm(future_to_task, desc="转换为PLINK", total=len(tasks)):
            try:
                if future.result():
                    success_count += 1
            except Exception as e:
                task_info = future_to_task[future]
                logger.error(f"PLINK转换任务 {task_info[0]} 的子进程发生异常: {e}")
    return success_count


def update_fam_file(fam_file, pheno_id, pheno_lookup_df):
    """更新单个FAM文件，支持浮点数表型"""
    try:
        fam_df = pd.read_csv(fam_file, sep='\s+', header=None, names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'])
        pheno_values = pheno_lookup_df.loc[pheno_id].values
        
        if len(fam_df) != len(pheno_values):
            logger.warning(f"样本数不匹配: {os.path.basename(fam_file)} ({len(fam_df)}) vs 表型文件 ({len(pheno_values)})。跳过。")
            return False
            
        fam_df['PHENO'] = pheno_values
        fam_df['PHENO'] = fam_df['PHENO'].fillna(-9)
        fam_df.to_csv(fam_file, sep=' ', header=False, index=False)
        return True
        
    except KeyError:
        logger.error(f"ID未找到: 在表型文件中未找到ID '{pheno_id}'，无法更新 {os.path.basename(fam_file)}")
        return False
    except Exception as e:
        logger.error(f"更新FAM文件 {os.path.basename(fam_file)} 时出错: {e}")
        return False


def main():
    """主函数，编排整个处理流程"""
    args = parse_arguments()
    
    logger.info("=" * 80)
    logger.info("启动CIWAS多分子表型数据提取流程 (Definitive Production Version)")
    logger.info(f"处理模式: {args.mode.upper()}")
    logger.info("=" * 80)
    
    check_dependencies()
    
    for f in [args.vcf_file, args.pos_file, args.pheno_file]:
        if not os.path.exists(f):
            logger.error(f"输入文件不存在: {f}")
            sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    if not check_vcf_file(args.vcf_file):
        sys.exit(1)

    master_df = load_data(args)
    if args.test_mode:
        master_df = master_df.head(20)
        logger.info(f"--- 测试模式已激活：仅处理前 {len(master_df)} 个分子表型 ---")
        
    logger.info("--- 步骤 1/3: 并行提取VCF区域 ---")
    extraction_results = run_parallel_extraction(master_df, args)
    
    successful_extractions = [r for r in extraction_results if r and r['success']]
    if not successful_extractions:
        logger.error("没有任何区域被成功提取，流程终止。")
        sys.exit(1)

    logger.info("--- 步骤 2/3: 并行转换VCF为PLINK格式 ---")
    plink_tasks = [(r['vcf_path'], os.path.join(args.output_dir, 'plink', r['id'])) for r in successful_extractions]
    successful_plinks = run_parallel_plink_conversion(plink_tasks, args.threads)
    if successful_plinks == 0:
        logger.error("没有任何VCF文件成功转换为PLINK格式，流程终止。")
        sys.exit(1)
    logger.info(f"PLINK转换完成: {successful_plinks}个成功, {len(plink_tasks) - successful_plinks}个失败。")

    logger.info("--- 步骤 3/3: 更新FAM文件中的表型值 ---")
    plink_dir = os.path.join(args.output_dir, 'plink')
    pheno_values_df = pd.read_csv(args.pheno_file, sep='\t')
    pheno_id_col = pheno_values_df.columns[0]
    pheno_lookup_df = pheno_values_df.set_index(pheno_id_col)
    
    fam_files = glob.glob(os.path.join(plink_dir, '*.fam'))
    if not fam_files:
        logger.warning("在PLINK目录中未找到任何 .fam 文件可供更新。")
    else:
        success_count = 0
        for fam_file in tqdm(fam_files, desc="更新FAM文件"):
            pheno_id = os.path.basename(fam_file).replace('.fam', '')
            if update_fam_file(fam_file, pheno_id, pheno_lookup_df):
                success_count += 1
        logger.info(f"FAM文件更新完成: {success_count}个成功, {len(fam_files) - success_count}个失败。")

    logger.info("=" * 80)
    logger.info("CIWAS多分子表型数据提取流程全部完成！")
    logger.info(f"所有输出文件位于: {args.output_dir}")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
