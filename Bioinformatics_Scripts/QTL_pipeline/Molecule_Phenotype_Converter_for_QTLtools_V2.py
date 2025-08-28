#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys
import argparse
import os
import subprocess
import logging
import shutil
import re
from pathlib import Path

# V2.0 新增：检查并导入scipy，用于数据转换
try:
    import scipy.stats as ss
except ImportError:
    print("错误: 本脚本的 'rint' 转换功能需要 'scipy' 库。", file=sys.stderr)
    print("请使用 'pip install scipy' 或 'conda install scipy' 进行安装。", file=sys.stderr)
    sys.exit(1)

# 定义一个合理的块大小，可以根据实际可用内存进行调整
CHUNK_SIZE = 100000

class MoleculePhenotypeConverterV2_1:
    """
    分子表型格式转换工具 (V2.1)
    
    设计目标 (V2.1 更新):
    - 新增 --tensortqtl-mode, 适配 tensorQTL 格式 (#chr, start, end, targetID)。
    - 为 QTLtools 提供标准格式的输入文件 (.bed.gz)。
    - 位置优先: 严格以位置信息文件为基准，过滤掉没有坐标的表型。
    - 高度灵活: 支持无标题行文件，允许用户通过列号指定。
    - 数据转换: 内置多种表型值转换方法 (log2, Z-score, RINT)。
    - 高精度: 使用 float64 保证表型数值的精度。
    - 高效率: 通过分块处理和优化的数据操作，实现较低的内存占用。
    - 鲁棒性: 自动检测分隔符，全面的异常捕获与日志报告。
    - 排序正确性: 使用自然排序 (natural sort) 保证染色体排序符合生物学常规 (chr1, chr2, ..., chr10)。
    - 安全性: 采用安全的外部命令调用方式，规避潜在的shell注入风险。
    """
    
    def __init__(self, output_file, default_strand='+', transform_method='none', rint_stochastic=False):
        self.logger = self._setup_logging(output_file)
        self.separators = ['\t', ',', ';', ' ']
        self.default_strand = default_strand
        self.abnormal_records = []
        self.transform_method = transform_method
        self.rint_stochastic = rint_stochastic
        self._check_external_tools()
        
    def _setup_logging(self, output_file):
        log_file = Path(output_file).with_suffix('.log')
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, mode='w', encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        return logging.getLogger(__name__)

    def write_abnormal_records(self, output_prefix):
        if not self.abnormal_records:
            self.logger.info("未发现异常或缺失数据行，不生成质控报告。")
            return
        abnormal_file_path = Path(f"{output_prefix}.abnormal.txt")
        self.logger.info(f"发现 {len(self.abnormal_records)} 条异常/缺失记录，正在写入质控报告: {abnormal_file_path}")
        try:
            with open(abnormal_file_path, 'w', encoding='utf-8') as f:
                f.write("OriginalLine\tReason\n")
                sorted_records = sorted(self.abnormal_records, key=lambda x: x[0])
                for line_num, reason in sorted_records:
                    f.write(f"{line_num}\t{reason}\n")
        except Exception as e:
            self.logger.error(f"写入异常报告文件失败: {e}")

    def _add_abnormal_record(self, line_num, reason):
        self.abnormal_records.append((line_num, reason))
    
    def _check_external_tools(self):
        self.logger.info("检查外部依赖工具...")
        for tool in ['bgzip', 'tabix']:
            if not shutil.which(tool):
                self.logger.error(f"关键外部工具缺失: {tool}")
                self.logger.error("请确保 htslib (提供 bgzip, tabix) 已安装并在系统的PATH环境变量中。")
                sys.exit(1)
            self.logger.info(f"  [✓] 找到工具: {tool}")
        return True
    
    def _detect_separator(self, filepath, manual_sep=None, no_header=False, max_lines=50):
        if manual_sep:
            self.logger.info(f"用户手动指定分隔符: '{manual_sep}'")
            return manual_sep
        self.logger.info(f"开始自动检测文件 '{Path(filepath).name}' 的分隔符...")
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                start_line = 0 if no_header else 1
                sample_lines = [line for i, line in enumerate(f) if i >= start_line and len(line.strip()) > 0][:max_lines]
            
            non_comment_lines = [line for line in sample_lines if not line.startswith('#')]
            if not non_comment_lines:
                self.logger.warning("文件头部没有可用于检测的非注释行，将使用默认制表符。")
                return '\t'
            
            results = []
            for sep in self.separators:
                counts = [len(line.split(sep)) for line in non_comment_lines]
                if not counts: continue
                mean_count = np.mean(counts)
                std_dev = np.std(counts)
                if mean_count > 1.5:
                    results.append({'sep': sep, 'std': std_dev, 'mean': mean_count})
            
            if not results:
                self.logger.warning("未检测到有效的多列分隔符，将使用默认制表符。")
                return '\t'
            
            best_result = sorted(results, key=lambda x: (x['std'], -x['mean']))[0]
            best_sep = best_result['sep']
            
            self.logger.info(f"检测到最可能的分隔符为: '{best_sep}' ({'制表符' if best_sep == chr(9) else best_sep}) - 列数稳定度(std): {best_result['std']:.2f}")
            return best_sep
        except Exception as e:
            self.logger.warning(f"分隔符检测失败: {e}, 将使用默认的制表符。")
            return '\t'

    def load_position_info(self, args, manual_sep=None):
        self.logger.info(f"开始加载位置信息文件: {args.position_file}")
        self.logger.info("注意：最终输出的分子表型将严格以此文件中的ID为基准。")
        try:
            sep = self._detect_separator(args.position_file, manual_sep, args.no_header_pos)
            header_arg = 'infer' if not args.no_header_pos else None
            use_cols_map, display_cols = self._get_col_map(args, 'pos', args.no_header_pos)
            self.logger.info(f"从位置文件读取列: {display_cols}")
            df = pd.read_csv(args.position_file, sep=sep, header=header_arg, usecols=use_cols_map.keys(), dtype=str, na_filter=False, comment='#')
            df.rename(columns=use_cols_map, inplace=True)
            position_info = {}
            for _, row in df.iterrows():
                gene_id = row['gene_id'].strip()
                if not gene_id: continue
                info = {'chr': row['chr'].strip()}
                if args.loop_mode:
                    coords = [self._safe_convert_single_coord(row[c]) for c in ['a1_start', 'a1_end', 'a2_start', 'a2_end']]
                    if any(v is None for v in coords): continue
                    info.update({'a1_start': coords[0], 'a1_end': coords[1], 'a2_start': coords[2], 'a2_end': coords[3]})
                else:
                    start_val = self._safe_convert_single_coord(row['start'])
                    if start_val is None: continue
                    info['start'] = start_val
                    if not args.single_base_mode or args.recalculate_pos:
                        end_val = self._safe_convert_single_coord(row.get('end'))
                        if end_val is None: continue
                        if start_val > end_val: start_val, end_val = end_val, start_val
                        info['end'] = end_val
                position_info[gene_id] = info
            self.logger.info(f"位置信息加载完成: 共 {len(position_info)} 个基准基因ID")
            if not position_info:
                self.logger.error("位置信息文件为空或未能解析出任何有效行，程序终止。")
                sys.exit(1)
            return position_info
        except Exception as e:
            self.logger.error(f"位置信息文件加载失败: {e}", exc_info=True)
            return None

    def load_strand_file(self, filepath, gene_col_arg, strand_col_arg, no_header, manual_sep=None):
        self.logger.info(f"开始加载链信息文件: {filepath}")
        try:
            sep = self._detect_separator(filepath, manual_sep, no_header)
            header_arg = 'infer' if not no_header else None
            gene_col = int(gene_col_arg) - 1 if no_header else gene_col_arg
            strand_col = int(strand_col_arg) - 1 if no_header else strand_col_arg
            self.logger.info(f"从链信息文件读取列: gene_id='{gene_col_arg}', strand='{strand_col_arg}'")
            df = pd.read_csv(filepath, sep=sep, header=header_arg, usecols=[gene_col, strand_col], dtype=str, na_filter=False, comment='#')
            df.rename(columns={gene_col: 'gene_id', strand_col: 'strand'}, inplace=True)
            df['normalized_strand'] = self._normalize_strand_vectorized(df['strand'])
            strand_dict = dict(zip(df['gene_id'].str.strip(), df['normalized_strand']))
            self.logger.info(f"链信息加载完成: 共 {len(strand_dict)} 个基因")
            return strand_dict
        except Exception as e:
            self.logger.error(f"链信息文件加载失败: {e}", exc_info=True)
            return None
    
    def load_gff3_strand(self, filepath, gene_attr='ID'):
        self.logger.info(f"开始从GFF3文件加载链信息: {filepath}")
        strand_info = {}
        attr_regex = re.compile(f'{re.escape(gene_attr)}=([^;]+)')
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('#') or not line.strip(): continue
                    parts = line.strip().split('\t')
                    if len(parts) != 9: continue
                    attributes_str = parts[8]
                    match = attr_regex.search(attributes_str)
                    if match:
                        gene_id = match.group(1).strip()
                        strand = parts[6]
                        if strand not in ['+', '-', '.']: strand = self.default_strand
                        strand_info[gene_id] = strand
            self.logger.info(f"GFF3加载完成: 获得 {len(strand_info)} 个基因的链信息")
            return strand_info
        except Exception as e:
            self.logger.error(f"GFF3文件加载失败: {e}", exc_info=True)
            return None

    def _chromosome_sort_key(self, chr_str):
        chr_str = str(chr_str).strip()
        
        cotton_match = re.match(r'^HC04_([AD])(\d+)$', chr_str)
        if cotton_match:
            letter = cotton_match.group(1)
            number = int(cotton_match.group(2))
            return (0 if letter == 'A' else 1, number)
        
        chr_match = re.match(r'^(?:chr)?(.+)$', chr_str, re.IGNORECASE)
        if chr_match:
            chr_content = chr_match.group(1)
            try:
                return (2, int(chr_content)) 
            except ValueError:
                return (3, chr_content) 
        
        return (4, chr_str)

    def _process_chunk(self, chunk_df, chunk_num, args, col_map, sample_cols, strand_info, position_info):
        work_df = chunk_df.copy()
        work_df.rename(columns=col_map, inplace=True)
        work_df['original_line'] = chunk_num * CHUNK_SIZE + work_df.index + (1 if args.no_header_main else 2)
        work_df['gene_id_clean'] = work_df['gene_id'].astype(str).str.strip()
        if args.position_file:
            master_id_set = set(position_info.keys())
            mask_in_master = work_df['gene_id_clean'].isin(master_id_set)
            work_df = work_df[mask_in_master]
            if work_df.empty: return None
        missing_pid_mask = work_df['gene_id_clean'] == ''
        for line_num in work_df.loc[missing_pid_mask, 'original_line']: self._add_abnormal_record(line_num, "缺失基因ID (gene_id)")
        work_df = work_df[~missing_pid_mask]
        if work_df.empty: return None
        if args.position_file:
            pos_map = work_df['gene_id_clean'].map(position_info)
            pos_df = pd.DataFrame(pos_map.dropna().tolist(), index=pos_map.dropna().index)
            work_df = work_df.join(pos_df)
        else:
            work_df['chr'] = work_df['chr'].astype(str).str.strip()
            if args.loop_mode:
                for col in ['a1_start', 'a1_end', 'a2_start', 'a2_end']: work_df[col] = self._safe_convert_coord(work_df[col])
            else:
                work_df['start'] = self._safe_convert_coord(work_df['start'])
                if not args.single_base_mode or args.recalculate_pos: work_df['end'] = self._safe_convert_coord(work_df['end'])
        coord_cols_to_check = {'chr': '缺失染色体信息'}
        if args.loop_mode: coord_cols_to_check.update({ 'a1_start': '无效anchor1_start', 'a1_end': '无效anchor1_end', 'a2_start': '无效anchor2_start', 'a2_end': '无效anchor2_end' })
        else:
            coord_cols_to_check['start'] = '无效start坐标'
            if not args.single_base_mode or args.recalculate_pos: coord_cols_to_check['end'] = '无效end坐标'
        for col, reason in coord_cols_to_check.items():
            if col not in work_df.columns: continue
            mask = work_df[col].isna()
            for line_num in work_df.loc[mask, 'original_line']: self._add_abnormal_record(line_num, reason)
            work_df.dropna(subset=[col], inplace=True)
        if work_df.empty: return None
        if args.loop_mode:
            if args.input_system == 'bed':
                a1s, a1e, a2s, a2e = work_df['a1_start'].astype(int)+1, work_df['a1_end'].astype(int), work_df['a2_start'].astype(int)+1, work_df['a2_end'].astype(int)
            else:
                a1s, a1e, a2s, a2e = work_df['a1_start'].astype(int), work_df['a1_end'].astype(int), work_df['a2_start'].astype(int), work_df['a2_end'].astype(int)
            mid1, mid2 = (a1s + a1e) // 2, (a2s + a2e) // 2
            final_mid = (mid1 + mid2) // 2
            bed_start, bed_end = final_mid - 1, final_mid
        else:
            start_raw = work_df['start'].astype(int)
            use_end_col = (not args.single_base_mode or args.recalculate_pos)
            if use_end_col: end_raw = work_df['end'].astype(int)
            if args.input_system == 'bed':
                start_1based = start_raw + 1
                if use_end_col: end_1based = end_raw
            else:
                start_1based = start_raw
                if use_end_col: end_1based = end_raw
            if use_end_col and (start_1based > end_1based).any():
                swap_mask = start_1based > end_1based
                start_1based[swap_mask], end_1based[swap_mask] = end_1based[swap_mask], start_1based[swap_mask]
            if args.recalculate_pos:
                midpoint = (start_1based + end_1based) // 2
                bed_start, bed_end = midpoint - 1, midpoint
            elif args.single_base_mode:
                bed_start, bed_end = start_1based - 1, start_1based
            else:
                bed_start, bed_end = start_1based - 1, end_1based
        work_df['bed_start'], work_df['bed_end'] = bed_start, bed_end
        
        # ### V2.1 变更: 仅在非tensorQTL模式下处理链信息 ###
        # 注意: 即使在tensorQTL模式下，此代码块也会运行，但其结果'strand'列不会被使用
        if args.use_strand:
            if strand_info: work_df['strand'] = work_df['gene_id_clean'].map(strand_info).fillna(self.default_strand)
            else: work_df['strand'] = self.default_strand
        else:
            work_df['strand'] = self.default_strand
            
        negative_mask = work_df['bed_start'] < 0
        for line_num in work_df.loc[negative_mask, 'original_line']: self._add_abnormal_record(line_num, "计算后的起始坐标为负数，已丢弃")
        work_df = work_df[~negative_mask]
        if work_df.empty: return None
        work_df = self._process_and_transform_phenotypes(work_df, sample_cols)
        
        ### V2.1 更新：根据 tensortqtl-mode 动态构建输出列 ###
        if args.tensortqtl_mode:
            # tensorQTL 模式: #chr, start, end, targetID
            base_cols_dict = {
                '#chr': work_df['chr'],
                'start': work_df['bed_start'].astype('Int64'),
                'end': work_df['bed_end'].astype('Int64'),
                'targetID': work_df['gene_id_clean']
            }
        else:
            # 默认 QTLtools 模式: #chr, start, end, pid, gid, strand
            base_cols_dict = {
                '#chr': work_df['chr'],
                'start': work_df['bed_start'].astype('Int64'),
                'end': work_df['bed_end'].astype('Int64'),
                'pid': work_df['gene_id_clean'],
                'gid': work_df['gene_id_clean'],
                'strand': work_df['strand']
            }

        base_cols_df = pd.DataFrame(base_cols_dict)
        phenotype_df = work_df[sample_cols]
        return pd.concat([base_cols_df, phenotype_df], axis=1)

    def process_file_in_chunks(self, args, strand_info, position_info, temp_output_file):
        self.logger.info("开始流式处理输入文件...")
        try:
            sep = self._detect_separator(args.input, args.separator, args.no_header_main)
            header_arg = 'infer' if not args.no_header_main else None
            col_map, _ = self._get_col_map(args, 'main', args.no_header_main)
            all_file_cols = pd.read_csv(args.input, sep=sep, header=header_arg, nrows=0, comment='#').columns
            meta_cols = list(col_map.keys())
            sample_cols = [c for c in all_file_cols if c not in meta_cols]
            if args.main_strand_col:
                if args.no_header_main:
                    strand_col_idx = int(args.main_strand_col) - 1
                    if strand_col_idx < len(all_file_cols):
                        strand_col_name = all_file_cols[strand_col_idx]
                        if strand_col_name in sample_cols:
                            sample_cols.remove(strand_col_name)
                            self.logger.info(f"已从样本列中排除指定的链信息列: (列号 {args.main_strand_col})")
                else:
                    if args.main_strand_col in sample_cols:
                        sample_cols.remove(args.main_strand_col)
                        self.logger.info(f"已从样本列中排除指定的链信息列: '{args.main_strand_col}'")
            self.logger.info(f"检测到 {len(sample_cols)} 个样本列。")
            if not sample_cols:
                self.logger.error("未能检测到任何样本列，请检查列名或列号指定是否正确。")
                return False
            reader = pd.read_csv(args.input, sep=sep, header=header_arg, dtype=str, na_filter=False, chunksize=CHUNK_SIZE, comment='#')
            is_first_chunk, total_rows_processed = True, 0
            for i, chunk_df in enumerate(reader):
                self.logger.info(f"正在处理数据块 {i+1}...")
                chunk_col_map = {**col_map, **{s: s for s in sample_cols}}
                processed_chunk = self._process_chunk(chunk_df, i, args, chunk_col_map, sample_cols, strand_info, position_info)
                if processed_chunk is not None and not processed_chunk.empty:
                    processed_chunk.to_csv(temp_output_file, sep='\t', index=False, header=is_first_chunk, mode='a' if not is_first_chunk else 'w', na_rep='NA')
                    total_rows_processed += len(processed_chunk)
                    if is_first_chunk: is_first_chunk = False
            self.logger.info(f"文件处理完成，共处理并写入 {total_rows_processed} 行有效数据到临时文件。")
            if args.position_file and total_rows_processed == 0: self.logger.warning("主数据文件中的基因ID与位置信息文件中的ID没有任何交集。")
            return True
        except Exception as e:
            self.logger.error(f"文件流式处理失败: {e}", exc_info=True)
            return False

    def sort_and_compress(self, unsorted_path, final_gz_path):
        if not Path(unsorted_path).exists() or os.path.getsize(unsorted_path) == 0:
            self.logger.warning("没有有效数据行输出，跳过排序和压缩步骤。")
            return True
        sorted_path = Path(f"{unsorted_path}.sorted")
        try:
            self.logger.info(f"开始对文件进行自定义排序: {unsorted_path}")
            
            df = pd.read_csv(unsorted_path, sep='\t', dtype={'start': 'Int64', 'end': 'Int64'})
            
            chr_col_name = df.columns[0]
            self.logger.info(f"使用列 '{chr_col_name}' 进行染色体排序。")
            df['_sort_key'] = df[chr_col_name].apply(self._chromosome_sort_key)
            
            df_sorted = df.sort_values(['_sort_key', 'start', 'end'])
            df_sorted = df_sorted.drop('_sort_key', axis=1)
            
            df_sorted.to_csv(sorted_path, sep='\t', index=False, na_rep='NA')
            
            self.logger.info("文件排序完成。")
            self.logger.info(f"开始使用bgzip压缩文件: {sorted_path} -> {final_gz_path}")
            with open(sorted_path, 'rb') as f_in, open(final_gz_path, 'wb') as f_out:
                subprocess.run(['bgzip', '-c'], stdin=f_in, stdout=f_out, check=True)
            self.logger.info(f"文件已压缩: {final_gz_path}")
            subprocess.run(['tabix', '-p', 'bed', str(final_gz_path)], check=True, capture_output=True, text=True)
            self.logger.info(f"索引已创建: {final_gz_path}.tbi")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"压缩或索引步骤失败. 命令: '{e.cmd}'. 标准错误: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"排序、压缩或索引过程中发生未知错误: {e}", exc_info=True)
            return False
        finally:
            if sorted_path.exists(): os.remove(sorted_path)
            if Path(unsorted_path).exists(): os.remove(unsorted_path)
            self.logger.info("已清理所有临时文件。")

    def _get_col_map(self, args, file_type, no_header):
        prefix = 'pos_' if file_type == 'pos' else ''
        cols = {}
        if file_type == 'main' and args.embedded_pos:
            cols = {'gene_id': args.gene_id, 'chr': args.chr}
            if args.loop_mode: cols.update({'a1_start': args.anchor1_start, 'a1_end': args.anchor1_end, 'a2_start': args.anchor2_start, 'a2_end': args.anchor2_end})
            else:
                cols['start'] = args.start
                if not args.single_base_mode or args.recalculate_pos: cols['end'] = args.end
        elif file_type == 'main':
            cols = {'gene_id': args.gene_id}
        elif file_type == 'pos':
            cols = {'gene_id': args.pos_gene_col, 'chr': args.pos_chr_col}
            if args.loop_mode: cols.update({'a1_start': args.pos_anchor1_start_col, 'a1_end': args.pos_anchor1_end_col, 'a2_start': args.pos_anchor2_start_col, 'a2_end': args.pos_anchor2_end_col})
            else:
                cols['start'] = args.pos_start_col
                if not args.single_base_mode or args.recalculate_pos: cols['end'] = args.pos_end_col
        cols = {k: v for k, v in cols.items() if v is not None}
        col_map = {int(v) - 1: k for k, v in cols.items()} if no_header else {v: k for k, v in cols.items()}
        display_cols = ", ".join([f"{k}='{v}'" for k, v in cols.items()])
        return col_map, display_cols

    def _safe_convert_coord(self, series):
        return pd.to_numeric(series, errors='coerce').astype('Int64')

    def _safe_convert_single_coord(self, value):
        try:
            if pd.isna(value) or str(value).strip().lower() in ['', 'nan', 'null', 'none', 'na']: return None
            coord = int(float(str(value)))
            return coord if coord >= 0 else None
        except (ValueError, TypeError): return None

    def _normalize_strand_vectorized(self, series):
        strand_map = { '+': '+', 'PLUS': '+', '1': '+', '-': '-', 'MINUS': '-', '-1': '-', '.': '.' }
        return series.astype(str).str.strip().str.upper().map(strand_map).fillna(self.default_strand)

    def _rank_to_normal(self, rank, c, n):
        return ss.norm.ppf((rank - c) / (n - 2 * c + 1))

    def _rank_INT(self, series, c=0.5, stochastic=False):
        assert isinstance(series, pd.Series)
        assert isinstance(c, float)
        assert isinstance(stochastic, bool)
        np.random.seed(123)
        orig_idx = series.index
        series = series.loc[~pd.isnull(series)]
        if series.empty: return pd.Series(np.nan, index=orig_idx)
        if stochastic:
            series = series.loc[np.random.permutation(series.index)]
            rank = ss.rankdata(series, method="ordinal")
        else:
            rank = ss.rankdata(series, method="average")
        rank = pd.Series(rank, index=series.index)
        transformed = rank.apply(self._rank_to_normal, c=c, n=len(rank))
        return transformed.reindex(orig_idx)

    def _process_and_transform_phenotypes(self, df, sample_cols):
        df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors='coerce', axis=1).astype('float64')
        for col in sample_cols:
            for line_num in df.loc[df[col].isna(), 'original_line']: self._add_abnormal_record(line_num, f"样本 '{col}' 的表型数据缺失或无效")
        if self.transform_method == 'log2':
            df[sample_cols] = np.log2(1 + df[sample_cols])
        elif self.transform_method == 'zscore':
            mean, std = df[sample_cols].mean(axis=1), df[sample_cols].std(axis=1)
            std[std == 0] = 1 
            df[sample_cols] = df[sample_cols].sub(mean, axis=0).div(std, axis=0)
        elif self.transform_method == 'rint':
            df[sample_cols] = df[sample_cols].apply(self._rank_INT, axis=1, stochastic=self.rint_stochastic)
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        return df

def parse_args():
    ### V2.1 更新: 修改描述信息 ###
    parser = argparse.ArgumentParser(description='分子表型格式转换为QTLtools/tensorQTL标准BED格式的工具 (V2.1)', formatter_class=argparse.RawTextHelpFormatter)
    req_group = parser.add_argument_group('输入/输出 (必选)')
    req_group.add_argument('-i', '--input', required=True, help='输入的分子表型数据文件')
    req_group.add_argument('-o', '--output', required=True, help='最终输出文件的前缀')
    req_group.add_argument('--gene-id', required=True, help='表型数据中的基因/分子ID列名或列号(1-based)')
    
    pos_group = parser.add_argument_group('位置信息来源 (必须选择一种)')
    pos_exclusive = pos_group.add_mutually_exclusive_group(required=True)
    pos_exclusive.add_argument('--embedded-pos', action='store_true', help='位置信息在表型数据文件中')
    pos_exclusive.add_argument('--position-file', help='单独的位置信息文件路径 (以此文件为基准)')
    
    col_group = parser.add_argument_group('列名/列号指定 (根据模式选择)')
    col_group.add_argument('--chr', help='染色体列名/列号 (embedded-pos 时必需)')
    col_group.add_argument('--start', help='起始位置列名/列号 (embedded-pos 时必需)')
    col_group.add_argument('--end', help='结束位置列名/列号 (embedded-pos, 非single-base模式时必需)')
    col_group.add_argument('--pos-gene-col', help='位置文件中基因ID列名/列号')
    col_group.add_argument('--pos-chr-col', help='位置文件中染色体列名/列号')
    col_group.add_argument('--pos-start-col', help='位置文件中起始位置列名/列号')
    col_group.add_argument('--pos-end-col', help='位置文件中结束位置列名/列号')
    col_group.add_argument('--anchor1-start', help='(Loop) 染色质环 anchor1 起始列名/列号')
    col_group.add_argument('--anchor1-end', help='(Loop) 染色质环 anchor1 结束列名/列号')
    col_group.add_argument('--anchor2-start', help='(Loop) 染色质环 anchor2 起始列名/列号')
    col_group.add_argument('--anchor2-end', help='(Loop) 染色质环 anchor2 结束列名/列号')
    col_group.add_argument('--pos-anchor1-start-col', help='(Loop) 位置文件中 anchor1 起始列名/列号')
    col_group.add_argument('--pos-anchor1-end-col', help='(Loop) 位置文件中 anchor1 结束列名/列号')
    col_group.add_argument('--pos-anchor2-start-col', help='(Loop) 位置文件中 anchor2 起始列名/列号')
    col_group.add_argument('--pos-anchor2-end-col', help='(Loop) 位置文件中 anchor2 结束列名/列号')
    
    header_group = parser.add_argument_group('标题行与特殊列')
    header_group.add_argument('--no-header-main', action='store_true', help='主数据文件无标题行。列名参数需提供列号(1-based)。')
    header_group.add_argument('--no-header-pos', action='store_true', help='位置信息文件无标题行。列名参数需提供列号(1-based)。')
    header_group.add_argument('--no-header-strand', action='store_true', help='链信息文件无标题行。列名参数需提供列号(1-based)。')
    header_group.add_argument('--main-strand-col', help='(可选) 主数据文件中包含的链信息列名/列号, 指定后此列将被忽略。')
    
    transform_group = parser.add_argument_group('表型数据转换')
    transform_group.add_argument('--transform', choices=['none', 'log2', 'zscore', 'rint'], default='none', help="对每个表型(行)的数据进行转换 (默认: none)")
    transform_group.add_argument('--rint-stochastic', action='store_true', help="当使用 'rint' 转换时, 启用随机排序处理ties (默认不启用)")
    
    ### V2.1 更新：添加新的输出格式控制参数组 ###
    output_format_group = parser.add_argument_group('输出格式控制 (V2.1 新增)')
    output_format_group.add_argument('--tensortqtl-mode', action='store_true', help='启用 tensorQTL 兼容模式 (输出列: #chr, start, end, targetID)')

    behavior_group = parser.add_argument_group('格式与行为控制 (通用)')
    behavior_group.add_argument('--input-system', choices=['bed', '1-based'], default='bed', help="输入坐标系统 (默认: bed)")
    behavior_group.add_argument('--default-strand', choices=['+', '-', '.'], default='+', help='默认链信息 (默认: +)')
    behavior_group.add_argument('--separator', help='手动指定输入文件的列分隔符')
    mode_group = behavior_group.add_mutually_exclusive_group()
    mode_group.add_argument('--loop-mode', action='store_true', help='染色质环模式')
    mode_group.add_argument('--recalculate-pos', action='store_true', help='中点模式')
    mode_group.add_argument('--single-base-mode', action='store_true', help='单碱基模式')
    
    strand_group = parser.add_argument_group('链信息 (可选, 非tensorQTL模式下生效)')
    strand_group.add_argument('--use-strand', action='store_true', help='启用链信息处理')
    strand_group.add_argument('--gff3', help='GFF3文件路径')
    strand_group.add_argument('--gff3-gene-attr', default='ID', help='GFF3中匹配基因ID的字段 (默认: ID)')
    strand_group.add_argument('--strand-file', help='包含链信息的独立文件路径')
    strand_group.add_argument('--strand-gene-col', help='链信息文件中基因ID列名/列号')
    strand_group.add_argument('--strand-col', help='链信息文件中链信息列名/列号')
    return parser.parse_args()

def validate_args(args):
    if not Path(args.input).exists(): return False, f"输入文件不存在: {args.input}"
    def is_int(s):
        if s is None: return True
        try: int(s); return True
        except (ValueError, TypeError): return False

    if args.no_header_main:
        main_cols = {'gene-id': args.gene_id, 'main-strand-col': args.main_strand_col}
        if args.embedded_pos:
            main_cols.update({'chr': args.chr})
            if args.loop_mode:
                main_cols.update({'anchor1-start': args.anchor1_start, 'anchor1-end': args.anchor1_end, 'anchor2-start': args.anchor2_start, 'anchor2-end': args.anchor2_end})
            else:
                main_cols.update({'start': args.start})
                if not args.single_base_mode or args.recalculate_pos: main_cols['end'] = args.end
        if not all(is_int(v) for v in main_cols.values()):
            return False, "使用 --no-header-main 时,所有主文件相关列参数必须是数字列号"

    if args.position_file:
        if not Path(args.position_file).exists(): return False, f"位置信息文件不存在: {args.position_file}"
        pos_cols = {'pos_gene_col': args.pos_gene_col, 'pos_chr_col': args.pos_chr_col}
        if args.loop_mode:
            pos_cols.update({'pos_anchor1_start_col': args.pos_anchor1_start_col, 'pos_anchor1_end_col': args.pos_anchor1_end_col, 'pos_anchor2_start_col': args.pos_anchor2_start_col, 'pos_anchor2_end_col': args.pos_anchor2_end_col})
        else:
            pos_cols.update({'pos_start_col': args.pos_start_col})
            if not args.single_base_mode or args.recalculate_pos:
                pos_cols['pos_end_col'] = args.pos_end_col
        if not all(pos_cols.values()):
            return False, "在 --position-file 模式下, 必须指定所有必需的 pos_* 列"
        if args.no_header_pos and not all(is_int(v) for v in pos_cols.values()):
            return False, "使用 --no-header-pos 时,所有位置文件相关列参数必须是数字列号"
    
    if args.use_strand and args.strand_file:
        if not Path(args.strand_file).exists(): return False, f"链信息文件不存在: {args.strand_file}"
        strand_cols = {'strand_gene_col': args.strand_gene_col, 'strand_col': args.strand_col}
        if not all(strand_cols.values()): return False, "使用 --strand-file 时, 必须同时指定 --strand-gene-col 和 --strand-col"
        if args.no_header_strand and not all(is_int(v) for v in strand_cols.values()):
            return False, "使用 --no-header-strand 时,所有链信息文件相关列参数必须是数字列号"
    
    return True, ""

def main():
    args = parse_args()
    valid, error_msg = validate_args(args)
    if not valid:
        print(f"参数错误: {error_msg}", file=sys.stderr)
        sys.exit(1)
    
    output_prefix, final_gz_path = args.output, Path(f"{args.output}.bed.gz")
    temp_unsorted_path = Path(f"{output_prefix}.unsorted.tmp.bed")
    temp_unsorted_path.parent.mkdir(parents=True, exist_ok=True)
    
    ### V2.1 更新: 实例化新类 ###
    converter = MoleculePhenotypeConverterV2_1(output_prefix, args.default_strand, args.transform, args.rint_stochastic)
    converter.logger.info("分子表型格式转换工具启动 (V2.1)")
    if args.tensortqtl_mode:
        converter.logger.info("已启用 tensorQTL 兼容模式。")
    converter.logger.info(f"运行参数: {' '.join(sys.argv[1:])}")
    
    try:
        position_info, strand_info = None, None
        if args.position_file:
            position_info = converter.load_position_info(args, args.separator)
            if position_info is None: sys.exit(1)
        
        # 在 tensorQTL 模式下，不加载任何链信息
        if not args.tensortqtl_mode and args.use_strand:
            if args.gff3: strand_info = converter.load_gff3_strand(args.gff3, args.gff3_gene_attr)
            elif args.strand_file: strand_info = converter.load_strand_file(args.strand_file, args.strand_gene_col, args.strand_col, args.no_header_strand, args.separator)
        
        if not converter.process_file_in_chunks(args, strand_info, position_info, temp_unsorted_path):
            converter.logger.error("文件处理失败，程序终止。")
        elif not converter.sort_and_compress(temp_unsorted_path, final_gz_path):
            converter.logger.error("转换完成但后续处理（排序/压缩/索引）失败。")
        else:
            converter.logger.info(f"转换成功完成！最终输出文件: {final_gz_path} and {final_gz_path}.tbi")
        converter.write_abnormal_records(output_prefix)
    except Exception as e:
        converter.logger.error(f"程序执行期间发生致命错误: {e}", exc_info=True)
        if temp_unsorted_path.exists(): os.remove(temp_unsorted_path)
        sys.exit(1)

if __name__ == "__main__":
    main()