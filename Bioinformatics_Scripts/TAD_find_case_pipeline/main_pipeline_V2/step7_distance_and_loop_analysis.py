#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP-TAD-Gene Analysis Core (v5.5) - Statistically Rigorous Quadruplet Strategy
Author: Based on user requirements
Date: 2025
Version: 5.5.0

严谨的统计分析流程：
1. 分支1：基因在TAD临界距离内 (三元组分析)。
2. 分支2：基因在TAD临界距离外 (四元组分析)。
   策略 (v5.2): 查找所有匹配的Loop。对每一个 (SNP, TAD, Gene, Loop) 组合（四元组）
                 进行独立和完整的统计分析，并为每一个通过的四元组输出单独一行。
   v5.3 修复: 修正了 check_overlap 的边界排斥问题 (<=)。
   v5.4 优化: 修复了 TypeError 并重写了日志报告。
    v5.5 修复: (关键Bug修复) 修正了 find_all_overlapping_loops... 函数中键名大小写不匹配的问题。
               此Bug曾导致 Loop_ID, Loop_Match_Type 等所有Loop物理信息在输出CSV中显示为NaN。
"""

import pandas as pd
import numpy as np
import sys
import os
from argparse import ArgumentParser
from multiprocessing import Pool
from pathlib import Path
from scipy.stats import mannwhitneyu, chi2_contingency, fisher_exact
from tqdm import tqdm
import warnings

# 抑制统计计算中的数值警告
warnings.filterwarnings("ignore", category=RuntimeWarning,
                        message=r"divide by zero encountered|invalid value encountered")
warnings.filterwarnings("ignore", category=UserWarning)

# --------------------------------------------------------------------------
# --- 命令行参数解析 ---
# --------------------------------------------------------------------------
def parse_args():
    parser = ArgumentParser(description="基于TAD-基因距离和方向性规则的严谨统计分析 (v5.5 - 四元组策略 + Bug修复)")
    
    # 必选参数
    parser.add_argument("triplets_with_positions", help="包含位置信息的三元组文件")
    parser.add_argument("snp_data", help="SNP基因型数据文件")
    parser.add_argument("gene_exp", help="基因表达数据文件")
    parser.add_argument("tad_is", help="TAD IS数据文件")
    parser.add_argument("tad_01", help="TAD 0-1矩阵文件")
    parser.add_argument("loop_pos", help="Loop位置信息文件")
    parser.add_argument("loop_strength", help="Loop强度数据文件")
    parser.add_argument("loop_01", help="Loop 0-1矩阵文件")
    parser.add_argument("num_cores", type=int, help="CPU核心数")
    
    # 可选参数
    parser.add_argument("--distance_threshold", type=int, default=80,
                        help="分支临界距离 (kb，默认80)")
    parser.add_argument("--max_distance", type=int, default=250,
                        help="分支2最远距离 (kb，默认250)")
    parser.add_argument("--gene_anchor_distance", type=int, default=40,
                        help="基因anchor扩展距离 (kb，默认40)")
    parser.add_argument("--loop_p_threshold", type=float, default=0.05,
                        help="Loop显著性p值阈值 (默认0.05, 独立应用于每个Loop)")
    parser.add_argument("--output_dir", type=str, default="analysis_output_v5_5_quadruplets",
                        help="输出文件夹路径 (默认: analysis_output_v5_5_quadruplets)")
    
    return parser.parse_args()

# --------------------------------------------------------------------------
# --- 全局配置 ---
# --------------------------------------------------------------------------
class Config:
    """全局配置参数"""
    EXPECTED_GENOTYPES = [0, 2]  # 期望的基因型值
    
    # 将在main函数中设置
    DISTANCE_THRESHOLD_KB = None
    MAX_DISTANCE_KB = None
    GENE_ANCHOR_DISTANCE_KB = None
    P_THRESHOLD_LOOP = None

# 多进程全局变量
_G = {}

# --------------------------------------------------------------------------
# --- 核心计算函数 (v5.3 修复版) ---
# --------------------------------------------------------------------------
def calculate_distance(start1, end1, start2, end2):
    """计算两个区间的最短距离"""
    if max(start1, start2) < min(end1, end2): 
        return 0  # 重叠
    return min(abs(start1 - end2), abs(end1 - start2))

def calculate_gene_tad_distance_and_relation(gene_start, gene_end, tad_start, tad_end):
    """计算基因与TAD的距离和位置关系"""
    distance = calculate_distance(gene_start, gene_end, tad_start, tad_end)
    
    if distance == 0:
        relation = "重叠"
    elif gene_end < tad_start:
        relation = "TAD上游"
    elif gene_start > tad_end:
        relation = "TAD下游"
    else:
        relation = "部分重叠" 
    
    return distance, relation

def check_overlap(start1, end1, start2, end2):
    """
    检查两个区间是否重叠 (v5.3 MODIFIED: 端点包容)
    """
    return max(start1, start2) <= min(end1, end2)

def determine_gene_tad_spatial_relationship(gene_start, gene_end, tad_start, tad_end):
    """确定基因相对于TAD的空间关系"""
    if gene_end < tad_start:
        return "upstream"
    elif gene_start > tad_end:
        return "downstream"
    else:
        return "overlap" 

# --------------------------------------------------------------------------
# --- 严谨的统计检验函数 ---
# --------------------------------------------------------------------------
def perform_mannwhitney_test(data, value_col, group_col):
    """
    严谨的Mann-Whitney U检验实现
    """
    try:
        sub = data[[value_col, group_col]].dropna()
        if len(sub) < 3:
            return {"p_value": np.nan, "stat": np.nan, "group0_size": 0, "group2_size": 0, "group0_median": np.nan, "group2_median": np.nan, "error": "Insufficient total sample size"}
        if sub[group_col].nunique() < 2:
            return {"p_value": np.nan, "stat": np.nan, "group0_size": 0, "group2_size": 0, "group0_median": np.nan, "group2_median": np.nan, "error": "Less than 2 groups available"}
        
        groups = sorted(sub[group_col].unique())
        if len(groups) != 2 or groups != [0, 2]:
            return {"p_value": np.nan, "stat": np.nan, "group0_size": 0, "group2_size": 0, "group0_median": np.nan, "group2_median": np.nan, "error": f"Unexpected groups: {groups}"}
        
        g0 = sub[sub[group_col] == 0][value_col]
        g2 = sub[sub[group_col] == 2][value_col]
        
        if len(g0) == 0 or len(g2) == 0:
            return {"p_value": np.nan, "stat": np.nan, "group0_size": len(g0), "group2_size": len(g2), "group0_median": float(g0.median()) if len(g0) > 0 else np.nan, "group2_median": float(g2.median()) if len(g2) > 0 else np.nan, "error": "Empty group(s)"}
        
        stat, p_value = mannwhitneyu(g0, g2, alternative='two-sided')
        
        return {
            "p_value": float(p_value),
            "stat": float(stat),
            "group0_size": len(g0),
            "group2_size": len(g2),
            "group0_median": float(g0.median()),
            "group2_median": float(g2.median()),
        }
    except Exception as e:
        return {"p_value": np.nan, "stat": np.nan, "group0_size": 0, "group2_size": 0, "group0_median": np.nan, "group2_median": np.nan, "error": f"Mann-Whitney test failed: {str(e)}"}

def _prepare_2x2_contingency_table(sub, categorical_col, group_col):
    """预处理2x2列联表并检查退化情况"""
    ct = pd.crosstab(sub[categorical_col], sub[group_col])
    ct = ct.reindex(index=[0, 1], columns=[0, 2], fill_value=0)
    if (ct.sum(axis=1) == 0).any() or (ct.sum(axis=0) == 0).any():
        return None
    return ct

def perform_independence_test(data, categorical_col, group_col):
    """
    严谨的2x2独立性检验
    """
    try:
        sub = data[[categorical_col, group_col]].dropna()
        if len(sub) < 4:
            return {"p_value": np.nan, "test_used": "none", "reason": "insufficient_data", "sample_size": len(sub), "expected_min": np.nan, "error": "Less than 4 observations"}
        
        ct = _prepare_2x2_contingency_table(sub, categorical_col, group_col)
        if ct is None:
            return {"p_value": np.nan, "test_used": "none", "reason": "degenerate_table", "sample_size": len(sub), "expected_min": np.nan, "error": "Empty row or column"}
        
        table = ct.values
        total_n = int(table.sum())
        _, _, _, expected = chi2_contingency(table, correction=True)
        expected = np.asarray(expected)
        expected_min = float(expected.min())
        n_cells_lt1 = int((expected < 1).sum())
        n_cells_lt5 = int((expected < 5).sum())
        
        if n_cells_lt1 > 0 or n_cells_lt5 > 0:
            _, p_fisher = fisher_exact(table, alternative="two-sided")
            return {
                "p_value": float(p_fisher),
                "test_used": "fisher_exact",
                "reason": f"expected_min={expected_min:.3f}_cells_lt5={n_cells_lt5}",
                "sample_size": total_n,
                "expected_min": expected_min,
                "contingency_table": table.tolist(),
            }
        else:
            chi2_stat, p_chi2, dof, _ = chi2_contingency(table, correction=False)
            return {
                "p_value": float(p_chi2),
                "test_used": "chi2_pearson",
                "reason": f"expected_min={expected_min:.3f}_all_ge5",
                "sample_size": total_n,
                "chi2_stat": float(chi2_stat),
                "expected_min": expected_min,
                "contingency_table": table.tolist(),
            }
    except Exception as e:
        return {"p_value": np.nan, "test_used": "error", "reason": "unexpected_error", "sample_size": 0, "expected_min": np.nan, "error": f"Independence test failed: {str(e)}"}

# --------------------------------------------------------------------------
# --- Loop搜索函数 (v5.5 KEY FIX) ---
# --------------------------------------------------------------------------
def find_all_overlapping_loops_for_branch2(gene_start, gene_end, tad_start, tad_end):
    """
    为分支2寻找 *所有* 合适的Loop (v5.5 修复: 字典键名必须与 initialize_stats_dict 一致)
    """
    found_loops_list = []
    
    try:
        gene_anchor_start = gene_start - Config.GENE_ANCHOR_DISTANCE_KB * 1000
        gene_anchor_end = gene_end + Config.GENE_ANCHOR_DISTANCE_KB * 1000
        gene_position = determine_gene_tad_spatial_relationship(gene_start, gene_end, tad_start, tad_end)
        
        if gene_position == "upstream":
            tad_anchor_start = tad_end
            tad_anchor_end = tad_end + Config.MAX_DISTANCE_KB * 1000
            tad_anchor_side = "downstream"
        elif gene_position == "downstream":
            tad_anchor_start = tad_start - Config.MAX_DISTANCE_KB * 1000
            tad_anchor_end = tad_end
            tad_anchor_side = "upstream_and_internal"
        else:
            return [] 

        for loop_id, loop_row in _G['loop_pos_df'].iterrows():
            anchor1_gene_overlap = check_overlap(gene_anchor_start, gene_anchor_end, loop_row['start1'], loop_row['end1'])
            anchor2_tad_overlap = check_overlap(tad_anchor_start, tad_anchor_end, loop_row['start2'], loop_row['end2'])
            anchor2_gene_overlap = check_overlap(gene_anchor_start, gene_anchor_end, loop_row['start2'], loop_row['end2'])
            anchor1_tad_overlap = check_overlap(tad_anchor_start, tad_anchor_end, loop_row['start1'], loop_row['end1'])

            # (v5.5 修复: 所有的键都必须与 initialize_stats_dict 中的键名(大写开头)完全匹配)
            loop_match_data = {
                "Loop_ID": loop_id,
                "Gene_Anchor_Start": gene_anchor_start,
                "Gene_Anchor_End": gene_anchor_end,
                "TAD_Anchor_Start": tad_anchor_start,
                "TAD_Anchor_End": tad_anchor_end,
                "Gene_Position": gene_position,
                "TAD_Anchor_Side": tad_anchor_side,
                "Loop_Anchor1_Start": loop_row['start1'],
                "Loop_Anchor1_End": loop_row['end1'],
                "Loop_Anchor2_Start": loop_row['start2'],
                "Loop_Anchor2_End": loop_row['end2'],
            }

            if anchor1_gene_overlap and anchor2_tad_overlap:
                loop_match_data["Loop_Match_Type"] = "anchor1_gene_anchor2_tad" # (v5.5 修复)
                found_loops_list.append(loop_match_data)
            elif anchor2_gene_overlap and anchor1_tad_overlap:
                loop_match_data["Loop_Match_Type"] = "anchor2_gene_anchor1_tad" # (v5.5 修复)
                found_loops_list.append(loop_match_data)
        
        return found_loops_list
    
    except Exception as e:
        print(f"Loop search warning: {e}", file=sys.stderr)
        return []

# --------------------------------------------------------------------------
# --- 方向性检验函数 ---
# --------------------------------------------------------------------------
def check_branch1_directionality(tad_is_median_0, tad_is_median_2, gene_exp_median_0, gene_exp_median_2):
    """
    分支1方向性检验：TAD IS与基因表达呈相反趋势
    """
    if any(pd.isna(x) for x in [tad_is_median_0, tad_is_median_2, gene_exp_median_0, gene_exp_median_2]):
        return False
    
    return ((tad_is_median_0 > tad_is_median_2 and gene_exp_median_0 < gene_exp_median_2) or
            (tad_is_median_0 < tad_is_median_2 and gene_exp_median_0 > gene_exp_median_2))

def check_branch2_directionality(tad_is_median_0, tad_is_median_2, 
                                 gene_exp_median_0, gene_exp_median_2,
                                 loop_strength_median_0, loop_strength_median_2,
                                 loop_count_0, loop_count_2):
    """
    分支2方向性检验：所有指标呈一致趋势
    """
    values = [tad_is_median_0, tad_is_median_2, gene_exp_median_0, gene_exp_median_2,
              loop_strength_median_0, loop_strength_median_2, loop_count_0, loop_count_2]
    if any(pd.isna(x) for x in values):
        return False
    
    condition1 = (tad_is_median_0 < tad_is_median_2 and
                  gene_exp_median_0 < gene_exp_median_2 and
                  (loop_strength_median_0 < loop_strength_median_2 or loop_count_0 < loop_count_2))
    
    condition2 = (tad_is_median_0 > tad_is_median_2 and
                  gene_exp_median_0 > gene_exp_median_2 and
                  (loop_strength_median_0 > loop_strength_median_2 or loop_count_0 > loop_count_2))
    
    return condition1 or condition2

# --------------------------------------------------------------------------
# --- 数据加载和验证 ---
# --------------------------------------------------------------------------
def load_and_validate_data(full_snp_id, gene_id, tad_id):
    """加载并验证三元组数据"""
    try:
        if full_snp_id not in _G["snp_data"].index:
            return None, f"SNP not found: {full_snp_id}"
        if gene_id not in _G["gene_exp"].index:
            return None, f"Gene not found: {gene_id}"
        if tad_id not in _G["tad_is"].index:
            return None, f"TAD not found in tad_is: {tad_id}"
        
        snp_values = _G["snp_data"].loc[full_snp_id, _G["common_samples"]].astype(float)
        gene_values = _G["gene_exp"].loc[gene_id, _G["common_samples"]].astype(float)
        tad_is_values = _G["tad_is"].loc[tad_id, _G["common_samples"]].astype(float)
        
        df = pd.DataFrame({
            "Genotype_num": snp_values,
            "TAD_IS_Score": tad_is_values,
            "Gene_Exp": gene_values
        })
        
        df = df[df["Genotype_num"].isin(Config.EXPECTED_GENOTYPES)]
        
        if df.empty:
            return None, "No valid genotype data after filtering"
        if df["Genotype_num"].nunique() < 2:
            return None, "Less than 2 genotype groups"
        if not set(df["Genotype_num"].unique()) == set(Config.EXPECTED_GENOTYPES):
            return None, f"Unexpected genotypes: {sorted(df['Genotype_num'].unique())}"
        
        return df, None
        
    except KeyError as e:
        return None, f"KeyError in data loading: {str(e)}"
    except Exception as e:
        return None, f"Unexpected error in data loading: {str(e)}"

# --------------------------------------------------------------------------
# --- 分析函数 (v5.2/v5.3 结构) ---
# --------------------------------------------------------------------------
def analyze_branch1(data_df, gene_start, gene_end, tad_start, tad_end):
    """分析分支1 (三元组级别)"""
    branch1_start = tad_start - Config.DISTANCE_THRESHOLD_KB * 1000
    branch1_end = tad_end + Config.DISTANCE_THRESHOLD_KB * 1000
    distance_pass = check_overlap(gene_start, gene_end, branch1_start, branch1_end)
    
    result = {
        "Branch1_Range_Start": branch1_start,
        "Branch1_Range_End": branch1_end,
        "Branch1_Distance_Pass": distance_pass
    }
    
    if not distance_pass:
        result["Failure_Reason"] = f"Branch1: Gene not within {Config.DISTANCE_THRESHOLD_KB}kb of TAD"
        return {"passed": False, "stats": result}
    
    median_is_0 = data_df[data_df["Genotype_num"] == 0]["TAD_IS_Score"].median()
    median_is_2 = data_df[data_df["Genotype_num"] == 2]["TAD_IS_Score"].median()
    median_exp_0 = data_df[data_df["Genotype_num"] == 0]["Gene_Exp"].median()
    median_exp_2 = data_df[data_df["Genotype_num"] == 2]["Gene_Exp"].median()
    
    result.update({
        "TAD_IS_Median_Group0": median_is_0,
        "TAD_IS_Median_Group2": median_is_2,
        "Gene_Exp_Median_Group0": median_exp_0,
        "Gene_Exp_Median_Group2": median_exp_2
    })
    
    direction_pass = check_branch1_directionality(median_is_0, median_is_2, median_exp_0, median_exp_2)
    result["Branch1_Direction_Pass"] = direction_pass
    
    if not direction_pass:
        result["Failure_Reason"] = "Branch1: Direction criteria not met"
        return {"passed": False, "stats": result}
    
    return {"passed": True, "stats": result}

def _test_single_loop_for_branch2(data_df_base, loop_info):
    """
    (v5.2 辅助函数): 独立测试单个Loop (四元组分析的核心)
    """
    loop_stats_results = {}
    current_loop_id = loop_info["Loop_ID"] # (v5.5 修复: 使用大写的键)
    
    try:
        if current_loop_id not in _G["loop_strength_df"].index or current_loop_id not in _G["loop_01_df"].index:
            return {"passed": False, "reason": "loop_data_not_found_in_matrix"}
            
        loop_strength_values = _G["loop_strength_df"].loc[current_loop_id, _G["common_samples"]].astype(float)
        loop_01_values = _G["loop_01_df"].loc[current_loop_id, _G["common_samples"]].astype(float)
        
        data_df_loop = data_df_base.copy()
        data_df_loop["Loop_Strength"] = loop_strength_values
        data_df_loop["Loop_Present"] = loop_01_values
        
        loop_strength_test = perform_mannwhitney_test(data_df_loop, "Loop_Strength", "Genotype_num")
        loop_01_test = perform_independence_test(data_df_loop, "Loop_Present", "Genotype_num")
        
        loop_strength_p = loop_strength_test.get("p_value", np.nan)
        loop_01_p = loop_01_test.get("p_value", np.nan)
        
        loop_significant = ((pd.notna(loop_strength_p) and loop_strength_p < Config.P_THRESHOLD_LOOP) or
                            (pd.notna(loop_01_p) and loop_01_p < Config.P_THRESHOLD_LOOP))

        if not loop_significant:
            return {"passed": False, "reason": f"not_significant(S_p={loop_strength_p:.2g}; P_p={loop_01_p:.2g})"}
        
        sig_tests = []
        if pd.notna(loop_strength_p) and loop_strength_p < Config.P_THRESHOLD_LOOP:
            sig_tests.append(f"strength_p={loop_strength_p:.4f}")
        if pd.notna(loop_01_p) and loop_01_p < Config.P_THRESHOLD_LOOP:
            sig_tests.append(f"presence_p={loop_01_p:.4f}")
        loop_stats_results["Loop_Significant_Tests"] = "; ".join(sig_tests)

        median_is_0 = data_df_loop[data_df_loop["Genotype_num"] == 0]["TAD_IS_Score"].median()
        median_is_2 = data_df_loop[data_df_loop["Genotype_num"] == 2]["TAD_IS_Score"].median()
        median_exp_0 = data_df_loop[data_df_loop["Genotype_num"] == 0]["Gene_Exp"].median()
        median_exp_2 = data_df_loop[data_df_loop["Genotype_num"] == 2]["Gene_Exp"].median()
        median_loop_0 = loop_strength_test.get("group0_median", np.nan)
        median_loop_2 = loop_strength_test.get("group2_median", np.nan)
        
        ct_loop_01 = _prepare_2x2_contingency_table(data_df_loop, "Loop_Present", "Genotype_num")
        count_loop_0 = ct_loop_01.loc[1, 0] if ct_loop_01 is not None and 1 in ct_loop_01.index else 0
        count_loop_2 = ct_loop_01.loc[1, 2] if ct_loop_01 is not None and 1 in ct_loop_01.index else 0

        direction_pass = check_branch2_directionality(median_is_0, median_is_2, median_exp_0, median_exp_2,
                                                      median_loop_0, median_loop_2, count_loop_0, count_loop_2)
        
        if not direction_pass:
            return {"passed": False, "reason": "direction_fail"}
        
        loop_stats_results.update({
            "Branch2_Loop_Significant": True,
            "Branch2_Direction_Pass": True,
            "Loop_Strength_P": loop_strength_p,
            "Loop_Strength_Stat": loop_strength_test.get("stat", np.nan),
            "Loop_Strength_N0": loop_strength_test.get("group0_size", 0),
            "Loop_Strength_N2": loop_strength_test.get("group2_size", 0),
            "Loop_01_P": loop_01_p,
            "Loop_01_Test_Method": loop_01_test.get("test_used", "unknown"),
            "Loop_01_Test_Reason": loop_01_test.get("reason", ""),
            "Loop_01_Sample_Size": loop_01_test.get("sample_size", 0),
            "Loop_01_Expected_Min": loop_01_test.get("expected_min", np.nan),
            "TAD_IS_Median_Group0": median_is_0,
            "TAD_IS_Median_Group2": median_is_2,
            "Gene_Exp_Median_Group0": median_exp_0,
            "Gene_Exp_Median_Group2": median_exp_2,
            "Loop_Strength_Median_Group0": median_loop_0,
            "Loop_Strength_Median_Group2": median_loop_2,
            "Loop_Count_Group0": count_loop_0,
            "Loop_Count_Group2": count_loop_2
        })
        
        return {"passed": True, "stats": loop_stats_results}

    except Exception as e:
        return {"passed": False, "reason": f"unexpected_error({str(e)})"}

def initialize_stats_dict(gene_id, tad_id, snp_id, full_snp_id,
                          gene_chrom, gene_start, gene_end,
                          tad_chrom, tad_start, tad_end,
                          gene_tad_distance, gene_tad_relation):
    """初始化统计结果字典 (v5.5: 修正了诊断列名称以减少混淆)"""
    return {
        # 基础信息
        "Gene_ID": gene_id, "TAD_ID": tad_id, "SNP_ID": snp_id, "Full_SNP_ID": full_snp_id,
        "Gene_chrom": gene_chrom, "Gene_start": gene_start, "Gene_end": gene_end,
        "TAD_chrom": tad_chrom, "TAD_start": tad_start, "TAD_end": tad_end,
        "Gene_TAD_Distance_bp": gene_tad_distance,
        "Gene_TAD_Distance_kb": round(gene_tad_distance / 1000, 2),
        "Gene_TAD_Relation": gene_tad_relation,
        
        # 分支判断
        "Branch_Pass": "None", 
        "Branch1_Distance_Pass": False,
        "Branch1_Direction_Pass": False,
        "Branch2_Distance_Pass": False,
        "Branch2_Loop_Found": False, 
        "Branch2_Loop_Significant": False, 
        "Branch2_Direction_Pass": False, 
        
        # 范围信息
        "Branch1_Range_Start": np.nan,
        "Branch1_Range_End": np.nan,
        "Branch2_Upstream_Range_Start": np.nan,
        "Branch2_Upstream_Range_End": np.nan,
        "Branch2_Downstream_Range_Start": np.nan,
        "Branch2_Downstream_Range_End": np.nan,
        
        # Loop信息
        "Loop_ID": np.nan, 
        "Loop_Match_Type": np.nan,
        "Triplet_B2_Candidate_Count": 0, # (v5.5 重命名: 澄清这是三元组级别的诊断)
        "Triplet_B2_Candidate_List_All": np.nan, # (v5.5 重命名)
        "Gene_Anchor_Start": np.nan,
        "Gene_Anchor_End": np.nan,
        "Gene_Anchor_Distance_kb": Config.GENE_ANCHOR_DISTANCE_KB,
        "TAD_Anchor_Start": np.nan,
        "TAD_Anchor_End": np.nan,
        "Gene_Position": np.nan,
        "TAD_Anchor_Side": np.nan,
        "Loop_Anchor1_Start": np.nan,
        "Loop_Anchor1_End": np.nan,
        "Loop_Anchor2_Start": np.nan,
        "Loop_Anchor2_End": np.nan,
        
        # 统计检验详细信息
        "Loop_Strength_P": np.nan,
        "Loop_Strength_Stat": np.nan,
        "Loop_Strength_N0": np.nan,
        "Loop_Strength_N2": np.nan,
        "Loop_Strength_Median0": np.nan,
        "Loop_Strength_Median2": np.nan,
        "Loop_01_P": np.nan,
        "Loop_01_Test_Method": np.nan,
        "Loop_01_Test_Reason": np.nan,
        "Loop_01_Sample_Size": np.nan,
        "Loop_01_Expected_Min": np.nan,
        "Loop_Significant_Tests": np.nan,
        
        # 中位数信息
        "TAD_IS_Median_Group0": np.nan,
        "TAD_IS_Median_Group2": np.nan,
        "Gene_Exp_Median_Group0": np.nan,
        "Gene_Exp_Median_Group2": np.nan,
        "Loop_Strength_Median_Group0": np.nan, 
        "Loop_Strength_Median_Group2": np.nan, 
        "Loop_Count_Group0": np.nan, 
        "Loop_Count_Group2": np.nan, 
        
        # 样本信息
        "Sample_Count_Group0": np.nan,
        "Sample_Count_Group2": np.nan,
        
        # 诊断信息
        "Failure_Reason": np.nan
    }

def analyze_single_triplet(triplet_dict):
    """
    分析单个三元组的主函数 (v5.5 结构)
    """
    
    stats_template = {} 
    try:
        gene_id = triplet_dict['Gene_ID']
        tad_id = triplet_dict['TAD_ID']
        snp_id = triplet_dict['SNP_ID']
        full_snp_id = triplet_dict['Full_SNP_ID']
        gene_chrom = triplet_dict['Gene_chrom']
        gene_start = triplet_dict['Gene_start']
        gene_end = triplet_dict['Gene_end']
        tad_chrom = triplet_dict['TAD_chrom']
        tad_start = triplet_dict['TAD_start']
        tad_end = triplet_dict['TAD_end']

        gene_tad_distance, gene_tad_relation = calculate_gene_tad_distance_and_relation(
            gene_start, gene_end, tad_start, tad_end
        )
        
        stats_template = initialize_stats_dict(gene_id, tad_id, snp_id, full_snp_id,
                                               gene_chrom, gene_start, gene_end,
                                               tad_chrom, tad_start, tad_end,
                                               gene_tad_distance, gene_tad_relation)
        
        data_df, failure_reason = load_and_validate_data(full_snp_id, gene_id, tad_id)
        if data_df is None:
            stats_template["Failure_Reason"] = failure_reason
            return [stats_template] 
        
        stats_template["Sample_Count_Group0"] = len(data_df[data_df["Genotype_num"] == 0])
        stats_template["Sample_Count_Group2"] = len(data_df[data_df["Genotype_num"] == 2])
        
        branch1_result = analyze_branch1(data_df, gene_start, gene_end, tad_start, tad_end)
        stats_template.update(branch1_result["stats"]) 

        if branch1_result["passed"]:
            stats_template["Branch_Pass"] = "Branch_1_Pass"
            return [stats_template] 
        
        upstream_start = tad_start - Config.MAX_DISTANCE_KB * 1000
        upstream_end = tad_start - Config.DISTANCE_THRESHOLD_KB * 1000
        downstream_start = tad_end + Config.DISTANCE_THRESHOLD_KB * 1000
        downstream_end = tad_end + Config.MAX_DISTANCE_KB * 1000
        distance_pass = (check_overlap(gene_start, gene_end, upstream_start, upstream_end) or
                         check_overlap(gene_start, gene_end, downstream_start, downstream_end))

        stats_template.update({
            "Branch2_Upstream_Range_Start": upstream_start,
            "Branch2_Upstream_Range_End": upstream_end,
            "Branch2_Downstream_Range_Start": downstream_start,
            "Branch2_Downstream_Range_End": downstream_end,
            "Branch2_Distance_Pass": distance_pass
        })
        
        if not distance_pass:
            stats_template["Failure_Reason"] = (f"{stats_template.get('Failure_Reason', 'B1_Fail')}; "
                                               f"Branch2: Gene not in {Config.DISTANCE_THRESHOLD_KB}-{Config.MAX_DISTANCE_KB}kb range")
            return [stats_template] 
        
        candidate_loops_list = find_all_overlapping_loops_for_branch2(gene_start, gene_end, tad_start, tad_end)
        # (v5.5 重命名)
        stats_template["Triplet_B2_Candidate_Count"] = len(candidate_loops_list)
        if candidate_loops_list:
            # (v5.5 重命名)
            stats_template["Triplet_B2_Candidate_List_All"] = ";".join([str(l['Loop_ID']) for l in candidate_loops_list]) 
            stats_template["Branch2_Loop_Found"] = True
        else:
            stats_template["Triplet_B2_Candidate_List_All"] = np.nan # (v5.5 重命名)
            stats_template["Failure_Reason"] = (f"{stats_template.get('Failure_Reason', 'B1_Fail')}; "
                                               f"Branch2: No overlapping loop found")
            return [stats_template] 

        passing_quadruplet_rows = []
        loop_failure_reasons = []

        for loop_info in candidate_loops_list:
            # loop_info 此时是一个包含 大写键 Loop_ID, Gene_Anchor_Start... 的字典
            loop_test_result = _test_single_loop_for_branch2(data_df, loop_info)
            
            if loop_test_result["passed"]:
                passing_row_stats = stats_template.copy()
                passing_row_stats.update(loop_info) # (v5.5 修复: 此处 loop_info 键已全部大写, 能够正确填充)
                passing_row_stats.update(loop_test_result["stats"]) 
                passing_row_stats["Branch_Pass"] = "Branch_2_Pass"
                passing_row_stats["Failure_Reason"] = np.nan 
                passing_quadruplet_rows.append(passing_row_stats)
            else:
                loop_failure_reasons.append(f"{loop_info['Loop_ID']}:{loop_test_result['reason']}") # (v5.5 修复: 使用大写 'Loop_ID')

        if passing_quadruplet_rows:
            return passing_quadruplet_rows
        else:
            stats_template["Failure_Reason"] = (f"{stats_template.get('Failure_Reason', 'B1_Fail')}; "
                                              f"Branch2: All {len(candidate_loops_list)} candidate loops failed tests. "
                                              f"Reasons: [{'; '.join(loop_failure_reasons)}]")
            return [stats_template]

    except Exception as e:
        if not stats_template:
            stats_template = {"Failure_Reason": f"Critical error: {str(e)}"}
        else:
            stats_template["Failure_Reason"] = f"Unexpected error: {str(e)}"
        return [stats_template] 

# --------------------------------------------------------------------------
# --- 多进程支持 ---
# --------------------------------------------------------------------------
def _init_globals(args_main):
    """多进程初始化函数"""
    try:
        data_files = {
            "snp_data": pd.read_csv(args_main.snp_data, sep='\t', index_col=0),
            "gene_exp": pd.read_csv(args_main.gene_exp, sep='\t', index_col=0),
            "tad_is": pd.read_csv(args_main.tad_is, sep='\t', index_col=0),
            "tad_01": pd.read_csv(args_main.tad_01, sep='\t', index_col=0),
            "loop_pos": pd.read_csv(args_main.loop_pos, sep='\t', index_col=0),
            "loop_strength": pd.read_csv(args_main.loop_strength, sep='\t', index_col=0),
            "loop_01": pd.read_csv(args_main.loop_01, sep='\t', index_col=0),
        }

        for df in data_files.values():
            if not df.columns.empty:
                df.columns = df.columns.str.strip()
        
        _G.update(data_files)
        _G["loop_strength_df"] = data_files["loop_strength"]
        _G["loop_01_df"] = data_files["loop_01"]
        _G["loop_pos_df"] = data_files["loop_pos"]
        
        if _G["loop_pos_df"].index.name is None and not _G["loop_pos_df"].index.is_unique:
             _G["loop_pos_df"] = _G["loop_pos_df"].reset_index(drop=True)
        elif _G["loop_pos_df"].index.name is not None and not _G["loop_pos_df"].index.is_unique:
             print("警告: loop_pos 索引列不唯一，但被保留。iterrows可能迭代重复的索引。", file=sys.stderr)

        _G["common_samples"] = data_files["snp_data"].columns.tolist()

    except Exception as e:
        print(f"子进程数据加载失败: {e}", file=sys.stderr)
        raise

def convert_triplets_to_dicts(triplets_df):
    """转换DataFrame为字典列表"""
    return [
        {
            'Gene_ID': row['Gene_ID'],
            'TAD_ID': row['TAD_ID'],
            'SNP_ID': row['SNP_ID'],
            'Full_SNP_ID': row['Full_SNP_ID'],
            'Gene_chrom': row['Gene_chrom'],
            'Gene_start': row['Gene_start'],
            'Gene_end': row['Gene_end'],
            'TAD_chrom': row['TAD_chrom'],
            'TAD_start': row['TAD_start'],
            'TAD_end': row['TAD_end']
        }
        for _, row in triplets_df.iterrows()
    ]

# --------------------------------------------------------------------------
# --- 结果保存和报告 (v5.5 MODIFIED: 匹配重命名) ---
# --------------------------------------------------------------------------
def save_results(all_output_rows, passed_output_rows, output_dir):
    """
    保存分析结果 (v5.5)
    """
    if not all_output_rows:
        print("没有结果数据", file=sys.stderr)
        return
    
    all_stats_df = pd.DataFrame(all_output_rows)
    
    # 列顺序
    output_cols = [
        # 基础信息
        "Gene_ID", "TAD_ID", "SNP_ID", "Full_SNP_ID",
        "Gene_chrom", "Gene_start", "Gene_end", "TAD_chrom", "TAD_start", "TAD_end",
        "Gene_TAD_Distance_bp", "Gene_TAD_Distance_kb", "Gene_TAD_Relation",
        
        # 分支判断
        "Branch_Pass", "Branch1_Distance_Pass", "Branch1_Direction_Pass", 
        "Branch2_Distance_Pass", "Branch2_Loop_Found", "Branch2_Loop_Significant", "Branch2_Direction_Pass",
        
        # 范围
        "Branch1_Range_Start", "Branch1_Range_End", 
        "Branch2_Upstream_Range_Start", "Branch2_Upstream_Range_End",
        "Branch2_Downstream_Range_Start", "Branch2_Downstream_Range_End",
        
        # Loop信息
        "Loop_ID", "Loop_Match_Type", 
        "Triplet_B2_Candidate_Count", "Triplet_B2_Candidate_List_All", # (v5.5 重命名)
        "Gene_Anchor_Start", "Gene_Anchor_End", "Gene_Anchor_Distance_kb",
        "TAD_Anchor_Start", "TAD_Anchor_End", "Gene_Position", "TAD_Anchor_Side",
        "Loop_Anchor1_Start", "Loop_Anchor1_End", "Loop_Anchor2_Start", "Loop_Anchor2_End",
        
        # 统计检验
        "Loop_Strength_P", "Loop_Strength_Stat", "Loop_Strength_N0", "Loop_Strength_N2",
        "Loop_Strength_Median0", "Loop_Strength_Median2",
        "Loop_01_P", "Loop_01_Test_Method", "Loop_01_Test_Reason", "Loop_01_Sample_Size", "Loop_01_Expected_Min",
        "Loop_Significant_Tests",
        
        # 中位数
        "TAD_IS_Median_Group0", "TAD_IS_Median_Group2", "Gene_Exp_Median_Group0", "Gene_Exp_Median_Group2",
        "Loop_Strength_Median_Group0", "Loop_Strength_Median_Group2", "Loop_Count_Group0", "Loop_Count_Group2",
        
        # 样本和诊断
        "Sample_Count_Group0", "Sample_Count_Group2", "Failure_Reason"
    ]
    
    output_cols = [col for col in output_cols if col in all_stats_df.columns]
    all_stats_df = all_stats_df.reindex(columns=output_cols)
    
    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 保存文件
    complete_file = output_path / "complete_analysis_results_quadruplets.csv"
    all_stats_df.to_csv(complete_file, sep=",", index=False, float_format='%.6g') 
    print(f"完整结果已保存: {complete_file.absolute()}", file=sys.stderr)
    
    if passed_output_rows:
        passed_df = pd.DataFrame(passed_output_rows).reindex(columns=output_cols)
        passed_file = output_path / "passed_quadruplets_results.csv"
        passed_df.to_csv(passed_file, sep=",", index=False, float_format='%.6g')
        print(f"通过筛选结果已保存: {passed_file.absolute()}", file=sys.stderr)

def print_summary_statistics(all_output_rows, passed_output_rows, total_triplets_input_count):
    """
    输出汇总统计 (v5.5, 匹配重命名的列)
    """
    
    print("-" * 60, file=sys.stderr)
    print("分析结果汇总 (v5.5 - 详细漏斗分析报告):", file=sys.stderr)
    print("=" * 60, file=sys.stderr)

    if not all_output_rows:
        print("警告: 未生成任何分析结果。", file=sys.stderr)
        return

    all_df = pd.DataFrame(all_output_rows)
    passed_df = pd.DataFrame(passed_output_rows) if passed_output_rows else pd.DataFrame(columns=all_df.columns)
    failed_df = all_df[all_df["Branch_Pass"] == "None"]

    unique_triplet_keys = ["SNP_ID", "TAD_ID", "Gene_ID"]
    total_triplets_processed = all_df[unique_triplet_keys].drop_duplicates().shape[0]

    # --- 1. 漏斗总览 ---
    print(f"--- 1. 分析漏斗总览 ---", file=sys.stderr)
    print(f"总共接收三元组: {total_triplets_input_count}", file=sys.stderr)
    if total_triplets_processed != total_triplets_input_count:
         print(f"警告: 已处理的独特三元组 ({total_triplets_processed}) 与输入数 ({total_triplets_input_count}) 不匹配!", file=sys.stderr)
    
    passed_b1_df = passed_df[passed_df["Branch_Pass"] == "Branch_1_Pass"]
    passed_b1_count = len(passed_b1_df)
    
    passed_b2_df = passed_df[passed_df["Branch_Pass"] == "Branch_2_Pass"]
    passed_b2_rows = len(passed_b2_df)
    passed_b2_unique_triplets = 0
    if passed_b2_rows > 0:
        passed_b2_unique_triplets = passed_b2_df[unique_triplet_keys].drop_duplicates().shape[0]

    failed_count = len(failed_df)

    print(f"\n[ 漏斗结果 ]:", file=sys.stderr)
    print(f"  > 成功 (分支1 通过): {passed_b1_count} 个三元组", file=sys.stderr)
    print(f"  > 成功 (分支2 通过): {passed_b2_unique_triplets} 个三元组 (产生了 {passed_b2_rows} 个成功的四元组)", file=sys.stderr)
    print(f"  > 失败 (所有分支):   {failed_count} 个三元组", file=sys.stderr)
    print(f"  ---------------------------------", file=sys.stderr)
    print(f"  总计:               {passed_b1_count + passed_b2_unique_triplets + failed_count} 个三元组被分类。", file=sys.stderr)

    # --- 2. 成功案例分析 ---
    print(f"\n--- 2. 成功案例分析 (共 {passed_b1_count + passed_b2_rows} 行) ---", file=sys.stderr)
    if passed_b2_unique_triplets > 0:
        print(f"分支2 成功详情 (来自 {passed_b2_unique_triplets} 个三元组):", file=sys.stderr)
        unique_passed_b2_triplets_df = passed_b2_df.drop_duplicates(subset=unique_triplet_keys)
        # (v5.5 重命名)
        avg_cands = unique_passed_b2_triplets_df['Triplet_B2_Candidate_Count'].mean()
        min_cands = unique_passed_b2_triplets_df['Triplet_B2_Candidate_Count'].min()
        max_cands = unique_passed_b2_triplets_df['Triplet_B2_Candidate_Count'].max()
        avg_passed_quads = passed_b2_rows / passed_b2_unique_triplets
        
        print(f"  - 平均每个成功的三元组产生了 {avg_passed_quads:.2f} 个成功的Loop通路。", file=sys.stderr)
        print(f"  - 这些成功的三元组平均匹配了 {avg_cands:.1f} 个候选Loop (范围从 {min_cands} 到 {max_cands})。", file=sys.stderr)

    # --- 3. 失败案例分析 ---
    print(f"\n--- 3. 失败案例分析 (共 {failed_count} 个三元组) ---", file=sys.stderr)
    print("失败原因分类统计 (所有类别):", file=sys.stderr)
    
    if failed_count > 0:
        failure_counts = failed_df["Failure_Reason"].value_counts().sort_values(ascending=False)
        for reason, count in failure_counts.items(): 
            reason_short = (reason[:130] + '...') if len(str(reason)) > 130 else str(reason)
            print(f"  - [{reason_short}]: {count} 个三元组", file=sys.stderr)

        failed_all_loops_df = failed_df[failed_df['Failure_Reason'].str.contains("All candidate loops failed tests", na=False)]
        failed_all_loops_count = len(failed_all_loops_df)
        if failed_all_loops_count > 0:
            print(f"\n  'All Loops Failed' 案例详情 (共 {failed_all_loops_count} 个):", file=sys.stderr)
            # (v5.5 重命名)
            avg_fail_cands = failed_all_loops_df['Triplet_B2_Candidate_Count'].mean()
            min_fail_cands = failed_all_loops_df['Triplet_B2_Candidate_Count'].min()
            max_fail_cands = failed_all_loops_df['Triplet_B2_Candidate_Count'].max()
            print(f"    - 这些失败的三元组平均测试了 {avg_fail_cands:.1f} 个Loop (范围从 {min_fail_cands} 到 {max_fail_cands})。", file=sys.stderr)
            print(f"    - ... 并且所有候选Loop均未通过检验 (p值过高 或 方向性错误)。", file=sys.stderr)
    else:
        print("  - 无失败记录。", file=sys.stderr)

    # --- 4. 统计方法使用情况 ---
    print(f"\n--- 4. 统计检验方法使用 (针对所有B2中测试过的Loop) ---", file=sys.stderr)
    loop_methods = all_df[all_df["Loop_01_Test_Method"].notna()]["Loop_01_Test_Method"].value_counts()
    if not loop_methods.empty:
        total_tests = loop_methods.sum()
        print(f"  总共执行了 {total_tests} 次 Loop 0-1 独立性检验:", file=sys.stderr)
        for method, count in loop_methods.items():
            pct = count / total_tests * 100
            print(f"    - {method}: {count} 次 ({pct:.1f}%)", file=sys.stderr)
    else:
        print("  - 未执行任何 Loop 0-1 检验。", file=sys.stderr)
        
    print("=" * 60, file=sys.stderr)

# --------------------------------------------------------------------------
# --- 主函数 ---
# --------------------------------------------------------------------------
def main():
    """主函数 (v5.4)""" # 逻辑与 v5.4 相同
    args = parse_args()
    
    # 设置配置
    Config.DISTANCE_THRESHOLD_KB = args.distance_threshold
    Config.MAX_DISTANCE_KB = args.max_distance
    Config.GENE_ANCHOR_DISTANCE_KB = args.gene_anchor_distance
    Config.P_THRESHOLD_LOOP = args.loop_p_threshold
    
    # 输出配置
    print("=" * 60, file=sys.stderr)
    print("严谨SNP-TAD-Gene分析开始 (v5.5 - 四元组策略 + 键名修复)", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"配置参数:", file=sys.stderr)
    print(f"   策略: 评估所有 (SNP,TAD,Gene,Loop) 四元组。", file=sys.stderr)
    print(f"   边界逻辑: 端点包容 (<=)。", file=sys.stderr)
    print(f"   Loop p值阈值: {Config.P_THRESHOLD_LOOP} (独立应用于每个候选Loop)", file=sys.stderr)
    print(f"   输出目录: {args.output_dir}", file=sys.stderr)
    print("-" * 60, file=sys.stderr)
    
    # 加载三元组
    triplets_df = pd.read_csv(args.triplets_with_positions)
    total_triplets_to_process = len(triplets_df)
    print(f"加载 {total_triplets_to_process} 个三元组进行分析", file=sys.stderr)
    
    triplets_list = convert_triplets_to_dicts(triplets_df)
    
    # 多进程分析
    all_output_rows_flat = [] 
    
    with Pool(processes=args.num_cores, initializer=_init_globals, initargs=(args,)) as pool:
        for result_list_per_triplet in tqdm(pool.imap_unordered(analyze_single_triplet, triplets_list, chunksize=5),
                                            total=total_triplets_to_process, desc="分析三元组", file=sys.stderr):
            
            all_output_rows_flat.extend(result_list_per_triplet)

    # 整理结果
    passed_output_rows = [row for row in all_output_rows_flat if row["Branch_Pass"] != "None"]
    
    print(f"\n分析完成: {total_triplets_to_process} 个三元组被处理。", file=sys.stderr)
    print(f"总共生成 {len(all_output_rows_flat)} 行输出 (包括 {len(passed_output_rows)} 行通过筛选的记录和 {len(all_output_rows_flat) - len(passed_output_rows)} 行失败记录)。", file=sys.stderr)
    
    # 保存和报告
    save_results(all_output_rows_flat, passed_output_rows, args.output_dir)
    print_summary_statistics(all_output_rows_flat, passed_output_rows, total_triplets_to_process)

if __name__ == "__main__":
    main()