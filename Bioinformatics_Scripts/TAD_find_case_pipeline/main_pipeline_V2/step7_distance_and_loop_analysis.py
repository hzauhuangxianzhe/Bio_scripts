#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP-TAD-Gene Analysis Core (v5) - Statistically Rigorous Version
Author: Based on user requirements  
Date: 2025
Version: 5.0.0

严谨的统计分析流程：
1. 分支1：基因在TAD临界距离内，检查TAD IS与基因表达的反向关联
2. 分支2：基因在TAD临界距离外，需要Loop介导，检查所有指标的一致性
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
    parser = ArgumentParser(description="基于TAD-基因距离和方向性规则的严谨统计分析")
    
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
                       help="Loop显著性p值阈值 (默认0.05)")
    parser.add_argument("--output_dir", type=str, default="analysis_output_v5",
                       help="输出文件夹路径 (默认: analysis_output_v5)")
    
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
# --- 核心计算函数 ---
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
    """检查两个区间是否重叠"""
    return max(start1, start2) < min(end1, end2)

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
    
    特点：
    - 成对删除缺失值
    - 严格的样本量检查
    - 双侧检验
    - 完整的诊断信息
    """
    try:
        # 成对删除策略
        sub = data[[value_col, group_col]].dropna()
        
        # 最小样本量和分组检查
        if len(sub) < 3:
            return {
                "p_value": np.nan,
                "stat": np.nan,
                "group0_size": 0,
                "group2_size": 0,
                "group0_median": np.nan,
                "group2_median": np.nan,
                "total_samples": len(sub),
                "error": "Insufficient total sample size"
            }
        
        if sub[group_col].nunique() < 2:
            return {
                "p_value": np.nan,
                "stat": np.nan,
                "group0_size": 0,
                "group2_size": 0,
                "group0_median": np.nan,
                "group2_median": np.nan,
                "total_samples": len(sub),
                "error": "Less than 2 groups available"
            }
        
        # 分组数据
        groups = sorted(sub[group_col].unique())
        if len(groups) != 2 or groups != [0, 2]:
            return {
                "p_value": np.nan,
                "stat": np.nan,
                "group0_size": 0,
                "group2_size": 0,
                "group0_median": np.nan,
                "group2_median": np.nan,
                "total_samples": len(sub),
                "error": f"Unexpected groups: {groups}"
            }
        
        g0 = sub[sub[group_col] == 0][value_col]
        g2 = sub[sub[group_col] == 2][value_col]
        
        # 每组至少需要1个样本
        if len(g0) == 0 or len(g2) == 0:
            return {
                "p_value": np.nan,
                "stat": np.nan,
                "group0_size": len(g0),
                "group2_size": len(g2),
                "group0_median": float(g0.median()) if len(g0) > 0 else np.nan,
                "group2_median": float(g2.median()) if len(g2) > 0 else np.nan,
                "total_samples": len(sub),
                "error": "Empty group(s)"
            }
        
        # 执行双侧Mann-Whitney U检验
        stat, p_value = mannwhitneyu(g0, g2, alternative='two-sided')
        
        return {
            "p_value": float(p_value),
            "stat": float(stat),
            "group0_size": len(g0),
            "group2_size": len(g2),
            "group0_median": float(g0.median()),
            "group2_median": float(g2.median()),
            "total_samples": len(sub)
        }
        
    except Exception as e:
        return {
            "p_value": np.nan,
            "stat": np.nan,
            "group0_size": 0,
            "group2_size": 0,
            "group0_median": np.nan,
            "group2_median": np.nan,
            "total_samples": 0,
            "error": f"Mann-Whitney test failed: {str(e)}"
        }

def _prepare_2x2_contingency_table(sub, categorical_col, group_col):
    """预处理2x2列联表并检查退化情况"""
    ct = pd.crosstab(sub[categorical_col], sub[group_col])
    # 强制标准2x2格式
    ct = ct.reindex(index=[0, 1], columns=[0, 2], fill_value=0)
    
    # 检查退化：任何行或列边际为0
    if (ct.sum(axis=1) == 0).any() or (ct.sum(axis=0) == 0).any():
        return None
    return ct

def perform_independence_test(data, categorical_col, group_col):
    """
    严谨的2x2独立性检验
    
    特点：
    - 科学的期望频数规则
    - Fisher vs 卡方的自适应选择
    - 详细的诊断信息
    - 完善的异常处理
    """
    try:
        # 成对删除缺失值
        sub = data[[categorical_col, group_col]].dropna()
        
        # 最小样本量检查
        if len(sub) < 4:
            return {
                "p_value": np.nan,
                "test_used": "none",
                "reason": "insufficient_data",
                "sample_size": len(sub),
                "expected_min": np.nan,
                "expected_max": np.nan,
                "error": "Less than 4 observations for 2x2 table"
            }
        
        # 准备2x2列联表
        ct = _prepare_2x2_contingency_table(sub, categorical_col, group_col)
        if ct is None:
            return {
                "p_value": np.nan,
                "test_used": "none",
                "reason": "degenerate_table",
                "sample_size": len(sub),
                "expected_min": np.nan,
                "expected_max": np.nan,
                "error": "Empty row or column in contingency table"
            }
        
        table = ct.values
        total_n = int(table.sum())
        
        # 计算期望频数
        _, _, _, expected = chi2_contingency(table, correction=True)
        expected = np.asarray(expected)
        
        expected_min = float(expected.min())
        expected_max = float(expected.max())
        n_cells_lt1 = int((expected < 1).sum())
        n_cells_lt5 = int((expected < 5).sum())
        
        # 严格的检验选择规则
        if n_cells_lt1 > 0 or n_cells_lt5 > 0:
            # 使用Fisher精确检验
            _, p_fisher = fisher_exact(table, alternative="two-sided")
            return {
                "p_value": float(p_fisher),
                "test_used": "fisher_exact",
                "reason": f"expected_min={expected_min:.3f}_cells_lt1={n_cells_lt1}_cells_lt5={n_cells_lt5}",
                "sample_size": total_n,
                "expected_min": expected_min,
                "expected_max": expected_max,
                "contingency_table": table.tolist(),
                "expected_frequencies": expected.tolist()
            }
        else:
            # 使用Pearson卡方检验（不使用Yates校正）
            chi2_stat, p_chi2, dof, _ = chi2_contingency(table, correction=False)
            return {
                "p_value": float(p_chi2),
                "test_used": "chi2_pearson",
                "reason": f"expected_min={expected_min:.3f}_all_ge5",
                "sample_size": total_n,
                "chi2_stat": float(chi2_stat),
                "degrees_of_freedom": int(dof),
                "expected_min": expected_min,
                "expected_max": expected_max,
                "contingency_table": table.tolist(),
                "expected_frequencies": expected.tolist()
            }
            
    except Exception as e:
        return {
            "p_value": np.nan,
            "test_used": "error",
            "reason": "unexpected_error",
            "sample_size": 0,
            "expected_min": np.nan,
            "expected_max": np.nan,
            "error": f"Independence test failed: {str(e)}"
        }

# --------------------------------------------------------------------------
# --- Loop搜索函数 ---
# --------------------------------------------------------------------------
def find_overlapping_loop_for_branch2(gene_start, gene_end, tad_start, tad_end):
    """
    为分支2寻找合适的Loop
    
    逻辑：
    1. 基因anchor：基因 ± gene_anchor_distance
    2. TAD anchor：根据基因位置确定对侧搜索范围
    3. 寻找两个anchor分别重叠的Loop
    """
    try:
        # 基因anchor范围
        gene_anchor_start = gene_start - Config.GENE_ANCHOR_DISTANCE_KB * 1000
        gene_anchor_end = gene_end + Config.GENE_ANCHOR_DISTANCE_KB * 1000
        
        # 确定基因位置
        gene_position = determine_gene_tad_spatial_relationship(gene_start, gene_end, tad_start, tad_end)
        
        if gene_position == "upstream":
            # 基因在TAD上游，TAD anchor在TAD下游搜索
            tad_anchor_start = tad_end
            tad_anchor_end = tad_end + Config.MAX_DISTANCE_KB * 1000
            tad_anchor_side = "downstream"
        elif gene_position == "downstream":
            # 基因在TAD下游，TAD anchor在TAD上游+内部搜索
            tad_anchor_start = tad_start - Config.MAX_DISTANCE_KB * 1000
            tad_anchor_end = tad_end
            tad_anchor_side = "upstream_and_internal"
        else:
            return {"error": "Gene overlaps with TAD in Branch 2"}

        # 搜索匹配的Loop
        for loop_id, loop_row in _G['loop_pos_df'].iterrows():
            # 情况1：anchor1与基因重叠，anchor2与TAD重叠
            anchor1_gene_overlap = check_overlap(gene_anchor_start, gene_anchor_end,
                                                  loop_row['start1'], loop_row['end1'])
            anchor2_tad_overlap = check_overlap(tad_anchor_start, tad_anchor_end,
                                               loop_row['start2'], loop_row['end2'])

            # 情况2：anchor2与基因重叠，anchor1与TAD重叠
            anchor2_gene_overlap = check_overlap(gene_anchor_start, gene_anchor_end,
                                                  loop_row['start2'], loop_row['end2'])
            anchor1_tad_overlap = check_overlap(tad_anchor_start, tad_anchor_end,
                                               loop_row['start1'], loop_row['end1'])

            if anchor1_gene_overlap and anchor2_tad_overlap:
                return {
                    "loop_id": loop_id,
                    "gene_anchor_start": gene_anchor_start,
                    "gene_anchor_end": gene_anchor_end,
                    "tad_anchor_start": tad_anchor_start,
                    "tad_anchor_end": tad_anchor_end,
                    "gene_position": gene_position,
                    "tad_anchor_side": tad_anchor_side,
                    "loop_anchor1_start": loop_row['start1'],
                    "loop_anchor1_end": loop_row['end1'],
                    "loop_anchor2_start": loop_row['start2'],
                    "loop_anchor2_end": loop_row['end2'],
                    "match_type": "anchor1_gene_anchor2_tad"
                }
            elif anchor2_gene_overlap and anchor1_tad_overlap:
                return {
                    "loop_id": loop_id,
                    "gene_anchor_start": gene_anchor_start,
                    "gene_anchor_end": gene_anchor_end,
                    "tad_anchor_start": tad_anchor_start,
                    "tad_anchor_end": tad_anchor_end,
                    "gene_position": gene_position,
                    "tad_anchor_side": tad_anchor_side,
                    "loop_anchor1_start": loop_row['start1'],
                    "loop_anchor1_end": loop_row['end1'],
                    "loop_anchor2_start": loop_row['start2'],
                    "loop_anchor2_end": loop_row['end2'],
                    "match_type": "anchor2_gene_anchor1_tad"
                }
        
        return None
    except Exception as e:
        return {"error": str(e)}

# --------------------------------------------------------------------------
# --- 方向性检验函数 ---
# --------------------------------------------------------------------------
def check_branch1_directionality(tad_is_median_0, tad_is_median_2, gene_exp_median_0, gene_exp_median_2):
    """
    分支1方向性检验：TAD IS与基因表达呈相反趋势
    """
    # 检查所有值是否有效
    if any(pd.isna(x) for x in [tad_is_median_0, tad_is_median_2, gene_exp_median_0, gene_exp_median_2]):
        return False
    
    # 相反趋势检验
    return ((tad_is_median_0 > tad_is_median_2 and gene_exp_median_0 < gene_exp_median_2) or
            (tad_is_median_0 < tad_is_median_2 and gene_exp_median_0 > gene_exp_median_2))

def check_branch2_directionality(tad_is_median_0, tad_is_median_2, 
                                gene_exp_median_0, gene_exp_median_2,
                                loop_strength_median_0, loop_strength_median_2,
                                loop_count_0, loop_count_2):
    """
    分支2方向性检验：所有指标呈一致趋势
    """
    # 检查所有值是否有效
    values = [tad_is_median_0, tad_is_median_2, gene_exp_median_0, gene_exp_median_2,
              loop_strength_median_0, loop_strength_median_2, loop_count_0, loop_count_2]
    if any(pd.isna(x) for x in values):
        return False
    
    # 一致趋势检验
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
        # 检查数据存在性
        if full_snp_id not in _G["snp_data"].index:
            return None, f"SNP not found: {full_snp_id}"
        if gene_id not in _G["gene_exp"].index:
            return None, f"Gene not found: {gene_id}"
        if tad_id not in _G["tad_is"].index:
            return None, f"TAD not found in tad_is: {tad_id}"
        
        # 加载数据
        snp_values = _G["snp_data"].loc[full_snp_id, _G["common_samples"]].astype(float)
        gene_values = _G["gene_exp"].loc[gene_id, _G["common_samples"]].astype(float)
        tad_is_values = _G["tad_is"].loc[tad_id, _G["common_samples"]].astype(float)
        
        df = pd.DataFrame({
            "Genotype_num": snp_values,
            "TAD_IS_Score": tad_is_values,
            "Gene_Exp": gene_values
        })
        
        # 过滤期望基因型
        df = df[df["Genotype_num"].isin(Config.EXPECTED_GENOTYPES)]
        
        # 验证数据质量
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
# --- 分析函数 ---
# --------------------------------------------------------------------------
def analyze_branch1(data_df, gene_start, gene_end, tad_start, tad_end):
    """分析分支1"""
    # 距离范围
    branch1_start = tad_start - Config.DISTANCE_THRESHOLD_KB * 1000
    branch1_end = tad_end + Config.DISTANCE_THRESHOLD_KB * 1000
    
    # 距离检查
    distance_pass = check_overlap(gene_start, gene_end, branch1_start, branch1_end)
    
    result = {
        "Branch1_Range_Start": branch1_start,
        "Branch1_Range_End": branch1_end,
        "Branch1_Distance_Pass": distance_pass
    }
    
    if not distance_pass:
        result["Failure_Reason"] = f"Branch1: Gene not within {Config.DISTANCE_THRESHOLD_KB}kb of TAD"
        return {"passed": False, "stats": result}
    
    # 计算中位数
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
    
    # 方向性检查
    direction_pass = check_branch1_directionality(median_is_0, median_is_2, median_exp_0, median_exp_2)
    result["Branch1_Direction_Pass"] = direction_pass
    
    if not direction_pass:
        result["Failure_Reason"] = "Branch1: Direction criteria not met"
        return {"passed": False, "stats": result}
    
    return {"passed": True, "stats": result}

def analyze_branch2(data_df, gene_start, gene_end, tad_start, tad_end):
    """分析分支2"""
    # 距离范围
    upstream_start = tad_start - Config.MAX_DISTANCE_KB * 1000
    upstream_end = tad_start - Config.DISTANCE_THRESHOLD_KB * 1000
    downstream_start = tad_end + Config.DISTANCE_THRESHOLD_KB * 1000
    downstream_end = tad_end + Config.MAX_DISTANCE_KB * 1000
    
    # 距离检查
    distance_pass = (check_overlap(gene_start, gene_end, upstream_start, upstream_end) or
                    check_overlap(gene_start, gene_end, downstream_start, downstream_end))
    
    result = {
        "Branch2_Upstream_Range_Start": upstream_start,
        "Branch2_Upstream_Range_End": upstream_end,
        "Branch2_Downstream_Range_Start": downstream_start,
        "Branch2_Downstream_Range_End": downstream_end,
        "Branch2_Distance_Pass": distance_pass
    }
    
    if not distance_pass:
        result["Failure_Reason"] = f"Branch2: Gene not in {Config.DISTANCE_THRESHOLD_KB}-{Config.MAX_DISTANCE_KB}kb range"
        return {"passed": False, "stats": result}
    
    # Loop搜索
    loop_info = find_overlapping_loop_for_branch2(gene_start, gene_end, tad_start, tad_end)
    
    if not loop_info or "loop_id" not in loop_info:
        result["Failure_Reason"] = "Branch2: No overlapping loop found"
        return {"passed": False, "stats": result}
    
    result["Branch2_Loop_Found"] = True
    result.update({
        "Loop_ID": loop_info["loop_id"],
        "Loop_Match_Type": loop_info["match_type"],
        "Gene_Anchor_Start": loop_info["gene_anchor_start"],
        "Gene_Anchor_End": loop_info["gene_anchor_end"],
        "TAD_Anchor_Start": loop_info["tad_anchor_start"],
        "TAD_Anchor_End": loop_info["tad_anchor_end"],
        "Gene_Position": loop_info["gene_position"],
        "TAD_Anchor_Side": loop_info["tad_anchor_side"],
        "Loop_Anchor1_Start": loop_info["loop_anchor1_start"],
        "Loop_Anchor1_End": loop_info["loop_anchor1_end"],
        "Loop_Anchor2_Start": loop_info["loop_anchor2_start"],
        "Loop_Anchor2_End": loop_info["loop_anchor2_end"]
    })
    
    # Loop显著性检验
    try:
        loop_strength_values = _G["loop_strength_df"].loc[loop_info["loop_id"], _G["common_samples"]].astype(float)
        loop_01_values = _G["loop_01_df"].loc[loop_info["loop_id"], _G["common_samples"]].astype(float)
        
        data_df = data_df.copy()
        data_df["Loop_Strength"] = loop_strength_values
        data_df["Loop_Present"] = loop_01_values
        
        # Mann-Whitney检验
        loop_strength_test = perform_mannwhitney_test(data_df, "Loop_Strength", "Genotype_num")
        result.update({
            "Loop_Strength_P": loop_strength_test["p_value"],
            "Loop_Strength_Stat": loop_strength_test.get("stat", np.nan),
            "Loop_Strength_N0": loop_strength_test.get("group0_size", 0),
            "Loop_Strength_N2": loop_strength_test.get("group2_size", 0),
            "Loop_Strength_Median0": loop_strength_test.get("group0_median", np.nan),
            "Loop_Strength_Median2": loop_strength_test.get("group2_median", np.nan)
        })
        
        # 独立性检验
        loop_01_test = perform_independence_test(data_df, "Loop_Present", "Genotype_num")
        result.update({
            "Loop_01_P": loop_01_test["p_value"],
            "Loop_01_Test_Method": loop_01_test.get("test_used", "unknown"),
            "Loop_01_Test_Reason": loop_01_test.get("reason", ""),
            "Loop_01_Sample_Size": loop_01_test.get("sample_size", 0),
            "Loop_01_Expected_Min": loop_01_test.get("expected_min", np.nan)
        })
        
        # 显著性检查
        loop_significant = ((pd.notna(result["Loop_Strength_P"]) and 
                            result["Loop_Strength_P"] < Config.P_THRESHOLD_LOOP) or
                           (pd.notna(result["Loop_01_P"]) and 
                            result["Loop_01_P"] < Config.P_THRESHOLD_LOOP))
        
        result["Branch2_Loop_Significant"] = loop_significant
        
        # 记录显著性信息
        sig_tests = []
        if pd.notna(result["Loop_Strength_P"]) and result["Loop_Strength_P"] < Config.P_THRESHOLD_LOOP:
            sig_tests.append(f"strength_p={result['Loop_Strength_P']:.4f}")
        if pd.notna(result["Loop_01_P"]) and result["Loop_01_P"] < Config.P_THRESHOLD_LOOP:
            sig_tests.append(f"presence_p={result['Loop_01_P']:.4f}")
        
        result["Loop_Significant_Tests"] = "; ".join(sig_tests) if sig_tests else "none"
        
        if not loop_significant:
            result["Failure_Reason"] = f"Branch2: Loop not significant ({result['Loop_Significant_Tests']})"
            return {"passed": False, "stats": result}
        
        # 计算所有中位数
        median_is_0 = data_df[data_df["Genotype_num"] == 0]["TAD_IS_Score"].median()
        median_is_2 = data_df[data_df["Genotype_num"] == 2]["TAD_IS_Score"].median()
        median_exp_0 = data_df[data_df["Genotype_num"] == 0]["Gene_Exp"].median()
        median_exp_2 = data_df[data_df["Genotype_num"] == 2]["Gene_Exp"].median()
        median_loop_0 = data_df[data_df["Genotype_num"] == 0]["Loop_Strength"].median()
        median_loop_2 = data_df[data_df["Genotype_num"] == 2]["Loop_Strength"].median()
        count_loop_0 = data_df[data_df["Genotype_num"] == 0]["Loop_Present"].sum()
        count_loop_2 = data_df[data_df["Genotype_num"] == 2]["Loop_Present"].sum()

        result.update({
            "TAD_IS_Median_Group0": median_is_0,
            "TAD_IS_Median_Group2": median_is_2,
            "Gene_Exp_Median_Group0": median_exp_0,
            "Gene_Exp_Median_Group2": median_exp_2,
            "Loop_Strength_Median_Group0": median_loop_0,
            "Loop_Strength_Median_Group2": median_loop_2,
            "Loop_Count_Group0": count_loop_0,
            "Loop_Count_Group2": count_loop_2
        })

        # 方向性检查
        direction_pass = check_branch2_directionality(median_is_0, median_is_2, median_exp_0, median_exp_2,
                                                     median_loop_0, median_loop_2, count_loop_0, count_loop_2)
        result["Branch2_Direction_Pass"] = direction_pass
        
        if not direction_pass:
            result["Failure_Reason"] = "Branch2: Direction criteria not met"
            return {"passed": False, "stats": result}
        
        return {"passed": True, "stats": result}
        
    except KeyError as e:
        result["Failure_Reason"] = f"Loop {loop_info['loop_id']} data not found: {str(e)}"
        return {"passed": False, "stats": result}
    except Exception as e:
        result["Failure_Reason"] = f"Loop analysis error: {str(e)}"
        return {"passed": False, "stats": result}

def initialize_stats_dict(gene_id, tad_id, snp_id, full_snp_id,
                         gene_chrom, gene_start, gene_end,
                         tad_chrom, tad_start, tad_end,
                         gene_tad_distance, gene_tad_relation):
    """初始化统计结果字典"""
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
    """分析单个三元组的主函数"""
    # 初始化stats以便在异常时也能返回
    stats = {}
    
    try:
        # 提取基础信息
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

        # 计算距离和关系
        gene_tad_distance, gene_tad_relation = calculate_gene_tad_distance_and_relation(
            gene_start, gene_end, tad_start, tad_end
        )
        
        # 初始化统计
        stats = initialize_stats_dict(gene_id, tad_id, snp_id, full_snp_id,
                                     gene_chrom, gene_start, gene_end,
                                     tad_chrom, tad_start, tad_end,
                                     gene_tad_distance, gene_tad_relation)
        
        # 加载数据
        data_df, failure_reason = load_and_validate_data(full_snp_id, gene_id, tad_id)
        if data_df is None:
            stats["Failure_Reason"] = failure_reason
            return {"pass": False, "stats": stats}
        
        # 记录样本量
        stats["Sample_Count_Group0"] = len(data_df[data_df["Genotype_num"] == 0])
        stats["Sample_Count_Group2"] = len(data_df[data_df["Genotype_num"] == 2])
        
        # 分支1分析
        branch1_result = analyze_branch1(data_df, gene_start, gene_end, tad_start, tad_end)
        if branch1_result["passed"]:
            stats.update(branch1_result["stats"])
            stats["Branch_Pass"] = "Branch_1_Pass"
            return {"pass": True, "stats": stats}
        else:
            stats.update(branch1_result["stats"])
        
        # 分支2分析
        branch2_result = analyze_branch2(data_df, gene_start, gene_end, tad_start, tad_end)
        if branch2_result["passed"]:
            stats.update(branch2_result["stats"])
            stats["Branch_Pass"] = "Branch_2_Pass"
            return {"pass": True, "stats": stats}
        else:
            stats.update(branch2_result["stats"])
            if not stats.get("Failure_Reason"):
                stats["Failure_Reason"] = "Did not pass any branch criteria"
            return {"pass": False, "stats": stats}
        
    except Exception as e:
        if not stats:
            stats = {"Failure_Reason": f"Critical error: {str(e)}"}
        else:
            stats["Failure_Reason"] = f"Unexpected error: {str(e)}"
        return {"pass": False, "stats": stats}

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

        # 清理列名
        for df in data_files.values():
            if not df.columns.empty:
                df.columns = df.columns.str.strip()
        
        _G.update(data_files)
        _G["loop_strength_df"] = data_files["loop_strength"]
        _G["loop_01_df"] = data_files["loop_01"]
        _G["loop_pos_df"] = data_files["loop_pos"]
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
# --- 结果保存和报告 ---
# --------------------------------------------------------------------------
def save_results(all_stats, passed_stats, output_dir):
    """保存分析结果"""
    if not all_stats:
        print("没有结果数据", file=sys.stderr)
        return
    
    all_stats_df = pd.DataFrame(all_stats)
    
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
        "Loop_ID", "Loop_Match_Type", "Gene_Anchor_Start", "Gene_Anchor_End", "Gene_Anchor_Distance_kb",
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
    complete_file = output_path / "complete_analysis_results.csv"
    all_stats_df.to_csv(complete_file, sep=",", index=False)
    print(f"完整结果已保存: {complete_file.absolute()}", file=sys.stderr)
    
    if passed_stats:
        passed_df = pd.DataFrame(passed_stats).reindex(columns=output_cols)
        passed_file = output_path / "passed_triplets_results.csv"
        passed_df.to_csv(passed_file, sep=",", index=False)
        print(f"通过筛选结果已保存: {passed_file.absolute()}", file=sys.stderr)

def print_summary_statistics(all_stats, passed_stats):
    """输出汇总统计"""
    print("-" * 60, file=sys.stderr)
    print("分析结果汇总:", file=sys.stderr)
    
    if passed_stats:
        passed_df = pd.DataFrame(passed_stats)
        print("\n通过各分支的数量:", file=sys.stderr)
        branch_counts = passed_df["Branch_Pass"].value_counts()
        for branch, count in branch_counts.items():
            print(f"  {branch}: {count}", file=sys.stderr)
    
    # 失败原因
    all_df = pd.DataFrame(all_stats)
    failed_df = all_df[all_df["Branch_Pass"] == "None"]
    if len(failed_df) > 0:
        print("\n失败原因统计:", file=sys.stderr)
        failure_counts = failed_df["Failure_Reason"].value_counts()
        for reason, count in failure_counts.head(10).items():
            print(f"  {reason}: {count}", file=sys.stderr)
    
    # 统计方法使用情况
    if len(all_df) > 0:
        print("\n统计检验方法使用:", file=sys.stderr)
        
        loop_methods = all_df[all_df["Loop_01_Test_Method"].notna()]["Loop_01_Test_Method"].value_counts()
        if not loop_methods.empty:
            total_tests = loop_methods.sum()
            for method, count in loop_methods.items():
                pct = count / total_tests * 100
                print(f"  {method}: {count}/{total_tests} ({pct:.1f}%)", file=sys.stderr)
    
    print("=" * 60, file=sys.stderr)

# --------------------------------------------------------------------------
# --- 主函数 ---
# --------------------------------------------------------------------------
def main():
    """主函数"""
    args = parse_args()
    
    # 设置配置
    Config.DISTANCE_THRESHOLD_KB = args.distance_threshold
    Config.MAX_DISTANCE_KB = args.max_distance
    Config.GENE_ANCHOR_DISTANCE_KB = args.gene_anchor_distance
    Config.P_THRESHOLD_LOOP = args.loop_p_threshold
    
    # 输出配置
    print("=" * 60, file=sys.stderr)
    print("严谨SNP-TAD-Gene分析开始", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"配置参数:", file=sys.stderr)
    print(f"  距离临界值: {Config.DISTANCE_THRESHOLD_KB} kb", file=sys.stderr)
    print(f"  最远距离: {Config.MAX_DISTANCE_KB} kb", file=sys.stderr)
    print(f"  基因anchor距离: {Config.GENE_ANCHOR_DISTANCE_KB} kb", file=sys.stderr)
    print(f"  Loop p值阈值: {Config.P_THRESHOLD_LOOP}", file=sys.stderr)
    print(f"  输出目录: {args.output_dir}", file=sys.stderr)
    print("-" * 60, file=sys.stderr)
    
    # 加载三元组
    triplets_df = pd.read_csv(args.triplets_with_positions)
    print(f"加载 {len(triplets_df)} 个三元组", file=sys.stderr)
    
    triplets_list = convert_triplets_to_dicts(triplets_df)
    
    # 多进程分析
    results = []
    with Pool(processes=args.num_cores, initializer=_init_globals, initargs=(args,)) as pool:
        for result in tqdm(pool.imap_unordered(analyze_single_triplet, triplets_list, chunksize=10),
                          total=len(triplets_list), desc="分析三元组", file=sys.stderr):
            results.append(result)

    # 整理结果
    all_stats = [r["stats"] for r in results]
    passed_stats = [r["stats"] for r in results if r["pass"]]
    
    print(f"\n分析完成: {len(all_stats)} 总数, {len(passed_stats)} 通过", file=sys.stderr)
    
    # 保存和报告
    save_results(all_stats, passed_stats, args.output_dir)
    print_summary_statistics(all_stats, passed_stats)

if __name__ == "__main__":
    main()
