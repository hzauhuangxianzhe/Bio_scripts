#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP-TAD-Gene Expression-Phenotype Analysis Core (v2) - Final Version
Author: Based on user requirements
Date: 2025
Version: 2.0.5 - Added comprehensive failure logging and reporting.
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
from collections import Counter # <-- 新增：用于统计失败原因

# ------------------------- Warning policy (D) -------------------------
warnings.filterwarnings(
    "ignore",
    message=r"divide by zero encountered in.*|invalid value encountered in.*",
    category=RuntimeWarning
)

# --------------------------------------------------------------------------
# --- CLI ---
# --------------------------------------------------------------------------
def parse_args():
    parser = ArgumentParser(description="SNP-TAD-Gene-Phenotype statistical analysis.")
    # Required
    parser.add_argument("triplets", help="Path to triplets file (Gene_ID, TAD_ID, SNP_ID)")
    parser.add_argument("snp_data", help="Path to SNP genotype data file")
    parser.add_argument("gene_exp", help="Path to gene expression data file")
    parser.add_argument("tad_is", help="Path to TAD IS data file")
    parser.add_argument("tad_01", help="Path to TAD 0-1 matrix file")
    parser.add_argument("phenotype", help="Path to phenotype data file")
    parser.add_argument("num_cores", type=int, help="Number of CPU cores for parallel processing")
    # Optional
    parser.add_argument("--output_dir", default="analysis_output", help="Output directory")
    parser.add_argument("--tad_is_p_threshold", type=float, default=0.05,
                        help="P-value threshold for TAD IS significance (default: 0.05)")
    parser.add_argument("--gene_exp_p_threshold", type=float, default=0.05,
                        help="P-value threshold for Gene Expression significance (default: 0.05)")
    return parser.parse_args()

# --------------------------------------------------------------------------
# --- Config ---
# --------------------------------------------------------------------------
class Config:
    EXPECTED_GENOTYPES = [0, 2]

# Global (for multiprocessing initializer)
_G = {}

# --------------------------------------------------------------------------
# --- Helpers ---
# --------------------------------------------------------------------------
def build_snp_index_map(snp_index: pd.Index):
    """
    Build O(1) map from 'SNP_ID' (prefix before first ':') to the full row index.
    Warn if duplicate prefixes exist; keep first occurrence.
    """
    s = snp_index.to_series().astype(str)
    base = s.str.split(":", n=1, expand=True).iloc[:, 0].astype(str)
    vc = base.value_counts()
    dups = vc[vc > 1].index.tolist()
    if dups:
        print(f"Warning: duplicated SNP index prefix found for {len(dups)} keys (showing up to 10): {dups[:10]}",
              file=sys.stderr)
    first_occ_mask = ~base.duplicated()
    return dict(zip(base[first_occ_mask], s[first_occ_mask]))

def _prep_2x2(sub: pd.DataFrame, categorical_col: str, group_col: str):
    """
    Force 2x2 crosstab to rows=[0,1], cols=[0,2]; return None if any row/col margin is zero.
    (A1)
    """
    ct = pd.crosstab(sub[categorical_col], sub[group_col])
    ct = ct.reindex(index=[0, 1], columns=[0, 2], fill_value=0)
    if (ct.sum(axis=1) == 0).any() or (ct.sum(axis=0) == 0).any():
        return None
    return ct

# --------------------------------------------------------------------------
# --- Statistical Tests ---
# --------------------------------------------------------------------------
def perform_mannwhitney_test(data, value_col, group_col):
    """
    Pairwise dropna per test; return p and sample sizes. (A2)
    """
    try:
        sub = data[[value_col, group_col]].dropna()
        if len(sub) < 3 or sub[group_col].nunique() < 2:
            return {"p_value": np.nan, "n0": 0, "n2": 0}
        groups = sorted(sub[group_col].unique())
        g0 = sub[sub[group_col] == groups[0]][value_col]
        g2 = sub[sub[group_col] == groups[1]][value_col]
        stat, p_value = mannwhitneyu(g0, g2, alternative='two-sided')
        return {"p_value": float(p_value), "n0": int(g0.size), "n2": int(g2.size)}
    except Exception:
        return {"p_value": np.nan, "n0": 0, "n2": 0}

def perform_independence_test(data, categorical_col, group_col):
    """
    2x2 independence test with expected-frequency-based routing. (A3)
    Rule: if any expected <1 or any expected <5 => Fisher (two-sided); else Pearson chi-square (no Yates).
    """
    try:
        sub = data[[categorical_col, group_col]].dropna()
        if len(sub) < 4:
            return {"test_used": "none", "p_value": np.nan, "reason": "insufficient_data"}

        ct = _prep_2x2(sub, categorical_col, group_col)
        if ct is None:
            return {"test_used": "none", "p_value": np.nan, "reason": "no_variation_in_table"}

        table = ct.values
        # expected from chi2 (Yates choice does not affect expected)
        chi2_y, p_y, dof, expected = chi2_contingency(table, correction=True)
        expected = np.asarray(expected)
        Emin = float(expected.min())
        n_lt5 = int((expected < 5).sum())
        any_lt1 = bool((expected < 1).any())

        if any_lt1 or n_lt5 > 0:
            _, p = fisher_exact(table, alternative="two-sided")
            return {
                "test_used": "fisher_exact",
                "p_value": float(p),
                "reason": f"expected_min={Emin:.2f}; n_cells_lt5={n_lt5}"
            }

        chi2, p_chi2, dof, _ = chi2_contingency(table, correction=False)  # Pearson, no Yates (power-friendly)
        return {
            "test_used": "chi2_pearson",
            "p_value": float(p_chi2),
            "reason": f"expected_min={Emin:.2f}; all_expected_ge_5"
        }

    except Exception as e:
        return {"test_used": "error", "p_value": np.nan, "reason": f"error_{str(e)}"}

# --------------------------------------------------------------------------
# --- Data IO & Validation ---
# --------------------------------------------------------------------------
def load_all_data(paths):
    print("Loading data files...", file=sys.stderr)
    data_list = {
        "triplets": pd.read_csv(paths["triplets"], sep='\t', header=0, names=["Gene_ID", "TAD_ID", "SNP_ID"]),
        "snp_data": pd.read_csv(paths["snp_data"], sep='\t', header=0, index_col=0),
        "gene_exp": pd.read_csv(paths["gene_exp"], sep='\t', header=0, index_col=0),
        "tad_is": pd.read_csv(paths["tad_is"], sep='\t', header=0, index_col=0),
        "tad_01": pd.read_csv(paths["tad_01"], sep='\t', header=0, index_col=0),
        "phenotype": pd.read_csv(paths["phenotype"], sep='\t', header=0, index_col=0)
    }
    print("Data loading complete.", file=sys.stderr)
    return data_list

def validate_sample_consistency(data_list):
    print("Validating sample consistency...", file=sys.stderr)
    file_samples = {
        "snp_data": set(data_list["snp_data"].columns),
        "gene_exp": set(data_list["gene_exp"].columns),
        "tad_is": set(data_list["tad_is"].columns),
        "tad_01": set(data_list["tad_01"].columns),
        "phenotype": set(data_list["phenotype"].columns)
    }
    reference_samples = file_samples["snp_data"]
    for file_name, samples in file_samples.items():
        if file_name == "snp_data":
            continue
        if samples != reference_samples:
            only_in_ref = reference_samples - samples
            only_in_current = samples - reference_samples
            error_msg = (f"Sample mismatch between snp_data and {file_name}:\n"
                         f"  snp_data: {len(reference_samples)} samples\n"
                         f"  {file_name}: {len(samples)} samples\n")
            if only_in_ref:
                error_msg += f"  Only in snp_data: {sorted(list(only_in_ref))[:10]}...\n"
            if only_in_current:
                error_msg += f"  Only in {file_name}: {sorted(list(only_in_current))[:10]}...\n"
            raise ValueError(error_msg)
    common_samples = sorted(list(reference_samples))
    print(f"Sample validation successful: {len(common_samples)} identical samples", file=sys.stderr)
    return common_samples

# --------------------------------------------------------------------------
# --- Worker & Globals (B1/B2 + A4) ---
# --------------------------------------------------------------------------
def extract_sample_data(snp_row_name, gene_id, tad_id):
    try:
        snp_values = _G["snp_data"].loc[snp_row_name, _G["common_samples"]].astype(float)
        gene_values = _G["gene_exp"].loc[gene_id, _G["common_samples"]].astype(float)
        tad_is_values = _G["tad_is"].loc[tad_id, _G["common_samples"]].astype(float)
        tad_01_values = _G["tad_01"].loc[tad_id, _G["common_samples"]].astype(float)
        fl_values = _G["phenotype"].loc["FL", _G["common_samples"]].astype(float)
        fs_values = _G["phenotype"].loc["FS", _G["common_samples"]].astype(float)
        return {
            "snp": snp_values,
            "gene": gene_values,
            "tad_is": tad_is_values,
            "tad_01": tad_01_values,
            "fl": fl_values,
            "fs": fs_values
        }
    except Exception:
        return None

def _process_triplet_worker(triplet_data):
    # (*** 已更新：添加了结构化的失败返回 ***)
    gene_id, tad_id, snp_id = triplet_data
    tad_thr = _G["tad_thr"]
    gene_thr = _G["gene_thr"]

    try:
        # 1. 检查 SNP
        snp_row_name = _G["SNP_MAP"].get(str(snp_id))  # O(1) lookup (B2)
        if snp_row_name is None:
            return {"success": False, "input_triplet": triplet_data, "message": f"SNP prefix not found in matrix index: {snp_id}"}

        # 2. 检查 Gene
        if gene_id not in _G["gene_exp"].index:
            return {"success": False, "input_triplet": triplet_data, "message": f"Gene_ID not found in gene_exp index: {gene_id}"}
        
        # 3. 检查 TAD
        if tad_id not in _G["tad_is"].index:
            return {"success": False, "input_triplet": triplet_data, "message": f"TAD_ID not found in tad_is index: {tad_id}"}
        if tad_id not in _G["tad_01"].index:
            return {"success": False, "input_triplet": triplet_data, "message": f"TAD_ID not found in tad_01 index: {tad_id}"}

        samples = extract_sample_data(snp_row_name, gene_id, tad_id)
        if samples is None:
            # 这通常是一个KeyError，如果extract_sample_data内部失败（不太可能发生，因为我们已经检查了所有ID）
            return {"success": False, "input_triplet": triplet_data, "message": "Data extraction failed (Internal KeyError)"}

        df = pd.DataFrame({
            "Genotype_num": samples["snp"],
            "TAD_IS_Score": samples["tad_is"],
            "TAD_Present": samples["tad_01"],
            "Gene_Exp": samples["gene"],
            "FL": samples["fl"],
            "FS": samples["fs"]
        })

        # Binary normalization & filter (A4)
        df = df[df["Genotype_num"].isin(Config.EXPECTED_GENOTYPES)]
        df["TAD_Present"] = (df["TAD_Present"] > 0).astype(int)

        # 4. 检查基因型分布
        uniq = sorted(df["Genotype_num"].dropna().unique())
        if len(uniq) != 2 or not all(g in Config.EXPECTED_GENOTYPES for g in uniq):
            return {"success": False, "input_triplet": triplet_data, "message": f"Data insufficient: requires both genotypes [0, 2], but found: {uniq}"}

        stats = {
            "Gene_ID": gene_id,
            "TAD_ID": tad_id,
            "SNP_ID": snp_id,
            "Full_SNP_ID": snp_row_name
        }

        # Mann-Whitney (A2)
        t_tad = perform_mannwhitney_test(df, "TAD_IS_Score", "Genotype_num")
        stats["MannWhitney_P_TAD"] = t_tad["p_value"]

        t_gene = perform_mannwhitney_test(df, "Gene_Exp", "Genotype_num")
        stats["MannWhitney_P_Gene"] = t_gene["p_value"]

        t_fl = perform_mannwhitney_test(df, "FL", "Genotype_num")
        stats["MannWhitney_P_FL"] = t_fl["p_value"]

        t_fs = perform_mannwhitney_test(df, "FS", "Genotype_num")
        stats["MannWhitney_P_FS"] = t_fs["p_value"]

        # 2x2 independence (A1/A3)
        indep = perform_independence_test(df, "TAD_Present", "Genotype_num")
        stats["Independence_P_TAD01_Genotype"] = indep["p_value"]
        stats["Independence_Test_Method"] = indep["test_used"]
        stats["Independence_Test_Reason"] = indep["reason"]

        # Significance flags (C2: thresholds respected)
        stats["TAD_IS_Significant"] = (pd.notna(stats["MannWhitney_P_TAD"])
                                      and stats["MannWhitney_P_TAD"] < tad_thr)
        stats["Gene_Exp_Significant"] = (pd.notna(stats["MannWhitney_P_Gene"])
                                        and stats["MannWhitney_P_Gene"] < gene_thr)
        stats["Both_Significant"] = bool(stats["TAD_IS_Significant"] and stats["Gene_Exp_Significant"])

        return {"success": True, "stats": stats}

    except Exception as e:
        # 5. 捕获任何其他意外错误
        return {"success": False, "input_triplet": triplet_data, "message": f"Unexpected processing error: {str(e)}"}

# Initializer (B1)
def _init_globals(snp_data, gene_exp, tad_is, tad_01, phenotype, common_samples,
                  tad_thr, gene_thr):
    _G["snp_data"] = snp_data
    _G["gene_exp"] = gene_exp
    _G["tad_is"] = tad_is
    _G["tad_01"] = tad_01
    _G["phenotype"] = phenotype
    _G["common_samples"] = common_samples
    _G["tad_thr"] = tad_thr
    _G["gene_thr"] = gene_thr
    _G["SNP_MAP"] = build_snp_index_map(snp_data.index)

# --------------------------------------------------------------------------
# --- Orchestrators ---
# --------------------------------------------------------------------------
def log_first_triplet_example(data_list, common_samples, args):
    print("\n" + "="*80, file=sys.stderr)
    print("STATISTICAL ANALYSIS METHODS DEMONSTRATION", file=sys.stderr)
    print("="*80, file=sys.stderr)
    first_triplet = data_list["triplets"].iloc[0]
    gene_id, tad_id, snp_id = first_triplet["Gene_ID"], first_triplet["TAD_ID"], first_triplet["SNP_ID"]
    
    # (*** 已更新：使用新的结构化返回 ***)
    r = _process_triplet_worker([gene_id, tad_id, snp_id])
    if r["success"]:
        s = r["stats"]
        print(f"Example triplet: {gene_id} - {tad_id} - {snp_id}", file=sys.stderr)
        print(f"1. Mann-Whitney (TAD IS vs Genotype): p={s['MannWhitney_P_TAD']}", file=sys.stderr)
        print(f"2. Mann-Whitney (Gene Exp vs Genotype): p={s['MannWhitney_P_Gene']}", file=sys.stderr)
        print(f"3. Mann-Whitney (FL vs Genotype): p={s['MannWhitney_P_FL']}", file=sys.stderr)
        print(f"4. Mann-Whitney (FS vs Genotype): p={s['MannWhitney_P_FS']}", file=sys.stderr)
        print(f"5. 2x2 Independence (TAD01 vs Genotype): method={s['Independence_Test_Method']}, "
              f"p={s['Independence_P_TAD01_Genotype']}", file=sys.stderr)
        print(f"   Reason: {s['Independence_Test_Reason']}", file=sys.stderr)
        print(f"Combined: Both significant = {s['Both_Significant']}", file=sys.stderr)
    else:
        # (*** 已更新：打印来自新返回结构的特定消息 ***)
        print(f"Failed to process example triplet ({gene_id} - {tad_id} - {snp_id}):", file=sys.stderr)
        print(f"  REASON: {r.get('message', 'Unknown failure')}", file=sys.stderr)
    print("="*80, file=sys.stderr)

def parallel_triplet_analysis(data_list, common_samples, args):
    triplets_list = data_list["triplets"].values.tolist()
    total = len(triplets_list)
    print(f"Starting analysis of {total} triplets using {args.num_cores} cores", file=sys.stderr)

    results = []
    with Pool(
        processes=args.num_cores,
        initializer=_init_globals,
        initargs=(data_list["snp_data"], data_list["gene_exp"], data_list["tad_is"],
                  data_list["tad_01"], data_list["phenotype"], common_samples,
                  args.tad_is_p_threshold, args.gene_exp_p_threshold)
    ) as pool:
        chunk = max(1, total // (args.num_cores * 8))  # scheduler-friendly
        for out in tqdm(pool.imap_unordered(_process_triplet_worker, triplets_list, chunksize=chunk),
                        total=total, desc="Processing triplets", file=sys.stderr):
            results.append(out)
    print("Parallel processing completed", file=sys.stderr)
    return results

# --------------------------------------------------------------------------
# --- Main ---
# --------------------------------------------------------------------------
def main():
    args = parse_args()
    print("SNP-TAD-Gene Expression-Phenotype Analysis Starting", file=sys.stderr)

    try:
        paths = {
            "triplets": args.triplets,
            "snp_data": args.snp_data,
            "gene_exp": args.gene_exp,
            "tad_is": args.tad_is,
            "tad_01": args.tad_01,
            "phenotype": args.phenotype
        }
        for name, path in paths.items():
            if not os.path.exists(path):
                raise FileNotFoundError(f"{name} file not found: {path}")

        data_list = load_all_data(paths)
        common_samples = validate_sample_consistency(data_list)

        # Init globals in main and demonstrate one triplet
        _init_globals(data_list["snp_data"], data_list["gene_exp"], data_list["tad_is"],
                      data_list["tad_01"], data_list["phenotype"], common_samples,
                      args.tad_is_p_threshold, args.gene_exp_p_threshold)
        log_first_triplet_example(data_list, common_samples, args)

        # Parallel
        results = parallel_triplet_analysis(data_list, common_samples, args)

        # (*** 核心修改：添加失败收集和报告 ***)
        
        # 1. 结果分类
        print("Consolidating results...", file=sys.stderr)
        success_count = 0
        all_stats = []
        all_failures = [] # <-- 新增：失败列表
        
        for r in results:
            if r["success"]:
                success_count += 1
                all_stats.append(r["stats"])
            else:
                all_failures.append(r) # <-- 新增：收集失败的返回对象

        print(f"Successfully processed {success_count}/{len(results)} triplets", file=sys.stderr)

        # 2. 报告失败（如果存在）
        if all_failures:
            print("\n" + "="*80, file=sys.stderr)
            print(f"FAILURE ANALYSIS: {len(all_failures)} out of {len(results)} triplets failed.", file=sys.stderr)
            print("Summary of failure reasons:", file=sys.stderr)
            
            # 提取所有标准化的错误消息
            failure_messages = [f.get("message", "Unknown reason") for f in all_failures]
            
            # 统计每种错误消息的频率
            reason_counts = Counter(failure_messages)
            
            # 按最常见的失败原因排序并打印
            for reason, count in reason_counts.most_common():
                print(f"  - [{count} triplets failed] Reason: {reason}", file=sys.stderr)
            
            print("="*80, file=sys.stderr)
            
            # (可选：如果您想保存详细的失败日志)
            # failure_log_path = Path(args.output_dir) / "failed_triplets_log.txt"
            # with open(failure_log_path, 'w') as f_log:
            #    for fail in all_failures:
            #        f_log.write(f"{fail.get('input_triplet')}\t{fail.get('message')}\n")
            # print(f"Detailed failure log saved to: {failure_log_path}", file=sys.stderr)

        # 3. 处理并保存成功的结果
        if all_stats:
            stats_df = pd.DataFrame(all_stats)
            column_order = [
                "Gene_ID", "TAD_ID", "SNP_ID", "Full_SNP_ID",
                "MannWhitney_P_TAD", "TAD_IS_Significant",
                "MannWhitney_P_Gene", "Gene_Exp_Significant",
                "Both_Significant",
                "MannWhitney_P_FL", "MannWhitney_P_FS",
                "Independence_P_TAD01_Genotype", "Independence_Test_Method", "Independence_Test_Reason"
            ]
            stats_df = stats_df.reindex(columns=column_order)

            outdir = Path(args.output_dir)
            outdir.mkdir(parents=True, exist_ok=True)

            csv_all = outdir / "all_triplet_statistics.csv"
            stats_df.to_csv(csv_all, sep=",", index=False)
            print(f"All results saved to: {csv_all}", file=sys.stderr)

            sig_df = stats_df[stats_df["Both_Significant"]]
            if not sig_df.empty:
                csv_sig = outdir / "significant_triplets.csv"
                sig_df.to_csv(csv_sig, sep=",", index=False)
                print(f"Significant results saved to: {csv_sig}", file=sys.stderr)
            else:
                print("No triplets met the significance criteria for 'Both_Significant'. No filtered file generated.", file=sys.stderr)

            # 打印成功结果的摘要
            print("\nAnalysis Summary (based on successfully processed triplets):", file=sys.stderr)
            print(f"  Total triplets processed successfully: {len(stats_df)}", file=sys.stderr)

            n_sig_tad = int((stats_df["MannWhitney_P_TAD"] < args.tad_is_p_threshold).sum())
            n_valid_tad = int(stats_df["MannWhitney_P_TAD"].notna().sum())
            pct_tad = (n_sig_tad / n_valid_tad * 100) if n_valid_tad else 0.0
            print(f"  MannWhitney_P_TAD @{args.tad_is_p_threshold}: {n_sig_tad}/{n_valid_tad} ({pct_tad:.1f}%)", file=sys.stderr)

            n_sig_gene = int((stats_df["MannWhitney_P_Gene"] < args.gene_exp_p_threshold).sum())
            n_valid_gene = int(stats_df["MannWhitney_P_Gene"].notna().sum())
            pct_gene = (n_sig_gene / n_valid_gene * 100) if n_valid_gene else 0.0
            print(f"  MannWhitney_P_Gene @{args.gene_exp_p_threshold}: {n_sig_gene}/{n_valid_gene} ({pct_gene:.1f}%)", file=sys.stderr)

            for col in ["MannWhitney_P_FL", "MannWhitney_P_FS", "Independence_P_TAD01_Genotype"]:
                n_sig = int((stats_df[col] < 0.05).sum())
                n_valid = int(stats_df[col].notna().sum())
                pct = (n_sig / n_valid * 100) if n_valid else 0.0
                print(f"  {col} @0.05: {n_sig}/{n_valid} ({pct:.1f}%)", file=sys.stderr)

            print(f"\n  TAD IS significant (p < {args.tad_is_p_threshold}): {int(stats_df['TAD_IS_Significant'].sum())}", file=sys.stderr)
            print(f"  Gene Expression significant (p < {args.gene_exp_p_threshold}): {int(stats_df['Gene_Exp_Significant'].sum())}", file=sys.stderr)
            print(f"  Both significant (saved to filtered file): {int(sig_df.shape[0])}", file=sys.stderr)

            method_counts = stats_df["Independence_Test_Method"].value_counts()
            print("Independence test methods used:", file=sys.stderr)
            for method, count in method_counts.items():
                pct = count / len(stats_df) * 100
                print(f"  {method}: {count}/{len(stats_df)} ({pct:.1f}%)", file=sys.stderr)
        else:
            print("No valid results to save (all processed triplets failed).", file=sys.stderr)

        print("\nAnalysis completed successfully", file=sys.stderr)

    except Exception as e:
        print(f"\nCRITICAL ERROR: Analysis failed unexpectedly: {str(e)}", file=sys.stderr)
        raise

if __name__ == "__main__":
    main()