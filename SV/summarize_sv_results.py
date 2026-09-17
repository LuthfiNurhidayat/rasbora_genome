#!/usr/bin/env python3

import gzip
import os
import pandas as pd


def parse_vcf_metrics(vcf_path, stage_name):
    """Parses SVTYPE and length metrics from a VCF file."""
    if not os.path.exists(vcf_path):
        print(f"[MISSING] {vcf_path}")
        return []

    records = []
    open_fn = gzip.open if vcf_path.endswith(".gz") else open

    with open_fn(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            info = parts[7]

            sv_type, sv_len = "UNKNOWN", 0
            for item in info.split(";"):
                if item.startswith("SVTYPE="):
                    sv_type = item.split("=")[1]
                elif item.startswith("SVLEN="):
                    try:
                        sv_len = abs(int(item.split("=")[1]))
                    except ValueError:
                        sv_len = 0

            records.append({
                "Stage": stage_name,
                "SVType": sv_type,
                "Length": sv_len,
            })

    return records


# ------------------------------------------------------------------------------
# 1. DEFINE VCF PATHS
# ------------------------------------------------------------------------------
vcf_files = {
    # Step 1: Intersected Platforms (Sniffles2 n Manta)
    "Male on Female (Intersected)": (
        "fixed_results/intersect_male_on_female.vcf"
    ),
    "Female on Male (Intersected)": (
        "fixed_results/intersect_female_on_male.vcf"
    ),
    "Male Baseline (Intersected)": (
        "fixed_results/intersect_male_baseline.vcf"
    ),
    "Female Baseline (Intersected)": (
        "fixed_results/intersect_female_baseline.vcf"
    ),
    # Step 2: Baseline-Subtracted (Sex-Specific High Confidence)
    "Male-Specific (Subtracted)": (
        "fixed_results/MALE_SPECIFIC_HIGH_CONFIDENCE.vcf"
    ),
    "Female-Specific (Subtracted)": (
        "fixed_results/FEMALE_SPECIFIC_HIGH_CONFIDENCE.vcf"
    ),
}

all_records = []
for stage, path in vcf_files.items():
    records = parse_vcf_metrics(path, stage)
    all_records.extend(records)

df = pd.DataFrame(all_records)

if df.empty:
    print(
        "Error: No records parsed. Check file paths in"
        " 'fixed_results/'."
    )
    exit(1)

# ------------------------------------------------------------------------------
# 2. GENERATE METRIC SUMMARY TABLES
# ------------------------------------------------------------------------------

# Table A: Overall Counts & Length Stats per Experimental Condition
summary_stats = (
    df.groupby("Stage")["Length"]
    .agg(
        Total_SVs="count",
        Mean_Size_bp="mean",
        Median_Size_bp="median",
        Min_Size_bp="min",
        Max_Size_bp="max",
    )
    .reset_index()
)

# Table B: SV Breakdown by Class Type (INS, DEL, INV, DUP, BND)
type_counts = (
    df.groupby(["Stage", "SVType"])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)

# Merge overall metrics with SV type breakdown
final_table = pd.merge(summary_stats, type_counts, on="Stage")

# Format output
os.makedirs("fixed_results/summary_tables", exist_ok=True)
tsv_out = "fixed_results/summary_tables/sv_matrix_summary.tsv"
csv_out = "fixed_results/summary_tables/sv_matrix_summary.csv"

final_table.to_csv(tsv_out, sep="\t", index=False)
final_table.to_csv(csv_out, index=False)

print("=== SV MATRIX SUMMARY TABLE ===")
print(final_table.to_string(index=False))
print(f"\n[SAVED] Output files:\n - {tsv_out}\n - {csv_out}")