#!/usr/bin/env python3

import gzip
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

def find_vcf_file(search_pattern):
    """Finds a VCF file matching a wildcard pattern, returning the first match."""
    matches = glob.glob(search_pattern)
    if matches:
        print(f"[FOUND] {matches[0]}")
        return matches[0]
    else:
        print(f"[NOT FOUND] Could not find file matching: {search_pattern}")
        return None

def parse_vcf(vcf_path, caller_name, direction_label, is_control):
    """Parses standard VCF files into a Pandas DataFrame."""
    if not vcf_path or not os.path.exists(vcf_path):
        return pd.DataFrame()

    records = []
    open_fn = gzip.open if vcf_path.endswith('.gz') else open

    with open_fn(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, sv_id, ref, alt, qual, filt, info = parts[:8]

            sv_type, sv_len = "UNKNOWN", 0
            for item in info.split(';'):
                if item.startswith('SVTYPE='):
                    sv_type = item.split('=')[1]
                elif item.startswith('SVLEN='):
                    try:
                        sv_len = abs(int(item.split('=')[1]))
                    except ValueError:
                        sv_len = 0

            records.append({
                'chrom': chrom,
                'pos': int(pos),
                'id': sv_id,
                'type': sv_type,
                'length': sv_len,
                'caller': caller_name,
                'direction': direction_label,
                'is_control': is_control
            })

    return pd.DataFrame(records)

print("=== Locating VCF Files Across Working Directory ===")

# Define flexible glob search patterns for all 8 VCFs
vcf_targets = [
    # Sniffles2 Files (Long Reads)
    ("*male*on*fem*.vcf*", "Sniffles2", "Male on Female Ref", False),
    ("*fem*on*male*.vcf*", "Sniffles2", "Female on Male Ref", False),
    ("*male*base*.vcf*", "Sniffles2", "Male Baseline Control", True),
    ("*fem*base*.vcf*", "Sniffles2", "Female Baseline Control", True),
    
    # Manta Files (Short Reads - searching subfolders for diploidSV.vcf.gz)
    ("*manta*male*on*fem*/results/variants/diploidSV.vcf.gz", "Manta", "Male on Female Ref", False),
    ("*manta*fem*on*male*/results/variants/diploidSV.vcf.gz", "Manta", "Female on Male Ref", False),
    ("*manta*male*on*male*/results/variants/diploidSV.vcf.gz", "Manta", "Male Baseline Control", True),
    ("*manta*fem*on*fem*/results/variants/diploidSV.vcf.gz", "Manta", "Female Baseline Control", True),
]

parsed_dfs = []

for pattern, caller, direction, is_control in vcf_targets:
    vcf_path = find_vcf_file(pattern)
    if vcf_path:
        df = parse_vcf(vcf_path, caller, direction, is_control)
        if not df.empty:
            parsed_dfs.append(df)

if not parsed_dfs:
    print("\nError: No valid VCF files could be parsed. Check your directory structure.")
    exit(1)

df_all = pd.concat(parsed_dfs, ignore_index=True)
df_all.to_csv("results/all_parsed_svs.tsv", sep='\t', index=False)
print(f"\nSuccessfully loaded {len(df_all)} total structural variants.")

# ------------------------------------------------------------------------------
# FIGURE 1: 2x2 Matrix Call Summary
# ------------------------------------------------------------------------------
plt.figure(figsize=(12, 6))
sns.set_theme(style="whitegrid")

ax = sns.barplot(
    data=df_all,
    x="direction",
    y="length",
    hue="caller",
    estimator=len,
    errorbar=None,
    palette="Dark2"
)

plt.title("Structural Variant Calls Across Reciprocal 2x2 Mapping Matrix", fontsize=12, fontweight="bold")
plt.xlabel("Mapping Direction & Reference Context", fontsize=10)
plt.ylabel("Variant Count", fontsize=10)
plt.xticks(rotation=15)
plt.legend(title="Caller / Platform")

for p in ax.patches:
    height = p.get_height()
    if not np.isnan(height) and height > 0:
        ax.annotate(f'{int(height)}',
                    (p.get_x() + p.get_width() / 2., height),
                    ha='center', va='bottom',
                    fontsize=9, xytext=(0, 3), textcoords='offset points')

plt.tight_layout()
plt.savefig("figures/fig1_2x2_matrix_sv_counts.png", dpi=300)
plt.close()

# ------------------------------------------------------------------------------
# FIGURE 2: Cross-Mapping Structural Variant Types
# ------------------------------------------------------------------------------
df_cross = df_all[~df_all['is_control']].copy()

if not df_cross.empty:
    plt.figure(figsize=(10, 5))
    ax2 = sns.countplot(
        data=df_cross,
        x="type",
        hue="direction",
        palette="Set1"
    )

    plt.title("Reciprocal Cross-Mapping Structural Variant Breakdown", fontsize=12, fontweight="bold")
    plt.xlabel("Structural Variant Class", fontsize=10)
    plt.ylabel("Variant Count", fontsize=10)
    plt.legend(title="Direction")

    for p in ax2.patches:
        height = p.get_height()
        if not np.isnan(height) and height > 0:
            ax2.annotate(f'{int(height)}',
                        (p.get_x() + p.get_width() / 2., height),
                        ha='center', va='bottom',
                        fontsize=9, xytext=(0, 3), textcoords='offset points')

    plt.tight_layout()
    plt.savefig("figures/fig2_reciprocal_sv_types.png", dpi=300)
    plt.close()

print("\nProcessing complete! Generated figures in 'figures/' directory.")