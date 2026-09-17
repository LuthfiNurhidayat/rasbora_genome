#!/usr/bin/env python3

import gzip
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_vcf_contigs(vcf_path, label):
    """Parses contig positions and SV types from a VCF file."""
    if not os.path.exists(vcf_path):
        print(f"[MISSING] File not found: {vcf_path}")
        return []

    records = []
    open_fn = gzip.open if vcf_path.endswith(".gz") else open

    with open_fn(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            info = parts[7]

            sv_type = "UNKNOWN"
            for item in info.split(";"):
                if item.startswith("SVTYPE="):
                    sv_type = item.split("=")[1]
                    break

            records.append({
                "Contig": chrom,
                "Position": pos,
                "SVType": sv_type,
                "Group": label,
            })

    return records


# ------------------------------------------------------------------------------
# 1. LOAD SUBTRACTED VCF DATA FROM FIXED_RESULTS
# ------------------------------------------------------------------------------
male_vcf = "fixed_results/MALE_SPECIFIC_HIGH_CONFIDENCE.vcf"
female_vcf = "fixed_results/FEMALE_SPECIFIC_HIGH_CONFIDENCE.vcf"

records_male = parse_vcf_contigs(male_vcf, "Male-Specific")
records_female = parse_vcf_contigs(female_vcf, "Female-Specific")

df = pd.DataFrame(records_male + records_female)

if df.empty:
    print(
        "Error: No SV records found in fixed_results/ VCF files. Check file"
        " paths."
    )
    exit(1)

# Sort contigs naturally (e.g., contig_1, contig_2... contig_10)
df["contig_num"] = df["Contig"].str.extract(r"(\d+)").astype(float)
df = df.sort_values("contig_num").drop(columns=["contig_num"])

# ------------------------------------------------------------------------------
# 2. AGGREGATE SV COUNTS PER CONTIG
# ------------------------------------------------------------------------------
counts_df = (
    df.groupby(["Contig", "Group"]).size().unstack(fill_value=0).reset_index()
)

# Filter top 25 most variant-rich contigs for clear visualization
top_contigs = df["Contig"].value_counts().head(25).index.tolist()
df_filtered = counts_df[counts_df["Contig"].isin(top_contigs)]

# ------------------------------------------------------------------------------
# 3. PLOT CANDIDATE SV DENSITY
# ------------------------------------------------------------------------------
os.makedirs("figures", exist_ok=True)
plt.figure(figsize=(14, 7))
sns.set_theme(style="whitegrid")

# Reshape for grouped barplot
df_melted = df_filtered.melt(
    id_vars="Contig",
    value_vars=["Male-Specific", "Female-Specific"],
    var_name="Candidate Class",
    value_name="SV Count",
)

ax = sns.barplot(
    data=df_melted,
    x="Contig",
    y="SV Count",
    hue="Candidate Class",
    palette=["#2b5c8f", "#d95f02"],
)

plt.title(
    "Genome-Wide Distribution of Candidate Sex-Differentiating SVs",
    fontsize=14,
    fontweight="bold",
    pad=15,
)
plt.xlabel("Contig / Chromosome", fontsize=11, labelpad=10)
plt.ylabel("High-Confidence SV Count", fontsize=11)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.legend(title="Candidate Dataset", title_fontsize="10", loc="upper right")

# Annotate counts on top of major bars
for p in ax.patches:
    height = p.get_height()
    if not pd.isna(height) and height > 15:
        ax.annotate(
            f"{int(height)}",
            (p.get_x() + p.get_width() / 2.0, height),
            ha="center",
            va="bottom",
            fontsize=8,
            xytext=(0, 2),
            textcoords="offset points",
        )

plt.tight_layout()
out_png = "figures/fig3_sex_specific_sv_contig_density.png"
plt.savefig(out_png, dpi=300)
plt.close()

print(f"[SUCCESS] Density plot generated and saved to: {out_png}")