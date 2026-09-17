#!/usr/bin/env bash
set -e

mkdir -p fixed_results

echo "=== STEP 1: Intersecting Sniffles2 & Manta via bedtools ==="

# Helper function to intersect Sniffles2 and Manta VCFs directly
intersect_tools() {
    local sniff_vcf="$1"
    local manta_vcf="$2"
    local out_vcf="$3"

    # Keeps Sniffles2 SVs that have at least 50% reciprocal overlap with Manta calls
    bedtools intersect -a "$sniff_vcf" -b "$manta_vcf" -f 0.5 -r -header > "$out_vcf"
    
    local count=$(grep -v "^#" "$out_vcf" | wc -l)
    echo "[INTERSECT] Generated $out_vcf ($count SVs)"
}

# Locate raw input files dynamically
SNIFF_MxF=$(find . -maxdepth 4 -path "*sv_male_on_fem.vcf*" | head -n 1)
SNIFF_FxM=$(find . -maxdepth 4 -path "*sv_fem_on_male.vcf*" | head -n 1)
SNIFF_MxM=$(find . -maxdepth 4 -path "*sv_male_baseline.vcf*" | head -n 1)
SNIFF_FxF=$(find . -maxdepth 4 -path "*sv_female_baseline.vcf*" | head -n 1)

MANTA_MxF=$(find . -maxdepth 4 -path "*manta*male*on*fem*/**/diploidSV.vcf.gz" | head -n 1)
MANTA_FxM=$(find . -maxdepth 4 -path "*manta*fem*on*male*/**/diploidSV.vcf.gz" | head -n 1)
MANTA_MxM=$(find . -maxdepth 4 -path "*manta*male*on*male*/**/diploidSV.vcf.gz" | head -n 1)
MANTA_FxF=$(find . -maxdepth 4 -path "*manta*fem*on*fem*/**/diploidSV.vcf.gz" | head -n 1)

# Intersect platform calls
intersect_tools "$SNIFF_MxF" "$MANTA_MxF" "fixed_results/intersect_male_on_female.vcf"
intersect_tools "$SNIFF_FxM" "$MANTA_FxM" "fixed_results/intersect_female_on_male.vcf"
intersect_tools "$SNIFF_MxM" "$MANTA_MxM" "fixed_results/intersect_male_baseline.vcf"
intersect_tools "$SNIFF_FxF" "$MANTA_FxF" "fixed_results/intersect_female_baseline.vcf"

echo ""
echo "=== STEP 2: Performing Background Subtractions ==="

# Male-Specific (Male-on-Female MINUS Female Baseline Control)
bedtools subtract \
  -A \
  -a fixed_results/intersect_male_on_female.vcf \
  -b fixed_results/intersect_female_baseline.vcf \
  -f 0.5 -r \
  > fixed_results/MALE_SPECIFIC_HIGH_CONFIDENCE.vcf

M_COUNT=$(grep -v "^#" fixed_results/MALE_SPECIFIC_HIGH_CONFIDENCE.vcf | wc -l)
echo "[FINAL] Male-Specific High-Confidence SVs: $M_COUNT"

# Female-Specific (Female-on-Male MINUS Male Baseline Control)
bedtools subtract \
  -A \
  -a fixed_results/intersect_female_on_male.vcf \
  -b fixed_results/intersect_male_baseline.vcf \
  -f 0.5 -r \
  > fixed_results/FEMALE_SPECIFIC_HIGH_CONFIDENCE.vcf

F_COUNT=$(grep -v "^#" fixed_results/FEMALE_SPECIFIC_HIGH_CONFIDENCE.vcf | wc -l)
echo "[FINAL] Female-Specific High-Confidence SVs: $F_COUNT"