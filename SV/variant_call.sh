#!/bin/bash
conda init bash
conda activate sniffles_env

#samtools index raw_male_to_male.bam
#samtools index raw_fem_to_fem.bam
#samtools index raw_male_to_fem.bam
#samtools index raw_fem_to_male.bam

# 1. Male Baseline
sniffles --input raw_male_to_male.bam --vcf variant_results/sv_male_baseline.vcf --threads 16

# 2. Female Baseline
sniffles --input raw_fem_to_fem.bam --vcf variant_results/sv_female_baseline.vcf --threads 16

# 3. Male Divergence relative to Female
sniffles --input raw_male_to_fem.bam --vcf variant_results/sv_male_on_fem.vcf --threads 16

# 4. Female Divergence relative to Male
sniffles --input raw_fem_to_male.bam --vcf variant_results/sv_fem_on_male.vcf --threads 16


conda deactivate



