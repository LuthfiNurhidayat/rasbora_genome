#!/bin/bash

bwa-mem2 index "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta"
bwa-mem2 index "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_fem_asm/Raslat_fem_polished.fasta"

bwa-mem2 mem -t 16 "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Male_R1_val_1.fq" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Male_R2_val_2.fq" | samtools sort -@ 8 -o male_illumina_to_male.bam

samtools index male_illumina_to_male.bam

bwa-mem2 mem -t 16 "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_fem_asm/Raslat_fem_polished.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Female_R1_val_1.fq" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Female_R2_val_2.fq" | samtools sort -@ 8 -o female_illumina_to_female.bam

samtools index female_illumina_to_female.bam

bwa-mem2 mem -t 16 "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Female_R1_val_1.fq" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Female_R2_val_2.fq" | samtools sort -@ 8 -o female_illumina_to_male.bam

samtools index female_illumina_to_male.bam

bwa-mem2 mem -t 16 "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_fem_asm/Raslat_fem_polished.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Male_R1_val_1.fq" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Male_R2_val_2.fq" | samtools sort -@ 8 -o male_illumina_to_female.bam

samtools index male_illumina_to_female.bam



