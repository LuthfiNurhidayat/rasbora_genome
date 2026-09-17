#!/bin/bash

##male to female mapping genome assembly contigs
minimap2 -ax asm5 -t 20 Raslat_fem_polished.fasta Raslat_male_polished.fasta | samtools sort -@ 8 -o male_vs_female.bam
samtools index male_vs_female.bam
minimap2 -ax asm5 -t 20 Raslat_male_polished.fasta Raslat_fem_polished.fasta | samtools sort -@ 8 -o female_vs_male.bam
samtools index female_vs_male.bam

samtools flagstat male_vs_female.bam
samtools flagstat female_vs_male.bam

samtools idxstats male_vs_female.bam
samtools idxstats female_vs_male.bam

samtools depth -aa male_vs_female.bam > male_vs_female.depth.txt
samtools depth -aa female_vs_male.bam > female_vs_male.depth.txt


##create paf file 
minimap2 -x asm5 -t 20 Raslat_fem_polished.fasta Raslat_male_polished.fasta > male_to_female.paf
minimap2 -x asm5 -t 20 Raslat_male_polished.fasta Raslat_fem_polished.fasta > female_to_male.paf