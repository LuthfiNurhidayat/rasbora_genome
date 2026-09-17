#!/bin/bash

##activate conda environment
conda init bash
conda activate asm

export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/bin:$PATH
export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/envs/asm/bin:$PATH

cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/
mkdir male_fem_map
cd male_fem_map
##male to female mapping genome assembly contigs
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/Raslat_male_polished.fasta" .
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_fem_asm/Raslat_fem_polished.fasta" .
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

##similarity check male to female
awk '($11 >= 10000 && $10/$11 >= 0.90)' male_to_female.paf > male_to_female.filtered.paf
#Calculate the identity distribution:
awk '{print ($10/$11)*100}' male_to_female.filtered.paf > male_to_female.identity.txt
awk '{sum += $10/$11;n++} END {print "Mean identity:", sum/n*100 }' male_to_female.filtered.paf

##similarity check female to male
awk '($11 >= 10000 && $10/$11 >= 0.90)' female_to_male.paf > female_to_male.filtered.paf
##Calculate the identity distribution:
awk '{print ($10/$11)*100}' female_to_male.filtered.paf > female_to_male.identity.txt
awk '{sum += $10/$11;n++} END {print "Mean identity:", sum/n*100 }' female_to_male.filtered.paf

## map male raw reads to male genome assembly
minimap2 -ax map-ont "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/ONT_reads_male.fastq" | samtools sort -o raw_male_to_male.bam

## map male raw reads to female genome assembly
minimap2 -ax map-ont "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_fem_asm/Raslat_fem_polished.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/ONT_reads_male.fastq" | samtools sort -o raw_male_to_fem.bam


##map female raw reads to female genome assembly

minimap2 -ax map-ont "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_fem_asm/Raslat_fem_polished.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/ONT_reads_female.fastq.gz" | samtools sort -o raw_fem_to_fem.bam 

##map female raw reads to male genome assembly

minimap2 -ax map-ont "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/ONT_reads_female.fastq.gz" | samtools sort -o raw_fem_to_male.bam

##map illumina raw reads to genome assembly
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



##map genome assembly to the zebrafish genome assembly
minimap2 -x asm5 -t 16 GCF_000002035.6_GRCz11_genomic.fna Raslat_male_polished.fasta > male_to_danio.paf
minimap2 -x asm5 -t 16 GCF_000002035.6_GRCz11_genomic.fna Raslat_fem_polished.fasta  > female_to_danio.paf

## prepare for the Ngenome sync file
##sort genome assembly

/data/6.PANEL_EXPERTS/luthfibio/apps/funannotate-docker sort -i Raslat_male_polished.fasta -o Raslat_male_polished_srt.fasta --minlen 200 -b contig
/data/6.PANEL_EXPERTS/luthfibio/apps/funannotate-docker sort -i Raslat_fem_polished.fasta -o Raslat_fem_polished_srt.fasta --minlen 200 -b contig

mkdir asm_synt
cd asm_synt
#male_fem
perl  "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/NGenomeSyn/bin/GetTwoGenomeSyn.pl"  -InGenomeA "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished_srt.fasta"  -InGenomeB "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_fem_polished_srt.fasta"  -OutPrefix male_fem   -MappingBin  minimap2    -BinDir    /home/luthfibio/miniforge3/envs/asm/bin/    -MinLenA  1000  -MinLenB 1000

#fem_male
perl  "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/NGenomeSyn/bin/GetTwoGenomeSyn.pl"  -InGenomeA "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_fem_polished_srt.fasta"  -InGenomeB "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished_srt.fasta"  -OutPrefix fem_male   -MappingBin  minimap2    -BinDir    /home/luthfibio/miniforge3/envs/asm/bin/    -MinLenA  1000  -MinLenB 1000

#fem_danio
perl  "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/NGenomeSyn/bin/GetTwoGenomeSyn.pl"  -InGenomeA "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_fem_polished_srt.fasta"  -InGenomeB "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/GCF_000002035.6_GRCz11_genomic.fna"  -OutPrefix fem_danio   -MappingBin  minimap2    -BinDir    /home/luthfibio/miniforge3/envs/asm/bin/    -MinLenA  1000  -MinLenB 1000

#danio_fem
perl  "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/NGenomeSyn/bin/GetTwoGenomeSyn.pl"  -InGenomeA "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/GCF_000002035.6_GRCz11_genomic.fna"  -InGenomeB "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_fem_polished_srt.fasta"  -OutPrefix danio_fem   -MappingBin  minimap2    -BinDir    /home/luthfibio/miniforge3/envs/asm/bin/    -MinLenA  1000  -MinLenB 1000

#male_danio
perl  "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/NGenomeSyn/bin/GetTwoGenomeSyn.pl"  -InGenomeA "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished_srt.fasta"  -InGenomeB "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/GCF_000002035.6_GRCz11_genomic.fna"  -OutPrefix male_danio   -MappingBin  minimap2    -BinDir    /home/luthfibio/miniforge3/envs/asm/bin/    -MinLenA  1000  -MinLenB 1000

#danio_male
perl  "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/NGenomeSyn/bin/GetTwoGenomeSyn.pl"  -InGenomeA "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/GCF_000002035.6_GRCz11_genomic.fna"  -InGenomeB "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished_srt.fasta"  -OutPrefix danio_male   -MappingBin  minimap2    -BinDir    /home/luthfibio/miniforge3/envs/asm/bin/    -MinLenA  1000  -MinLenB 1000



"/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/NGenomeSyn/bin/NGenomeSyn"










#######map reads to genome 
##copy genome assembly from BIO6000 (run  from BIO6000
scp luthfibio@10.4.100.106:/datadrive/drive_a/luthfibio/RasLat/Raslat_fem/Raslat_fem_annotate/annotate_results/Rasbora_lateristriata_fem.scaffolds.fa luthfibio@10.4.100.103:/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/malefemalemap/Rasbora_lateristriata_fem.scaffolds.fa 

##female reads ro female genome
minimap2 -ax map-ont Rasbora_lateristriata_fem.scaffolds.fa "/data/0.RAW/WGS/2026-30/20260203_1831_3F_PAW84497_644f1157/fastq_pass/2026-30.fastq.gz" | samtools sort -o raslat_female.bam
samtools index raslat_female.bam

### male read to male genome
minimap2 -ax map-ont Rasbora_lateristriata_male.scaffolds.fa "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/malefemalemap/SRR38014054.fastq" | samtools sort -o raslat_male-male.bam
samtools index raslat_male-male.bam

#### male reads to female genome
minimap2 -ax map-ont Rasbora_lateristriata_fem.scaffolds.fa "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/malefemalemap/SRR38014054.fastq" | samtools sort -o raslat_male.bam
samtools index raslat_male.bam

python sex_chr_detection.py -m raslat_male.bam -f raslat_female.bam -r Rasbora_lateristriata_fem.scaffolds.fa

#####check assembly

quast.py Rasbora_lateristriata_fem.scaffolds.fa Rasbora_lateristriata_male.scaffolds.fa -r GCF_049306965.1_GRCz12tu_genomic.fna -g GCF_049306965.1_GRCz12tu_genomic.gff --split-scaffolds --labels Rlatfemale,Rlatmale --eukaryote --large --circos -o raslat

### the results of coverage analysis conclude that there is no sex-link chromosome

##try with SNP analysis
bcftools mpileup -Ou -f Rasbora_lateristriata_fem.scaffolds.fa raslat_male.bam raslat_female.bam | bcftools call -mv -Oz -o variants.vcf.gz

bcftools index variants.vcf.gz
###filtering
bcftools filter -e 'QUAL<30 || DP<10' variants.vcf.gz -Oz -o filtered.vcf.gz
bcftools index filtered.vcf.gz


####explore SNP only in DMRT1 and DMRT3A region
bcftools mpileup -r scaffold_2:23300000-23600000 -f Rasbora_lateristriata_fem.scaffolds.fa raslat_male.bam raslat_female.bam -d 500 -q 20 -Q 20 -Ou | bcftools call -m -v -Oz -o dmrt_region.vcf.gz

bcftools index dmrt_region.vcf.gz

###Remove low-confidence variants:
bcftools filter -i 'QUAL>30 && DP>20' dmrt_region.vcf.gz -Oz -o dmrt_region.filtered.vcf.gz
bcftools index dmrt_region.filtered.vcf.gz

####Export genotypes
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' dmrt_region.filtered.vcf.gz > dmrt_genotypes.tsv

python find_dmrt_candidates.py

###Generate coverage profiles
#male
samtools depth -r scaffold_2:23300000-23600000 raslat_male.bam > male.dmrt.depth
#female
samtools depth -r scaffold_2:23300000-23600000 raslat_female.bam > female.dmrt.depth

##Plot coverage and SNP density
python plot_dmrt_region.py




