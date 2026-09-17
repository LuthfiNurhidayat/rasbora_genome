#!/bin/bash

##activate conda environment
conda init bash
conda activate asm

export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/bin:$PATH
export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/envs/genann/bin:$PATH

cd /home/luthfibio/raslat_fem/
mkdir raslat_fem_asm
flye --nano-raw "/data/0.RAW/WGS/2026-30/20260203_1831_3F_PAW84497_644f1157/fastq_pass/2026-30.fastq.gz" --threads 64 --out-dir raslat_fem_asm

flye --nano-hq "/data/0.RAW/WGS/2026-30/20260203_1831_3F_PAW84497_644f1157/fastq_pass/2026-30.fastq.gz" --threads 64 --out-dir raslat_fem_asm -g 1.2G --asm-coverage 40

cd  /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/raslat_fem_asm/
assembly-stats assembly.fasta > raslat_fem_flye_stats1.txt
assembly_stats assembly.fasta > raslat_fem_flye_stats2.txt
busco -i assembly.fasta -m genome -l actinopterygii_odb10 -o BUSCO_Raslat_contig -c 32

conda activate asm
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/raslat_fem_asm/
mkdir adapter_outputdir
"/data/6.PANEL_EXPERTS/luthfibio/apps/run_fcsadaptor.sh" --fasta-input assembly.fasta --output-dir ./adapter_outputdir --euk --container-engine singularity --image "/data/6.PANEL_EXPERTS/luthfibio/apps/fcs-adaptor.sif"

cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/adapter_outputdir/cleaned_sequences/
cat "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/assembly.fasta" | python3 "/data/6.PANEL_EXPERTS/luthfibio/apps/fcs.py" --no-report-analytics --image "/data/6.PANEL_EXPERTS/luthfibio/apps/fcs-gx.sif" clean genome --action-report /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/adapter_outputdir/fcs_adaptor_report.txt --output raslat_male_clean.fasta --contam-fasta-out contam.fasta 

##database FCS-GX: /home/prom/Datadrive/database/data/gxdb/
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm
mkdir gx_out
python3 "/data/6.PANEL_EXPERTS/luthfibio/apps/fcs.py" screen genome --fasta "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/adapter_outputdir/cleaned_sequences/raslat_male_clean.fasta" --out-dir gx_out --gx-db "/home/prom/Datadrive/database/data/gxdb" --tax-id 590941 

cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/gx_out/raslat_fem_clean.590941.fcs_gx_report.txt" .
cat "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/adapter_outputdir/cleaned_sequences/raslat_male_clean.fasta" | python3 "/data/6.PANEL_EXPERTS/luthfibio/apps/fcs.py" --no-report-analytics --image "/data/6.PANEL_EXPERTS/luthfibio/apps/fcs-gx.sif" clean genome --action-report "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/raslat_male_clean.590941.fcs_gx_report.txt" --output raslat_male_ctg_clean.fasta --contam-fasta-out contam.fasta 


###purge dup contigs
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
mkdir purge_ctg
cd purge_ctg
conda activate asm
minimap2 -x map-ont -t 20 "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/raslat_male_ctg_clean.fasta" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/ONT_reads_male.fastq" > aln_male_ctg.paf


# 1. Generate coverage stats
pbcstat aln_male_ctg.paf

# 2. Calculate cutoffs
#calcuts PB.stat > cutoffs
echo "20 75 200" > cutoffs

# 3. Split assembly
split_fa "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/raslat_male_ctg_clean.fasta" > raslat_male_ctg_clean.split.fasta

# 4. Self-alignment
minimap2 -x asm5 -DP raslat_male_ctg_clean.split.fasta raslat_male_ctg_clean.split.fasta > asm_ctg.paf

# 5. Identify duplicates
purge_dups -2 -T cutoffs -c PB.base.cov asm_ctg.paf > dups_ctg.bed

# 6. Purge
get_seqs -p raslat_male_ctg_clean dups_ctg.bed /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/raslat_male_ctg_clean.fasta

##RUN BUSCO
conda activate asm
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/purge_ctg/raslat_male_ctg_clean.purged.fa" .
busco -i raslat_male_ctg_clean.purged.fa -m genome -l actinopterygii_odb10 -o BUSCO_Raslat_purge_ctg -c 32
assembly-stats raslat_male_ctg_clean.purged.fa > raslat_male_ctg_purge1.txt

busco -i raslat_male_ctg_clean.fasta -m genome -l actinopterygii_odb10 -o BUSCO_Raslat_ctg -c 32
assembly-stats raslat_male_ctg_clean.fasta > raslat_male_ctg1.txt
assembly_stats raslat_male_ctg_clean.fasta > raslat_male_ctg2.txt

busco -i raslat_male_ctg_clean.hap.fa -m genome -l actinopterygii_odb10 -o BUSCO_Raslat_hap_ctg -c 32

##align raw read to the purged assembly
cd purge_ctg 
/data/6.PANEL_EXPERTS/luthfibio/apps/funannotate-docker sort -i raslat_male_ctg_clean.purged.fa -b contig -o raslat_male_ctg_clean.purged.srt.fa --minlen 200

minimap2 -ax map-ont raslat_male_ctg_clean.purged.srt.fa "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/ONT_reads_male.fastq" | samtools sort -o male_purged.bam
samtools index male_purged.bam
samtools depth male_purged.bam > male_purged.depth
samtools coverage male_purged.bam > male_purged.coverage.tsv

##polishing using shortread
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/
prefetch SRR39758492 --max-size 30G
fasterq-dump SRR39758492 --split-files -e 8 -p
mv SRR39758492_1.fastq Fish_Male_R1.fastq
mv SRR39758492_2.fastq Fish_Male_R2.fastq
#adapter trimming
trim_galore --cores 8 --paired Fish_Male_R1.fastq Fish_Male_R2.fastq
mkdir shortreads_trimmed
cd shortreads_trimmed
mv /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/*.fq .
mv /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/*.json .
mv /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/*.txt .

#QC
fastqc *.fq

##polishing
bash polishing_QC.sh

### copy the polished genome to the female genome folder
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/polishing_results/male/male.polished.fasta" Raslat_male_polished.fasta
##check busco 
busco -i Raslat_male_polished.fasta -m genome -l actinopterygii_odb10 -o BUSCO_Raslat_plsh -c 32
assembly-stats Raslat_male_polished.fasta > raslat_male_plsh.txt

###Scaffolding male genome to Danio
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
ragtag.py scaffold GCF_000002035.6_GRCz11_genomic.fna Raslat_male_polished.fasta

##Clean sort and rename scaffold using funannotate
##Clean sort and rename scaffold
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/ragtag_output/
/data/6.PANEL_EXPERTS/luthfibio/apps/funannotate-docker sort -i ragtag.scaffold.fasta -o raslat_male_scaff_srt.fasta --minlen 200
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
mv /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/ragtag_output/raslat_male_scaff_srt.fasta .

##calculating assembly stat (scaffolds)0

assembly-stats raslat_male_scaff_srt.fasta > raslat_male_scaffold_sort_stat.txt
assembly_stats raslat_male_scaff_srt.fasta > raslat_male_scaffold_sort_stat2.txt

##RUN BUSCO

busco -i raslat_male_scaff_srt.fasta -m genome -l actinopterygii_odb10 -o BUSCO_Raslat_scaff -c 32

##create visualization for genome stat n busco
conda activate btk
mkdir Rasbora_lateristriata_male
blobtools create --fasta raslat_male_scaff_srt.fasta --busco "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/BUSCO_Raslat_scaff/run_actinopterygii_odb10/full_table.tsv" Rasbora_lateristriata_male
blobtools view --plot --view snail Rasbora_lateristriata_male


##copy to BIO6000 computer
scp luthfibio@10.4.100.103:"/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/raslat_male_scaff_srt.fasta" luthfibio@10.4.100.106:/datadrive/drive_a/luthfibio/RasLat/male_annot/raslat_male_scaff_srt.fasta
scp -r luthfibio@10.4.100.103:/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/BUSCO_Raslat_scaff luthfibio@10.4.100.106:/datadrive/drive_a/luthfibio/RasLat/male_annot/

scp -r luthfibio@10.4.100.103:/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/Rasbora_lateristriata_male luthfibio@10.4.100.106:/datadrive/drive_a/luthfibio/RasLat/male_annot/

scp -r luthfibio@10.4.100.103:/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/raslat_fem_asm/ragtag_output luthfibio@10.4.100.106:/datadrive/drive_a/luthfibio/RasLat/Raslat_fem/

scp luthfibio@10.4.100.103:/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/raslat_fem_asm/assembly.fasta luthfibio@10.4.100.106:/datadrive/drive_a/luthfibio/RasLat/Raslat_fem/assembly.fasta


###################################################

##estimates genome size
conda activate genann
mkdir gen_size
cd gen_size
gunzip -c /data/0.RAW/WGS/2026-30/20260203_1831_3F_PAW84497_644f1157/fastq_pass/2026-30.fastq.gz > 2026-30.fastq
jellyfish count -t 16 -C -m 21 -s 5G -o RasLat_21mer_out 2026-30.fastq
jellyfish histo -o RasLat_21mer_out.histo RasLat_21mer_out
mkdir RasLat_fem_21mer
Rscript /data/6.PANEL_EXPERTS/luthfibio/apps/genomescope2.0-2.0.1/genomescope.R -i RasLat_21mer_out.histo -o RasLat_fem_21mer -k 21
rm 2026-30.fastq

