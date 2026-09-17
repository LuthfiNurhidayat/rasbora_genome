#!/bin/bash

##activate conda environment
conda init bash
conda init
conda activate asm

export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/bin:$PATH
export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/envs/asm/bin:$PATH
export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/envs/btk/bin:$PATH

conda activate asm
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/
mkdir mcscanx
cd mcscanx
##download protein and annotation database
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/903/798/145/GCF_903798145.1_fDanAes4.1/GCF_903798145.1_fDanAes4.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/903/798/145/GCF_903798145.1_fDanAes4.1/GCF_903798145.1_fDanAes4.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/340/385/GCF_018340385.1_ASM1834038v1/GCF_018340385.1_ASM1834038v1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/340/385/GCF_018340385.1_ASM1834038v1/GCF_018340385.1_ASM1834038v1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_protein.faa.gz

gunzip *.gz
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/Rlat_male_Annotate/annotate_results/Rasbora_lateristriata_male.proteins.fa" .
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/Rlat_male_Annotate/annotate_results/Rasbora_lateristriata_male.gff3" .
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_fem_asm/Rlat_fem_Annotate/annotate_results/Rasbora_lateristriata_fem.gff3" .
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_fem_asm/Rlat_fem_Annotate/annotate_results/Rasbora_lateristriata_fem.proteins.fa" .

mkdir data

##prepare mcscan data
python "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/prepare_mcscanx_unified.py"  --type refseq --gff GCF_903798145.1_fDanAes4.1_genomic.gff --fasta GCF_903798145.1_fDanAes4.1_protein.faa --outdir data --prefix Daes

python "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/prepare_mcscanx_unified.py"  --type refseq --gff GCF_000002035.6_GRCz11_genomic.gff --fasta GCF_000002035.6_GRCz11_protein.faa --outdir data --prefix Drer

python "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/prepare_mcscanx_unified.py"  --type refseq --gff GCF_018340385.1_ASM1834038v1_genomic.gff --fasta GCF_018340385.1_ASM1834038v1_protein.faa --outdir data --prefix Ccar

python "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/prepare_mcscanx_unified.py"  --type custom --gff Rasbora_lateristriata_male.gff3 --fasta Rasbora_lateristriata_male.proteins.fa --outdir data --prefix Rlatm

python "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/prepare_mcscanx_unified.py"  --type custom --gff Rasbora_lateristriata_fem.gff3 --fasta Rasbora_lateristriata_fem.proteins.fa --outdir data --prefix Rlatf

cd data
makeblastdb -in Drer.fa -dbtype prot
makeblastdb -in Daes.fa -dbtype prot
makeblastdb -in Ccar.fa -dbtype prot
makeblastdb -in Rlatm.fa -dbtype prot
makeblastdb -in Rlatf.fa -dbtype prot

cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/
mkdir duplication


mkdir blast_all
bash blast_all.sh

cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/data/

python "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/merge_all_bed.py" --beds Rlatm.bed Rlatf.bed Drer.bed Daes.bed Ccar.bed --prefixes Rlm Rlf Dr Da Cc --output /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/cyprin/cyprin2.gff --sort

cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/mcscanx/
"/data/6.PANEL_EXPERTS/luthfibio/apps/MCScanX-1.0.0/MCScanX" -s 5 -e 1e-10 -m 20 -b 0  cyprin/cyprin ##unsorted bed merge
"/data/6.PANEL_EXPERTS/luthfibio/apps/MCScanX-1.0.0/MCScanX" -s 5 -e 1e-10 -m 20 -b 0  cyprin/cyprin2






