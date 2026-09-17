#!/bin/bash

##activate conda environment
conda init bash
conda activate genann






## trimming adapter
cd /datadrive/drive_a/luthfibio/rasbora_mito/
porechop -i ASF004_Raslat_mito.fastq.gz -b ASF004_Raslat_mito_trim.fastq.gz -t 32

mkdir /datadrive/drive_a/luthfibio/rasbora_mito/
flye --nano-hq ASF004_Raslat_mito.fastq.gz --threads 32 --out-dir RasLat_Mito
/datadrive/drive_a/luthfibio/script/Flye/bin/flye --nano-hq ASF004_Raslat_mito.fastq.gz --genome-size 17k --threads 20 --meta --scaffold --iterations 5 --read-error 0.05 --min-overlap 100 --out-dir RasLat_Mito

/datadrive/drive_a/luthfibio/script/Flye/bin/flye --nano-hq ASF004_Raslat_mito.fastq.gz --genome-size 17k --threads 20 --plasmid --scaffold --iterations 5 --min-overlap 100 --out-dir RasLat_Mito1

minimap2 -d danrer.mmi NC_002333.2[1..16596].fa                     # indexing
minimap2 -a danrer.mmi /data/0.RAW/WGS/2026-30/20260203_1831_3F_PAW84497_644f1157/fastq_pass/2026-30.fastq.gz > RasLat2Danrer.sam   # alignment

minimap2 -ax asm10 -k10 NC_002333.2[1..16596].fa ASF004_Raslat_mito.fastq -t 20 > RasLat2Danrer.sam

samtools view -S -b RasLat2Danrer.sam > RasLat2Danrer.bam

samtools view -b -F 4 RasLat2Danrer.bam > mapped_RasLat2Danrer.bam
samtools view -bq 1 mapped_RasLat2Danrer.bam > mapped_RasLat2Danrer_Unique.bam

samtools fastq mapped_RasLat2Danrer_Unique.bam > mapped_RasLat2Danrer_Unique.fastq
rm RasLat2Danrer.sam
rm RasLat2Danrer.bam


/datadrive/drive_a/luthfibio/script/Flye/bin/flye --nano-hq mapped_RasLat2Danrer_Unique.fastq --threads 20 --meta --scaffold --iterations 5 --min-overlap 100 --out-dir RasLat_Mito_1