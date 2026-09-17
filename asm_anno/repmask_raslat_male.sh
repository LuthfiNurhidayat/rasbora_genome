#!/bin/bash

##activate conda environment
conda init bash
conda init
conda activate asm

export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/bin:$PATH
export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/envs/asm/bin:$PATH
export PATH=/data/6.PANEL_EXPERTS/luthfibio/apps/miniforge3/envs/btk/bin:$PATH


cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
mkdir repeat_annot
cd repeat_annot

####REPEAT MODELLER
mkdir raslat_repmodeller
cd raslat_repmodeller
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/raslat_male_scaff_srt.fasta" .
# Create RepeatModeler database
singularity exec /data/6.PANEL_EXPERTS/luthfibio/apps/dfam-tetools-latest.sif BuildDatabase -name raslat_male raslat_male_scaff_srt.fasta
## Run RepeatModeler
singularity exec /data/6.PANEL_EXPERTS/luthfibio/apps/dfam-tetools-latest.sif RepeatModeler -threads 32 -database raslat_male
rm raslat_male_scaff_srt.fasta

###### Repeat masker using repeat modeller results
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/repeat_annot/
mkdir raslat_repmask
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/raslat_male_scaff_srt.fasta" .

singularity exec -B /data/6.PANEL_EXPERTS/luthfibio/apps/Libraries:/opt/RepeatMasker/Libraries /data/6.PANEL_EXPERTS/luthfibio/apps/dfam-tetools-latest.sif RepeatMasker -pa 16 -species "Actinopterygii" -a -dir raslat_repmask -xsmall raslat_male_scaff_srt.fasta

mkdir raslat_repmask_denovo
singularity exec -B /data/6.PANEL_EXPERTS/luthfibio/apps/Libraries:/opt/RepeatMasker/Libraries /data/6.PANEL_EXPERTS/luthfibio/apps/dfam-tetools-latest.sif RepeatMasker -pa 16 -a -dir raslat_repmask_denovo -xsmall -lib /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/repeat_annot/raslat_repmodeller/raslat_male-families.fa raslat_male_scaff_srt.fasta
