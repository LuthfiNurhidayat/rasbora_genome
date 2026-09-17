#!/bin/bash

##activate conda environment
conda init bash
conda init
conda activate genann

export PATH=/home/luthfibio/miniforge3/bin:$PATH
export PATH=/home/luthfibio/miniforge3/envs/genann/bin:$PATH
export PATH=/home/luthfibio/miniforge3/envs/buscoenv/bin:$PATH


##copy genome assembly from repeat masker de novo
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/repeat_annot/raslat_repmask_denovo/raslat_male_scaff_srt.fasta.masked" raslat_male_scaff_srt_mask.fasta

##Structural annotation using galba (Braker-derived script for large genome size)
mkdir struct_annot
conda activate asm


## download and concatenate protein from closely related species 
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/cyprinoidei.faa" .
cp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_anno/galba.sif" .
singularity exec galba.sif galba.pl --genome=raslat_male_scaff_srt_mask.fasta --prot_seq=cyprinoidei.faa --workingdir=struct_annot --threads 16 --gff3 --AUGUSTUS_CONFIG_PATH=/home/luthfibio/miniforge3/envs/asm/config/

#check galba busco since it did not use compleasm
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/struct_annot
conda activate asm
busco -i galba.aa -m proteins -l actinopterygii_odb10 -o BUSCO_Raslatmale_prot -c 32

cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
rm galba.sif
#####################################################################################

##Functional Annotation
#run interprosca
conda activate asm
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
mkdir funct_annot

##remove * in the protein sequence

sed -i "s/\*//g" "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/struct_annot/galba.aa"

/data/6.PANEL_EXPERTS/luthfibio/apps/interproscan-5.61-93.0/interproscan.sh --input /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/struct_annot/galba.aa --formats xml --cpu 16 --output-dir funct_annot -dp

#"/data/6.PANEL_EXPERTS/luthfibio/apps/funannotate-docker" iprscan --input /data/6.PANEL_EXPERTS/luthfibio/Results/Ular_fajar/Nsput/struct_annot/galba.aa --method local --out /data/6.PANEL_EXPERTS/luthfibio/Results/Ular_fajar/Nsput/funct_annot/Nsput_iprscan.xml --iprscan_path /data/6.PANEL_EXPERTS/luthfibio/apps/interproscan-5.61-93.0/interproscan.sh --cpus 16

#run signalp

cd funct_annot
signalp -fasta /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/struct_annot/galba.aa -org euk -format short -prefix Rlat_male


### run funannotate annotate 
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/
mkdir Rlat_male_Annotate
"/data/6.PANEL_EXPERTS/luthfibio/apps/funannotate-docker" annotate --gff "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/struct_annot/galba.gff3" --fasta "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/raslat_male_scaff_srt.fasta" --species "Rasbora lateristriata male" --out Rlat_male_Annotate --cpus 16 --iprscan "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/funct_annot/galba.aa.xml" --signalp "/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/new_male_asm/funct_annot/Rlat_male_summary.signalp5" --busco vertebrata --force 

##run busco for final annotation results
cd /data/6.PANEL_EXPERTS/luthfibio/Results/Ular_fajar/Nsput/Nsput_Annotate/annotate_results/
conda init
conda activate busco
busco -i Naja_sputatrix.proteins.fa -m proteins -l sauropsida_odb10 -o BUSCO_Nsput_annot -c 32