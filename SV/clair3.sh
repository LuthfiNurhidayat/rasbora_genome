#!/bin/bash
conda init bash

conda activate variant_env

MODEL_NAME=r1041_e82_400bps_sup_v520 # Adjust to match your ONT model


# Direction 1: Male mapped to Female Reference
singularity exec clair3_v2.0.2.sif /opt/bin/run_clair3.sh \
  --bam_fn=raw_male_to_fem.bam \
  --ref_fn="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_fem_polished.fasta" \
  --threads=16 \
  --platform="ont" \
  --model_path=/opt/models/r1041_e82_400bps_sup_v520 \
  --output=variant_results/clair3_male_on_fem_snps --include_all_ctgs --sample_name=male_on_fem_snps
  

# Direction 2: Female mapped to Male Reference
singularity exec clair3_v2.0.2.sif /opt/bin/run_clair3.sh \
  --bam_fn=raw_fem_to_male.bam \
  --ref_fn="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta" \
  --threads=16 \
  --platform="ont" \
  --model_path=/opt/models/r1041_e82_400bps_sup_v520 \
  --output=variant_results/clair3_fem_on_male_snps --include_all_ctgs --sample_name=fem_on_male_snps
  
# Direction 3: Male mapped to Male Reference
singularity exec clair3_v2.0.2.sif /opt/bin/run_clair3.sh \
  --bam_fn=raw_male_to_male.bam \
  --ref_fn="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta" \
  --threads=16 \
  --platform="ont" \
  --model_path=/opt/models/r1041_e82_400bps_sup_v520 \
  --output=variant_results/clair3_male_on_male_snps --include_all_ctgs --sample_name=male_on_male_snps

# Direction 4: Female mapped to Female Reference
singularity exec clair3_v2.0.2.sif /opt/bin/run_clair3.sh \
  --bam_fn=raw_fem_to_fem.bam \
  --ref_fn="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_fem_polished.fasta" \
  --threads=16 \
  --platform="ont" \
  --model_path=/opt/models/r1041_e82_400bps_sup_v520 \
  --output=variant_results/clair3_fem_on_fem_snps --include_all_ctgs --sample_name=fem_on_fem_snps
  
conda deactivate

