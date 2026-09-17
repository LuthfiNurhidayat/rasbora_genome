#!/bin/bash
conda init bash
conda activate manta

# Direction 1: Male mapped to Female Reference
configManta.py \
  --bam male_illumina_to_female.bam \
  --reference /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_fem_polished.fasta \
  --runDir variant_results/manta_male_on_fem_out

python2 variant_results/manta_male_on_fem_out/runWorkflow.py -j 16

# Direction 2: Female mapped to Male Reference
configManta.py \
  --bam female_illumina_to_male.bam \
  --reference /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta \
  --runDir variant_results/manta_fem_on_male_out

python2 variant_results/manta_fem_on_male_out/runWorkflow.py -j 16

# Direction 3: Female mapped to Female Reference
configManta.py \
  --bam female_illumina_to_female.bam \
  --reference /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_fem_polished.fasta \
  --runDir variant_results/manta_fem_on_fem_out

python2 variant_results/manta_fem_on_fem_out/runWorkflow.py -j 16

# Direction 4: male mapped to Male Reference
configManta.py \
  --bam male_illumina_to_male.bam \
  --reference /data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/male_fem_map/Raslat_male_polished.fasta \
  --runDir variant_results/manta_male_on_male_out

python2 variant_results/manta_male_on_male_out/runWorkflow.py -j 16

conda deactivate


