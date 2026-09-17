#!/bin/bash

#!/usr/bin/env bash

#!/usr/bin/env bash

###############################################################################
# Polypolish + Merqury polishing workflow
#
# Species: Rasbora lateristriata
#
# Workflow:
#
#   PURGED ASSEMBLY
#        |
#        +---- R1 -- BWA-MEM2 -a --> R1.sam ----+
#        |                                       |
#        +---- R2 -- BWA-MEM2 -a --> R2.sam ----+--> Polypolish
#                                                |
#                                                v
#                                      polished assembly
#
# Additional:
#   - BAM generation for alignment QC
#   - samtools statistics
#   - QUAST
#   - Merqury
#
# IMPORTANT:
#   R1 and R2 MUST be aligned separately.
#   Do NOT combine R1 and R2 in one bwa-mem2 command.
###############################################################################

set -euo pipefail

###############################################################################
# USER SETTINGS
###############################################################################

# Working directory
WORKDIR="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano"

# Purged assemblies
FEMALE_ASM="${WORKDIR}/new_fem_asm/purge_ctg/raslat_fem_ctg_clean.purged.srt.fa"
MALE_ASM="${WORKDIR}/new_male_asm/purge_ctg/raslat_male_ctg_clean.purged.srt.fa"

# Illumina reads
FEMALE_R1="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Female_R1_val_1.fq"
FEMALE_R2="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Female_R2_val_2.fq"

MALE_R1="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Male_R1_val_1.fq"
MALE_R2="/data/6.PANEL_EXPERTS/luthfibio/Results/Raslat_nano/Rawreads/shortreads_trimmed/Fish_Male_R2_val_2.fq"

# Number of CPU threads
THREADS=16

# Threads for samtools compression/sorting
SAMTOOLS_THREADS=8

# Output directory
OUTDIR="${WORKDIR}/polishing_results"

###############################################################################
# PROGRAMS
###############################################################################

BWA="bwa-mem2"
SAMTOOLS="samtools"
POLYPOLISH="polypolish"
QUAST="quast.py"
MERQURY="merqury.sh"
MERYL="meryl"

###############################################################################
# CREATE DIRECTORY STRUCTURE
###############################################################################

mkdir -p "${OUTDIR}"

mkdir -p "${OUTDIR}/female"
mkdir -p "${OUTDIR}/female/sam"
mkdir -p "${OUTDIR}/female/bam"
mkdir -p "${OUTDIR}/female/qc"
mkdir -p "${OUTDIR}/female/logs"

mkdir -p "${OUTDIR}/male"
mkdir -p "${OUTDIR}/male/sam"
mkdir -p "${OUTDIR}/male/bam"
mkdir -p "${OUTDIR}/male/qc"
mkdir -p "${OUTDIR}/male/logs"

mkdir -p "${OUTDIR}/merqury"
mkdir -p "${OUTDIR}/quast"

###############################################################################
# LOGGING
###############################################################################

LOG="${OUTDIR}/pipeline.log"

exec > >(tee -a "${LOG}") 2>&1

echo "============================================================"
echo "Rasbora lateristriata genome polishing workflow"
echo "Started: $(date)"
echo "============================================================"

###############################################################################
# CHECK INPUT FILES
###############################################################################

echo ""
echo "Checking input files..."

for FILE in \
    "${FEMALE_ASM}" \
    "${MALE_ASM}" \
    "${FEMALE_R1}" \
    "${FEMALE_R2}" \
    "${MALE_R1}" \
    "${MALE_R2}"
do
    if [[ ! -f "${FILE}" ]]; then
        echo "ERROR: File not found:"
        echo "${FILE}"
        exit 1
    fi
done

echo "All input files found."

###############################################################################
# CHECK PROGRAMS
###############################################################################

echo ""
echo "Checking required programs..."

for PROG in \
    "${BWA}" \
    "${SAMTOOLS}" \
    "${POLYPOLISH}" \
    "${QUAST}" \
    "${MERQURY}" \
    "${MERYL}"
do
    if ! command -v "${PROG}" >/dev/null 2>&1; then
        echo "ERROR: Program not found: ${PROG}"
        exit 1
    fi
done

echo "All required programs found."

###############################################################################
# VERSION INFORMATION
###############################################################################

echo ""
echo "============================================================"
echo "Software versions"
echo "============================================================"

echo "BWA-MEM2:"
${BWA} version || true

echo ""
echo "Samtools:"
${SAMTOOLS} --version | head -1

echo ""
echo "Polypolish:"
${POLYPOLISH} --version

echo ""
echo "Meryl:"
${MERYL} --version

###############################################################################
# INDEX FEMALE ASSEMBLY
###############################################################################

echo ""
echo "============================================================"
echo "Indexing female assembly"
echo "============================================================"

${BWA} index \
    "${FEMALE_ASM}" \
    > "${OUTDIR}/female/logs/bwa_index.log" \
    2>&1

###############################################################################
# FEMALE R1 ALIGNMENT
###############################################################################

echo ""
echo "============================================================"
echo "Female: BWA-MEM2 alignment of R1"
echo "============================================================"

${BWA} mem \
    -t "${THREADS}" \
    -a \
    "${FEMALE_ASM}" \
    "${FEMALE_R1}" \
    > "${OUTDIR}/female/sam/female_R1.sam" \
    2> "${OUTDIR}/female/logs/bwa_R1.log"

###############################################################################
# FEMALE R2 ALIGNMENT
###############################################################################

echo ""
echo "============================================================"
echo "Female: BWA-MEM2 alignment of R2"
echo "============================================================"

${BWA} mem \
    -t "${THREADS}" \
    -a \
    "${FEMALE_ASM}" \
    "${FEMALE_R2}" \
    > "${OUTDIR}/female/sam/female_R2.sam" \
    2> "${OUTDIR}/female/logs/bwa_R2.log"

###############################################################################
# FEMALE BAM FOR QC
###############################################################################

echo ""
echo "Creating female BAM files for QC..."

${SAMTOOLS} view \
    -@ "${SAMTOOLS_THREADS}" \
    -b \
    "${OUTDIR}/female/sam/female_R1.sam" \
    > "${OUTDIR}/female/bam/female_R1.bam"

${SAMTOOLS} view \
    -@ "${SAMTOOLS_THREADS}" \
    -b \
    "${OUTDIR}/female/sam/female_R2.sam" \
    > "${OUTDIR}/female/bam/female_R2.bam"

${SAMTOOLS} sort \
    -@ "${SAMTOOLS_THREADS}" \
    -o "${OUTDIR}/female/bam/female_R1.sorted.bam" \
    "${OUTDIR}/female/bam/female_R1.bam"

${SAMTOOLS} sort \
    -@ "${SAMTOOLS_THREADS}" \
    -o "${OUTDIR}/female/bam/female_R2.sorted.bam" \
    "${OUTDIR}/female/bam/female_R2.bam"

${SAMTOOLS} index \
    "${OUTDIR}/female/bam/female_R1.sorted.bam"

${SAMTOOLS} index \
    "${OUTDIR}/female/bam/female_R2.sorted.bam"

###############################################################################
# FEMALE BAM QC
###############################################################################

echo ""
echo "Female alignment QC..."

${SAMTOOLS} flagstat \
    "${OUTDIR}/female/bam/female_R1.sorted.bam" \
    > "${OUTDIR}/female/qc/female_R1.flagstat.txt"

${SAMTOOLS} flagstat \
    "${OUTDIR}/female/bam/female_R2.sorted.bam" \
    > "${OUTDIR}/female/qc/female_R2.flagstat.txt"

${SAMTOOLS} stats \
    "${OUTDIR}/female/bam/female_R1.sorted.bam" \
    > "${OUTDIR}/female/qc/female_R1.stats.txt"

${SAMTOOLS} stats \
    "${OUTDIR}/female/bam/female_R2.sorted.bam" \
    > "${OUTDIR}/female/qc/female_R2.stats.txt"

###############################################################################
# FEMALE POLYPOLISH
###############################################################################

echo ""
echo "============================================================"
echo "Female: Polypolish"
echo "============================================================"

${POLYPOLISH} polish \
    "${FEMALE_ASM}" \
    "${OUTDIR}/female/sam/female_R1.sam" \
    "${OUTDIR}/female/sam/female_R2.sam" \
    > "${OUTDIR}/female/female.polished.fasta" \
    2> "${OUTDIR}/female/logs/polypolish.log"

echo "Female polishing completed."

###############################################################################
# INDEX MALE ASSEMBLY
###############################################################################

echo ""
echo "============================================================"
echo "Indexing male assembly"
echo "============================================================"

${BWA} index \
    "${MALE_ASM}" \
    > "${OUTDIR}/male/logs/bwa_index.log" \
    2>&1

###############################################################################
# MALE R1 ALIGNMENT
###############################################################################

echo ""
echo "============================================================"
echo "Male: BWA-MEM2 alignment of R1"
echo "============================================================"

${BWA} mem \
    -t "${THREADS}" \
    -a \
    "${MALE_ASM}" \
    "${MALE_R1}" \
    > "${OUTDIR}/male/sam/male_R1.sam" \
    2> "${OUTDIR}/male/logs/bwa_R1.log"

###############################################################################
# MALE R2 ALIGNMENT
###############################################################################

echo ""
echo "============================================================"
echo "Male: BWA-MEM2 alignment of R2"
echo "============================================================"

${BWA} mem \
    -t "${THREADS}" \
    -a \
    "${MALE_ASM}" \
    "${MALE_R2}" \
    > "${OUTDIR}/male/sam/male_R2.sam" \
    2> "${OUTDIR}/male/logs/bwa_R2.log"

###############################################################################
# MALE BAM FOR QC
###############################################################################

echo ""
echo "Creating male BAM files for QC..."

${SAMTOOLS} view \
    -@ "${SAMTOOLS_THREADS}" \
    -b \
    "${OUTDIR}/male/sam/male_R1.sam" \
    > "${OUTDIR}/male/bam/male_R1.bam"

${SAMTOOLS} view \
    -@ "${SAMTOOLS_THREADS}" \
    -b \
    "${OUTDIR}/male/sam/male_R2.sam" \
    > "${OUTDIR}/male/bam/male_R2.bam"

${SAMTOOLS} sort \
    -@ "${SAMTOOLS_THREADS}" \
    -o "${OUTDIR}/male/bam/male_R1.sorted.bam" \
    "${OUTDIR}/male/bam/male_R1.bam"

${SAMTOOLS} sort \
    -@ "${SAMTOOLS_THREADS}" \
    -o "${OUTDIR}/male/bam/male_R2.sorted.bam" \
    "${OUTDIR}/male/bam/male_R2.bam"

${SAMTOOLS} index \
    "${OUTDIR}/male/bam/male_R1.sorted.bam"

${SAMTOOLS} index \
    "${OUTDIR}/male/bam/male_R2.sorted.bam"

###############################################################################
# MALE BAM QC
###############################################################################

echo ""
echo "Male alignment QC..."

${SAMTOOLS} flagstat \
    "${OUTDIR}/male/bam/male_R1.sorted.bam" \
    > "${OUTDIR}/male/qc/male_R1.flagstat.txt"

${SAMTOOLS} flagstat \
    "${OUTDIR}/male/bam/male_R2.sorted.bam" \
    > "${OUTDIR}/male/qc/male_R2.flagstat.txt"

${SAMTOOLS} stats \
    "${OUTDIR}/male/bam/male_R1.sorted.bam" \
    > "${OUTDIR}/male/qc/male_R1.stats.txt"

${SAMTOOLS} stats \
    "${OUTDIR}/male/bam/male_R2.sorted.bam" \
    > "${OUTDIR}/male/qc/male_R2.stats.txt"

###############################################################################
# MALE POLYPOLISH
###############################################################################

echo ""
echo "============================================================"
echo "Male: Polypolish"
echo "============================================================"

${POLYPOLISH} polish \
    "${MALE_ASM}" \
    "${OUTDIR}/male/sam/male_R1.sam" \
    "${OUTDIR}/male/sam/male_R2.sam" \
    > "${OUTDIR}/male/male.polished.fasta" \
    2> "${OUTDIR}/male/logs/polypolish.log"

echo "Male polishing completed."

###############################################################################
# QUAST
###############################################################################

echo ""
echo "============================================================"
echo "Running QUAST"
echo "============================================================"

${QUAST} \
    "${FEMALE_ASM}" \
    "${OUTDIR}/female/female.polished.fasta" \
    "${MALE_ASM}" \
    "${OUTDIR}/male/male.polished.fasta" \
    -o "${OUTDIR}/quast" \
    -t "${THREADS}" \
    --eukaryote \
    > "${OUTDIR}/quast/quast.stdout.log" \
    2> "${OUTDIR}/quast/quast.stderr.log"

###############################################################################
# MERQURY
###############################################################################

echo ""
echo "============================================================"
echo "Preparing Merqury k-mer database"
echo "============================================================"

MERQURY_DIR="${OUTDIR}/merqury"
MERYL_DIR="${MERQURY_DIR}/meryl"

mkdir -p "${MERYL_DIR}"

###############################################################################
# MERYL: FEMALE ILLUMINA
###############################################################################

echo ""
echo "Building female Illumina k-mer database..."

${MERYL} k=21 count \
    "${FEMALE_R1}" \
    output "${MERYL_DIR}/female_R1.meryl"

${MERYL} k=21 count \
    "${FEMALE_R2}" \
    output "${MERYL_DIR}/female_R2.meryl"

${MERYL} union-sum \
    "${MERYL_DIR}/female_R1.meryl" \
    "${MERYL_DIR}/female_R2.meryl" \
    output "${MERYL_DIR}/female_illumina.meryl"

###############################################################################
# MERYL: MALE ILLUMINA
###############################################################################

echo ""
echo "Building male Illumina k-mer database..."

${MERYL} k=21 count \
    "${MALE_R1}" \
    output "${MERYL_DIR}/male_R1.meryl"

${MERYL} k=21 count \
    "${MALE_R2}" \
    output "${MERYL_DIR}/male_R2.meryl"

${MERYL} union-sum \
    "${MERYL_DIR}/male_R1.meryl" \
    "${MERYL_DIR}/male_R2.meryl" \
    output "${MERYL_DIR}/male_illumina.meryl"

###############################################################################
# MERQURY: FEMALE ASSEMBLY
###############################################################################

echo ""
echo "============================================================"
echo "Merqury: female assembly"
echo "============================================================"

${MERQURY} \
    "${MERYL_DIR}/female_illumina.meryl" \
    "${FEMALE_ASM}" \
    "${OUTDIR}/merqury/female_unpolished"

${MERQURY} \
    "${MERYL_DIR}/female_illumina.meryl" \
    "${OUTDIR}/female/female.polished.fasta" \
    "${OUTDIR}/merqury/female_polished"

###############################################################################
# MERQURY: MALE ASSEMBLY
###############################################################################

echo ""
echo "============================================================"
echo "Merqury: male assembly"
echo "============================================================"

${MERQURY} \
    "${MERYL_DIR}/male_illumina.meryl" \
    "${MALE_ASM}" \
    "${OUTDIR}/merqury/male_unpolished"

${MERQURY} \
    "${MERYL_DIR}/male_illumina.meryl" \
    "${OUTDIR}/male/male.polished.fasta" \
    "${OUTDIR}/merqury/male_polished"

###############################################################################
# BASIC FASTA STATISTICS
###############################################################################

echo ""
echo "============================================================"
echo "Generating final FASTA statistics"
echo "============================================================"

if command -v seqkit >/dev/null 2>&1; then

    seqkit stats \
        "${FEMALE_ASM}" \
        "${OUTDIR}/female/female.polished.fasta" \
        "${MALE_ASM}" \
        "${OUTDIR}/male/male.polished.fasta" \
        > "${OUTDIR}/assembly_seqkit_stats.txt"

else

    echo "seqkit not installed; skipping seqkit statistics."

fi

###############################################################################
# FINAL SUMMARY
###############################################################################

echo ""
echo "============================================================"
echo "WORKFLOW COMPLETED"
echo "============================================================"

echo ""
echo "Female polished assembly:"
echo "${OUTDIR}/female/female.polished.fasta"

echo ""
echo "Male polished assembly:"
echo "${OUTDIR}/male/male.polished.fasta"

echo ""
echo "QUAST:"
echo "${OUTDIR}/quast"

echo ""
echo "Merqury:"
echo "${OUTDIR}/merqury"

echo ""
echo "Full pipeline log:"
echo "${LOG}"

echo ""
echo "Finished: $(date)"
echo "============================================================"