#!/usr/bin/env python3
# # extract_pcg.py
#
# This script is a pipeline build for mitochondial genes aligment to build phylogenetic tree.
#? PCG extraction
#? Gene-wise alignment (MAFFT)
#? Concatenation
#? Partition file
#? Robust logging (pcg_log.tsv)
#? Safe overwrite (no duplication issues)
#
# Usage:
#   python extract_pcg.py -i genbank_files -o results
#
# The report file should be tab-delimited with a header including at least:
#   sequence_name    start    end
# Coordinates are 1-based inclusive.


#!/usr/bin/env python3

import os
import argparse
import subprocess
from Bio import SeqIO

# =========================
# CONFIG
# =========================

TARGET_GENES = [
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "COXI", "COXII", "COXIII",
    "ATP6", "ATP8",
    "CYTB"
]

# =========================
# UTILITIES
# =========================

def prepare_output_dir(path):
    if os.path.exists(path):
        for f in os.listdir(path):
            os.remove(os.path.join(path, f))
    else:
        os.makedirs(path)

def clean_species_name(name):
    return name.replace(" ", "_")

# =========================
# ? ROBUST GENE DETECTION
# =========================

def detect_gene(feature):
    qualifiers = feature.qualifiers

    text = ""
    if "gene" in qualifiers:
        text += qualifiers["gene"][0] + " "
    if "product" in qualifiers:
        text += qualifiers["product"][0]

    text = text.upper().replace("-", "").replace("_", "").strip()

    # -------- ND genes --------
    if "ND1" in text:
        return "ND1"
    elif "ND2" in text:
        return "ND2"
    elif "ND3" in text:
        return "ND3"
    elif "ND4L" in text:
        return "ND4L"
    elif "ND4" in text:
        return "ND4"
    elif "ND5" in text:
        return "ND5"
    elif "ND6" in text:
        return "ND6"

    # -------- COX genes (FIXED properly) --------
    elif "COX" in text or "CYTOCHROME C OXIDASE" in text:

        # IMPORTANT: order matters (III before II before I)
        if "III" in text:
            return "COXIII"
        elif "II" in text:
            return "COXII"
        elif "I" in text:
            return "COXI"

    # -------- ATP genes --------
    elif "ATP6" in text or "ATPASE6" in text:
        return "ATP6"
    elif "ATP8" in text or "ATPASE8" in text:
        return "ATP8"

    # -------- CYTB --------
    elif "CYTB" in text or "CYTOCHROME B" in text or "CYTB" in text:
        return "CYTB"

    return None

# =========================
# STEP 1: EXTRACT + LOG
# =========================

def extract_pcgs(input_dir, output_dir, log_file):
    print("Step 1: Extracting PCGs...")
    prepare_output_dir(output_dir)

    species_gene_matrix = {}

    for file in os.listdir(input_dir):
        if file.endswith(".gb") or file.endswith(".gbk"):

            species = clean_species_name(file.split('.')[0])
            filepath = os.path.join(input_dir, file)

            species_gene_matrix[species] = {g: 0 for g in TARGET_GENES}

            for record in SeqIO.parse(filepath, "genbank"):
                for feature in record.features:

                    if feature.type != "CDS":
                        continue

                    gene = detect_gene(feature)

                    if gene:
                        seq = feature.extract(record.seq)
                        out_file = os.path.join(output_dir, f"{gene}.fasta")

                        with open(out_file, "a") as f:
                            f.write(f">{species}\n{seq}\n")

                        species_gene_matrix[species][gene] = 1

    # Write log
    with open(log_file, "w") as f:
        header = ["Species"] + TARGET_GENES
        f.write("\t".join(header) + "\n")

        for species, genes in species_gene_matrix.items():
            row = [species] + [str(genes[g]) for g in TARGET_GENES]
            f.write("\t".join(row) + "\n")

    print(f"? Log file: {log_file}")

# =========================
# STEP 2: ALIGNMENT
# =========================

def align_genes(pcg_dir, aligned_dir):
    print("Step 2: Aligning genes with MAFFT...")
    prepare_output_dir(aligned_dir)

    for gene_file in os.listdir(pcg_dir):
        if not gene_file.endswith(".fasta"):
            continue

        infile = os.path.join(pcg_dir, gene_file)
        outfile = os.path.join(aligned_dir, gene_file.replace(".fasta", "_aln.fasta"))

        cmd = f"mafft --auto {infile} > {outfile}"
        subprocess.run(cmd, shell=True)

    print("? Alignment complete")

# =========================
# STEP 3: CONCATENATION + CODON PARTITION
# =========================

def concatenate_alignments(aligned_dir, output_prefix):
    print("Step 3: Concatenating alignments...")

    gene_files = sorted([f for f in os.listdir(aligned_dir) if f.endswith("_aln.fasta")])

    sequences = {}
    partitions = []
    current_pos = 1

    for gene_file in gene_files:
        gene_name = gene_file.split("_")[0]
        filepath = os.path.join(aligned_dir, gene_file)

        alignment = list(SeqIO.parse(filepath, "fasta"))
        if not alignment:
            continue

        aln_length = len(alignment[0].seq)

        for record in alignment:
            if record.id not in sequences:
                sequences[record.id] = ""
            sequences[record.id] += str(record.seq)

        start = current_pos
        end = current_pos + aln_length - 1

        # codon partition
        partitions.append(f"DNA, {gene_name}_pos1 = {start}-{end}\\3")
        partitions.append(f"DNA, {gene_name}_pos2 = {start+1}-{end}\\3")
        partitions.append(f"DNA, {gene_name}_pos3 = {start+2}-{end}\\3")

        current_pos = end + 1

    with open(f"{output_prefix}.fasta", "w") as f:
        for sp, seq in sequences.items():
            f.write(f">{sp}\n{seq}\n")

    with open(f"{output_prefix}.partition", "w") as f:
        for p in partitions:
            f.write(p + "\n")

    print("? Concatenation complete")

# =========================
# MAIN
# =========================

def main():
    parser = argparse.ArgumentParser(description="Mitochondrial PCG pipeline")

    parser.add_argument("-i", "--input", required=True, help="GenBank folder")
    parser.add_argument("-o", "--output", required=True, help="Output folder")

    args = parser.parse_args()

    pcg_dir = os.path.join(args.output, "pcg_raw")
    aln_dir = os.path.join(args.output, "pcg_aligned")
    concat_prefix = os.path.join(args.output, "concatenated")
    log_file = os.path.join(args.output, "pcg_log.tsv")

    os.makedirs(args.output, exist_ok=True)

    extract_pcgs(args.input, pcg_dir, log_file)
    align_genes(pcg_dir, aln_dir)
    concatenate_alignments(aln_dir, concat_prefix)

    print("\n?? Pipeline completed successfully!")

if __name__ == "__main__":
    main()