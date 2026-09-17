# prepare_mcscanx_unified.py
#
# This script prepares renamed protein FASTA files and BED files with gene symbols for MCScanX analysis.
# It supports annotations from RefSeq, Ensembl, and custom GFF3 files.
#
# Usage:
#   python prepare_mcscanx_unified.py \
#       --type refseq|ensemble|custom \
#       --gff path/to/annotations.gff3 \
#       --fasta path/to/proteins.fasta \
#       --outdir path/to/output_dir \
#       --prefix speciesprefix \
#       [--force]
#
# Options:
#   --type     Annotation type: 'refseq', 'ensemble', or 'custom'
#   --gff      Path to the GFF3 file
#   --fasta    Path to the protein FASTA file
#   --outdir   Output directory
#   --prefix   Prefix for gene symbols and output files (e.g. Ggeck)
#   --force    Overwrite existing files and regenerate GFF database

import argparse
import os
import gffutils
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare protein FASTA and BED files with gene symbols for MCScanX." )
    parser.add_argument('--type', required=True, choices=['refseq', 'ensemble', 'custom'],
                        help='Annotation type')
    parser.add_argument('--gff', required=True, help='Path to GFF3 file')
    parser.add_argument('--fasta', required=True, help='Path to protein FASTA file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--prefix', required=True,
                        help='Prefix for gene symbols and output filenames (e.g. Ggeck)')
    parser.add_argument('--force', action='store_true', help='Overwrite existing outputs and rebuild DB')
    return parser.parse_args()


def create_or_load_db(gff_path, outdir, force):
    db_path = os.path.join(outdir, os.path.basename(gff_path) + '.db')
    if force or not os.path.exists(db_path):
        print(f"Creating GFF database at {db_path}...")
        return gffutils.create_db(
            gff_path, dbfn=db_path, force=True,
            keep_order=True, merge_strategy='create_unique', sort_attribute_values=True
        )
    return gffutils.FeatureDB(db_path, keep_order=True)


def build_gene_map(gff_db, annotation_type, prefix):
    gene_map = {}
    count = defaultdict(int)
    for gene in gff_db.features_of_type('gene'):
        gid = gene.attributes.get('ID', [None])[0]
        if not gid:
            continue
        if annotation_type == 'refseq':
            symbol = gene.attributes.get('gene', gene.attributes.get('Name', [gid]))[0]
        else:
            symbol = gene.attributes.get('Name', [gid])[0]
        base = f"{prefix}_{symbol}"
        count[base] += 1
        name = f"{base}_{count[base]}" if count[base] > 1 else base
        gene_map[gid.split('.')[0]] = name
    return gene_map


def build_protein_map_refseq(gff_db, gene_map):
    prot_map = {}
    for cds in gff_db.features_of_type('CDS'):
        pids = cds.attributes.get('protein_id', [])
        # find gene via transcript parent
        tid = cds.attributes.get('Parent', [None])[0]
        if not tid:
            continue
        try:
            transcript = gff_db[tid]
            gid = transcript.attributes.get('Parent', [None])[0]
        except Exception:
            gid = None
        symbol = gene_map.get(gid.split('.')[0]) if gid else None
        for pid in pids:
            key = pid.split('.')[0]
            prot_map[key] = symbol if symbol else key
    return prot_map


def build_protein_map_ensembl(gff_db, gene_map):
    prot_map = {}
    # map transcript to gene
    tr2gene = {}
    for mrna in gff_db.features_of_type('mRNA'):
        tid = mrna.attributes.get('ID', [None])[0]
        gid = mrna.attributes.get('Parent', [None])[0]
        if tid and gid:
            tr2gene[tid.split('.')[0]] = gid.split('.')[0]
    # map protein via CDS
    for cds in gff_db.features_of_type('CDS'):
        pids = cds.attributes.get('protein_id', [])
        tid = cds.attributes.get('Parent', [None])[0]
        gene_id = tr2gene.get(tid.split('.')[0]) if tid else None
        symbol = gene_map.get(gene_id) if gene_id else None
        for pid in pids:
            key = pid.split('.')[0]
            prot_map[key] = symbol if symbol else key
    return prot_map


def build_protein_map_custom(gff_db, gene_map):
    prot_map = {}
    # mRNA -> gene symbol
    for mrna in gff_db.features_of_type('mRNA'):
        tid_full = mrna.attributes.get('ID', [None])[0]
        gid_full = mrna.attributes.get('Parent', [None])[0]
        if tid_full and gid_full:
            tid = tid_full.split('.')[0]
            gid = gid_full.split('.')[0]
            symbol = gene_map.get(gid)
            prot_map[tid] = symbol if symbol else tid
    return prot_map


def rename_fasta_headers(fasta_path, out_path, prot_map):
    records = []
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        key = rec.id.split('.')[0]
        name = prot_map.get(key, key)
        rec.id = name
        rec.description = ''
        records.append(rec)
    SeqIO.write(records, out_path, 'fasta')


def write_bed(gff_db, out_path, gene_map):
    with open(out_path, 'w') as fo:
        for gene in gff_db.features_of_type('gene'):
            gid = gene.attributes.get('ID', [None])[0]
            key = gid.split('.')[0] if gid else None
            name = gene_map.get(key, key)
            fo.write(f"{gene.chrom}\t{name}\t{gene.start}\t{gene.stop}\n")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    gff_db = create_or_load_db(args.gff, args.outdir, args.force)

    # build gene symbol map
    gene_map = build_gene_map(gff_db, args.type, args.prefix)

    # build protein map
    if args.type == 'refseq':
        prot_map = build_protein_map_refseq(gff_db, gene_map)
    elif args.type == 'ensemble':
        prot_map = build_protein_map_ensembl(gff_db, gene_map)
    else:
        prot_map = build_protein_map_custom(gff_db, gene_map)

    # write FASTA
    fasta_out = os.path.join(args.outdir, f"{args.prefix}.fa")
    rename_fasta_headers(args.fasta, fasta_out, prot_map)

    # write BED
    bed_out = os.path.join(args.outdir, f"{args.prefix}.bed")
    write_bed(gff_db, bed_out, gene_map)

if __name__ == '__main__':
    main()
