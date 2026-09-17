#!/usr/bin/env python3
# merge_all_beds.py
#
# This script merges multiple BED files for MCScanX and assigns unique chromosome names
# per input file by applying a per-file prefix to chromosome codes.
#
# Usage:
#   python merge_all_beds.py \
#       --beds path/to/bed1.bed path/to/bed2.bed ... \
#       --prefixes P1 P2 ... \
#       --output path/to/merged.bed \
#       [--force] [--sort] [--filter-unknown]
#
# Arguments:
#   --beds          Space-separated list of BED file paths to merge
#   --prefixes      Same number of prefixes as BED files; each prefix will be prepended to chromosome names in that file
#   --output        Path to output merged BED file
#   --force         Overwrite output if it exists
#   --sort          Sort merged entries by chromosome then start position
#   --filter-unknown Remove entries with gene names starting with 'UNKNOWN_'

import argparse
import os
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge multiple BED files and apply per-file chromosome prefixes for MCScanX."
    )
    parser.add_argument('--beds', nargs='+', required=True,
                        help='Paths to BED files to merge')
    parser.add_argument('--prefixes', nargs='+', required=True,
                        help='Prefixes corresponding to each BED file')
    parser.add_argument('--output', required=True,
                        help='Path to output merged BED file')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite the output file if it already exists')
    parser.add_argument('--sort', action='store_true',
                        help='Sort by chromosome and start position')
    parser.add_argument('--filter-unknown', action='store_true',
                        help="Filter out entries with gene names starting with 'UNKNOWN_'")
    return parser.parse_args()


def rename_chromosomes(df):
    """
    Rename raw chromosome/scaffold IDs to a normalized form.
    """
    def rename(x):
        if x.startswith('NC_'):
            # Zebrafish refseq example
            if '.2' in x:
                return 'drMT'
            num = int(x.split('_')[-1].split('.')[0][-2:]) - 11
            return f'{num}'
        if x.startswith('NW_'):
            num = x.split('_')[-1][:-2]
            return f'{num}'
        if x.startswith('JAPD'):
            return f'{x.split("_")[-1]}'
        if x.startswith('scaffold'):
            return f'{x.split("_")[-1]}'
        return x

    df[0] = df[0].apply(rename)
    return df


def merge_beds(bed_paths, prefixes, output_bed, force, sort, filter_unknown):
    if len(bed_paths) != len(prefixes):
        raise ValueError('Number of prefixes must match number of BED files')
    if os.path.exists(output_bed) and not force:
        raise IOError(f"File {output_bed} already exists. Use --force to overwrite.")

    merged_frames = []
    for bed, prefix in zip(bed_paths, prefixes):
        df = pd.read_csv(bed, sep='\t', header=None)
        df = rename_chromosomes(df)
        # apply per-file prefix to chromosome column
        df[0] = df[0].apply(lambda x: f"{prefix}_{x}")
        merged_frames.append(df)

    merged_df = pd.concat(merged_frames, ignore_index=True)

    if filter_unknown:
        merged_df = merged_df[~merged_df[1].str.startswith('UNKNOWN_')]

    if sort:
        merged_df.sort_values(by=[0, 2], inplace=True)

    merged_df.to_csv(output_bed, sep='\t', header=False, index=False)


def main():
    args = parse_args()
    merge_beds(
        args.beds,
        args.prefixes,
        args.output,
        args.force,
        args.sort,
        args.filter_unknown
    )

if __name__ == '__main__':
    main()
