#!/bin/bash

set -euo pipefail

#############################################
# USER SETTINGS
#############################################

MALE_PAF="male_to_female.paf"
FEMALE_PAF="female_to_male.paf"

OUTDIR="assembly_comparison_paf"

#############################################
# CHECK INPUTS
#############################################

if [[ ! -f "$MALE_PAF" ]]; then
    echo "ERROR: Cannot find $MALE_PAF"
    exit 1
fi

if [[ ! -f "$FEMALE_PAF" ]]; then
    echo "ERROR: Cannot find $FEMALE_PAF"
    exit 1
fi

mkdir -p "$OUTDIR"

#############################################
# FUNCTIONS
#############################################

analyze_paf () {

    PAF="$1"
    PREFIX="$2"

    echo "=========================================="
    echo "Analyzing: $PAF"
    echo "=========================================="

    #########################################
    # 1. Basic alignment statistics
    #########################################

    awk '
    BEGIN {
        total=0;
        matches=0;
        aln=0;
        n=0;
    }
    {
        qlen=$2;
        qstart=$3;
        qend=$4;
        strand=$5;
        tlen=$6;
        tstart=$7;
        tend=$8;
        nmatch=$10;
        block=$11;

        total++;
        matches += nmatch;
        aln += block;

        if (block > maxblock) {
            maxblock=block;
            maxq=$1;
            maxt=$6;
        }
    }
    END {
        print "Total PAF alignments:", total;
        print "Total aligned bases:", aln;
        print "Total matching bases:", matches;

        if (aln > 0)
            print "Overall sequence identity:", matches/aln;

        print "Longest alignment:", maxblock;
    }
    ' "$PAF" > "$OUTDIR/${PREFIX}.basic_stats.txt"


    #########################################
    # 2. Alignment length distribution
    #########################################

    awk '
    {
        print $11
    }
    ' "$PAF" |
    sort -n > "$OUTDIR/${PREFIX}.alignment_lengths.txt"


    #########################################
    # 3. Identity distribution
    #########################################

    awk '
    {
        if ($11 > 0)
            print $10/$11
    }
    ' "$PAF" |
    sort -n > "$OUTDIR/${PREFIX}.identity_values.txt"


    #########################################
    # 4. High-confidence alignments
    #########################################
    # Criteria:
    # alignment >= 10 kb
    # identity >= 95%
    #########################################

    awk '
    ($11 >= 10000 && ($10/$11) >= 0.95)
    ' "$PAF" > "$OUTDIR/${PREFIX}.highconf_10kb_95id.paf"


    #########################################
    # 5. More stringent alignments
    #########################################

    awk '
    ($11 >= 100000 && ($10/$11) >= 0.98)
    ' "$PAF" > "$OUTDIR/${PREFIX}.highconf_100kb_98id.paf"


    #########################################
    # 6. Summary of high-confidence alignment
    #########################################

    awk '
    BEGIN {
        n=0;
        bases=0;
        matches=0;
    }
    {
        n++;
        bases += $11;
        matches += $10;
    }
    END {
        print "Number of high-confidence alignments:", n;
        print "Aligned bases:", bases;
        if (bases > 0)
            print "Mean identity:", matches/bases;
    }
    ' "$OUTDIR/${PREFIX}.highconf_10kb_95id.paf" \
    > "$OUTDIR/${PREFIX}.highconf_stats.txt"


    #########################################
    # 7. Query coverage
    #########################################
    # NOTE:
    # This calculates coverage as the sum of aligned
    # blocks per query sequence. It is NOT corrected
    # for overlapping alignments.
    #########################################

    awk '
    {
        q=$1;
        len=$2;
        cov[q] += $4-$3;
        qlen[q] = len;
    }
    END {
        for (q in cov)
            print q, qlen[q], cov[q], cov[q]/qlen[q]
    }
    ' "$OUTDIR/${PREFIX}.highconf_10kb_95id.paf" |
    sort -k4,4nr \
    > "$OUTDIR/${PREFIX}.query_coverage.txt"


    #########################################
    # 8. Target coverage
    #########################################

    awk '
    {
        t=$6;
        len=$7;
        cov[t] += $8-$7;
        tlen[t] = $7;
    }
    END {
        for (t in cov)
            print t, cov[t], cov[t]/tlen[t]
    }
    ' "$OUTDIR/${PREFIX}.highconf_10kb_95id.paf" |
    sort -k3,3nr \
    > "$OUTDIR/${PREFIX}.target_coverage.txt"


    #########################################
    # 9. Strand information
    #########################################

    awk '
    {
        strand[$5]++;
    }
    END {
        print "+", strand["+"];
        print "-", strand["-"];
    }
    ' "$PAF" > "$OUTDIR/${PREFIX}.strand_stats.txt"


    #########################################
    # 10. Alignment identity bins
    #########################################

    awk '
    {
        id=$10/$11;

        if (id >= 0.99) bin=">=99%";
        else if (id >= 0.98) bin="98-99%";
        else if (id >= 0.95) bin="95-98%";
        else if (id >= 0.90) bin="90-95%";
        else bin="<90%";

        count[bin]++;
    }
    END {
        print ">=99%", count[">=99%"];
        print "98-99%", count["98-99%"];
        print "95-98%", count["95-98%"];
        print "90-95%", count["90-95%"];
        print "<90%", count["<90%"];
    }
    ' "$PAF" > "$OUTDIR/${PREFIX}.identity_bins.txt"


    #########################################
    # 11. Generate BED files
    #########################################

    awk '
    ($11 >= 10000 && ($10/$11) >= 0.95) {
        print $1, $3, $4
    }
    ' OFS="\t" \
    "$PAF" > "$OUTDIR/${PREFIX}.query.highconf.bed"


    awk '
    ($11 >= 10000 && ($10/$11) >= 0.95) {
        print $6, $7, $8
    }
    ' OFS="\t" \
    "$PAF" > "$OUTDIR/${PREFIX}.target.highconf.bed"


    echo "Finished: $PREFIX"
    echo

}


#############################################
# RUN BOTH DIRECTIONS
#############################################

analyze_paf "$MALE_PAF" "male_vs_female"

analyze_paf "$FEMALE_PAF" "female_vs_male"


#############################################
# RECIPROCAL SUMMARY
#############################################

echo "=========================================="
echo "Generating reciprocal summary"
echo "=========================================="

{
    echo -e "Direction\tTotal_alignments\tAligned_bases\tMatching_bases\tMean_identity"

    awk '
    BEGIN {n=0; a=0; m=0}
    {
        n++;
        a += $11;
        m += $10;
    }
    END {
        if(a>0)
            printf "Male_to_Female\t%d\t%d\t%d\t%.6f\n",n,a,m,m/a;
    }
    ' "$MALE_PAF"

    awk '
    BEGIN {n=0; a=0; m=0}
    {
        n++;
        a += $11;
        m += $10;
    }
    END {
        if(a>0)
            printf "Female_to_Male\t%d\t%d\t%d\t%.6f\n",n,a,m,m/a;
    }
    ' "$FEMALE_PAF"

} > "$OUTDIR/reciprocal_summary.tsv"


#############################################
# FINISHED
#############################################

echo
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo
echo "Results:"
echo "  $OUTDIR/"
echo
ls -lh "$OUTDIR"