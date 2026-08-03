#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# INPUTS
# ============================================================

SHARED_SOURCE="COSMIC_at_common_shared_WG.bed"

CPG_FILE="cpgIsland_hg38_noScaffolds_sorted.bed.gz"
MICROSAT_FILE="microsat_merged_1to22.bed.gz"

# Adjust this path if the file is in the current directory
GENOME_FILE="../chr_lengths_WG_.txt"

OUTPUT_DIR="hotspot_overlap_final_WG"

# ============================================================
# CHECKS
# ============================================================

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools was not found."
    exit 1
}

for input_file in \
    "$SHARED_SOURCE" \
    "$CPG_FILE" \
    "$MICROSAT_FILE" \
    "$GENOME_FILE"
do
    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: file not found or empty: $input_file"
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

TASK_TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TASK_TMP_DIR"' EXIT

AUTOSOME_GENOME="${TASK_TMP_DIR}/autosomes.genome"
SHARED_BED="${OUTPUT_DIR}/common_COSMIC_shared_WG.unique.bed"
CPG_BED="${OUTPUT_DIR}/CpG_islands.autosomes.merged.bed"
MICROSAT_BED="${OUTPUT_DIR}/microsatellites.autosomes.merged.bed"

SUMMARY_FILE="${OUTPUT_DIR}/hotspot_overlap_summary.tsv"

# ============================================================
# FUNCTIONS
# ============================================================

read_bed_file() {
    local input_file="$1"

    if [[ "$input_file" == *.gz ]]; then
        gzip -cd "$input_file"
    else
        command cat "$input_file"
    fi
}

prepare_region_file() {
    local input_file="$1"
    local output_file="$2"

    read_bed_file "$input_file" |
    awk 'BEGIN{OFS="\t"}
        $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ &&
        $2 ~ /^[0-9]+$/ &&
        $3 ~ /^[0-9]+$/ &&
        $3 > $2 {
            print $1,$2,$3
        }
    ' |
    bedtools sort -faidx "$AUTOSOME_GENOME" |
    bedtools merge \
    > "$output_file"
}

analyze_hotspot() {
    local label="$1"
    local region_bed="$2"
    local raw_output="$3"

    local overlap_count
    local percent
    local region_count
    local covered_bases
    local fisher_line
    local fisher_left
    local fisher_right
    local fisher_two_tail
    local fisher_or

    region_count=$(wc -l < "$region_bed" | tr -d ' ')

    covered_bases=$(
        awk '{sum += $3-$2} END{print sum+0}' "$region_bed"
    )

    overlap_count=$(
        bedtools intersect \
            -sorted \
            -g "$AUTOSOME_GENOME" \
            -u \
            -a "$SHARED_BED" \
            -b "$region_bed" |
        wc -l |
        tr -d ' '
    )

    percent=$(
        awk -v overlap="$overlap_count" -v total="$N_SHARED" \
            'BEGIN{printf "%.6f", (overlap/total)*100}'
    )

    bedtools fisher \
        -g "$AUTOSOME_GENOME" \
        -a "$SHARED_BED" \
        -b "$region_bed" \
        > "$raw_output"

    fisher_line=$(tail -n 1 "$raw_output")

    read -r \
        fisher_left \
        fisher_right \
        fisher_two_tail \
        fisher_or \
        <<< "$fisher_line"

    printf "%s\t%d\t%d\t%.6f\t%d\t%d\t%s\t%s\t%s\n" \
        "$label" \
        "$N_SHARED" \
        "$overlap_count" \
        "$percent" \
        "$region_count" \
        "$covered_bases" \
        "$fisher_or" \
        "$fisher_right" \
        "$fisher_two_tail" \
        >> "$SUMMARY_FILE"

    echo
    echo "$label"
    echo "  Shared loci:               $N_SHARED"
    echo "  Loci overlapping context:  $overlap_count"
    echo "  Percentage:                ${percent}%"
    echo "  Merged context intervals:  $region_count"
    echo "  Context coverage:          $covered_bases bp"
    echo "  Fisher OR:                 $fisher_or"
    echo "  Fisher right-tail P:       $fisher_right"
    echo "  Fisher two-tail P:         $fisher_two_tail"
}

# ============================================================
# AUTOSOMAL GENOME FILE
# ============================================================

awk 'BEGIN{OFS="\t"}
    $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {
        print $1,$2
    }
' "$GENOME_FILE" > "$AUTOSOME_GENOME"

N_CHROMOSOMES=$(wc -l < "$AUTOSOME_GENOME" | tr -d ' ')

if [[ "$N_CHROMOSOMES" -ne 22 ]]; then
    echo "ERROR: expected 22 autosomes in $GENOME_FILE, found $N_CHROMOSOMES."
    echo "Check chromosome names and genome-file path."
    exit 1
fi

# ============================================================
# SHARED COMMON-COSMIC COORDINATES
# ============================================================

cut -f1-3 "$SHARED_SOURCE" |
awk 'BEGIN{OFS="\t"}
    $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ &&
    $2 ~ /^[0-9]+$/ &&
    $3 ~ /^[0-9]+$/ &&
    $3 > $2 {
        print $1,$2,$3
    }
' |
sort -k1,1V -k2,2n -k3,3n -u |
bedtools sort -faidx "$AUTOSOME_GENOME" \
> "$SHARED_BED"

N_SHARED=$(wc -l < "$SHARED_BED" | tr -d ' ')

echo "Unique common-COSMIC shared coordinates: $N_SHARED"

if [[ "$N_SHARED" -ne 59569 ]]; then
    echo "ERROR: expected 59,569 shared coordinates, but obtained $N_SHARED."
    echo "Do not continue until the shared-locus input has been verified."
    exit 1
fi

# ============================================================
# PREPARE ANNOTATION DATASETS
# ============================================================

echo "Preparing CpG-island intervals..."
prepare_region_file "$CPG_FILE" "$CPG_BED"

echo "Preparing microsatellite intervals..."
prepare_region_file "$MICROSAT_FILE" "$MICROSAT_BED"

# ============================================================
# ANALYSES
# ============================================================

printf \
"context\tshared_loci\toverlapping_loci\tpercentage\tmerged_intervals\tcovered_bases\tOR\tP_right\tP_two_tail\n" \
> "$SUMMARY_FILE"

analyze_hotspot \
    "CpG_islands" \
    "$CPG_BED" \
    "${OUTPUT_DIR}/CpG_islands.bedtools_fisher.txt"

analyze_hotspot \
    "Microsatellites" \
    "$MICROSAT_BED" \
    "${OUTPUT_DIR}/microsatellites.bedtools_fisher.txt"

echo
echo "Analysis completed."
echo "Summary: $SUMMARY_FILE"
echo "Processed CpG BED: $CPG_BED"
echo "Processed microsatellite BED: $MICROSAT_BED"
