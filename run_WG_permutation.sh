#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# INPUT FILES
# ============================================================

COMMON_BED="SNPs_WG.sorted.bed"
COSMIC_BED="cosmic_WG.sorted.bed"
GENOME_FILE="../chr_lengths_WG_.txt"

# Number of permutations and initial seed
N_PERM="${N_PERM:-1000}"
BASE_SEED="${BASE_SEED:-123}"

# Output directory
OUTPUT_DIR="permutation_HGDP_common_vs_COSMIC"
COUNTS_FILE="${OUTPUT_DIR}/permutation_overlap_counts.tsv"
SUMMARY_FILE="${OUTPUT_DIR}/permutation_summary.txt"

# ============================================================
# CHECK DEPENDENCIES AND INPUTS
# ============================================================

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools was not found."
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found."
    exit 1
}

for input_file in "$COMMON_BED" "$COSMIC_BED" "$GENOME_FILE"; do
    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: File not found or empty: $input_file"
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

TASK_TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TASK_TMP_DIR"' EXIT

COMMON_UNIQUE="${TASK_TMP_DIR}/common_unique.sorted.bed"
COSMIC_UNIQUE="${TASK_TMP_DIR}/cosmic_unique.sorted.bed"
SHUFFLED_BED="${TASK_TMP_DIR}/common_shuffled.sorted.bed"

# ============================================================
# PREPARE UNIQUE GENOMIC COORDINATES
# ============================================================

echo "Preparing unique common SNV coordinates..."

cut -f1-3 "$COMMON_BED" |
sort -k1,1V -k2,2n -k3,3n -u |
bedtools sort -faidx "$GENOME_FILE" \
> "$COMMON_UNIQUE"

echo "Preparing unique COSMIC coordinates..."

cut -f1-3 "$COSMIC_BED" |
sort -k1,1V -k2,2n -k3,3n -u |
bedtools sort -faidx "$GENOME_FILE" \
> "$COSMIC_UNIQUE"

N_COMMON=$(wc -l < "$COMMON_UNIQUE" | tr -d ' ')
N_COSMIC=$(wc -l < "$COSMIC_UNIQUE" | tr -d ' ')

echo "Unique HGDP common SNV coordinates: $N_COMMON"
echo "Unique COSMIC coordinates: $N_COSMIC"

# ============================================================
# OBSERVED OVERLAP
# ============================================================

OBSERVED=$(
    bedtools intersect \
        -sorted \
        -g "$GENOME_FILE" \
        -u \
        -a "$COMMON_UNIQUE" \
        -b "$COSMIC_UNIQUE" |
    wc -l |
    tr -d ' '
)

echo "Observed overlap: $OBSERVED"

if [[ "$OBSERVED" -ne 59569 ]]; then
    echo "WARNING: expected 59,569 observed overlaps, but obtained $OBSERVED."
    echo "Check whether the correct final BED files were provided."
fi

# ============================================================
# PERMUTATIONS
# ============================================================

printf "iteration\tseed\toverlap\n" > "$COUNTS_FILE"

echo "Starting $N_PERM chromosome-constrained permutations..."

for ((iteration=1; iteration<=N_PERM; iteration++)); do

    current_seed=$((BASE_SEED + iteration - 1))

    bedtools shuffle \
        -i "$COMMON_UNIQUE" \
        -g "$GENOME_FILE" \
        -chrom \
        -noOverlapping \
        -seed "$current_seed" |
    bedtools sort -faidx "$GENOME_FILE" \
    > "$SHUFFLED_BED"

    permuted_overlap=$(
        bedtools intersect \
            -sorted \
            -g "$GENOME_FILE" \
            -u \
            -a "$SHUFFLED_BED" \
            -b "$COSMIC_UNIQUE" |
        wc -l |
        tr -d ' '
    )

    printf "%d\t%d\t%d\n" \
        "$iteration" \
        "$current_seed" \
        "$permuted_overlap" \
        >> "$COUNTS_FILE"

    if (( iteration % 25 == 0 || iteration == N_PERM )); then
        echo "Completed $iteration/$N_PERM permutations"
    fi

done

# ============================================================
# EMPIRICAL STATISTICS
# ============================================================

python3 - \
    "$COUNTS_FILE" \
    "$SUMMARY_FILE" \
    "$OBSERVED" \
    "$N_COMMON" \
    "$N_COSMIC" \
    "$BASE_SEED" <<'PY'
import csv
import math
import statistics
import sys

counts_file = sys.argv[1]
summary_file = sys.argv[2]
observed = int(sys.argv[3])
n_common = int(sys.argv[4])
n_cosmic = int(sys.argv[5])
base_seed = int(sys.argv[6])

values = []

with open(counts_file, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        values.append(int(row["overlap"]))

if not values:
    raise RuntimeError("No permutation results were found.")


def percentile(data, probability):
    ordered = sorted(data)

    if len(ordered) == 1:
        return float(ordered[0])

    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)

    if lower == upper:
        return float(ordered[lower])

    fraction = position - lower

    return (
        ordered[lower] * (1 - fraction)
        + ordered[upper] * fraction
    )


n_permutations = len(values)
null_mean = statistics.mean(values)
null_median = statistics.median(values)
null_sd = statistics.stdev(values) if n_permutations > 1 else 0.0

lower_95 = percentile(values, 0.025)
upper_95 = percentile(values, 0.975)

b = sum(value >= observed for value in values)
empirical_p = (b + 1) / (n_permutations + 1)

fold_over_null = observed / null_mean

summary = f"""Chromosome-constrained permutation analysis
================================================
Initial random seed: {base_seed}
Number of permutations: {n_permutations}

Unique HGDP common SNVs: {n_common:,}
Unique COSMIC coordinates: {n_cosmic:,}
Observed overlap: {observed:,}

Null mean: {null_mean:,.3f}
Null median: {null_median:,.3f}
Null standard deviation: {null_sd:,.3f}
Empirical 95% interval: {lower_95:,.3f}–{upper_95:,.3f}

Permutations with overlap >= observed (b): {b}
Empirical P value: ({b} + 1) / ({n_permutations} + 1)
Empirical P value: {empirical_p:.6g}

Observed/null mean ratio: {fold_over_null:.4f}
"""

print()
print(summary)

with open(summary_file, "w") as handle:
    handle.write(summary)
PY

echo "Permutation counts: $COUNTS_FILE"
echo "Summary: $SUMMARY_FILecho "Summary: $SUMMARY_FILE"
