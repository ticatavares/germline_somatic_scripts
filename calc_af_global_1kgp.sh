#!/bin/bash
set -e

# ==========================================
# Calculate GLOBAL allele frequencies
# using LD-pruned 1KGP variants
# ==========================================

mkdir -p AF_global

for chr in {1..22}; do

    echo "====================================="
    echo "Calculating GLOBAL AF - chr${chr}"
    echo "====================================="

    plink2 \
        --bfile ld_pruned/chr${chr}.LDpruned \
        --freq \
        --out AF_global/chr${chr}

    echo "Finished chr${chr}"

done

echo "====================================="
echo "ALL GLOBAL AF CALCULATIONS FINISHED"
echo "====================================="
