#!/bin/bash
set -e

# ==========================================
# Calculate allele frequencies by superpopulation
# for LD-pruned 1KGP chromosome files
# ==========================================

mkdir -p AF_superpop
mkdir -p logs

POPS=("AFR" "AMR" "EAS" "EUR" "SAS")

for chr in {1..22}; do

    echo "====================================="
    echo "Chromosome ${chr}"
    echo "====================================="

    BFILE="ld_pruned/chr${chr}.LDpruned"

    if [ ! -f "${BFILE}.bed" ]; then
        echo "ERROR: Missing ${BFILE}.bed"
        exit 1
    fi

    for pop in "${POPS[@]}"; do

        echo "Calculating AF for chr${chr} - ${pop}"

        KEEP="keep_superpop_plink/${pop}.keep"
        OUT="AF_superpop/chr${chr}.${pop}"

        if [ ! -f "$KEEP" ]; then
            echo "ERROR: Missing keep file: $KEEP"
            exit 1
        fi

        plink2 \
            --bfile "$BFILE" \
            --keep "$KEEP" \
            --freq \
            --out "$OUT"

        echo "Finished chr${chr} - ${pop}"

    done

done

echo "====================================="
echo "ALL AF CALCULATIONS FINISHED"
echo "====================================="
