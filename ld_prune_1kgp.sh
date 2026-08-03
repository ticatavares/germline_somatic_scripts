#!/bin/bash
set -e

mkdir -p ld_pruned

for chr in {1..22}; do

    echo "====================================="
    echo "LD pruning chromosome ${chr}"
    echo "====================================="

    plink2 \
        --bfile qc_plink/chr${chr}.qc \
        --indep-pairwise 200kb 1 0.4 \
        --out ld_pruned/chr${chr}

    plink2 \
        --bfile qc_plink/chr${chr}.qc \
        --extract ld_pruned/chr${chr}.prune.in \
        --make-bed \
        --out ld_pruned/chr${chr}.LDpruned

    echo "Finished chr${chr}"
    echo "Variants retained after LD pruning:"
    wc -l ld_pruned/chr${chr}.LDpruned.bim

done

echo "ALL LD PRUNING FINISHED"
