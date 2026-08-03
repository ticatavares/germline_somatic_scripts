#!/bin/bash
set -e

mkdir -p qc_plink

for chr in {1..22}; do

    echo "====================================="
    echo "QC chromosome ${chr}"
    echo "====================================="

    plink2 \
        --bfile plink_files/chr${chr} \
        --mind 0.1 \
        --geno 0.1 \
        --make-bed \
        --out qc_plink/chr${chr}.qc

    echo "Finished chr${chr}"
    wc -l qc_plink/chr${chr}.qc.bim
    wc -l qc_plink/chr${chr}.qc.fam

done

echo "ALL QC FINISHED"
