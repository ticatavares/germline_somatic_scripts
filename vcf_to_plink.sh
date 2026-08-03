#!/bin/bash
set -e

mkdir -p plink_files

for chr in {1..22}; do

    echo "====================================="
    echo "Converting chromosome ${chr}"
    echo "====================================="

    INPUT="filtered_vcfs/chr${chr}.biallelic.snps.vcf.gz"
    OUTPUT="plink_files/chr${chr}"

    plink2 \
        --vcf "$INPUT" \
        --make-bed \
        --out "$OUTPUT"

    echo "Finished chromosome ${chr}"

done

echo "====================================="
echo "ALL CONVERSIONS FINISHED"
echo "====================================="
