#!/bin/bash

# ==========================================
# FILTER 1KGP HIGH-COVERAGE VCF FILES
# Keep only biallelic SNVs
# No FILTER / ExcHet filtering at this stage
# ==========================================

mkdir -p filtered_vcfs
mkdir -p logs

for chr in {1..22}; do

    echo "====================================="
    echo "Processing chromosome ${chr}"
    echo "====================================="

    INPUT="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    OUTPUT="filtered_vcfs/chr${chr}.biallelic.snps.vcf.gz"

    if [ ! -f "$INPUT" ]; then
        echo "ERROR: Input file not found: $INPUT"
        exit 1
    fi

    bcftools view \
        --threads 4 \
        -m2 -M2 \
        -v snps \
        "$INPUT" \
        -Oz \
        -o "$OUTPUT"

    tabix -f -p vcf "$OUTPUT"

    echo "Finished chromosome ${chr}"
    echo "Variant count:"
    bcftools index -n "$OUTPUT"

done

echo "====================================="
echo "ALL CHROMOSOMES FINISHED"
echo "====================================="
