#!/bin/bash

BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

for chr in {1..22}; do

    FILE="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

    echo "Downloading chr${chr}..."

    curl -L --fail --retry 5 --retry-delay 10 -C - -O "${BASE}/${FILE}"

    curl -L --fail --retry 5 --retry-delay 10 -C - -O "${BASE}/${FILE}.tbi"

done

echo "DONE"

