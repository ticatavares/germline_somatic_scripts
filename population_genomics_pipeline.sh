#!/bin/bash

# ======================================================================
# Population Genomics VCF Processing Pipeline
# Author: Thais Tavares
# Date: 09/04/2025
# Description: Processes raw VCF files through quality filtering, LD pruning,
#              and format conversion for downstream population genetics analyses
# ======================================================================

# ----------------------------------------------------------------------
# Configuration Section (USER MUST MODIFY)
# ----------------------------------------------------------------------

# Define input and output directories
input_dir=""       # Directory containing input VCF files
output_dir=""      # Directory for processed output files

# Path to bcftools (modify if not in PATH)
BCFTOOLS="bcftools"

# ----------------------------------------------------------------------
# Input Validation
# ----------------------------------------------------------------------

# Check if directories are provided
if [[ -z "$input_dir" || -z "$output_dir" ]]; then
echo "ERROR: Please set input_dir and output_dir variables" >&2
exit 1
fi

# Check if bcftools is available
if ! command -v "$BCFTOOLS" &> /dev/null; then
echo "ERROR: bcftools not found. Please install or specify path." >&2
exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# ----------------------------------------------------------------------
# Main Processing Loop (Chromosomes 1-21)
# ----------------------------------------------------------------------

for i in {1..21}; do
echo "Processing chromosome ${i}..."

# Input and output file names
input_vcf="${input_dir}/${i}file1_hgdp_snps.vcf.gz"
base_out="${output_dir}/${i}chrm_hgdp"

# --------------------------------------------------------------
# Step 1: Initial VCF Processing
# --------------------------------------------------------------

# 1.1 Index the input VCF file
echo "Indexing input VCF..."
$BCFTOOLS index -t "$input_vcf"

# 1.2 Filter for biallelic SNPs only
echo "Filtering biallelic SNPs..."
$BCFTOOLS view -v snps --max-alleles 2 --threads 50 "$input_vcf" \
-Oz -o "${base_out}_biallelic.vcf.gz"

# 1.3 Filter for PASS variants only
echo "Filtering PASS variants..."
$BCFTOOLS view -f PASS --threads 50 "${base_out}_biallelic.vcf.gz" \
-Oz -o "${base_out}_bial_PASS.vcf.gz"

# --------------------------------------------------------------
# Step 2: Sample Filtering
# --------------------------------------------------------------

# 2.1 Remove problematic transcriptome samples
echo "Removing transcriptome samples..."
$BCFTOOLS view --force-samples \
-s "^HGDP00213,HGDP00721,HGDP00237,HGDP00877,HGDP00992,HGDP00959,HGDP00720,HGDP00858,HGDP01262,HGDP01259,HGDP00462,HGDP00467,HGDP00239,HGDP00258,HGDP00222,HGDP00715,HGDP00471,HGDP00854,HGDP00948,HGDP00950,HGDP00955,HGDP00967,HGDP00712,HGDP01275,HGDP00856,HGDP00868,HGDP00964,HGDP01277,HGDP01258,HGDP00247,HGDP00711,HGDP00716,HGDP01029,HGDP01081,HGDP01264" \
"${base_out}_bial_PASS.vcf.gz" --threads 50 \
-Oz -o "${base_out}_bial_PASS_noTrans.vcf.gz"

# 2.2 Filter samples based on population list
echo "Filtering populations..."
$BCFTOOLS view --force-samples \
-S "${input_dir}/ids_popsKeepUmap_line.txt" \
"${base_out}_bial_PASS_noTrans.vcf.gz" --threads 50 \
-Oz -o "${base_out}_bial_PASS_noTrans_filterUmapPops.vcf.gz"

# --------------------------------------------------------------
# Step 3: PLINK2 Conversion and LD Pruning
# --------------------------------------------------------------

# 3.1 Convert to PLINK2 format with unique variant IDs
echo "Converting to PLINK2 format..."
plink2 --threads 50 \
--vcf "${base_out}_bial_PASS_noTrans_filterUmapPops.vcf.gz" \
--set-all-var-ids '@:#[b37]\$r,\$a' \
--new-id-max-allele-len 227 \
--mind 0.1 \
--geno 0.1 \
--make-bed \
--out "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind"

# 3.2 Perform LD pruning (200 SNP window, 25 SNP step, r² = 0.4)
echo "Performing LD pruning..."
plink2 --threads 50 \
--bfile "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind" \
--indep-pairwise 200 25 0.4 \
--out "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind_LD"

# 3.3 Create pruned dataset (keeping LD-independent variants)
echo "Creating pruned dataset..."
plink2 --threads 50 \
--bfile "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind" \
--extract "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind_LD.prune.in" \
--make-bed \
--out "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind_LDpruned"

# --------------------------------------------------------------
# Step 4: Final VCF Export
# --------------------------------------------------------------

# 4.1 Export back to VCF format
echo "Exporting to VCF..."
plink2 --threads 60 \
--bfile "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind_LDpruned" \
--export vcf \
--out "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind_LDpruned_vcf"

# 4.2 Compress and index final VCF
echo "Compressing and indexing final VCF..."
pigz -p 50 "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind_LDpruned_vcf.vcf"
$BCFTOOLS index -t "${base_out}_bial_PASS_noTrans_filterUmapPops_genoMind_LDpruned_vcf.vcf.gz"

echo "Completed processing for chromosome ${i}"
echo "--------------------------------------------------"
done

echo "Pipeline completed successfully for all chromosomes!"