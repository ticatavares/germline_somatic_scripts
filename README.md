# Germline-somatic positional recurrence analysis

This repository contains scripts used to preprocess population variation data and to analyze exact-coordinate recurrence between germline single-nucleotide variants (SNVs) and somatic cancer mutations. The workflow supports the analyses described in the associated manuscript, including population-frequency classification, linkage-disequilibrium (LD) pruning, nucleotide-substitution concordance, chromosome-constrained permutation testing, and overlap with canonical mutational contexts.

## Study design and reproducibility scope

The Human Genome Diversity Project (HGDP) was the primary population reference panel, and the high-coverage 1000 Genomes Project (1KGP) was used as an independent replication panel. Both panels were analyzed using the same overall framework and the same common/rare frequency thresholds:

- only autosomal, biallelic SNVs on GRCh38 were retained;
- population SNVs were quality controlled and LD pruned;
- global and group-specific allele frequencies were calculated;
- a variant was classified as **common** when its global allele frequency (AF) was at least 1% and its AF was at least 1% in at least two dataset-defined population groups;
- a variant was classified as **rare** when its global AF was below 1%; and
- variants with global AF at least 1% that did not meet the two-group criterion were retained as an intermediate category but were not included in the common-versus-rare comparisons.

The repository provides the organized 1KGP implementation as the worked population-data workflow. The HGDP analysis followed the same analytical logic and classification thresholds, but required dataset-specific inputs, quality annotations, sample metadata, and population-group labels. HGDP used seven source-defined geographic regions, whereas 1KGP used the five official superpopulations (AFR, AMR, EAS, EUR, and SAS). Thus, the 1KGP scripts document the shared computational workflow without implying that the two source datasets have identical file structures or metadata.

The exact analysis-ready BED files used for the HGDP discovery analysis and the 1KGP replication analysis are available from Zenodo:

**Zenodo:** `[ZENODO_RECORD_TITLE]` - `[ZENODO_RECORD_URL]`  
**DOI:** `<ZENODO_DOI>`

## Repository contents

| Script | Purpose | Main output |
|---|---|---|
| `download_1kgp.sh` | Downloads phased high-coverage 1KGP chromosome VCFs and indexes | Per-chromosome VCF and TBI files |
| `filter_1kgp.sh` | Retains biallelic SNVs | `filtered_vcfs/*.biallelic.snps.vcf.gz` |
| `vcf_to_plink.sh` | Converts filtered VCFs to PLINK binary format | `plink_files/chr*` |
| `qc_plink.sh` | Removes samples and variants with more than 10% missingness | `qc_plink/*.qc.{bed,bim,fam}` |
| `ld_prune_1kgp.sh` | Performs LD pruning and writes pruned PLINK datasets | `ld_pruned/*.LDpruned.{bed,bim,fam}` |
| `calc_af_global_1kgp.sh` | Calculates global AFs in the LD-pruned panel | `AF_global/*.afreq` |
| `calc_af_superpop_1kgp.sh` | Calculates AFs within the five 1KGP superpopulations | `AF_superpop/*.afreq` |
| `classify_1kgp_common_rare_to_bed.py` | Applies the study definitions of common, rare, and intermediate variants | Classification BED files and summary TSV |
| `compare_shared_alleles.py` | Compares REF>ALT substitutions at exact-coordinate germline-somatic overlaps | Per-locus and summary TSV files |
| `run_WG_permutation.sh` | Runs chromosome-constrained permutations of HGDP common-SNV coordinates | Permutation counts and empirical summary |
| `check_hotspot_overlap.sh` | Tests shared loci for overlap with CpG islands and microsatellites | BEDTools Fisher outputs and summary TSV |

## Software requirements

- Bash
- `curl`
- BCFtools and HTSlib/Tabix
- PLINK 2
- BEDTools
- Python 3
- Python package: `pandas`

Record the exact software versions used in the final archived release. They can be reported with, for example:

```bash
bcftools --version
plink2 --version
bedtools --version
python3 --version
python3 -c "import pandas; print(pandas.__version__)"
```

## Population-data workflow: high-coverage 1KGP

Run the scripts from the repository root in the following order:

```bash
bash download_1kgp.sh
bash filter_1kgp.sh
bash vcf_to_plink.sh
bash qc_plink.sh
bash ld_prune_1kgp.sh
bash calc_af_global_1kgp.sh
bash calc_af_superpop_1kgp.sh
python3 classify_1kgp_common_rare_to_bed.py
```

### Required superpopulation sample files

`calc_af_superpop_1kgp.sh` expects the following PLINK keep files:

```text
keep_superpop_plink/AFR.keep
keep_superpop_plink/AMR.keep
keep_superpop_plink/EAS.keep
keep_superpop_plink/EUR.keep
keep_superpop_plink/SAS.keep
```

These files must contain the sample identifiers assigned to each official 1KGP superpopulation and must match the identifiers in the PLINK sample files. The sample-to-superpopulation metadata source and the commands used to generate these keep files should be included with the archived workflow.

### Adapting the population workflow to HGDP

The same sequence of operations was used for HGDP, with the following dataset-specific substitutions:

| Component | 1KGP implementation | HGDP adaptation |
|---|---|---|
| Source input | High-coverage 1KGP chromosome VCFs | HGDP GRCh38 VCF input used in the study |
| Population grouping | Five official superpopulations | Seven source-defined geographic regions |
| Source-specific filtering | Biallelic SNV selection | `PASS` SNVs, excluding records annotated as `LOW_VQSLOD` or `ExcHet`, followed by biallelic SNV selection |
| Missingness QC | `--mind 0.1 --geno 0.1` | Same 10% sample- and variant-missingness thresholds |
| AF classification | Global AF and five superpopulation AFs | Global AF and seven regional AFs |
| Common/rare thresholds | Common: global AF >= 0.01 and AF >= 0.01 in at least two groups; rare: global AF < 0.01 | Same thresholds |

Because HGDP is the primary discovery panel, the final archived release should also include either the exact HGDP command script or a dataset-specific configuration/command log sufficient to regenerate the deposited HGDP BED files.

## ClinVar and COSMIC preprocessing

ClinVar and COSMIC preprocessing should be documented as executable scripts even when the operations consist primarily of BCFtools commands. Small scripts make the exact database release, input file, genome build, filter expression, normalization, deduplication, and output schema auditable.

The recommended additions are:

- `prepare_clinvar.sh`: records the ClinVar release date; retains autosomal biallelic SNVs on GRCh38; separates benign/likely benign and pathogenic/likely pathogenic records; applies the study's review-status criterion to pathogenic records; and produces the MANE-CDS analysis inputs.
- `prepare_cosmic.sh`: records the COSMIC release and downloaded file; retains autosomal somatic SNVs on GRCh38; normalizes and deduplicates records by `chromosome:position:REF:ALT`; and produces the BED inputs used in the study.

Do not substitute prose descriptions for the exact BCFtools expressions. The scripts should contain the commands that were actually run, including the database release identifiers.

### COSMIC access and redistribution

COSMIC source data require a registered account and are governed by the COSMIC licence. This repository should provide processing code and instructions for obtaining the specified release directly from COSMIC, but it should not publicly redistribute COSMIC source data or a reconstruction of the database unless written permission permits it. Confirm that every COSMIC-derived file planned for Zenodo is compatible with the applicable licence before deposition.

## BED schemas

The scripts use more than one BED-like schema. These formats must not be interchanged without conversion.

### Coordinate analyses

Fisher overlap tests, permutation testing, and hotspot overlap require at least three zero-based, half-open BED columns:

```text
chromosome    start    end
```

### Population-classification output

`classify_1kgp_common_rare_to_bed.py` writes six columns:

```text
chromosome    start    end    variant_ID    global_AF    number_of_groups_with_AF_ge_1pct
```

### Allele-concordance input

`compare_shared_alleles.py` expects exactly six tab-delimited columns, with no header:

```text
chromosome    start    end    consequence    REF    ALT
```

Therefore, the direct output of `classify_1kgp_common_rare_to_bed.py` is not an allele-concordance input. REF and ALT must be retained or restored from the source/annotated VCF when generating the six-column files used by `compare_shared_alleles.py`.

## Exact-substitution concordance

Example for HGDP common SNVs and COSMIC mutations:

```bash
python3 compare_shared_alleles.py \
  --germline SNPs_WG.sorted.bed \
  --somatic cosmic_WG.sorted.bed \
  --prefix SNPs_vs_COSMIC_WG
```

A shared coordinate is counted once. A locus is classified as concordant when at least one identical REF>ALT substitution is present in both datasets; all substitutions are retained at multiallelic loci.

Both input files must be coordinate sorted and must use the same genome build and chromosome naming convention.

## Chromosome-constrained permutation analysis

`run_WG_permutation.sh` randomizes unique HGDP common-SNV coordinates within their original autosomes while preserving the number of loci and preventing overlap among shuffled intervals. Exact-coordinate overlap with unique COSMIC coordinates is recalculated for each permutation.

Defaults:

- permutations: 1,000;
- initial seed: 123;
- empirical one-sided P value: `(b + 1) / (B + 1)`, where `b` is the number of permutations with overlap at least as large as observed and `B` is the total number of permutations.

Run with defaults:

```bash
bash run_WG_permutation.sh
```

Override the defaults if required:

```bash
N_PERM=10000 BASE_SEED=123 bash run_WG_permutation.sh
```

## CpG-island and microsatellite overlap

`check_hotspot_overlap.sh` prepares autosomal annotation intervals, merges overlapping intervals, intersects them with the unique common-COSMIC shared coordinates, and performs BEDTools Fisher tests.

The paths at the beginning of the script must point to:

- the common-COSMIC shared-coordinate BED;
- the GRCh38 CpG-island BED;
- the GRCh38 microsatellite BED; and
- an autosomal GRCh38 chromosome-length file.

Run:

```bash
bash check_hotspot_overlap.sh
```

## Data provenance

| Resource | Role in the study | Access |
|---|---|---|
| HGDP | Primary population discovery panel | Bergström et al. (2020); GRCh38 WGS release `hgdp_wgs.20190516`: https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/ |
| High-coverage 1KGP | Independent population replication panel | The download URL is defined in `download_1kgp.sh` |
| ClinVar | Benign/likely benign and pathogenic/likely pathogenic germline SNVs in MANE-defined CDS | Release used: `<CLINVAR_RELEASE>`; source: <https://www.ncbi.nlm.nih.gov/clinvar/> |
| COSMIC | Somatic cancer SNVs | Release used: `<COSMIC_RELEASE>`; registered access: <https://cancer.sanger.ac.uk/cosmic/> |
| Zenodo | Analysis-ready files that may legally be redistributed | **CDS** `<https://doi.org/10.5281/zenodo.21779127>` |
| Zenodo | Analysis-ready files that may legally be redistributed | **WG** `<https://doi.org/10.5281/zenodo.21779579>` |

## Citation

If you use this workflow, please cite:

> `<FULL MANUSCRIPT CITATION OR BIORXIV DOI>`


Please also cite the original databases and software tools used in the corresponding analysis.

## Contact

Thais Tavares  
GitHub: [@ticatavares](https://github.com/ticatavares)

## Licence

The source code in this repository is available under the MIT License. See the LICENSE file for details.
This license applies only to the original source code provided in this repository. Source data and database-derived files are not covered by the MIT License and remain subject to the respective terms of use and licensing conditions of HGDP, the 1000 Genomes Project, ClinVar, COSMIC, and other third-party resources.
