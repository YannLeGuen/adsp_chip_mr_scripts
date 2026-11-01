# ADSP CHIP and Alzheimer disease analyses
Creative Commons NonCommercial (CC BY-NC)

This repository identifies CHIP carriers in ADSP whole-genome sequencing (WGS)
and whole-exome sequencing (WES), fits CHIP–Alzheimer disease association models,
and performs bidirectional Mendelian randomization (MR) using CHIP and Alzheimer
disease GWAS summary statistics.

Place the repository at `/main`. All default input and output paths use that
directory. Each analysis has one fixed implementation and produces auditable
sample counts, variant filters, and analysis settings.

## Scripts

| Script | Function |
| --- | --- |
| `scripts/chip/select_chip_variants.py` | Select candidate driver variants using the Bouzid reference and annotated chromosome, position, REF, ALT, gene, and protein consequence |
| `scripts/chip/call_chip_status.py` | Extract allelic read depths, apply WGS/WES filters, and assign individual CHIP status |
| `scripts/associate_chip_ad.py` | Fit age-adjusted pooled and APOE-stratified association models |
| `scripts/mr/prepare_mr_inputs.py` | Prepare matched GWAS data with independently selected instruments for each MR direction |
| `scripts/mr/run_mr.R` | Harmonize alleles and run unpruned and 1 Mb-filtered MR analyses, with optional scatterplots |

## Installation

Requirements are Python 3.10 or later, the packages in `requirements.txt`, R with
TwoSampleMR, and bcftools for indexed VCF/BCF access.

```bash
python3 -m venv /main/.venv
/main/.venv/bin/python -m pip install -r /main/requirements.txt

mkdir -p /main/software/R
export R_LIBS_USER=/main/software/R
Rscript -e 'install.packages("TwoSampleMR", repos=c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))'
```

Keep `R_LIBS_USER` set when running R with this installation. Install bcftools
using the package manager for the analysis environment. Compressed VCFs must
be BGZF-compressed and indexed for indexed access. `--reader scan` supports
sequential reading of plain or gzip VCFs.

## Required data

| Input | Path |
| --- | --- |
| Formatted Bouzid driver reference | `/main/reference/bouzid_chip_variants.xlsx` |
| WGS annotations | `/main/data/adsp_wgs/variant_annotations.tsv` |
| WES annotations | `/main/data/adsp_wes/variant_annotations.tsv` |
| Cohort VCF/BCF shards | `/main/data/adsp_wgs/vcf/` and `/main/data/adsp_wes/vcf/` |
| Eligible blood samples | `/main/data/adsp_wgs/samples.tsv` and `/main/data/adsp_wes/samples.tsv` |
| Harmonized phenotypes | `/main/data/adsp_wgs/phenotypes.tsv` and `/main/data/adsp_wes/phenotypes.tsv` |
| CHIP GWAS | `/main/data/mr/chip_gwas.csv` |
| Alzheimer disease GWAS | `/main/data/mr/ad_gwas.tsv.gz` |

### Driver annotations and allele depths

The reference workbook format is specified in [reference/README.md](reference/README.md).
Annotation tables require `SNP`, `SYMBOL`, `Consequence`, `HGVSp`, and
`Amino_acids`; `Protein_position` is optional. `SNP` must encode
`CHROM:POS:REF:ALT`, with one alternate allele per row. Supply one aligned
variant/transcript annotation per row. A variant qualifies when a supplied
transcript annotation matches a driver rule. Duplicate variant hits are merged.

Annotations and VCFs must use the same reference assembly and normalized allele
representation. VCF `FORMAT/AD` must contain REF read depth followed by one
depth for each ALT allele. Depth is `sum(AD)`, and the VAF for each ALT is its
read depth divided by that total. Other ALT reads remain in the denominator at
multiallelic sites. Variant matching always uses chromosome, position, REF,
and ALT. `candidate_regions.tsv` uses 1-based positions.

### Sample manifest and QC population

Both sample manifests must contain:

```text
sample_id	participant_id	body_site	technical_replicate
```

Supply one sample per unique participant, with `body_site=Blood` and
`technical_replicate=0`. Sample selection must be independent of CHIP and AD
status. The manifest includes all eligible blood participants in the sequencing
cohort, before exclusions for clinical status or missing regression covariates.
These participants define the variant-QC population and the WES prevalence
denominator.

VCF sample names must match `sample_id`. Shards must have the same sample set
and disjoint variant records; sample-column order may differ. Every candidate
allele must be present in the supplied shards. Missing samples, overlapping
records, and absent candidate alleles produce errors.

### CHIP definition

Variant-level statistics are computed among manifest participants with at least
one read supporting the relevant ALT allele.

| Criterion | WGS | WES |
| --- | --- | --- |
| Driver definition | Match to Bouzid reference | Match to Bouzid reference |
| Median carrier VAF | <0.40 | <0.40 |
| Site FILTER | PASS | No additional site-FILTER restriction |
| Site read support | Median depth × median VAF >=3 | Individual depth requirement below |
| Carrier frequency | At most 36 participants | Strictly <0.1% of manifest participants |
| Individual VAF | >0.02 | >0.02 |
| Individual allelic depth | >15 | >30 |

The output `CHIP` column is 1 when a participant has a qualifying call and 0
when no qualifying call is detected. It is NA when no positive call exists
and a passing site has missing or malformed allele depths. Status 0 denotes
no detected qualifying call; it does not establish the absence of clones below
the assay's detection limit.

## Association models

Phenotype files require the following columns. They may contain the clinically
annotated subset of the manifest; unmatched rows are counted in the analysis
summary. Sample-to-participant mappings must agree with the CHIP status table.

| Column | Coding |
| --- | --- |
| `sample_id`, `participant_id` | Sample and participant identifiers |
| `AD` | 1 = AD case; 0 = control |
| `sex` | 1 = female; 0 = male |
| `age` | Age at onset for cases; age at last examination for controls; `90+` becomes 90 |
| `APOE` | 22, 23, 24, 33, 34, or 44; E-prefixed forms are accepted |
| `PC1`–`PC4` | WGS genetic ancestry components |
| `PC1`, `PC2` | WES European genetic ancestry components |
| `ancestry` | Required for WES; EUR participants are retained |

Missing status or model covariates are excluded. The pooled and stratified
models use the same complete-case cohort. APOE epsilon2 and epsilon4 dosages
are calculated from the genotype.

```text
WGS: AD ~ CHIP + sex + age + APOE2 + APOE4 + PC1 + PC2 + PC3 + PC4
WES: AD ~ CHIP + sex + age + APOE2 + APOE4 + PC1 + PC2
```

Models use a binomial GLM with a logit link. APOE-stratified analyses report
E23, E33, E34, and E44, with APOE dosage terms omitted within strata. Outputs
include odds ratios, two-sided Wald 95% confidence intervals, p-values,
participant counts, and demographics. Nonestimable models receive an explicit
status instead of a finite effect estimate.

## Run CHIP and association analyses

```bash
/main/.venv/bin/python /main/scripts/chip/select_chip_variants.py --cohort WGS
/main/.venv/bin/python /main/scripts/chip/call_chip_status.py --cohort WGS
/main/.venv/bin/python /main/scripts/associate_chip_ad.py --cohort WGS

/main/.venv/bin/python /main/scripts/chip/select_chip_variants.py --cohort WES
/main/.venv/bin/python /main/scripts/chip/call_chip_status.py --cohort WES
/main/.venv/bin/python /main/scripts/associate_chip_ad.py --cohort WES
```

Each cohort's output directory contains `candidate_variants.tsv`,
`candidate_regions.tsv`, `variant_qc.tsv`, `allele_observations.tsv`,
`status.tsv`, and `run_settings.json`. The `association/` subdirectory contains
`chip_ad_association.tsv`, `demographics.tsv`, and `analysis_counts.json`.

## Bidirectional MR

The CHIP GWAS CSV requires `SNP`, `Beta_unadj`, `SE_unadj`, `P_unadj`, and `AAF`.
`SNP` encodes `CHROM:POS:REF:ALT`; beta and allele frequency refer to ALT.

The AD GWAS TSV requires `chromosome`, `base_pair_location`, `effect_allele`,
`other_allele`, `beta`, `standard_error`, `p_value`, and
`effect_allele_frequency`.

Both inputs must use the same genome assembly, normalized forward-strand
alleles, and log-odds betas. The default assembly is GRCh38. Build declarations
are checked for agreement; the scripts do not perform coordinate liftover.
Reverse AD → CHIP analysis requires genome-wide CHIP outcome statistics,
including nonsignificant variants. Confirm this property before setting
`--chip-genome-wide`.

```bash
/main/.venv/bin/python /main/scripts/mr/prepare_mr_inputs.py --chip-genome-wide
Rscript /main/scripts/mr/run_mr.R
```

For each direction, instruments are selected using exposure p <5e-8, matched
to available outcomes by coordinate and allele pair, and analyzed both
unpruned and with greedy 1 Mb thinning. Thinning retains variants in ascending
exposure p-value order and excludes a candidate within <=1,000,000 bp of any
retained variant on the same chromosome. Outcome p-values do not select
instruments. Distance thinning does not impose an LD r-squared threshold.

`TwoSampleMR::harmonise_data(action=1)` aligns effects under the required
forward-strand input convention. Retained harmonized data are saved and used
for both estimates and plots. IVW is reported with weighted-median and MR-Egger
estimates when enough instruments are available. Single-instrument analyses
use the Wald ratio.

To add p <5e-6 sensitivity analyses, use `--include-relaxed` on both commands:

```bash
/main/.venv/bin/python /main/scripts/mr/prepare_mr_inputs.py --chip-genome-wide --include-relaxed
Rscript /main/scripts/mr/run_mr.R --include-relaxed
```

Add `--plots` to the R command for scatterplots. Results under
`/main/results/mr/` include `mr_results.tsv`, `ivw_summary.tsv`,
`analysis_summary.tsv`, `instrument_retention.tsv`, `harmonised/`,
`run_settings.tsv`, and `session_info.txt`. Instrument counts refer to retained
harmonized variants. Approximate instrument F statistics are `(beta/SE)^2`.

## Synthetic tests

```bash
/main/.venv/bin/python /main/tests/smoke_test.py
```

The tests generate artificial reference tables, VCFs, phenotypes, and GWAS
records in temporary directories. They check allele matching, threshold
boundaries, missing data, participant identity, and exposure-specific MR input
selection.

## References

- [Bouzid et al., Nature Medicine 2023](https://pubmed.ncbi.nlm.nih.gov/37322115/)
- [Kessler et al., Nature 2022](https://pubmed.ncbi.nlm.nih.gov/36450978/)
- [Bellenguez et al., Nature Genetics 2022](https://pubmed.ncbi.nlm.nih.gov/35379992/)
- [TwoSampleMR documentation](https://mrcieu.github.io/TwoSampleMR/)
- [bcftools documentation](https://samtools.github.io/bcftools/bcftools.html)
- [statsmodels GLM documentation](https://www.statsmodels.org/stable/glm.html)

Creative Commons NonCommercial (CC BY-NC)
