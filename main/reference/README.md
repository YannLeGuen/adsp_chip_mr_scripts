# Bouzid CHIP reference

Place the curated, formatted Bouzid reference workbook here as:

`/main/reference/bouzid_chip_variants.xlsx`

This curated driver catalog is a required input to variant selection.

The first worksheet must contain these columns:

| Column | Meaning |
| --- | --- |
| `Gene name` | Gene symbol |
| `Variants_category` | Curated broad classes, including Frameshift, nonsense, or splice, where applicable |
| `Exact_Missense` | Comma- or semicolon-separated exact protein changes from the curated list; blank if none |
| `Limit` | Optional inclusive protein-position interval, such as `p.1-100`; blank if unrestricted |

The selector also accepts a tab-delimited equivalent supplied through
`--reference`. It validates unsupported or malformed rules instead of guessing.
Protein changes may use one- or three-letter amino-acid notation. Preserve the
curated entries and their provenance when adding the workbook to the repository.

Reference: Bouzid H, et al. *Clonal hematopoiesis is associated with protection
from Alzheimer's disease*. Nature Medicine 29, 1662–1670 (2023).
[PubMed record](https://pubmed.ncbi.nlm.nih.gov/37322115/)
