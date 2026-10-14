# Source data and verification scope

Large public expression matrices are download-only inputs. Users obtain them
from the official sources below and keep them locally; uploading them to this
chat or committing them to the repository is not required. Selected processed
inputs and reference outputs are retained for the documented downstream runs.

| Source | Official download location | Local input / workflow |
|---|---|---|
| Adult GTEx v11 bulk gene TPM | https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression | Tissue `.gct`/`.gct.gz` files in `Expression/GTEx/`; step 06B uses subject sex metadata from `Expression/GTEx/Annotations/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt` |
| HPA transcript TPM and sample metadata | https://www.proteinatlas.org/about/download | `transcript_rna_tissue.tsv.zip` and `rna_tissue_hpa_samples.tsv` (or ZIP) in `Expression/HPA/`; steps 06 and 06C |
| GTEx single-nucleus atlas | https://gtexportal.org/home/downloads/adult-gtex/single_cell | `GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad` in `Expression/GTEx/`; step 07B |
| FlyCellAtlas 10x datasets | https://flycellatlas.org/ | Fifteen H5AD inputs listed in `Expression/FlyCellAtlas/README.md`; step 07 |

Input environment variables and interactive selectors are documented in the
README and source-specific guides. Download the correct release and matrix
variant; filenames and portal links alone do not prove historical equivalence.
Record the release, exact file URL, retrieval date and checksum when source
files are obtained. Where these are unavailable, provenance remains incomplete.

## Verified in this environment

- Bulk comparisons and multi-metric consensus were rebuilt from the author's
  selected latest sex-specific tissue summaries using steps 08B and 08C.
- Step 09 reproduced all sixteen supplied aging output tables from processed
  fly and human inputs; human sample-level aging summaries were also checked.
- Steps 10 and 11 ran with validated 25-position RCL identities, rebuilt bulk
  consensus, supplied cellular correlations and reproduced aging evidence.

HPA transcript aggregation and sample-sex joins were verified against the
retained compact target-transcript extract and uploaded sample metadata.

## Source checks and remaining independent-verification limits

- GTEx bulk TPM-to-subject-sex joins were additionally verified for the three
  available tissues using unchanged six-gene extracts: all 3,624 gene/sample
  rows matched sex metadata and all 30 resulting tissue/gene/sex summary medians
  match the selected tissue inputs. Full coverage remains unverified because
  the other original matrices were unavailable.
- HPA sample-to-sex joins are now verified: all 1,116 gene/sample rows joined
  to the uploaded 186-sample metadata with no missing sex or tissue mismatch.
  All 318 selected tissue/gene/sex medians and counts reproduce.
- FlyCellAtlas H5AD extraction and original cell metadata correspondence:
  supplied processed tables passed consistency and aggregation checks.
- GTEx single-nucleus original H5AD expression extraction remains unverified.
  The supplied donor-system, donor-tissue, sex-summary and cell-metadata tables
  now pass aggregation, count and coverage checks. Step 08D independently
  reproduces all 72 cellular correlation output rows from these processed inputs.
- Original AFCA cell-level extraction: the large single-nucleus export could
  not be transferred; processed aging input/output agreement was verified.

These are declared verification limits, not pending requirements for the
author to upload large files. A user with local source data can rerun the
upstream scripts. Successful downstream runs do not validate missing upstream
joins, dataset releases or biological annotations.
