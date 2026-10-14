# GTEx input and validation status

Eighteen GTEx v11 bulk tissue files were supplied. Three transferred within the upload tool's 32 MiB limit: spinal cord (204 samples), ovary (193 samples), and terminal ileum (207 samples). The other 15 could not be transferred. Original compressed files are retained locally and ignored by Git. `Validation_Reports/bulk_input_manifest.csv` records availability and hashes for transferred originals.

All three available files passed gzip integrity, declared GCT row/column dimensions, unique sample IDs, six unique selected SERPIN rows, nonnegative target TPM checks, and complete subject-phenotype ID matching. All 18 gene-by-tissue summaries match independent TPM calculations. This is validation of three tissues, not complete GTEx coverage.

`SERPIN_Subsets/` contains compact GCT extracts retaining the six original target-gene rows, original sample headers, and original TPM text values. Only declared row counts were adjusted to six. Running step 06 against these extracts reproduces the full-file target summaries. Annotations are preserved unchanged under `Annotations/`; their DD spreadsheets are dictionaries, while SubjectPhenotypesDS is the actual subject metadata table.

From the repository root:

```sh
Rscript Scripts/06_GTEx_HPA_SERPINE1_TissueExpression.R --gtex-only
SERPINE1_GTEX_INPUT=Expression/GTEx/SERPIN_Subsets Rscript Scripts/06_GTEx_HPA_SERPINE1_TissueExpression.R --gtex-only
```

The unchanged GTEx calculations run independently of HPA in this CLI mode. The official HPA archive has now been downloaded and validated; use `SERPINE1_HPA_INPUT=Expression/HPA/SERPIN_Subsets` for the preserved compact target subset. The current partial summary contains only three tissues and should not be used as a complete cross-species input.

## Single-nucleus extraction

`Scripts/07B_GTEx_snRNA_SERPIN_Extraction.R` is the supplied extraction code with portable paths. It requires `GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad`, or an explicit `SERPINE1_GTEX_SNRNA_INPUT` path. `SERPINE1_GTEX_SNRNA_OUTPUT` overrides outputs. The script's categorical decoding, chunked CSR extraction, donor-level expression percentages, system mapping, exclusions, parameters, and output names are preserved. Its rhdf5 dependency is installed. An isolated synthetic HDF5 fixture passed categorical decoding, CSR extraction, empty target counts, and donor-level detection percentage checks. The original atlas has not been supplied or validated.

## Download access and retained metadata

Official bulk and single-cell pages were supplied by the author. Following publication, GTEx and HPA pages are accessible. The official GTEx API identifies `storage.googleapis.com` as the data host. That host remains blocked with HTTP 403. Saving its addition failed because this chat draft has a stale base after publication; a fresh setup chat must reconcile and save the proposed domain addition. Exact download URLs and file sizes are preserved in `Validation_Reports/official_download_manifest.json`.

Original annotation spreadsheets retain source creator/editor metadata: François Aguet, Jared Nedzel, and the profile label `Admnistrator`. These are flagged for publication review and were not silently removed from original source documents. Study sample/subject identifiers are retained for valid metadata joins and donor-level analyses; they are study identifiers, not local computer paths.

Large source files are download-only inputs. See [source data and verification scope](../../DATA_AVAILABILITY.md) for official portals and upstream steps not independently rerun here. Uploading large matrices to this chat is not required.

## Partial subject-sex join verification

Step 06B now runs with installed data.table. The portable
SERPINE1_GTEX_SUBJECT_INPUT override selects subject metadata independently
of the expression input folder. With the three unchanged compact GCT extracts,
all 3,624 gene/sample rows match sex metadata; all 30 tissue/gene/sex summary
medians exactly match the supplied selected GTEx table. This verifies three
tissues, not all original tissue matrices. Large data remain download-only.
