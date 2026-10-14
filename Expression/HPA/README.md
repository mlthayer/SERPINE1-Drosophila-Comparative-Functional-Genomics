# HPA tissue transcript data

Downloaded the official archive from `https://www.proteinatlas.org/download/tsv/transcript_rna_tissue.tsv.zip` on 2026-10-13. The archive contains 205,541 transcript rows and 186 TPM sample columns. Its ZIP stream passed integrity validation. All six target SERPIN genes were recovered across 24 transcripts. This is the current official download; equivalence to an earlier author-used HPA release is not claimed. Exact archive URL, byte count, and SHA256 are recorded in `download_validation.json`.

`SERPIN_Subsets/transcript_rna_tissue.tsv.zip` retains the original target transcript rows, original header, and every sample column unchanged. This compact subset allows the existing workflow to run without loading unrelated transcript rows. Original gene-level aggregation and tissue statistics are preserved.

From the repository root:

```sh
SERPINE1_GTEX_INPUT=Expression/GTEx/SERPIN_Subsets SERPINE1_HPA_INPUT=Expression/HPA/SERPIN_Subsets Rscript Scripts/06_GTEx_HPA_SERPINE1_TissueExpression.R
```

Step 06 completed on these inputs. Independent checks matched all 1,116 gene/sample transcript sums and all 240 gene/tissue median and mean summaries. GTEx coverage remains limited to three tissues, so the resulting cross-resource comparisons are partial and are not a complete publication result.

Large source files are download-only inputs. See [source data and verification scope](../../DATA_AVAILABILITY.md) for official portals and upstream steps not independently rerun here. Uploading large matrices to this chat is not required.

## Sample-sex join validation

The uploaded rna_tissue_hpa_samples.tsv.zip provides 186 unique sample records.
Step 06C joined all 1,116 gene/sample TPM records with no missing sex or tissue
mismatch. It generated 384 tissue/gene/sex summaries; every median and sample
count in the 318-row selected HPA tissue table matches after case normalization
of tissue labels. The additional rows correspond to broader tissue coverage
in the full input, not inserted data. sex_join_validation.json records checks.
