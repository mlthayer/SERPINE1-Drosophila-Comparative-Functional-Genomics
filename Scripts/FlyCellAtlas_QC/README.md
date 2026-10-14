# FlyCellAtlas supporting checks

These four helpers were adapted from the author's uploaded extraction, inspection,
metadata-summary, and matrix-verification documents. Run from the repository root.
The primary analysis remains `Scripts/07_FlyCellAtlas_CellExpression.R`.
`SERPINE1_FCA_INPUT` selects a downloaded input folder; default:
`Expression/FlyCellAtlas`. Source datasets are available through
<https://flycellatlas.org/>. Preserve the exact release and matrix variants used
by the study; similarly named H5AD files need not contain equivalent matrices.

- `07_FCA_Antenna_Extraction.R`: extracts six genes from the older flattened
  `raw.X` CSC matrix in `s_fca_biohub_antenna_10x.h5ad`. Keeps the supplied
  expression values and detection calculations. Now rejects non-CSC storage.
- `07_FCA_H5AD_Inspection.R`: inspects the first discovered file using
  zellkonverter, SingleCellExperiment and SummarizedExperiment. The first file
  is not necessarily antenna. Requires those additional Bioconductor packages;
  their availability and this helper's runtime have not been verified.
- `07_FCA_Verify_X_Matrices.R`: checks seven named newer files. Corrected the
  orientation diagnostic: for AnnData cells-by-genes storage, cell-length
  pointers correspond to CSR and gene-length pointers to CSC. Encoding and
  shape must agree; pointer length alone is inconclusive for square matrices.
  This helper reads the entire stored data vector and may require substantial RAM.
- `07_FCA_Antenna_Metadata_Summary.R`: continuation after antenna extraction.
  Requires `sce` from the same file, `serpin_expression`, `n_cells`, and
  `fca_raw_cell_barcodes` from the source obs index in raw.X row order. Rejects
  unverified barcode order before the original index-based join. Retains the
  supplied gene/sex/annotation grouping and all statistics. It does not group
  by sample despite the prose in the supplied document; no sample grouping was
  invented. The main step 07 remains the multi-tissue analysis.

All helpers pass R parsing. Real-data extraction, metadata correspondence and
matrix-scale checks remain unverified until the original H5AD files are available.
The attached documents included prose; only their R code was incorporated.
Source annotations named the matrix "raw", but this is not proof of raw UMI
counts. The primary workflow treats stored values as normalized/transformed
and uses positive-expression percentages for comparisons.
