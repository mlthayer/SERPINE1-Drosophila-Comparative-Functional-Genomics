# Requires the extraction objects and sce read from the SAME antenna file.
# Supply fca_raw_cell_barcodes in raw.X row order from that file's obs index.
if (!exists("fca_raw_cell_barcodes")) {
  stop("Supply raw.X cell barcodes to verify metadata ordering before joining.")
}
# ------------------------------------------------------------
# 10. Attach FCA cell metadata
# ------------------------------------------------------------

cell_metadata <- as.data.frame(
  SummarizedExperiment::colData(sce)
) %>%
  tibble::rownames_to_column("Cell_Barcode") %>%
  dplyr::mutate(
    Cell_Index = dplyr::row_number()
  )

cat("\nNumber of metadata rows:\n")
print(nrow(cell_metadata))

cat("\nNumber of raw.X cells:\n")
print(n_cells)

if (nrow(cell_metadata) != n_cells) {
  stop("Cell metadata does not match raw.X cell count.")
}

# Join expression and metadata.

if (anyDuplicated(cell_metadata$Cell_Barcode) ||
    !identical(as.character(cell_metadata$Cell_Barcode), as.character(fca_raw_cell_barcodes))) {
  stop("Metadata barcodes do not match raw.X row order; do not join by Cell_Index.")
}
serpin_cell_data <- serpin_expression %>%
  dplyr::left_join(
    cell_metadata,
    by = "Cell_Index"
  )

# Summarize by gene, sex and annotation.

antenna_celltype_summary <- serpin_cell_data %>%
  dplyr::group_by(
    Gene,
    sex,
    annotation
  ) %>%
  dplyr::summarise(
    N_Cells = dplyr::n(),
    N_Expressing = sum(Expression > 0, na.rm = TRUE),
    Percent_Expressing = 100 * N_Expressing / N_Cells,
    Mean_Expression_All = mean(Expression, na.rm = TRUE),
    Median_Expression_All = median(Expression, na.rm = TRUE),
    Mean_Expression_Positive = ifelse(
      N_Expressing > 0,
      mean(Expression[Expression > 0], na.rm = TRUE),
      0
    ),
    Max_Expression = max(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    Gene,
    dplyr::desc(Percent_Expressing),
    dplyr::desc(Mean_Expression_Positive)
  )

print(
  antenna_celltype_summary,
  n = 100
)
