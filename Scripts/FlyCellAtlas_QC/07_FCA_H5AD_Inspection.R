# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Supporting FCA quality-control helper; run from the repository root.
# ============================================================
# 07_FCA_H5AD_Inspection.R
#
# Purpose:
# Inspect raw Fly Cell Atlas .h5ad files before re-analysis
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)
library(stringr)
library(tibble)


# ------------------------------------------------------------
# 2. File path
# ------------------------------------------------------------

flycell_dir <- Sys.getenv("SERPINE1_FCA_INPUT", unset = "Expression/FlyCellAtlas")


# ------------------------------------------------------------
# 3. Identify .h5ad files
# ------------------------------------------------------------

h5ad_files <- list.files(
  flycell_dir,
  pattern = "\\.h5ad$",
  full.names = TRUE
)


cat("\n========================================\n")
cat("FLY CELL ATLAS FILES\n")
cat("========================================\n")

print(
  basename(h5ad_files)
)

cat(
  "\nNumber of files:",
  length(h5ad_files),
  "\n"
)


# ------------------------------------------------------------
# 4. Check expected tissues
# ------------------------------------------------------------

expected_files <- c(
  "fat_body.h5ad",
  "gut.h5ad",
  "heart.h5ad",
  "leg.h5ad",
  "oenocyte.h5ad",
  "ovary.h5ad",
  "testis.h5ad",
  "trachea.h5ad"
)


missing_files <- setdiff(
  expected_files,
  basename(h5ad_files)
)


if (length(missing_files) > 0) {

  cat(
    "\nMissing expected files:\n"
  )

  print(missing_files)

} else {

  cat(
    "\nAll expected Fly Cell Atlas files found.\n"
  )
}


# ============================================================
# 5. Inspect the first file
# ============================================================

if (!length(h5ad_files)) stop("No FCA .h5ad inputs found in ", flycell_dir)
first_file <- h5ad_files[1]


cat("\n========================================\n")
cat("READING FIRST H5AD FILE\n")
cat("========================================\n")

cat(
  "\nFile:",
  basename(first_file),
  "\n"
)


sce <- readH5AD(
  first_file
)


# ------------------------------------------------------------
# 6. Basic dimensions
# ------------------------------------------------------------

cat("\nObject dimensions:\n")

cat(
  "Genes/features:",
  nrow(sce),
  "\n"
)

cat(
  "Cells/nuclei:",
  ncol(sce),
  "\n"
)


# ------------------------------------------------------------
# 7. Available assays
# ------------------------------------------------------------

cat("\nAvailable assays:\n")

print(
  assayNames(sce)
)


# ------------------------------------------------------------
# 8. Gene metadata
# ------------------------------------------------------------

cat("\n========================================\n")
cat("GENE METADATA COLUMNS\n")
cat("========================================\n")

print(
  colnames(
    rowData(sce)
  )
)


cat("\nFirst gene metadata rows:\n")

print(
  as.data.frame(
    rowData(sce)
  ) %>%
    head(10)
)


# ------------------------------------------------------------
# 9. Cell metadata
# ------------------------------------------------------------

cat("\n========================================\n")
cat("CELL METADATA COLUMNS\n")
cat("========================================\n")

print(
  colnames(
    colData(sce)
  )
)


cat("\nFirst cell metadata rows:\n")

print(
  as.data.frame(
    colData(sce)
  ) %>%
    head(10)
)


# ------------------------------------------------------------
# 10. Row names / gene identifiers
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FIRST 25 FEATURE NAMES\n")
cat("========================================\n")

print(
  head(
    rownames(sce),
    25
  )
)


# ------------------------------------------------------------
# 11. Check candidate serpins
# ------------------------------------------------------------

selected_serpins <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)


cat("\n========================================\n")
cat("SERPIN ROW-NAME CHECK\n")
cat("========================================\n")


serpin_check <- tibble(
  Gene = selected_serpins,
  Present_in_rownames =
    selected_serpins %in% rownames(sce)
)


print(
  serpin_check
)


# ------------------------------------------------------------
# 12. Search all rowData columns for serpin names
# ------------------------------------------------------------

rowdata_df <- as.data.frame(
  rowData(sce)
)


if (ncol(rowdata_df) > 0) {

  cat("\nSearching rowData for candidate serpins...\n")

  for (column_name in names(rowdata_df)) {

    values <-
      as.character(
        rowdata_df[[column_name]]
      )

    matches <-
      values[
        str_detect(
          values,
          regex(
            paste(
              selected_serpins,
              collapse = "|"
            ),
            ignore_case = TRUE
          )
        )
      ]

    if (length(matches) > 0) {

      cat(
        "\nColumn:",
        column_name,
        "\n"
      )

      print(
        unique(matches)
      )
    }
  }
}


# ------------------------------------------------------------
# 13. Inspect likely annotation columns
# ------------------------------------------------------------

metadata_names <- colnames(
  colData(sce)
)


annotation_candidates <- metadata_names[
  str_detect(
    metadata_names,
    regex(
      "cell|type|cluster|annotation|tissue|sex|gender|organ|sample",
      ignore_case = TRUE
    )
  )
]


cat("\n========================================\n")
cat("POSSIBLE ANNOTATION COLUMNS\n")
cat("========================================\n")

print(
  annotation_candidates
)


# ------------------------------------------------------------
# 14. Print unique values from likely annotation columns
# ------------------------------------------------------------

for (column_name in annotation_candidates) {

  cat(
    "\n----------------------------------------\n"
  )

  cat(
    "Column:",
    column_name,
    "\n"
  )

  values <-
    unique(
      as.character(
        colData(sce)[[column_name]]
      )
    )

  print(
    head(
      values,
      50
    )
  )
}


# ------------------------------------------------------------
# 15. Inspect assay values
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ASSAY VALUE INSPECTION\n")
cat("========================================\n")


for (assay_name in assayNames(sce)) {

  cat(
    "\nAssay:",
    assay_name,
    "\n"
  )

  assay_matrix <-
    assay(
      sce,
      assay_name
    )

  sample_values <-
    as.numeric(
      assay_matrix[
        seq_len(
          min(
            50,
            nrow(assay_matrix)
          )
        ),
        seq_len(
          min(
            50,
            ncol(assay_matrix)
          )
        )
      ]
    )

  cat(
    "Minimum:",
    min(
      sample_values,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Maximum:",
    max(
      sample_values,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Mean:",
    mean(
      sample_values,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Fraction integer-like:",
    mean(
      abs(
        sample_values -
          round(sample_values)
      ) < 1e-8,
      na.rm = TRUE
    ),
    "\n"
  )
}


# ------------------------------------------------------------
# 16. Completion
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FLY CELL ATLAS INSPECTION COMPLETE\n")
cat("========================================\n")