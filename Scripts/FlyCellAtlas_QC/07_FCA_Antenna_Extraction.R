# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Supporting FCA quality-control helper; run from the repository root.
# ============================================================
# 07_FCA_Antenna_Extraction.R
#
# Extract six SERPIN genes from the full raw FCA matrix
# without loading the entire matrix into memory
# ============================================================


library(hdf5r)
library(dplyr)
library(tidyr)
library(tibble)

# ------------------------------------------------------------
# Open antenna H5AD
# ------------------------------------------------------------

flycell_dir <- Sys.getenv("SERPINE1_FCA_INPUT", unset = "Expression/FlyCellAtlas")

antenna_file <- file.path(
  flycell_dir,
  "s_fca_biohub_antenna_10x.h5ad"
)

h5 <- hdf5r::H5File$new(
  antenna_file,
  mode = "r"
)

# ------------------------------------------------------------
# 1. Target genes
# ------------------------------------------------------------

selected_serpins <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)


# ------------------------------------------------------------
# 2. Confirm raw matrix dimensions
# ------------------------------------------------------------

raw_shape <- h5[["raw.X"]]$
  attr_open("h5sparse_shape")$
  read()

n_cells <- raw_shape[1]
n_genes <- raw_shape[2]

cat("\nRaw matrix dimensions:\n")
cat("Cells:", n_cells, "\n")
cat("Genes:", n_genes, "\n")


# ------------------------------------------------------------
# 3. Check sparse matrix format
# ------------------------------------------------------------

cat("\nSparse matrix format:\n")

raw_format <- h5[["raw.X"]]$
  attr_open("h5sparse_format")$
  read()

print(raw_format)
if (!identical(tolower(as.character(raw_format)), "csc")) stop("This helper requires CSC raw.X storage.")


# ------------------------------------------------------------
# 4. Read raw gene annotation
# ------------------------------------------------------------

raw_var <- h5[["raw.var"]][]


# Gene names are stored in the index column
gene_names <- as.character(
  raw_var$index
)


# ------------------------------------------------------------
# 5. Locate target SERPIN columns
# ------------------------------------------------------------

gene_positions <- match(
  selected_serpins,
  gene_names
)


serpin_positions <- tibble(
  Gene = selected_serpins,
  Gene_Column = gene_positions,
  Present = !is.na(gene_positions)
)


cat("\n========================================\n")
cat("SERPIN GENE POSITIONS\n")
cat("========================================\n")

print(
  serpin_positions,
  n = Inf
)


if (any(is.na(gene_positions))) {

  stop(
    "One or more selected SERPIN genes are missing from raw.var."
  )
}


# ------------------------------------------------------------
# 6. Read CSC pointer vector
# ------------------------------------------------------------

indptr <- h5[["raw.X/indptr"]][]


cat("\nLength of indptr:\n")
print(length(indptr))

cat("\nExpected length:\n")
print(n_genes + 1)


if (length(indptr) != n_genes + 1) {

  stop(
    "indptr length does not match a CSC matrix with genes as columns."
  )
}


# ------------------------------------------------------------
# 7. Function to extract one gene column
# ------------------------------------------------------------

extract_gene <- function(
    h5_object,
    gene_position,
    gene_name,
    n_cells,
    indptr
) {

  # indptr values are zero-based offsets into data/indices
  start_zero <- indptr[gene_position]
  end_zero   <- indptr[gene_position + 1] - 1

  # Initialize all cells as zero expression
  expression <- numeric(n_cells)

  # Only read if this gene has detected values
  if (end_zero >= start_zero) {

    read_positions <- seq.int(
      start_zero + 1,
      end_zero + 1
    )

    gene_values <- h5_object[["raw.X/data"]][
      read_positions
    ]

    cell_indices_zero <- h5_object[["raw.X/indices"]][
      read_positions
    ]

    # Convert zero-based cell indices to one-based R indices
    cell_indices <- cell_indices_zero + 1

    expression[cell_indices] <- gene_values
  }

  tibble::tibble(
    Cell_Index = seq_len(n_cells),
    Gene = gene_name,
    Expression = expression
  )
}

# ------------------------------------------------------------
# 8. Extract all six genes
# ------------------------------------------------------------

serpin_expression <- dplyr::bind_rows(
  lapply(
    seq_along(selected_serpins),
    function(i) {

      extract_gene(
        h5_object = h5,
        gene_position = gene_positions[i],
        gene_name = selected_serpins[i],
        n_cells = n_cells,
        indptr = indptr
      )
    }
  )
)


cat("\n========================================\n")
cat("SERPIN EXPRESSION TABLE\n")
cat("========================================\n")

print(
  head(
    serpin_expression,
    20
  )
)

# ------------------------------------------------------------
# 9. Verify detection counts
# ------------------------------------------------------------

detection_check <- serpin_expression %>%
  group_by(Gene) %>%
  summarize(
    N_Cells_Total = n(),
    N_Cells_Expressing = sum(
      Expression > 0,
      na.rm = TRUE
    ),
    Percent_Expressing =
      100 *
      N_Cells_Expressing /
      N_Cells_Total,
    Min_Expression = min(
      Expression,
      na.rm = TRUE
    ),
    Max_Expression = max(
      Expression,
      na.rm = TRUE
    ),
    Mean_Expression = mean(
      Expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


cat("\n========================================\n")
cat("DETECTION CHECK\n")
cat("========================================\n")

print(
  detection_check,
  n = Inf
)

# ------------------------------------------------------------
# Close H5AD
# ------------------------------------------------------------

h5$close_all()
rm(h5)
gc()