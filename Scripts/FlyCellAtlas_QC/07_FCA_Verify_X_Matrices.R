# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Supporting FCA quality-control helper; run from the repository root.
flycell_dir <- Sys.getenv("SERPINE1_FCA_INPUT", unset = "Expression/FlyCellAtlas")
# ============================================================
# VERIFY NEW FCA X MATRICES
# ============================================================

library(hdf5r)

new_files <- c(
  "fat_body.h5ad",
  "gut.h5ad",
  "heart.h5ad",
  "oenocyte.h5ad",
  "ovary.h5ad",
  "testis.h5ad",
  "trachea.h5ad"
)

target_genes <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)

for (f in new_files) {

  cat("\n\n========================================\n")
  cat("FILE:", f, "\n")
  cat("========================================\n")

  h5 <- hdf5r::H5File$new(
    filename = file.path(flycell_dir, f),
    mode = "r"
  )

  x <- h5[["X"]]

  # ----------------------------------------------------------
  # MATRIX ATTRIBUTES
  # ----------------------------------------------------------

  cat("\nEncoding type:\n")
  print(
    hdf5r::h5attr(
      x,
      "encoding-type"
    )
  )

  cat("\nShape:\n")

  x_shape <- hdf5r::h5attr(
    x,
    "shape"
  )

  print(x_shape)

  # ----------------------------------------------------------
  # CELLS AND GENES
  # ----------------------------------------------------------

  cell_ids <-
    h5[["obs"]][["_index"]][]

  gene_names <-
    h5[["var"]][["_index"]][]

  n_cells <- length(cell_ids)
  n_genes <- length(gene_names)

  indptr <-
    x[["indptr"]][]

  cat("\nNumber of cells:", n_cells, "\n")
  cat("Number of genes:", n_genes, "\n")
  cat("indptr length:", length(indptr), "\n")

  # ----------------------------------------------------------
  # DETERMINE SPARSE MATRIX ORIENTATION
  # ----------------------------------------------------------

  cat("\nMatrix orientation check:\n")

  if (length(indptr) == n_cells + 1) {

    cat(
      "indptr = n_cells + 1\n",
      "Pointer length is compatible with CSR cells x genes; confirm encoding-type and shape.\n"
    )

  } else if (length(indptr) == n_genes + 1) {

    cat(
      "indptr = n_genes + 1\n",
      "Pointer length is compatible with CSC cells x genes; confirm encoding-type and shape.\n"
    )

  } else {

    cat(
      "WARNING: indptr length does not match cells or genes.\n"
    )
  }

  # ----------------------------------------------------------
  # TARGET-GENE PRESENCE
  # ----------------------------------------------------------

  cat("\nTarget gene presence:\n")

  gene_check <- data.frame(
    Gene = target_genes,
    Present = target_genes %in% gene_names,
    Position = match(
      target_genes,
      gene_names
    )
  )

  print(
    gene_check,
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # EXPRESSION VALUES
  # ----------------------------------------------------------

  x_data <-
    x[["data"]][]

  cat("\nExpression-value summary:\n")

  print(
    summary(x_data)
  )

  cat(
    "\nMinimum:",
    min(
      x_data,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Maximum:",
    max(
      x_data,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Any negative values:",
    any(
      x_data < 0,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Fraction approximately integer:",
    mean(
      abs(
        x_data -
          round(x_data)
      ) < 1e-8,
      na.rm = TRUE
    ),
    "\n"
  )

  # ----------------------------------------------------------
  # SEX CATEGORIES
  # ----------------------------------------------------------

  cat("\nSex categories:\n")

  obs_categories <-
    h5[["obs"]][["__categories"]]

  if (
    "sex" %in%
      names(obs_categories)
  ) {

    print(
      obs_categories[["sex"]][]
    )
  }

  # ----------------------------------------------------------
  # ANNOTATION CATEGORIES
  # ----------------------------------------------------------

  cat("\nAnnotation categories:\n")

  if (
    "annotation" %in%
      names(obs_categories)
  ) {

    print(
      obs_categories[["annotation"]][]
    )
  }

  h5$close_all()
}
