# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 07_FlyCellAtlas_CellExpression.R
# Run from the repository root.
# ============================================================
# 07C_FlyCellAtlas_MultiTissue_SERPIN_Reanalysis.R
#
# Multi-tissue Fly Cell Atlas re-analysis of six serpin genes
#
# Target genes:
#   Spn42Da
#   Spn42Db
#   Spn42Dc
#   Spn42Dd
#   Spn28B
#   Spn28F
#
# SUPPORTED FCA H5AD STRUCTURES
#
# 1. Older flattened FCA format
#       raw.X
#       raw.var
#       obs = compound H5D dataset
#       sparse format = CSC
#
# 2. Newer BioHub nested format
#       raw/X
#       raw/var
#       obs = H5Group
#       sparse format = CSC
#
# 3. Testis dataset (BioHub 10x raw version unavailable)
#       X
#       var/_index
#       obs = H5Group
#       sparse format = CSR
#
# IMPORTANT
#
# - Expression detection is defined as Expression > 0.
#
# - The stored expression values are non-integer and are
#   treated as normalized/transformed expression values,
#   not raw UMI counts.
#
# - Expression-value scales differ between FCA file groups.
#   Therefore, Percent_Expressing is the primary metric used
#   for cross-tissue and cross-file comparisons.
#
# - Mean expression values are retained in output tables,
#   but should not be interpreted as directly comparable
#   between different FCA expression-matrix formats.
#
# - Pooled cell-type summaries are descriptive.
#
# - A minimum of 20 pooled cells is used for plot-ready
#   descriptive summaries.
#
# ============================================================



# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(hdf5r)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)



# ------------------------------------------------------------
# 1. Working directory
# ------------------------------------------------------------

# Download source H5AD datasets from https://flycellatlas.org/ (10x data).
# See Expression/FlyCellAtlas/README.md for the exact input manifest.
flycell_dir <- Sys.getenv("SERPINE1_FCA_INPUT", unset = "Expression/FlyCellAtlas")


output_dir <- file.path(
  flycell_dir,
  "Reanalysis_Output"
)


table_dir <- file.path(
  output_dir,
  "Tables"
)


figure_dir <- file.path(
  output_dir,
  "Figures"
)


dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


dir.create(
  table_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)



# ------------------------------------------------------------
# 2. Define and validate the exact 15 H5AD input files
#
# Using an explicit manifest prevents deleted/backup/legacy H5AD
# files from being included accidentally.
# ------------------------------------------------------------

expected_h5ad_basenames <- c(
  "s_fca_biohub_antenna_10x.h5ad",
  "s_fca_biohub_body_wall_10x.h5ad",
  "s_fca_biohub_fat_body_10x.h5ad",
  "s_fca_biohub_gut_10x.h5ad",
  "s_fca_biohub_haltere_10x.h5ad",
  "s_fca_biohub_heart_10x.h5ad",
  "s_fca_biohub_leg_10x.h5ad",
  "s_fca_biohub_male_reproductive_glands_10x.h5ad",
  "s_fca_biohub_malpighian_tubule_10x.h5ad",
  "s_fca_biohub_oenocyte_10x.h5ad",
  "s_fca_biohub_ovary_10x.h5ad",
  "s_fca_biohub_proboscis_and_maxillary_palps_10x.h5ad",
  "testis.h5ad",
  "s_fca_biohub_trachea_10x.h5ad",
  "s_fca_biohub_wing_10x.h5ad"
)

h5ad_files <- file.path(
  flycell_dir,
  expected_h5ad_basenames
)

missing_input_files <-
  expected_h5ad_basenames[
    !file.exists(h5ad_files)
  ]

if (
  length(missing_input_files) > 0
) {
  stop(
    paste0(
      "Download the required FCA datasets from https://flycellatlas.org/ and place them in ",
      flycell_dir,
      " (or set SERPINE1_FCA_INPUT). Missing required H5AD file(s): ",
      paste(
        missing_input_files,
        collapse = ", "
      )
    )
  )
}

cat("\n========================================\n")
cat("FCA H5AD FILES VALIDATED\n")
cat("========================================\n")

print(
  basename(h5ad_files)
)

cat(
  "\nNumber of validated H5AD files:",
  length(h5ad_files),
  "\n"
)

if (
  length(h5ad_files) != 15L
) {
  stop(
    "The validated FCA input manifest must contain exactly 15 H5AD files."
  )
}


# ------------------------------------------------------------
# 3. Target serpin genes
# ------------------------------------------------------------

selected_serpins <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)



# ============================================================
# HELPER FUNCTIONS
# ============================================================



# ------------------------------------------------------------
# 4. Safely retrieve metadata column
# ------------------------------------------------------------

get_metadata_column <- function(
    metadata,
    column_name,
    default = NA_character_
) {

  if (
    column_name %in%
      colnames(metadata)
  ) {

    return(
      as.character(
        metadata[[column_name]]
      )
    )
  }


  rep(
    default,
    nrow(metadata)
  )
}



# ------------------------------------------------------------
# 5. Decode zero-based AnnData categorical codes
# ------------------------------------------------------------

decode_annData_category <- function(
    codes,
    categories
) {

  codes <- as.integer(
    codes
  )


  categories <- as.character(
    categories
  )


  result <- rep(
    NA_character_,
    length(codes)
  )


  valid <-
    !is.na(codes) &
    codes >= 0 &
    codes < length(categories)


  result[valid] <-
    categories[
      codes[valid] + 1L
    ]


  result
}



# ------------------------------------------------------------
# 6. Read FCA metadata from original/BioHub files
#
# Supports:
#   A. H5D compound obs
#   B. H5Group obs
# ------------------------------------------------------------

read_fca_obs <- function(
    h5_object
) {

  obs_object <-
    h5_object[["obs"]]


  # ==========================================================
  # A. Older flattened compound obs
  # ==========================================================

  if (
    inherits(
      obs_object,
      "H5D"
    )
  ) {

    obs_df <- as.data.frame(
      obs_object[]
    )


    categorical_fields <- c(
      "fca_id",
      "sample_id",
      "sex",
      "tissue",
      "annotation",
      "annotation_broad",
      "annotation_broad_extrapolated",
      "R_annotation",
      "R_annotation_broad"
    )


    if (
      !"uns" %in%
        names(h5_object)
    ) {

      stop(
        "Older H5AD obs detected but uns group is missing."
      )
    }


    uns_names <- names(
      h5_object[["uns"]]
    )


    for (
      field in categorical_fields
    ) {

      category_name <- paste0(
        field,
        "_categories"
      )


      if (
        field %in%
          colnames(obs_df) &&
        category_name %in%
          uns_names
      ) {

        categories <-
          h5_object[[
            paste0(
              "uns/",
              category_name
            )
          ]][]


        obs_df[[field]] <-
          decode_annData_category(
            codes =
              obs_df[[field]],
            categories =
              categories
          )
      }
    }


    return(
      tibble::as_tibble(
        obs_df
      )
    )
  }


  # ==========================================================
  # B. Newer H5Group obs
  # ==========================================================

  if (
    inherits(
      obs_object,
      "H5Group"
    )
  ) {

    obs_names <- names(
      obs_object
    )


    wanted_fields <- c(
      "index",
      "_index",
      "age",
      "fca_id",
      "sample_id",
      "sex",
      "tissue",
      "annotation",
      "annotation_broad",
      "annotation_broad_extrapolated",
      "R_annotation",
      "R_annotation_broad",
      "technology"
    )


    wanted_fields <- intersect(
      wanted_fields,
      obs_names
    )


    if (
      length(wanted_fields) == 0
    ) {

      stop(
        "No usable metadata fields were found in obs."
      )
    }


    obs_list <- vector(
      mode = "list",
      length = length(
        wanted_fields
      )
    )


    names(
      obs_list
    ) <- wanted_fields


    category_group_exists <-
      "__categories" %in%
      obs_names


    category_names <- character(0)


    if (
      category_group_exists
    ) {

      category_names <- names(
        obs_object[["__categories"]]
      )
    }


    for (
      field in wanted_fields
    ) {

      raw_values <-
        obs_object[[field]][]


      if (
        category_group_exists &&
        field %in%
          category_names
      ) {

        categories <-
          obs_object[[
            "__categories"
          ]][[
            field
          ]][]


        obs_list[[field]] <-
          decode_annData_category(
            codes =
              raw_values,
            categories =
              categories
          )

      } else {

        obs_list[[field]] <-
          raw_values
      }
    }


    return(
      tibble::as_tibble(
        obs_list
      )
    )
  }


  stop(
    paste(
      "Unsupported obs object type:",
      paste(
        class(obs_object),
        collapse = ", "
      )
    )
  )
}



# ------------------------------------------------------------
# 7. Read metadata from new X-only tissue datasets
#
# These files contain:
#   _index
#   annotation
#   sex
#   technology
#
# Tissue is derived from the filename.
# ------------------------------------------------------------

read_x_only_obs <- function(
    h5_object,
    file_path
) {

  obs <-
    h5_object[["obs"]]


  if (
    !inherits(
      obs,
      "H5Group"
    )
  ) {

    stop(
      paste(
        "X-only file does not contain H5Group obs:",
        basename(file_path)
      )
    )
  }


  obs_names <- names(
    obs
  )


  if (
    !"_index" %in%
      obs_names
  ) {

    stop(
      paste(
        "X-only dataset lacks obs/_index:",
        basename(file_path)
      )
    )
  }


  category_names <- character(0)


  if (
    "__categories" %in%
      obs_names
  ) {

    category_names <- names(
      obs[["__categories"]]
    )
  }


  decode_field <- function(
      field
  ) {

    if (
      !field %in%
        obs_names
    ) {

      return(
        rep(
          NA_character_,
          length(
            obs[["_index"]][]
          )
        )
      )
    }


    raw_values <-
      obs[[field]][]


    if (
      field %in%
        category_names
    ) {

      categories <-
        obs[[
          "__categories"
        ]][[
          field
        ]][]


      return(
        decode_annData_category(
          codes =
            raw_values,
          categories =
            categories
        )
      )
    }


    as.character(
      raw_values
    )
  }


  cell_ids <-
    as.character(
      obs[["_index"]][]
    )


  n_cells <-
    length(
      cell_ids
    )


  tissue_name <-
    tools::file_path_sans_ext(
      basename(
        file_path
      )
    )


  tibble::tibble(

    Cell_Index =
      seq_len(
        n_cells
      ),

    Cell_Barcode =
      cell_ids,

    tissue =
      tissue_name,

    sex =
      decode_field(
        "sex"
      ),

    age =
      NA_character_,

    sample_id =
      NA_character_,

    fca_id =
      NA_character_,

    annotation =
      decode_field(
        "annotation"
      ),

    annotation_broad =
      NA_character_,

    annotation_broad_extrapolated =
      NA_character_,

    R_annotation =
      NA_character_,

    R_annotation_broad =
      NA_character_,

    technology =
      decode_field(
        "technology"
      )
  )
}



# ------------------------------------------------------------
# 8. Standardize metadata from original/BioHub files
# ------------------------------------------------------------

standardize_fca_metadata <- function(
    metadata
) {

  # ----------------------------------------------------------
  # Determine barcode/index column
  # ----------------------------------------------------------

  if (
    "Cell_Barcode" %in%
      colnames(metadata)
  ) {

    cell_barcode <-
      as.character(
        metadata$Cell_Barcode
      )

  } else if (
    "index" %in%
      colnames(metadata)
  ) {

    cell_barcode <-
      as.character(
        metadata$index
      )

  } else if (
    "_index" %in%
      colnames(metadata)
  ) {

    cell_barcode <-
      as.character(
        metadata$`_index`
      )

  } else {

    cell_barcode <-
      as.character(
        seq_len(
          nrow(metadata)
        )
      )
  }


  tibble::tibble(

    Cell_Index =
      seq_len(
        nrow(metadata)
      ),

    Cell_Barcode =
      cell_barcode,

    tissue =
      get_metadata_column(
        metadata,
        "tissue"
      ),

    sex =
      get_metadata_column(
        metadata,
        "sex"
      ),

    age =
      get_metadata_column(
        metadata,
        "age"
      ),

    sample_id =
      get_metadata_column(
        metadata,
        "sample_id"
      ),

    fca_id =
      get_metadata_column(
        metadata,
        "fca_id"
      ),

    annotation =
      get_metadata_column(
        metadata,
        "annotation"
      ),

    annotation_broad =
      get_metadata_column(
        metadata,
        "annotation_broad"
      ),

    annotation_broad_extrapolated =
      get_metadata_column(
        metadata,
        "annotation_broad_extrapolated"
      ),

    R_annotation =
      get_metadata_column(
        metadata,
        "R_annotation"
      ),

    R_annotation_broad =
      get_metadata_column(
        metadata,
        "R_annotation_broad"
      ),

    technology =
      get_metadata_column(
        metadata,
        "technology"
      )
  )
}



# ------------------------------------------------------------
# 9. Read gene names from nested var group
#
# Supports index or _index.
# ------------------------------------------------------------

read_gene_names_group <- function(
    var_group
) {

  var_names <- names(
    var_group
  )


  if (
    "index" %in%
      var_names
  ) {

    return(
      as.character(
        var_group[["index"]][]
      )
    )
  }


  if (
    "_index" %in%
      var_names
  ) {

    return(
      as.character(
        var_group[["_index"]][]
      )
    )
  }


  stop(
    "Could not identify gene-name index in var group."
  )
}



# ------------------------------------------------------------
# 10. Extract one gene from CSC matrix
#
# Matrix orientation:
#   cells x genes
#
# CSC indptr:
#   one pointer per gene column
# ------------------------------------------------------------

extract_gene_csc <- function(
    x_group,
    gene_position,
    gene_name,
    n_cells,
    indptr
) {

  # gene_position is 1-based R position

  start_zero <-
    indptr[
      gene_position
    ]


  end_zero <-
    indptr[
      gene_position + 1L
    ] - 1L


  expression <- numeric(
    n_cells
  )


  if (
    end_zero >=
      start_zero
  ) {

    read_positions <-
      seq.int(
        from =
          start_zero + 1L,
        to =
          end_zero + 1L
      )


    gene_values <-
      x_group[["data"]][
        read_positions
      ]


    cell_indices_zero <-
      x_group[["indices"]][
        read_positions
      ]


    cell_indices <-
      as.integer(
        cell_indices_zero
      ) + 1L


    expression[
      cell_indices
    ] <-
      as.numeric(
        gene_values
      )
  }


  tibble::tibble(

    Cell_Index =
      seq_len(
        n_cells
      ),

    Gene =
      gene_name,

    Expression =
      expression
  )
}



# ------------------------------------------------------------
# 11. Extract target genes from CSR matrix
#
# Matrix orientation:
#   cells x genes
#
# CSR indptr:
#   one pointer per cell row
#
# This is processed in cell chunks so that very large testis
# matrices do not require loading all sparse indices at once.
# ------------------------------------------------------------

extract_target_genes_csr <- function(
    x_group,
    gene_names,
    selected_genes,
    n_cells,
    chunk_size = 2000L
) {

  gene_positions <-
    match(
      selected_genes,
      gene_names
    )


  available_genes <-
    selected_genes[
      !is.na(
        gene_positions
      )
    ]


  available_positions <-
    gene_positions[
      !is.na(
        gene_positions
      )
    ]


  if (
    length(
      available_genes
    ) == 0
  ) {

    stop(
      "None of the selected genes were present in CSR matrix."
    )
  }


  # AnnData column indices are zero-based

  target_columns_zero <-
    available_positions -
    1L


  indptr <-
    as.integer(
      x_group[["indptr"]][]
    )


  if (
    length(
      indptr
    ) !=
      n_cells + 1L
  ) {

    stop(
      paste0(
        "CSR indptr length = ",
        length(indptr),
        "; expected ",
        n_cells + 1L
      )
    )
  }


  hit_list <-
    list()


  hit_counter <-
    0L


  chunk_starts <-
    seq.int(
      from = 1L,
      to = n_cells,
      by = chunk_size
    )


  for (
    cell_start in chunk_starts
  ) {

    cell_end <-
      min(
        cell_start +
          chunk_size -
          1L,
        n_cells
      )


    # Boundaries for rows cell_start:cell_end

    row_boundaries <-
      indptr[
        cell_start:
          (cell_end + 1L)
      ]


    start_zero <-
      row_boundaries[1]


    end_zero <-
      row_boundaries[
        length(
          row_boundaries
        )
      ] - 1L


    # No stored entries in this chunk

    if (
      end_zero <
        start_zero
    ) {

      next
    }


    read_positions <-
      seq.int(
        from =
          start_zero + 1L,
        to =
          end_zero + 1L
      )


    chunk_indices <-
      as.integer(
        x_group[["indices"]][
          read_positions
        ]
      )


    target_hit <-
      chunk_indices %in%
      target_columns_zero


    if (
      !any(
        target_hit
      )
    ) {

      next
    }


    local_hit_positions <-
      which(
        target_hit
      )


    global_offsets_zero <-
      start_zero +
      local_hit_positions -
      1L


    # Determine the cell row robustly from the CSR row lengths.
    # This handles empty cells (repeated indptr boundaries) correctly.

    row_lengths <-
      diff(
        row_boundaries
      )

    row_ids <-
      rep(
        seq.int(
          cell_start,
          cell_end
        ),
        times =
          row_lengths
      )

    local_offsets <-
      local_hit_positions

    hit_cells <-
      row_ids[
        local_offsets
      ]


    hit_gene_number <-
      match(
        chunk_indices[
          local_hit_positions
        ],
        target_columns_zero
      )


    hit_genes <-
      available_genes[
        hit_gene_number
      ]


    hit_values <-
      as.numeric(
        x_group[["data"]][
          read_positions[
            local_hit_positions
          ]
        ]
      )


    hit_counter <-
      hit_counter + 1L


    hit_list[[
      hit_counter
    ]] <-
      tibble::tibble(

        Cell_Index =
          hit_cells,

        Gene =
          hit_genes,

        Expression =
          hit_values
      )
  }


  # ----------------------------------------------------------
  # Build complete cell x selected-gene table
  # ----------------------------------------------------------

  complete_grid <-
    tidyr::crossing(

      Cell_Index =
        seq_len(
          n_cells
        ),

      Gene =
        available_genes
    )


  if (
    length(
      hit_list
    ) == 0
  ) {

    return(
      complete_grid %>%
        dplyr::mutate(
          Expression = 0
        )
    )
  }


  expressed <-
    dplyr::bind_rows(
      hit_list
    )


  # Sparse matrix should have max one value per cell x gene.
  # Summarise defensively in case duplicate entries occur.

  expressed <-
    expressed %>%

    dplyr::group_by(
      Cell_Index,
      Gene
    ) %>%

    dplyr::summarise(
      Expression =
        sum(
          Expression,
          na.rm = TRUE
        ),
      .groups = "drop"
    )


  complete_grid %>%

    dplyr::left_join(
      expressed,
      by = c(
        "Cell_Index",
        "Gene"
      )
    ) %>%

    dplyr::mutate(
      Expression =
        dplyr::coalesce(
          Expression,
          0
        )
    )
}



# ============================================================
# 12. PROCESS ONE FCA H5AD FILE
# ============================================================

process_fca_file <- function(
    file_path,
    selected_genes
) {

  cat("\n\n========================================\n")
  cat("PROCESSING FILE\n")
  cat("========================================\n")


  cat(
    basename(
      file_path
    ),
    "\n"
  )


  # ----------------------------------------------------------
  # A. Open H5AD
  # ----------------------------------------------------------

  h5 <-
    hdf5r::H5File$new(
      filename =
        file_path,
      mode =
        "r"
    )


  on.exit(
    try(
      h5$close_all(),
      silent = TRUE
    ),
    add = TRUE
  )


  top_level <-
    names(
      h5
    )


  # ----------------------------------------------------------
  # B. Determine matrix layout FIRST
  # ----------------------------------------------------------

  matrix_format <-
    NA_character_


  expression_scale_group <-
    NA_character_


  x_group <-
    NULL


  gene_names <-
    NULL


  cell_metadata <-
    NULL


  matrix_shape <-
    NULL


  # ==========================================================
  # FORMAT 1:
  # Older flattened raw.X / raw.var
  # ==========================================================

  if (
    "raw.X" %in%
      top_level &&
    "raw.var" %in%
      top_level
  ) {

    matrix_format <-
      "raw_flattened_csc"


    expression_scale_group <-
      "raw_full_gene_transformed"


    x_group <-
      h5[["raw.X"]]


    matrix_shape <-
      as.integer(
        x_group$
          attr_open(
            "h5sparse_shape"
          )$
          read()
      )


    raw_var <-
      h5[["raw.var"]][]


    if (
      !"index" %in%
        names(
          raw_var
        )
    ) {

      stop(
        paste(
          "raw.var does not contain index in",
          basename(file_path)
        )
      )
    }


    gene_names <-
      as.character(
        raw_var$index
      )


    raw_metadata <-
      read_fca_obs(
        h5
      )


    cell_metadata <-
      standardize_fca_metadata(
        raw_metadata
      )


    cat(
      "Matrix format: flattened raw.X CSC\n"
    )


  # ==========================================================
  # FORMAT 2:
  # Nested raw/X / raw/var
  # ==========================================================

  } else if (
    "raw" %in%
      top_level &&
    "X" %in%
      names(
        h5[["raw"]]
      ) &&
    "var" %in%
      names(
        h5[["raw"]]
      )
  ) {

    matrix_format <-
      "raw_nested_csc"


    expression_scale_group <-
      "raw_full_gene_transformed"


    x_group <-
      h5[["raw"]][["X"]]


    matrix_shape <-
      as.integer(
        hdf5r::h5attr(
          x_group,
          "shape"
        )
      )


    gene_names <-
      read_gene_names_group(
        h5[["raw"]][["var"]]
      )


    raw_metadata <-
      read_fca_obs(
        h5
      )


    cell_metadata <-
      standardize_fca_metadata(
        raw_metadata
      )


    cat(
      "Matrix format: nested raw/X CSC\n"
    )


  # ==========================================================
  # FORMAT 3:
  # X + var/_index CSR
  # ==========================================================

  } else if (
    "X" %in%
      top_level &&
    "var" %in%
      top_level &&
    inherits(
      h5[["X"]],
      "H5Group"
    ) &&
    "_index" %in%
      names(
        h5[["var"]]
      )
  ) {

    matrix_format <-
      "X_csr"


    expression_scale_group <-
      "X_full_gene_transformed"


    x_group <-
      h5[["X"]]


    matrix_shape <-
      as.integer(
        hdf5r::h5attr(
          x_group,
          "shape"
        )
      )


    gene_names <-
      as.character(
        h5[["var"]][["_index"]][]
      )


    cell_metadata <-
      read_x_only_obs(
        h5_object =
          h5,
        file_path =
          file_path
      )


    cat(
      "Matrix format: X CSR\n"
    )


  } else {

    stop(
      paste(
        "No supported FCA expression matrix was found in",
        basename(
          file_path
        )
      )
    )
  }


  # ----------------------------------------------------------
  # C. Validate dimensions
  # ----------------------------------------------------------

  if (
    length(
      matrix_shape
    ) !=
      2
  ) {

    stop(
      paste(
        "Unexpected matrix shape in",
        basename(
          file_path
        )
      )
    )
  }


  n_cells <-
    as.integer(
      matrix_shape[1]
    )


  n_genes <-
    as.integer(
      matrix_shape[2]
    )


  cat(
    "Cells:",
    n_cells,
    "\n"
  )


  cat(
    "Genes:",
    n_genes,
    "\n"
  )


  cat(
    "Metadata rows:",
    nrow(
      cell_metadata
    ),
    "\n"
  )


  # ----------------------------------------------------------
  # D. Validate metadata alignment
  # ----------------------------------------------------------

  if (
    nrow(
      cell_metadata
    ) !=
      n_cells
  ) {

    stop(
      paste0(
        "Metadata/matrix mismatch in ",
        basename(file_path),
        ": metadata rows = ",
        nrow(cell_metadata),
        "; matrix cells = ",
        n_cells
      )
    )
  }


  # ----------------------------------------------------------
  # E. Validate gene-name count
  # ----------------------------------------------------------

  if (
    length(
      gene_names
    ) !=
      n_genes
  ) {

    stop(
      paste0(
        "Gene-name count mismatch in ",
        basename(file_path),
        ": gene names = ",
        length(gene_names),
        "; matrix genes = ",
        n_genes
      )
    )
  }


  # ----------------------------------------------------------
  # F. Validate sparse format and indptr
  # ----------------------------------------------------------

  indptr <-
    as.integer(
      x_group[["indptr"]][]
    )


  if (
    matrix_format %in%
      c(
        "raw_flattened_csc",
        "raw_nested_csc"
      )
  ) {

    expected_indptr <-
      n_genes + 1L


    if (
      length(
        indptr
      ) !=
        expected_indptr
    ) {

      stop(
        paste0(
          "CSC indptr mismatch in ",
          basename(file_path),
          ": observed = ",
          length(indptr),
          "; expected = ",
          expected_indptr
        )
      )
    }

  } else if (
    matrix_format ==
      "X_csr"
  ) {

    expected_indptr <-
      n_cells + 1L


    if (
      length(
        indptr
      ) !=
        expected_indptr
    ) {

      stop(
        paste0(
          "CSR indptr mismatch in ",
          basename(file_path),
          ": observed = ",
          length(indptr),
          "; expected = ",
          expected_indptr
        )
      )
    }
  }


  # ----------------------------------------------------------
  # G. Metadata console checks
  # ----------------------------------------------------------

  cat(
    "Tissue:",
    paste(
      unique(
        cell_metadata$tissue
      ),
      collapse = ", "
    ),
    "\n"
  )


  cat(
    "Sex:",
    paste(
      unique(
        cell_metadata$sex
      ),
      collapse = ", "
    ),
    "\n"
  )


  # ----------------------------------------------------------
  # H. Gene presence
  # ----------------------------------------------------------

  gene_positions <-
    match(
      selected_genes,
      gene_names
    )


  gene_presence <-
    tibble::tibble(

      File =
        basename(
          file_path
        ),

      Matrix_Format =
        matrix_format,

      Expression_Scale_Group =
        expression_scale_group,

      Gene =
        selected_genes,

      Gene_Position =
        gene_positions,

      Present =
        !is.na(
          gene_positions
        )
    )


  cat(
    "\nTarget gene presence:\n"
  )


  print(
    gene_presence,
    n = Inf
  )


  if (
    !any(
      gene_presence$Present
    )
  ) {

    stop(
      paste(
        "None of the target serpins were present in",
        basename(file_path)
      )
    )
  }


  # ----------------------------------------------------------
  # I. Extract expression
  # ----------------------------------------------------------

  if (
    matrix_format ==
      "X_csr"
  ) {

    serpin_expression <-
      extract_target_genes_csr(

        x_group =
          x_group,

        gene_names =
          gene_names,

        selected_genes =
          selected_genes,

        n_cells =
          n_cells
      )


  } else {

    available_genes <-
      selected_genes[
        selected_genes %in%
          gene_names
      ]


    expression_list <-
      lapply(
        available_genes,
        function(
            gene_name
        ) {

          gene_position <-
            match(
              gene_name,
              gene_names
            )


          extract_gene_csc(

            x_group =
              x_group,

            gene_position =
              gene_position,

            gene_name =
              gene_name,

            n_cells =
              n_cells,

            indptr =
              indptr
          )
        }
      )


    serpin_expression <-
      dplyr::bind_rows(
        expression_list
      )
  }


  # ----------------------------------------------------------
  # J. Attach metadata
  # ----------------------------------------------------------

  serpin_cell_data <-
    serpin_expression %>%

    dplyr::left_join(
      cell_metadata,
      by =
        "Cell_Index"
    ) %>%

    dplyr::mutate(

      Source_File =
        basename(
          file_path
        ),

      Matrix_Format =
        matrix_format,

      Expression_Scale_Group =
        expression_scale_group
    )


  # ----------------------------------------------------------
  # K. Overall gene detection by file
  # ----------------------------------------------------------

  overall_detection <-
    serpin_cell_data %>%

    dplyr::group_by(
      Gene
    ) %>%

    dplyr::summarise(

      N_Cells_Total =
        dplyr::n(),

      N_Cells_Expressing =
        sum(
          Expression > 0,
          na.rm = TRUE
        ),

      Percent_Expressing =
        100 *
        N_Cells_Expressing /
        N_Cells_Total,

      Mean_Expression_All =
        mean(
          Expression,
          na.rm = TRUE
        ),

      Median_Expression_All =
        median(
          Expression,
          na.rm = TRUE
        ),

      Mean_Expression_Positive =
        ifelse(
          N_Cells_Expressing > 0,
          mean(
            Expression[
              Expression > 0
            ],
            na.rm = TRUE
          ),
          0
        ),

      Median_Expression_Positive =
        ifelse(
          N_Cells_Expressing > 0,
          median(
            Expression[
              Expression > 0
            ],
            na.rm = TRUE
          ),
          0
        ),

      .groups =
        "drop"
    ) %>%

    dplyr::mutate(

      Source_File =
        basename(
          file_path
        ),

      Matrix_Format =
        matrix_format,

      Expression_Scale_Group =
        expression_scale_group
    )


  # ----------------------------------------------------------
  # L. Cell-type summary within file
  #
  # The testis X-only dataset does not provide the full BioHub
  # metadata set, so unavailable metadata are retained as NA.
  #
  # Source_File is explicitly included in the grouping so that
  # datasets with missing Sample_ID can never be pooled
  # accidentally at this stage.
  # ----------------------------------------------------------

  sample_summary <-
    serpin_cell_data %>%

    dplyr::group_by(

      Source_File,

      Matrix_Format,

      Expression_Scale_Group,

      Gene,

      tissue,

      sex,

      age,

      sample_id,

      fca_id,

      annotation,

      annotation_broad,

      annotation_broad_extrapolated,

      R_annotation,

      R_annotation_broad

    ) %>%

    dplyr::summarise(

      N_Cells =
        dplyr::n(),

      N_Expressing =
        sum(
          Expression > 0,
          na.rm = TRUE
        ),

      Percent_Expressing =
        100 *
        N_Expressing /
        N_Cells,

      Mean_Expression_All =
        mean(
          Expression,
          na.rm = TRUE
        ),

      Median_Expression_All =
        median(
          Expression,
          na.rm = TRUE
        ),

      Mean_Expression_Positive =
        ifelse(
          N_Expressing > 0,
          mean(
            Expression[
              Expression > 0
            ],
            na.rm = TRUE
          ),
          0
        ),

      Median_Expression_Positive =
        ifelse(
          N_Expressing > 0,
          median(
            Expression[
              Expression > 0
            ],
            na.rm = TRUE
          ),
          0
        ),

      Max_Expression =
        max(
          Expression,
          na.rm = TRUE
        ),

      .groups =
        "drop"
    )


  # ----------------------------------------------------------
  # M. Standardize output names
  # ----------------------------------------------------------

  sample_summary <-
    sample_summary %>%

    dplyr::rename(

      Tissue =
        tissue,

      Sex =
        sex,

      Age =
        age,

      Sample_ID =
        sample_id,

      FCA_ID =
        fca_id,

      Annotation =
        annotation,

      Annotation_Broad =
        annotation_broad,

      Annotation_Broad_Extrapolated =
        annotation_broad_extrapolated,

      R_Annotation =
        R_annotation,

      R_Annotation_Broad =
        R_annotation_broad
    )


  # ----------------------------------------------------------
  # N. Return
  # ----------------------------------------------------------

  list(

    presence =
      gene_presence,

    detection =
      overall_detection,

    sample_summary =
      sample_summary
  )
}



# ============================================================
# 13. PROCESS ALL H5AD FILES
# ============================================================

fca_results <- vector(
  mode =
    "list",
  length =
    length(
      h5ad_files
    )
)


names(
  fca_results
) <-
  basename(
    h5ad_files
  )


failed_files <- list()


for (
  i in seq_along(
    h5ad_files
  )
) {

  cat(
    "\n\n########################################\n"
  )


  cat(
    "FILE",
    i,
    "OF",
    length(
      h5ad_files
    ),
    "\n"
  )


  cat(
    "########################################\n"
  )


  current_file <-
    h5ad_files[i]


  result <-
    tryCatch(

      {

        process_fca_file(

          file_path =
            current_file,

          selected_genes =
            selected_serpins
        )
      },

      error =
        function(
            e
        ) {

          cat(
            "\nERROR PROCESSING:\n",
            basename(
              current_file
            ),
            "\n"
          )


          cat(
            "ERROR MESSAGE:\n",
            conditionMessage(
              e
            ),
            "\n"
          )


          failed_files[[
            basename(
              current_file
            )
          ]] <<-
            conditionMessage(
              e
            )


          NULL
        }
    )


  fca_results[[i]] <-
    result


  gc()
}



# ------------------------------------------------------------
# 14. Remove failed files
# ------------------------------------------------------------

successful_results <-
  fca_results[
    !vapply(
      fca_results,
      is.null,
      logical(1)
    )
  ]


if (
  length(
    successful_results
  ) == 0
) {

  stop(
    "No FCA H5AD files were processed successfully."
  )
}



# ------------------------------------------------------------
# 15. Combine results
# ------------------------------------------------------------

combined_gene_presence <-
  dplyr::bind_rows(
    lapply(
      successful_results,
      `[[`,
      "presence"
    )
  )


combined_overall_detection <-
  dplyr::bind_rows(
    lapply(
      successful_results,
      `[[`,
      "detection"
    )
  )


combined_sample_summary <-
  dplyr::bind_rows(
    lapply(
      successful_results,
      `[[`,
      "sample_summary"
    )
  )



# ------------------------------------------------------------
# 16. Save primary tables
# ------------------------------------------------------------

readr::write_csv(

  combined_gene_presence,

  file.path(
    table_dir,
    "FCA_SERPIN_Gene_Presence_By_File.csv"
  )
)


readr::write_csv(

  combined_overall_detection,

  file.path(
    table_dir,
    "FCA_SERPIN_Overall_Detection_By_File.csv"
  )
)


readr::write_csv(

  combined_sample_summary,

  file.path(
    table_dir,
    "FCA_SERPIN_SampleLevel_CellType_Expression.csv"
  )
)



# ============================================================
# 17. POOLED DESCRIPTIVE CELL-TYPE SUMMARY
# ============================================================

pooled_celltype_summary <-
  combined_sample_summary %>%

  dplyr::group_by(

    Gene,

    Tissue,

    Sex,

    Annotation,

    Annotation_Broad

  ) %>%

  dplyr::summarise(

    N_Source_Files =
      dplyr::n_distinct(
        Source_File
      ),

    N_Samples =
      if (
        all(
          is.na(
            Sample_ID
          ) |
          Sample_ID == ""
        )
      ) {

        NA_integer_

      } else {

        dplyr::n_distinct(
          Sample_ID[
            !is.na(
              Sample_ID
            ) &
            Sample_ID != ""
          ]
        )
      },

    N_Expression_Scale_Groups =
      dplyr::n_distinct(
        Expression_Scale_Group
      ),

    Weighted_Mean_Expression_All =
      if (
        dplyr::n_distinct(
          Expression_Scale_Group
        ) == 1L
      ) {

        weighted.mean(
          Mean_Expression_All,
          w =
            N_Cells,
          na.rm =
            TRUE
        )

      } else {

        NA_real_
      },

    Weighted_Mean_Expression_Positive =
      if (
        dplyr::n_distinct(
          Expression_Scale_Group
        ) != 1L
      ) {

        NA_real_

      } else if (
        sum(
          N_Expressing,
          na.rm = TRUE
        ) > 0
      ) {

        weighted.mean(
          Mean_Expression_Positive,
          w =
            N_Expressing,
          na.rm =
            TRUE
        )

      } else {

        0
      },

    Total_Cells =
      sum(
        N_Cells,
        na.rm = TRUE
      ),

    Total_Expressing =
      sum(
        N_Expressing,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%

  dplyr::mutate(

    Percent_Expressing =
      100 *
      Total_Expressing /
      Total_Cells

  ) %>%

  dplyr::rename(

    N_Cells =
      Total_Cells,

    N_Expressing =
      Total_Expressing,

    Mean_Expression_All =
      Weighted_Mean_Expression_All,

    Mean_Expression_Positive =
      Weighted_Mean_Expression_Positive
  )


readr::write_csv(

  pooled_celltype_summary,

  file.path(
    table_dir,
    "FCA_SERPIN_Pooled_CellType_Expression.csv"
  )
)



# ============================================================
# 18. PLOT-READY TABLE
#
# Minimum 20 cells and at least one expressing cell.
# ============================================================

plot_ready_celltypes <-
  pooled_celltype_summary %>%

  dplyr::filter(
    N_Cells >= 20,
    N_Expressing > 0
  )


readr::write_csv(

  plot_ready_celltypes,

  file.path(
    table_dir,
    "FCA_SERPIN_PlotReady_CellType_Expression_N20.csv"
  )
)



# ------------------------------------------------------------
# 19. Male/female-only summary
# ------------------------------------------------------------

sex_specific_summary <-
  combined_sample_summary %>%

  dplyr::filter(
    tolower(
      Sex
    ) %in%
      c(
        "female",
        "male"
      )
  )


readr::write_csv(

  sex_specific_summary,

  file.path(
    table_dir,
    "FCA_SERPIN_SampleLevel_MaleFemale_Only.csv"
  )
)



# ------------------------------------------------------------
# 20. Spn42Dd-focused table
# ------------------------------------------------------------

spn42dd_summary <-
  plot_ready_celltypes %>%

  dplyr::filter(
    Gene ==
      "Spn42Dd"
  ) %>%

  dplyr::arrange(
    Tissue,
    Sex,
    dplyr::desc(
      Percent_Expressing
    )
  )


readr::write_csv(

  spn42dd_summary,

  file.path(
    table_dir,
    "FCA_Spn42Dd_CellType_Expression_N20.csv"
  )
)



# ------------------------------------------------------------
# 21. Failed-file log
# ------------------------------------------------------------

if (
  length(
    failed_files
  ) > 0
) {

  failed_file_table <-
    tibble::tibble(

      File =
        names(
          failed_files
        ),

      Error =
        unlist(
          failed_files,
          use.names =
            FALSE
        )
    )


  readr::write_csv(

    failed_file_table,

    file.path(
      table_dir,
      "FCA_SERPIN_Failed_Files_Log.csv"
    )
  )

} else {

  failed_file_table <-
    tibble::tibble(
      File = character(),
      Error = character()
    )
}



# ============================================================
# 22. SANITY CHECKS
# ============================================================

cat("\n\n========================================\n")
cat("ANALYSIS SANITY CHECKS\n")
cat("========================================\n")


cat(
  "\nSuccessfully processed files:",
  length(
    successful_results
  ),
  "of",
  length(
    h5ad_files
  ),
  "\n"
)


cat(
  "\nSuccessful files:\n"
)


print(
  names(
    successful_results
  )
)


if (
  length(
    failed_files
  ) > 0
) {

  cat(
    "\nFAILED FILES:\n"
  )


  print(
    failed_file_table,
    n = Inf
  )

} else {

  cat(
    "\nNo files failed.\n"
  )
}


cat(
  "\nTissues represented in combined analysis:\n"
)


analysis_tissues <-
  sort(
    unique(
      combined_sample_summary$Tissue
    )
  )


print(
  analysis_tissues
)


cat(
  "\nNumber of tissues:",
  length(
    analysis_tissues
  ),
  "\n"
)


cat(
  "\nMatrix formats represented:\n"
)


print(
  combined_sample_summary %>%

    dplyr::distinct(
      Source_File,
      Matrix_Format,
      Expression_Scale_Group
    )
)


cat(
  "\nGene presence:\n"
)


print(
  combined_gene_presence,
  n = Inf
)



# ============================================================
# EXPECTED 15-TISSUE VALIDATION
# ============================================================

expected_tissues <- c(
  "antenna",
  "body_wall",
  "fat_body",
  "gut",
  "haltere",
  "heart",
  "leg",
  "male_reproductive_glands",
  "malpighian_tubule",
  "oenocyte",
  "ovary",
  "proboscis_and_maxpalp",
  "testis",
  "trachea",
  "wing"
)


missing_expected_tissues <-
  setdiff(
    expected_tissues,
    analysis_tissues
  )


if (
  length(
    missing_expected_tissues
  ) > 0
) {

  warning(
    paste(
      "Expected tissues missing from analysis:",
      paste(
        missing_expected_tissues,
        collapse = ", "
      )
    )
  )

} else {

  cat(
    "\nAll 15 expected tissues are present.\n"
  )
}

# ============================================================
# FINAL PRE-PLOT VALIDATION
# ============================================================

plot_tissues <-
  sort(
    unique(
      plot_ready_celltypes$Tissue
    )
  )

cat(
  "\nTissues entering figures:\n"
)

print(
  plot_tissues
)

cat(
  "\nNumber of tissues entering figures:",
  length(plot_tissues),
  "\n"
)


if (
  length(plot_tissues) != 15
) {

  stop(
    paste0(
      "STOPPED: plotting dataset contains ",
      length(plot_tissues),
      " tissues instead of the expected 15. ",
      "Do not generate figures until the upstream analysis is rerun."
    )
  )
}

# ============================================================
# FIGURE-DATA VALIDATION
# ============================================================

if (
  length(
    unique(
      plot_ready_celltypes$Tissue
    )
  ) != 15L
) {
  stop(
    "Plot-ready dataset does not contain all 15 tissues."
  )
}

cat(
  "\nPlot-ready rows:",
  nrow(plot_ready_celltypes),
  "\n"
)

cat(
  "Plot-ready tissues:",
  paste(
    sort(
      unique(
        plot_ready_celltypes$Tissue
      )
    ),
    collapse = ", "
  ),
  "\n"
)


# ============================================================
# 23. ALL-SERPIN DOT PLOT
#
# IMPORTANT:
# Mean expression is NOT used as a common color scale because
# expression scales differ between FCA file groups.
#
# Dot size = Percent_Expressing
# ============================================================

dotplot_data <-
  plot_ready_celltypes %>%

  dplyr::mutate(

    Tissue_CellType =
      paste(
        Tissue,
        Annotation,
        sep = " | "
      ),

    Gene =
      factor(
        Gene,
        levels =
          selected_serpins
      )
  )


# Use ALL plot-ready cell populations from all 15 tissues.
# No additional top-N filtering is applied here.

dotplot_final <-
  dotplot_data


readr::write_csv(

  dotplot_final,

  file.path(
    table_dir,
    "FCA_SERPIN_AllGenes_DotPlot_Data.csv"
  )
)


serpin_dotplot <-
  ggplot2::ggplot(

    dotplot_final,

    ggplot2::aes(
      x =
        Gene,
      y =
        Tissue_CellType
    )
  ) +

  ggplot2::geom_point(

    ggplot2::aes(
      size =
        Percent_Expressing
    ),

    alpha =
      0.8
  ) +

  ggplot2::scale_size_continuous(

    name =
      "% expressing",

    range =
      c(
        1,
        9
      )
  ) +

  ggplot2::labs(

    title =
      "Fly Cell Atlas serpin expression by tissue and cell type",

    subtitle =
      "Dot size represents percentage of cells expressing each gene",

    x =
      NULL,

    y =
      NULL
  ) +

  ggplot2::theme_bw(
    base_size = 11
  ) +

  ggplot2::theme(

    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 11,
        face = "bold"
      ),

    axis.text.y =
      ggplot2::element_text(
        size = 9
      ),

    legend.title =
      ggplot2::element_text(
        size = 11,
        face = "bold"
      ),

    legend.text =
      ggplot2::element_text(
        size = 9
      ),

    plot.title =
      ggplot2::element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      ),

    plot.subtitle =
      ggplot2::element_text(
        hjust = 0.5,
        size = 11
      )
  )


print(
  serpin_dotplot
)


ggplot2::ggsave(

  filename =
    file.path(
      figure_dir,
      "FCA_SERPIN_AllGenes_CellType_DotPlot.png"
    ),

  plot =
    serpin_dotplot,

  width =
    11,

  height =
    16,

  units =
    "in",

  dpi =
    600
)



# ============================================================
# 24. SEX-SPECIFIC DOT PLOTS
# ============================================================

sex_dotplot_data <-
  plot_ready_celltypes %>%

  dplyr::filter(
    tolower(
      Sex
    ) %in%
      c(
        "female",
        "male"
      )
  ) %>%

  dplyr::mutate(

    Sex =
      tolower(
        Sex
      ),

    Tissue_CellType =
      paste(
        Tissue,
        Annotation,
        sep = " | "
      ),

    Gene =
      factor(
        Gene,
        levels =
          selected_serpins
      )
  )


# Use ALL plot-ready female and male populations.
# Sex-specific tissue counts are expected to differ biologically:
# female excludes male-only testis/male reproductive gland data;
# male excludes ovary.

female_dotplot_data <-
  sex_dotplot_data %>%

  dplyr::filter(
    Sex ==
      "female"
  )


male_dotplot_data <-
  sex_dotplot_data %>%

  dplyr::filter(
    Sex ==
      "male"
  )


max_percent <-
  max(
    sex_dotplot_data$Percent_Expressing,
    na.rm = TRUE
  )



# ------------------------------------------------------------
# Female
# ------------------------------------------------------------

female_dotplot <-
  ggplot2::ggplot(

    female_dotplot_data,

    ggplot2::aes(
      x =
        Gene,
      y =
        Tissue_CellType
    )
  ) +

  ggplot2::geom_point(

    ggplot2::aes(
      size =
        Percent_Expressing
    ),

    alpha =
      0.8
  ) +

  ggplot2::scale_size_continuous(

    name =
      "% expressing",

    limits =
      c(
        0,
        max_percent
      ),

    range =
      c(
        1,
        9
      )
  ) +

  ggplot2::labs(

    title =
      "Female Fly Cell Atlas serpin expression",

    subtitle =
      "Dot size represents percentage of cells expressing each gene",

    x =
      NULL,

    y =
      NULL
  ) +

  ggplot2::theme_bw(
    base_size = 11
  ) +

  ggplot2::theme(

    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 11,
        face = "bold"
      ),

    axis.text.y =
      ggplot2::element_text(
        size = 9
      ),

    legend.title =
      ggplot2::element_text(
        size = 11,
        face = "bold"
      ),

    legend.text =
      ggplot2::element_text(
        size = 9
      ),

    plot.title =
      ggplot2::element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      ),

    plot.subtitle =
      ggplot2::element_text(
        hjust = 0.5,
        size = 11
      )
  )


print(
  female_dotplot
)


ggplot2::ggsave(

  filename =
    file.path(
      figure_dir,
      "FCA_SERPIN_Female_CellType_DotPlot.png"
    ),

  plot =
    female_dotplot,

  width =
    11,

  height =
    16,

  units =
    "in",

  dpi =
    600
)



# ------------------------------------------------------------
# Male
# ------------------------------------------------------------

male_dotplot <-
  ggplot2::ggplot(

    male_dotplot_data,

    ggplot2::aes(
      x =
        Gene,
      y =
        Tissue_CellType
    )
  ) +

  ggplot2::geom_point(

    ggplot2::aes(
      size =
        Percent_Expressing
    ),

    alpha =
      0.8
  ) +

  ggplot2::scale_size_continuous(

    name =
      "% expressing",

    limits =
      c(
        0,
        max_percent
      ),

    range =
      c(
        1,
        9
      )
  ) +

  ggplot2::labs(

    title =
      "Male Fly Cell Atlas serpin expression",

    subtitle =
      "Dot size represents percentage of cells expressing each gene",

    x =
      NULL,

    y =
      NULL
  ) +

  ggplot2::theme_bw(
    base_size = 11
  ) +

  ggplot2::theme(

    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 11,
        face = "bold"
      ),

    axis.text.y =
      ggplot2::element_text(
        size = 9
      ),

    legend.title =
      ggplot2::element_text(
        size = 11,
        face = "bold"
      ),

    legend.text =
      ggplot2::element_text(
        size = 9
      ),

    plot.title =
      ggplot2::element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      ),

    plot.subtitle =
      ggplot2::element_text(
        hjust = 0.5,
        size = 11
      )
  )


print(
  male_dotplot
)


ggplot2::ggsave(

  filename =
    file.path(
      figure_dir,
      "FCA_SERPIN_Male_CellType_DotPlot.png"
    ),

  plot =
    male_dotplot,

  width =
    11,

  height =
    16,

  units =
    "in",

  dpi =
    600
)



# ------------------------------------------------------------
# Save sex-specific plotting tables
# ------------------------------------------------------------

readr::write_csv(

  female_dotplot_data,

  file.path(
    table_dir,
    "FCA_SERPIN_Female_DotPlot_Data.csv"
  )
)


readr::write_csv(

  male_dotplot_data,

  file.path(
    table_dir,
    "FCA_SERPIN_Male_DotPlot_Data.csv"
  )
)



# ============================================================
# 25. Spn42Dd-ONLY TABLES AND PLOTS
# ============================================================

spn42dd_plot_data <-
  plot_ready_celltypes %>%

  dplyr::filter(
    Gene ==
      "Spn42Dd"
  ) %>%

  dplyr::mutate(

    Tissue_CellType =
      paste(
        Tissue,
        Annotation,
        sep = " | "
      ),

    Sex =
      tolower(
        Sex
      )
  )


readr::write_csv(

  spn42dd_plot_data,

  file.path(
    table_dir,
    "FCA_Spn42Dd_All_CellType_DotPlot_Data.csv"
  )
)


spn42dd_female_data <-
  spn42dd_plot_data %>%

  dplyr::filter(
    Sex ==
      "female"
  )


spn42dd_male_data <-
  spn42dd_plot_data %>%

  dplyr::filter(
    Sex ==
      "male"
  )


spn42dd_max_percent <-
  max(
    spn42dd_plot_data$Percent_Expressing,
    na.rm = TRUE
  )



# ------------------------------------------------------------
# Female Spn42Dd
# ------------------------------------------------------------

spn42dd_female_plot <-
  ggplot2::ggplot(

    spn42dd_female_data,

    ggplot2::aes(
      x =
        "Spn42Dd",
      y =
        Tissue_CellType
    )
  ) +

  ggplot2::geom_point(

    ggplot2::aes(
      size =
        Percent_Expressing
    ),

    alpha =
      0.8
  ) +

  ggplot2::scale_size_continuous(

    name =
      "% expressing",

    limits =
      c(
        0,
        spn42dd_max_percent
      ),

    range =
      c(
        1,
        9
      )
  ) +

  ggplot2::labs(

    title =
      "Female Fly Cell Atlas Spn42Dd expression",

    x =
      NULL,

    y =
      NULL
  ) +

  ggplot2::theme_bw(
    base_size = 11
  ) +

  ggplot2::theme(

    axis.text.x =
      ggplot2::element_text(
        size = 11,
        face = "bold"
      ),

    axis.text.y =
      ggplot2::element_text(
        size = 9
      ),

    legend.title =
      ggplot2::element_text(
        size = 11,
        face = "bold"
      ),

    legend.text =
      ggplot2::element_text(
        size = 9
      ),

    plot.title =
      ggplot2::element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      )
  )


print(
  spn42dd_female_plot
)


ggplot2::ggsave(

  filename =
    file.path(
      figure_dir,
      "FCA_Spn42Dd_Female_CellType_DotPlot.png"
    ),

  plot =
    spn42dd_female_plot,

  width =
    8,

  height =
    16,

  units =
    "in",

  dpi =
    600
)



# ------------------------------------------------------------
# Male Spn42Dd
# ------------------------------------------------------------

spn42dd_male_plot <-
  ggplot2::ggplot(

    spn42dd_male_data,

    ggplot2::aes(
      x =
        "Spn42Dd",
      y =
        Tissue_CellType
    )
  ) +

  ggplot2::geom_point(

    ggplot2::aes(
      size =
        Percent_Expressing
    ),

    alpha =
      0.8
  ) +

  ggplot2::scale_size_continuous(

    name =
      "% expressing",

    limits =
      c(
        0,
        spn42dd_max_percent
      ),

    range =
      c(
        1,
        9
      )
  ) +

  ggplot2::labs(

    title =
      "Male Fly Cell Atlas Spn42Dd expression",

    x =
      NULL,

    y =
      NULL
  ) +

  ggplot2::theme_bw(
    base_size = 11
  ) +

  ggplot2::theme(

    axis.text.x =
      ggplot2::element_text(
        size = 11,
        face = "bold"
      ),

    axis.text.y =
      ggplot2::element_text(
        size = 9
      ),

    legend.title =
      ggplot2::element_text(
        size = 11,
        face = "bold"
      ),

    legend.text =
      ggplot2::element_text(
        size = 9
      ),

    plot.title =
      ggplot2::element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      )
  )


print(
  spn42dd_male_plot
)


ggplot2::ggsave(

  filename =
    file.path(
      figure_dir,
      "FCA_Spn42Dd_Male_CellType_DotPlot.png"
    ),

  plot =
    spn42dd_male_plot,

  width =
    8,

  height =
    16,

  units =
    "in",

  dpi =
    600
)


readr::write_csv(

  spn42dd_female_data,

  file.path(
    table_dir,
    "FCA_Spn42Dd_Female_DotPlot_Data.csv"
  )
)


readr::write_csv(

  spn42dd_male_data,

  file.path(
    table_dir,
    "FCA_Spn42Dd_Male_DotPlot_Data.csv"
  )
)



# ============================================================
# 26. TOP FIVE CELL TYPES FOR EACH GENE AND SEX
#
# Ranking metric:
#   Percent_Expressing
#
# This is the primary cross-dataset metric.
# ============================================================

serpin_sex_data <-
  plot_ready_celltypes %>%

  dplyr::filter(
    tolower(
      Sex
    ) %in%
      c(
        "female",
        "male"
      )
  ) %>%

  dplyr::mutate(

    Sex =
      stringr::str_to_title(
        tolower(
          Sex
        )
      ),

    Gene =
      factor(
        Gene,
        levels =
          selected_serpins
      ),

    Tissue_CellType =
      paste(
        Tissue,
        Annotation,
        sep = " | "
      )
  )



# ------------------------------------------------------------
# Rank top 5 by percent expressing
#
# Ties:
#   1. Percent_Expressing
#   2. N_Expressing
#   3. N_Cells
# ------------------------------------------------------------

top5_gene_sex <-
  serpin_sex_data %>%

  dplyr::group_by(
    Sex,
    Gene
  ) %>%

  dplyr::arrange(

    dplyr::desc(
      Percent_Expressing
    ),

    dplyr::desc(
      N_Expressing
    ),

    dplyr::desc(
      N_Cells
    ),

    .by_group = TRUE
  ) %>%

  dplyr::slice_head(
    n = 5
  ) %>%

  dplyr::mutate(
    Rank =
      dplyr::row_number()
  ) %>%

  dplyr::ungroup() %>%

  dplyr::mutate(

    Plot_Label =
      paste0(
        Rank,
        ". ",
        Tissue,
        " | ",
        Annotation
      ),

    Panel =
      paste(
        Sex,
        Gene,
        sep = " | "
      )
  )



# ------------------------------------------------------------
# Save top-five table
# ------------------------------------------------------------

readr::write_csv(

  top5_gene_sex,

  file.path(
    table_dir,
    "FCA_SERPIN_Top5_CellTypes_By_Gene_and_Sex.csv"
  )
)



# ------------------------------------------------------------
# Concise publication-ready table
# ------------------------------------------------------------

top5_table <-
  top5_gene_sex %>%

  dplyr::select(

    Sex,

    Gene,

    Rank,

    Tissue,

    Annotation,

    Annotation_Broad,

    N_Samples,

    N_Cells,

    N_Expressing,

    Percent_Expressing,

    Mean_Expression_All,

    Mean_Expression_Positive
  )


readr::write_csv(

  top5_table,

  file.path(
    table_dir,
    "FCA_SERPIN_Top5_CellTypes_Concise_Table.csv"
  )
)


print(
  top5_table,
  n = Inf,
  width = Inf
)



# ============================================================
# 27. TOP-FIVE FIGURE — REVISED
#
# Separate panel for each Gene x Sex combination.
#
# x-axis = Percent expressing
# point size = Number of expressing cells
#
# Figure dimensions retained at 15 x 14 inches.
#
# IMPORTANT:
# Two facet columns are used because the cell-type labels are
# too long for a three-column layout.
# ============================================================


# ------------------------------------------------------------
# A. Prepare cleaner plotting labels
# ------------------------------------------------------------

top5_plot_data <-
  top5_gene_sex %>%

  dplyr::arrange(
    Sex,
    Gene,
    Rank
  ) %>%

  dplyr::mutate(

    # Make tissue names publication-friendly
    Tissue_Clean =
      dplyr::recode(
        as.character(
          Tissue
        ),

        "body_wall" =
          "body wall",

        "fat_body" =
          "fat body",

        "male_reproductive_glands" =
          "male reproductive glands",

        "malpighian_tubule" =
          "Malpighian tubule",

        "proboscis_and_maxpalp" =
          "proboscis/maxillary palp",

        .default =
          stringr::str_replace_all(
            as.character(
              Tissue
            ),
            "_",
            " "
          )
      ),

    # Label shown beside each point
    Display_Label =
      paste0(
        Rank,
        ". ",
        Tissue_Clean,
        " | ",
        Annotation
      ),

    # Wrap long labels so they do not destroy panel width
    Display_Label =
      stringr::str_wrap(
        Display_Label,
        width = 42
      ),

    # Cleaner facet heading
    Panel =
      paste0(
        Sex,
        "\n",
        Gene
      ),

    # Unique internal factor label
    Unique_Label =
      paste(
        Panel,
        Rank,
        Display_Label,
        sep = "___"
      )
  )


# ------------------------------------------------------------
# B. Preserve rank order independently within every panel
# ------------------------------------------------------------

top5_plot_data <-
  top5_plot_data %>%

  dplyr::mutate(

    Unique_Label =
      factor(
        Unique_Label,
        levels =
          rev(
            unique(
              Unique_Label
            )
          )
      )
  )


# ------------------------------------------------------------
# C. Function for displaying only the clean label
# ------------------------------------------------------------

clean_top5_label <- function(
    x
) {

  sub(
    "^.*___[0-9]+___",
    "",
    x
  )
}


# ------------------------------------------------------------
# D. Create figure
# ------------------------------------------------------------

top5_serpin_plot <-
  ggplot2::ggplot(

    top5_plot_data,

    ggplot2::aes(
      x =
        Percent_Expressing,
      y =
        Unique_Label
    )
  ) +

  # Horizontal lollipop line
  ggplot2::geom_segment(

    ggplot2::aes(
      x =
        0,
      xend =
        Percent_Expressing,
      yend =
        Unique_Label
    ),

    linewidth =
      0.5,

    alpha =
      0.5
  ) +

  # Point
  ggplot2::geom_point(

    ggplot2::aes(
      size =
        N_Expressing
    ),

    alpha =
      0.85
  ) +

  # ----------------------------------------------------------
  # TWO columns instead of three
  # ----------------------------------------------------------

  ggplot2::facet_wrap(

    ~ Panel,

    scales =
      "free_y",

    ncol =
      2
  ) +

  # ----------------------------------------------------------
  # Y labels
  # ----------------------------------------------------------

  ggplot2::scale_y_discrete(

    labels =
      clean_top5_label
  ) +

  # ----------------------------------------------------------
  # X-axis
  # ----------------------------------------------------------

  ggplot2::scale_x_continuous(

    limits =
      c(
        0,
        100
      ),

    breaks =
      c(
        0,
        25,
        50,
        75,
        100
      ),

    labels =
      function(
          x
      ) {

        paste0(
          x,
          "%"
        )
      },

    expand =
      ggplot2::expansion(
        mult =
          c(
            0,
            0.04
          )
      )
  ) +

  # ----------------------------------------------------------
  # Point-size legend
  # ----------------------------------------------------------

  ggplot2::scale_size_continuous(

    name =
      "Expressing cells",

    range =
      c(
        2,
        7
      )
  ) +

  # ----------------------------------------------------------
  # Labels
  # ----------------------------------------------------------

  ggplot2::labs(

    title =
      "Top five serpin-expressing cell types by gene and sex",

    subtitle =
      paste0(
        "15-tissue Fly Cell Atlas reanalysis; ",
        "cell types ranked by percentage of cells expressing each gene"
      ),

    x =
      "% expressing",

    y =
      NULL
  ) +

  # ----------------------------------------------------------
  # Theme
  #
  # All plotted text is >= 9 pt.
  # ----------------------------------------------------------

  ggplot2::theme_bw(
    base_size = 10
  ) +

  ggplot2::theme(

    axis.title.x =
      ggplot2::element_text(
        size = 11,
        face = "bold",
        margin =
          ggplot2::margin(
            t = 8
          )
      ),

    axis.text.x =
      ggplot2::element_text(
        size = 9
      ),

    axis.text.y =
      ggplot2::element_text(
        size = 9,
        lineheight = 0.95
      ),

    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 10,
        lineheight = 1.0
      ),

    strip.background =
      ggplot2::element_rect(
        linewidth = 0.5
      ),

    legend.title =
      ggplot2::element_text(
        size = 10,
        face = "bold"
      ),

    legend.text =
      ggplot2::element_text(
        size = 9
      ),

    plot.title =
      ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 15
      ),

    plot.subtitle =
      ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        margin =
          ggplot2::margin(
            b = 10
          )
      ),

    panel.grid.major.y =
      ggplot2::element_blank(),

    panel.grid.minor =
      ggplot2::element_blank(),

    panel.spacing =
      grid::unit(
        0.8,
        "lines"
      ),

    plot.margin =
      ggplot2::margin(
        t = 10,
        r = 15,
        b = 10,
        l = 10
      )
  )


# ------------------------------------------------------------
# E. Display
# ------------------------------------------------------------

print(
  top5_serpin_plot
)


# ------------------------------------------------------------
# F. Save
# ------------------------------------------------------------

ggplot2::ggsave(

  filename =
    file.path(
      figure_dir,
      "FCA_SERPIN_Top5_CellTypes_By_Gene_and_Sex.png"
    ),

  plot =
    top5_serpin_plot,

  width =
    15,

  height =
    14,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)

# ============================================================
# 28. FINAL FIGURE OUTPUT VALIDATION
# ============================================================

expected_figure_files <- c(
  "FCA_SERPIN_AllGenes_CellType_DotPlot.png",
  "FCA_SERPIN_Female_CellType_DotPlot.png",
  "FCA_SERPIN_Male_CellType_DotPlot.png",
  "FCA_Spn42Dd_Female_CellType_DotPlot.png",
  "FCA_Spn42Dd_Male_CellType_DotPlot.png",
  "FCA_SERPIN_Top5_CellTypes_By_Gene_and_Sex.png"
)

expected_figure_paths <-
  file.path(
    figure_dir,
    expected_figure_files
  )

if (
  !all(
    file.exists(
      expected_figure_paths
    )
  )
) {
  stop(
    paste0(
      "One or more expected figure files were not created: ",
      paste(
        expected_figure_files[
          !file.exists(
            expected_figure_paths
          )
        ],
        collapse = ", "
      )
    )
  )
}

cat(
  "\nAll expected figures were created successfully.\n"
)

print(
  file.info(
    expected_figure_paths
  )[
    ,
    c(
      "size",
      "mtime"
    )
  ]
)
