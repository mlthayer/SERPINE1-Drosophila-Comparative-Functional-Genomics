# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 08_CrossSpecies_TissueSimilarity.R
# Run from the repository root.
# ==============================================================================
# 11_CrossSpecies_SERPIN_4Dataset_TissueSimilarity_REVISED.R
#
# Cross-species SERPIN tissue-expression similarity using FOUR independent
# expression resources:
#
#   1. FlyAtlas2 bulk Drosophila tissue expression
#   2. Fly Cell Atlas (FCA) single-cell tissue prevalence
#   3. Human GTEx tissue expression
#   4. Human Protein Atlas (HPA) tissue expression
#
# Four comparisons are kept INDEPENDENT:
#
#   FlyAtlas2 vs GTEx
#   FlyAtlas2 vs HPA
#   FCA       vs GTEx
#   FCA       vs HPA
#
# The datasets are NOT averaged together.
#
# Similarity statistic:
#   Spearman correlation across matched cross-species tissue systems.
#
# Final cross-species systems:
#   Heart
#   Skeletal muscle
#   Liver-metabolic
#   Digestive
#   Renal/excretory
#   Respiratory
#   Ovary
#   Male reproductive
#   Sensory/neural
#
# FCA tissue grouping:
#   heart                              -> Heart
#   body_wall + leg + wing + haltere  -> Skeletal muscle
#   fat_body + oenocyte               -> Liver-metabolic
#   gut                                -> Digestive
#   malpighian_tubule                  -> Renal/excretory
#   trachea                            -> Respiratory
#   ovary                              -> Ovary
#   male_reproductive_glands           -> Male reproductive
#   antenna + proboscis_and_maxpalp    -> Sensory/neural
#
# FCA testis is excluded.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)


# ------------------------------------------------------------------------------
# 1. EXACT PROJECT DIRECTORIES
# ------------------------------------------------------------------------------

# Fly Cell Atlas output directory
fca_output_dir <- "Expression/FlyCellAtlas/Reanalysis_Output"

# FlyAtlas2 analysis directory
flyatlas2_dir <- "Expression/FlyAtlas2_Analysis"

# Human GTEx + HPA reanalysis directory
human_dir <- "Expression/Human_SERPIN_reanalysis"


# ------------------------------------------------------------------------------
# 2. OUTPUT DIRECTORY
# ------------------------------------------------------------------------------

output_dir <- file.path(
  human_dir,
  "CrossSpecies_4Dataset_TissueSimilarity"
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
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. FCA INPUT FILE
# ------------------------------------------------------------------------------

fca_file <- file.path(
  fca_output_dir,
  "Tables",
  "FCA_SERPIN_PlotReady_CellType_Expression_N20.csv"
)

if (!file.exists(fca_file)) {
  stop(
    paste0(
      "FCA input file not found:\n",
      fca_file
    )
  )
}


# ------------------------------------------------------------------------------
# 4. INPUT FILE DISCOVERY FOR FLYATLAS2 / GTEx / HPA
#
# Because the exact filenames may differ from earlier working names, this script
# searches the confirmed directories recursively for likely CSV files.
#
# If more than one strong candidate is found, the script stops and prints them
# so you can explicitly choose the correct one.
# ------------------------------------------------------------------------------

find_candidate_file <- function(
    search_dir,
    filename_regex,
    dataset_name
) {

  if (!dir.exists(search_dir)) {

    stop(
      paste0(
        dataset_name,
        " directory does not exist:\n",
        search_dir
      )
    )
  }

  files <- list.files(
    search_dir,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0) {

    stop(
      paste0(
        "No CSV files found under ",
        dataset_name,
        " directory:\n",
        search_dir
      )
    )
  }

  matches <- files[
    stringr::str_detect(
      basename(files),
      stringr::regex(
        filename_regex,
        ignore_case = TRUE
      )
    )
  ]

  if (length(matches) == 0) {

    cat(
      "\nCSV files found under ",
      dataset_name,
      " directory:\n",
      sep = ""
    )

    print(files)

    stop(
      paste0(
        "\nNo CSV filename matched the expected ",
        dataset_name,
        " pattern.\n",
        "Edit the corresponding *_file object manually near the top of this script."
      )
    )
  }

  # Prefer summary/output files over raw/intermediate files
  preferred <- matches[
    stringr::str_detect(
      basename(matches),
      stringr::regex(
        "summary|tissue|expression|serpin|reanaly",
        ignore_case = TRUE
      )
    )
  ]

  if (length(preferred) == 1) {
    return(preferred)
  }

  if (length(preferred) > 1) {

    cat(
      "\nMultiple likely ",
      dataset_name,
      " files were found:\n",
      sep = ""
    )

    print(preferred)

    stop(
      paste0(
        "\nPlease set ",
        dataset_name,
        " file explicitly near the top of the script."
      )
    )
  }

  if (length(matches) == 1) {
    return(matches)
  }

  cat(
    "\nMultiple ",
    dataset_name,
    " CSV files matched:\n",
    sep = ""
  )

  print(matches)

  stop(
    paste0(
      "\nPlease set ",
      dataset_name,
      " file explicitly near the top of the script."
    )
  )
}


# ------------------------------------------------------------------------------
# 5. AUTO-DISCOVER INPUT FILES
#
# If you already know the exact filenames, you may replace these three blocks
# with direct file.path(...) assignments.
# ------------------------------------------------------------------------------

flyatlas2_file <- find_candidate_file(
  search_dir = flyatlas2_dir,
  filename_regex = "flyatlas|fly_atlas|fly atlas",
  dataset_name = "FlyAtlas2"
)

gtex_file <- find_candidate_file(
  search_dir = human_dir,
  filename_regex = "gtex",
  dataset_name = "GTEx"
)

hpa_file <- find_candidate_file(
  search_dir = human_dir,
  filename_regex = "hpa|protein.?atlas",
  dataset_name = "HPA"
)

cat("\nInput files selected:\n")
cat("\nFCA:\n", fca_file, "\n")
cat("\nFlyAtlas2:\n", flyatlas2_file, "\n")
cat("\nGTEx:\n", gtex_file, "\n")
cat("\nHPA:\n", hpa_file, "\n")


# ------------------------------------------------------------------------------
# 6. GENES
# ------------------------------------------------------------------------------

fly_serpins <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)

human_serpins <- c(
  "SERPINE1",
  "SERPINI1",
  "SERPINB3",
  "SERPINB4",
  "SERPINB11",
  "SERPINB9"
)

system_order <- c(
  "Heart",
  "Skeletal muscle",
  "Liver-metabolic",
  "Digestive",
  "Renal/excretory",
  "Respiratory",
  "Ovary",
  "Male reproductive",
  "Sensory/neural"
)


# ------------------------------------------------------------------------------
# 7. FINAL FCA TISSUE-TO-SYSTEM MAPPING
# ------------------------------------------------------------------------------

fca_mapping <- tibble::tribble(
  ~Tissue,                       ~System,

  "heart",                       "Heart",

  "body_wall",                   "Skeletal muscle",
  "leg",                         "Skeletal muscle",
  "wing",                        "Skeletal muscle",
  "haltere",                     "Skeletal muscle",

  "fat_body",                    "Liver-metabolic",
  "oenocyte",                    "Liver-metabolic",

  "gut",                         "Digestive",

  "malpighian_tubule",           "Renal/excretory",

  "trachea",                     "Respiratory",

  "ovary",                       "Ovary",

  "male_reproductive_glands",    "Male reproductive",

  "antenna",                     "Sensory/neural",
  "proboscis_and_maxpalp",       "Sensory/neural"
)

readr::write_csv(
  fca_mapping,
  file.path(
    table_dir,
    "FCA_Final_Tissue_to_System_Mapping.csv"
  )
)


# ------------------------------------------------------------------------------
# 8. FLYATLAS2 TISSUE CLASSIFIER
# ------------------------------------------------------------------------------

classify_flyatlas2_tissue <- function(x) {

  x2 <- x %>%
    as.character() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[_-]", " ") %>%
    stringr::str_squish()

  dplyr::case_when(

    stringr::str_detect(
      x2,
      "\\bheart\\b|cardiac"
    ) ~ "Heart",

    stringr::str_detect(
      x2,
      "carcass|body wall|skeletal muscle|thorax|leg|wing|haltere"
    ) ~ "Skeletal muscle",

    stringr::str_detect(
      x2,
      "fat body|oenocyte"
    ) ~ "Liver-metabolic",

    stringr::str_detect(
      x2,
      "midgut|hindgut|foregut|crop|gut|digest"
    ) ~ "Digestive",

    stringr::str_detect(
      x2,
      "malpigh"
    ) ~ "Renal/excretory",

    stringr::str_detect(
      x2,
      "trache"
    ) ~ "Respiratory",

    stringr::str_detect(
      x2,
      "\\bovary\\b|ovaries"
    ) ~ "Ovary",

    stringr::str_detect(
      x2,
      "accessory gland|male reproductive|seminal vesicle"
    ) ~ "Male reproductive",

    stringr::str_detect(
      x2,
      "antenna|brain|eye|head|sensory|proboscis|maxillary"
    ) ~ "Sensory/neural",

    TRUE ~ NA_character_
  )
}


# ------------------------------------------------------------------------------
# 9. HUMAN TISSUE CLASSIFIER
# ------------------------------------------------------------------------------

classify_human_tissue <- function(x) {

  x2 <- x %>%
    as.character() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[_-]", " ") %>%
    stringr::str_squish()

  dplyr::case_when(

    stringr::str_detect(
      x2,
      "\\bheart\\b|cardiac"
    ) ~ "Heart",

    stringr::str_detect(
      x2,
      "skeletal muscle|muscle skeletal|muscle - skeletal"
    ) ~ "Skeletal muscle",

    stringr::str_detect(
      x2,
      "\\bliver\\b|adipose"
    ) ~ "Liver-metabolic",

    stringr::str_detect(
      x2,
      paste0(
        "stomach|esophagus|oesophagus|small intestine|duodenum|",
        "jejunum|ileum|colon|large intestine|rectum|appendix|",
        "gallbladder|gall bladder|digestive"
      )
    ) ~ "Digestive",

    stringr::str_detect(
      x2,
      "kidney|renal"
    ) ~ "Renal/excretory",

    stringr::str_detect(
      x2,
      "\\blung\\b|respiratory"
    ) ~ "Respiratory",

    stringr::str_detect(
      x2,
      "\\bovary\\b"
    ) ~ "Ovary",

    stringr::str_detect(
      x2,
      "prostate|seminal|epididym|vas deferens|male reproductive"
    ) ~ "Male reproductive",

    stringr::str_detect(
      x2,
      "brain|nerve|neural|nervous|sensory"
    ) ~ "Sensory/neural",

    TRUE ~ NA_character_
  )
}


# ------------------------------------------------------------------------------
# 10. FLEXIBLE COLUMN HELPERS
# ------------------------------------------------------------------------------

first_existing_column <- function(df, candidates) {

  hits <- candidates[
    candidates %in% colnames(df)
  ]

  if (length(hits) == 0) {
    return(NA_character_)
  }

  hits[1]
}


standardize_flyatlas2 <- function(df) {

  gene_col <- first_existing_column(
    df,
    c(
      "Gene",
      "gene",
      "gene_symbol",
      "Gene_Symbol",
      "symbol",
      "Symbol"
    )
  )

  tissue_col <- first_existing_column(
    df,
    c(
      "Tissue",
      "tissue",
      "tissue_name",
      "Tissue_Name"
    )
  )

  expression_col <- first_existing_column(
    df,
    c(
      "MedianFPKM",
      "median_FPKM",
      "median_fpkm",
      "MeanFPKM",
      "mean_FPKM",
      "mean_fpkm",
      "FPKM",
      "fpkm"
    )
  )

  if (
    any(
      is.na(
        c(
          gene_col,
          tissue_col,
          expression_col
        )
      )
    )
  ) {

    stop(
      paste0(
        "FlyAtlas2 table needs recognizable Gene, Tissue, and FPKM columns.\n",
        "Columns found:\n",
        paste(
          colnames(df),
          collapse = ", "
        )
      )
    )
  }

  df %>%
    dplyr::transmute(
      Gene = as.character(
        .data[[gene_col]]
      ),
      Tissue = as.character(
        .data[[tissue_col]]
      ),
      FlyAtlas2_Expression = as.numeric(
        .data[[expression_col]]
      )
    )
}


standardize_human <- function(df, source_name) {

  gene_col <- first_existing_column(
    df,
    c(
      "Gene",
      "gene",
      "gene_symbol",
      "Gene_Symbol",
      "symbol",
      "Symbol"
    )
  )

  tissue_col <- first_existing_column(
    df,
    c(
      "Tissue",
      "tissue",
      "tissue_name",
      "Tissue_Name"
    )
  )

  expression_col <- first_existing_column(
    df,
    c(
      "median_TPM",
      "Median_TPM",
      "median_tpm",
      "median",
      "Median",
      "TPM",
      "tpm"
    )
  )

  if (
    any(
      is.na(
        c(
          gene_col,
          tissue_col,
          expression_col
        )
      )
    )
  ) {

    stop(
      paste0(
        source_name,
        " table needs recognizable Gene, Tissue, and median TPM columns.\n",
        "Columns found:\n",
        paste(
          colnames(df),
          collapse = ", "
        )
      )
    )
  }

  df %>%
    dplyr::transmute(
      Gene = as.character(
        .data[[gene_col]]
      ),
      Tissue = as.character(
        .data[[tissue_col]]
      ),
      Human_Expression = as.numeric(
        .data[[expression_col]]
      )
    )
}


# ------------------------------------------------------------------------------
# 11. READ FCA
# ------------------------------------------------------------------------------

fca_raw <- readr::read_csv(
  fca_file,
  show_col_types = FALSE
)

fca_required <- c(
  "Gene",
  "Tissue",
  "N_Cells",
  "N_Expressing"
)

fca_missing <- setdiff(
  fca_required,
  colnames(fca_raw)
)

if (length(fca_missing) > 0) {

  stop(
    paste0(
      "FCA table is missing: ",
      paste(
        fca_missing,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------------------------
# 12. FCA: CELL TYPES -> TISSUES -> SYSTEMS
#
# Testis excluded.
# System-level percent expressing is cell-count weighted:
#
#   100 * sum(N_Expressing) / sum(N_Cells)
# ------------------------------------------------------------------------------

fca_tissue <- fca_raw %>%

  dplyr::filter(
    Gene %in% fly_serpins,
    Tissue != "testis"
  ) %>%

  dplyr::group_by(
    Gene,
    Tissue
  ) %>%

  dplyr::summarise(
    N_Cells = sum(
      N_Cells,
      na.rm = TRUE
    ),
    N_Expressing = sum(
      N_Expressing,
      na.rm = TRUE
    ),
    Percent_Expressing =
      100 *
      N_Expressing /
      N_Cells,
    .groups = "drop"
  )


# Diagnostic: FCA tissues not in mapping
fca_unmapped <- fca_tissue %>%

  dplyr::anti_join(
    fca_mapping,
    by = "Tissue"
  ) %>%

  dplyr::distinct(
    Tissue
  ) %>%

  dplyr::arrange(
    Tissue
  )

readr::write_csv(
  fca_unmapped,
  file.path(
    table_dir,
    "FCA_Unmapped_Tissues_Diagnostic.csv"
  )
)


fca_system <- fca_tissue %>%

  dplyr::inner_join(
    fca_mapping,
    by = "Tissue"
  ) %>%

  dplyr::group_by(
    Gene,
    System
  ) %>%

  dplyr::summarise(
    N_Cells = sum(
      N_Cells,
      na.rm = TRUE
    ),
    N_Expressing = sum(
      N_Expressing,
      na.rm = TRUE
    ),
    Expression =
      100 *
      N_Expressing /
      N_Cells,
    .groups = "drop"
  ) %>%

  dplyr::mutate(
    Fly_Dataset = "FCA",
    Fly_Metric = "Percent expressing"
  )

readr::write_csv(
  fca_system,
  file.path(
    table_dir,
    "FCA_SERPIN_SystemExpression.csv"
  )
)


# ------------------------------------------------------------------------------
# 13. READ + SUMMARIZE FLYATLAS2
# ------------------------------------------------------------------------------

flyatlas2_raw <- readr::read_csv(
  flyatlas2_file,
  show_col_types = FALSE
)

flyatlas2_std <- standardize_flyatlas2(
  flyatlas2_raw
) %>%

  dplyr::filter(
    Gene %in% fly_serpins
  ) %>%

  dplyr::mutate(
    System =
      classify_flyatlas2_tissue(
        Tissue
      )
  )


flyatlas2_unmapped <- flyatlas2_std %>%

  dplyr::filter(
    is.na(System)
  ) %>%

  dplyr::distinct(
    Tissue
  ) %>%

  dplyr::arrange(
    Tissue
  )

readr::write_csv(
  flyatlas2_unmapped,
  file.path(
    table_dir,
    "FlyAtlas2_Unmapped_Tissues_Diagnostic.csv"
  )
)


# Several source tissues can map to one system.
# Median bulk expression is used across source tissues.
flyatlas2_system <- flyatlas2_std %>%

  dplyr::filter(
    !is.na(System)
  ) %>%

  dplyr::group_by(
    Gene,
    System
  ) %>%

  dplyr::summarise(
    Expression =
      stats::median(
        FlyAtlas2_Expression,
        na.rm = TRUE
      ),
    N_Source_Tissues =
      dplyr::n(),
    .groups = "drop"
  ) %>%

  dplyr::mutate(
    Fly_Dataset = "FlyAtlas2",
    Fly_Metric = "Median FPKM"
  )

readr::write_csv(
  flyatlas2_system,
  file.path(
    table_dir,
    "FlyAtlas2_SERPIN_SystemExpression.csv"
  )
)


# ------------------------------------------------------------------------------
# 14. READ GTEx + HPA
# ------------------------------------------------------------------------------

gtex_raw <- readr::read_csv(
  gtex_file,
  show_col_types = FALSE
)

hpa_raw <- readr::read_csv(
  hpa_file,
  show_col_types = FALSE
)

gtex_std <- standardize_human(
  gtex_raw,
  "GTEx"
)

hpa_std <- standardize_human(
  hpa_raw,
  "HPA"
)


# ------------------------------------------------------------------------------
# 15. SUMMARIZE HUMAN TISSUES -> SYSTEMS
#
# Multiple human tissues assigned to the same system are summarized by
# median TPM. GTEx and HPA remain separate.
# ------------------------------------------------------------------------------

summarize_human_systems <- function(
    df,
    human_dataset
) {

  mapped <- df %>%

    dplyr::filter(
      Gene %in% human_serpins
    ) %>%

    dplyr::mutate(
      System =
        classify_human_tissue(
          Tissue
        )
    )


  unmapped <- mapped %>%

    dplyr::filter(
      is.na(System)
    ) %>%

    dplyr::distinct(
      Tissue
    ) %>%

    dplyr::arrange(
      Tissue
    )


  readr::write_csv(
    unmapped,
    file.path(
      table_dir,
      paste0(
        human_dataset,
        "_Unmapped_Tissues_Diagnostic.csv"
      )
    )
  )


  mapped %>%

    dplyr::filter(
      !is.na(System)
    ) %>%

    dplyr::group_by(
      Gene,
      System
    ) %>%

    dplyr::summarise(
      Expression =
        stats::median(
          Human_Expression,
          na.rm = TRUE
        ),
      N_Source_Tissues =
        dplyr::n(),
      .groups = "drop"
    ) %>%

    dplyr::mutate(
      Human_Dataset =
        human_dataset,
      Human_Metric =
        "Median TPM"
    )
}


gtex_system <- summarize_human_systems(
  gtex_std,
  "GTEx"
)

hpa_system <- summarize_human_systems(
  hpa_std,
  "HPA"
)

readr::write_csv(
  gtex_system,
  file.path(
    table_dir,
    "GTEx_SERPIN_SystemExpression.csv"
  )
)

readr::write_csv(
  hpa_system,
  file.path(
    table_dir,
    "HPA_SERPIN_SystemExpression.csv"
  )
)


# ------------------------------------------------------------------------------
# 16. SYSTEM COVERAGE DIAGNOSTIC
# ------------------------------------------------------------------------------

coverage_table <- dplyr::bind_rows(

  flyatlas2_system %>%
    dplyr::transmute(
      Dataset = "FlyAtlas2",
      Gene,
      System
    ),

  fca_system %>%
    dplyr::transmute(
      Dataset = "FCA",
      Gene,
      System
    ),

  gtex_system %>%
    dplyr::transmute(
      Dataset = "GTEx",
      Gene,
      System
    ),

  hpa_system %>%
    dplyr::transmute(
      Dataset = "HPA",
      Gene,
      System
    )

) %>%

  dplyr::distinct() %>%

  dplyr::mutate(
    System =
      factor(
        System,
        levels = system_order
      )
  ) %>%

  dplyr::arrange(
    Dataset,
    Gene,
    System
  )


readr::write_csv(
  coverage_table,
  file.path(
    table_dir,
    "FourDataset_System_Coverage.csv"
  )
)


# ------------------------------------------------------------------------------
# 17. PAIRWISE COMPARISON FUNCTION
#
# Each fly gene is compared with each human gene using systems available in both.
# Minimum matched systems = 3.
# ------------------------------------------------------------------------------

compare_pair <- function(
    fly_df,
    human_df,
    fly_dataset_name,
    human_dataset_name
) {

  pairs <- tidyr::crossing(
    Fly_Gene = fly_serpins,
    Human_Gene = human_serpins
  )

  results <- vector(
    "list",
    nrow(pairs)
  )

  matched_rows <- vector(
    "list",
    nrow(pairs)
  )

  for (i in seq_len(nrow(pairs))) {

    fg <- pairs$Fly_Gene[i]
    hg <- pairs$Human_Gene[i]


    fly_one <- fly_df %>%

      dplyr::filter(
        Gene == fg
      ) %>%

      dplyr::select(
        System,
        Fly_Expression = Expression
      )


    human_one <- human_df %>%

      dplyr::filter(
        Gene == hg
      ) %>%

      dplyr::select(
        System,
        Human_Expression = Expression
      )


    matched <- fly_one %>%

      dplyr::inner_join(
        human_one,
        by = "System"
      ) %>%

      dplyr::filter(
        is.finite(Fly_Expression),
        is.finite(Human_Expression)
      ) %>%

      dplyr::mutate(
        Fly_Gene = fg,
        Human_Gene = hg,
        Fly_Dataset = fly_dataset_name,
        Human_Dataset = human_dataset_name,
        Comparison = paste0(
          fly_dataset_name,
          " vs ",
          human_dataset_name
        )
      )


    n_matched <- nrow(matched)


    if (
      n_matched >= 3 &&
      dplyr::n_distinct(
        matched$Fly_Expression
      ) > 1 &&
      dplyr::n_distinct(
        matched$Human_Expression
      ) > 1
    ) {

      rho <- suppressWarnings(
        stats::cor(
          matched$Fly_Expression,
          matched$Human_Expression,
          method = "spearman",
          use = "complete.obs"
        )
      )

      p_value <- suppressWarnings(
        stats::cor.test(
          matched$Fly_Expression,
          matched$Human_Expression,
          method = "spearman",
          exact = FALSE
        )$p.value
      )

    } else {

      rho <- NA_real_
      p_value <- NA_real_
    }


    results[[i]] <- tibble::tibble(
      Fly_Dataset = fly_dataset_name,
      Human_Dataset = human_dataset_name,
      Comparison = paste0(
        fly_dataset_name,
        " vs ",
        human_dataset_name
      ),
      Fly_Gene = fg,
      Human_Gene = hg,
      N_Matched_Systems = n_matched,
      Matched_Systems = paste(
        matched$System,
        collapse = "; "
      ),
      Spearman_rho = rho,
      P_Value = p_value
    )


    matched_rows[[i]] <- matched
  }


  list(
    correlations =
      dplyr::bind_rows(results),
    matched_data =
      dplyr::bind_rows(matched_rows)
  )
}


# ------------------------------------------------------------------------------
# 18. RUN FOUR INDEPENDENT COMPARISONS
# ------------------------------------------------------------------------------

fa2_gtex <- compare_pair(
  fly_df = flyatlas2_system,
  human_df = gtex_system,
  fly_dataset_name = "FlyAtlas2",
  human_dataset_name = "GTEx"
)

fa2_hpa <- compare_pair(
  fly_df = flyatlas2_system,
  human_df = hpa_system,
  fly_dataset_name = "FlyAtlas2",
  human_dataset_name = "HPA"
)

fca_gtex <- compare_pair(
  fly_df = fca_system,
  human_df = gtex_system,
  fly_dataset_name = "FCA",
  human_dataset_name = "GTEx"
)

fca_hpa <- compare_pair(
  fly_df = fca_system,
  human_df = hpa_system,
  fly_dataset_name = "FCA",
  human_dataset_name = "HPA"
)


# ------------------------------------------------------------------------------
# 19. COMBINE CORRELATION RESULTS
#
# Expression values are NOT merged.
# FDR is calculated separately within each of the four comparisons.
# ------------------------------------------------------------------------------

all_correlations <- dplyr::bind_rows(
  fa2_gtex$correlations,
  fa2_hpa$correlations,
  fca_gtex$correlations,
  fca_hpa$correlations
) %>%

  dplyr::group_by(
    Comparison
  ) %>%

  dplyr::mutate(
    FDR =
      stats::p.adjust(
        P_Value,
        method = "BH"
      )
  ) %>%

  dplyr::ungroup()


readr::write_csv(
  all_correlations,
  file.path(
    table_dir,
    "FourDataset_CrossSpecies_SERPIN_Spearman_Correlations.csv"
  )
)


# ------------------------------------------------------------------------------
# 20. SAVE EXACT MATCHED-SYSTEM DATA USED
# ------------------------------------------------------------------------------

all_matched_data <- dplyr::bind_rows(
  fa2_gtex$matched_data,
  fa2_hpa$matched_data,
  fca_gtex$matched_data,
  fca_hpa$matched_data
)

readr::write_csv(
  all_matched_data,
  file.path(
    table_dir,
    "FourDataset_CrossSpecies_MatchedSystem_Data.csv"
  )
)


# ------------------------------------------------------------------------------
# 21. RANK SIMILARITIES WITHIN EACH COMPARISON
# ------------------------------------------------------------------------------

ranked_similarity <- all_correlations %>%

  dplyr::arrange(
    Comparison,
    dplyr::desc(
      Spearman_rho
    )
  ) %>%

  dplyr::group_by(
    Comparison
  ) %>%

  dplyr::mutate(
    Similarity_Rank =
      dplyr::row_number()
  ) %>%

  dplyr::ungroup()


readr::write_csv(
  ranked_similarity,
  file.path(
    table_dir,
    "FourDataset_CrossSpecies_SERPIN_Ranked_Similarities.csv"
  )
)


# ------------------------------------------------------------------------------
# 22. CONSENSUS SUMMARY
#
# This summarizes agreement among the FOUR independent correlations.
# Raw expression values remain separate.
# ------------------------------------------------------------------------------

consensus_summary <- all_correlations %>%

  dplyr::group_by(
    Fly_Gene,
    Human_Gene
  ) %>%

  dplyr::summarise(

    N_Comparisons_Available =
      sum(
        is.finite(
          Spearman_rho
        )
      ),

    Median_Rho =
      ifelse(
        any(
          is.finite(
            Spearman_rho
          )
        ),
        stats::median(
          Spearman_rho,
          na.rm = TRUE
        ),
        NA_real_
      ),

    Mean_Rho =
      ifelse(
        any(
          is.finite(
            Spearman_rho
          )
        ),
        mean(
          Spearman_rho,
          na.rm = TRUE
        ),
        NA_real_
      ),

    Minimum_Rho =
      ifelse(
        any(
          is.finite(
            Spearman_rho
          )
        ),
        min(
          Spearman_rho,
          na.rm = TRUE
        ),
        NA_real_
      ),

    Maximum_Rho =
      ifelse(
        any(
          is.finite(
            Spearman_rho
          )
        ),
        max(
          Spearman_rho,
          na.rm = TRUE
        ),
        NA_real_
      ),

    N_Positive =
      sum(
        Spearman_rho > 0,
        na.rm = TRUE
      ),

    N_Rho_GE_0_5 =
      sum(
        Spearman_rho >= 0.5,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%

  dplyr::arrange(
    dplyr::desc(
      Median_Rho
    ),
    dplyr::desc(
      N_Comparisons_Available
    )
  )


readr::write_csv(
  consensus_summary,
  file.path(
    table_dir,
    "FourDataset_CrossSpecies_SERPIN_Consensus_Summary.csv"
  )
)


# ------------------------------------------------------------------------------
# 23. HEATMAP FUNCTION
# ------------------------------------------------------------------------------

make_heatmap <- function(
    df,
    comparison_name
) {

  plot_df <- df %>%

    dplyr::filter(
      Comparison == comparison_name
    ) %>%

    dplyr::mutate(

      Fly_Gene =
        factor(
          Fly_Gene,
          levels = rev(
            fly_serpins
          )
        ),

      Human_Gene =
        factor(
          Human_Gene,
          levels =
            human_serpins
        )
    )


  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = Human_Gene,
      y = Fly_Gene,
      fill = Spearman_rho
    )
  ) +

    ggplot2::geom_tile(
      linewidth = 0.5
    ) +

    ggplot2::geom_text(
      ggplot2::aes(
        label =
          ifelse(
            is.na(
              Spearman_rho
            ),
            "",
            sprintf(
              "%.2f",
              Spearman_rho
            )
          )
      ),
      size = 3.5
    ) +

    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      midpoint = 0,
      name = "Spearman\nrho"
    ) +

    ggplot2::labs(
      title =
        paste0(
          comparison_name,
          ": cross-species tissue-expression similarity"
        ),
      subtitle =
        "Nine-system framework; systems available in both datasets are used",
      x = "Human serpin",
      y = "Drosophila serpin"
    ) +

    ggplot2::theme_bw(
      base_size = 11
    ) +

    ggplot2::theme(

      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          size = 10
        ),

      axis.text.y =
        ggplot2::element_text(
          size = 10
        ),

      axis.title =
        ggplot2::element_text(
          size = 11,
          face = "bold"
        ),

      plot.title =
        ggplot2::element_text(
          size = 14,
          face = "bold"
        ),

      plot.subtitle =
        ggplot2::element_text(
          size = 10
        ),

      legend.title =
        ggplot2::element_text(
          size = 10
        ),

      legend.text =
        ggplot2::element_text(
          size = 9
        )
    )
}


# ------------------------------------------------------------------------------
# 24. CREATE FOUR HEATMAPS
# ------------------------------------------------------------------------------

comparison_names <- c(
  "FlyAtlas2 vs GTEx",
  "FlyAtlas2 vs HPA",
  "FCA vs GTEx",
  "FCA vs HPA"
)


for (
  comparison_name in
    comparison_names
) {

  p <- make_heatmap(
    all_correlations,
    comparison_name
  )

  print(p)

  safe_name <- comparison_name %>%
    stringr::str_replace_all(
      " ",
      "_"
    )

  ggplot2::ggsave(
    filename =
      file.path(
        figure_dir,
        paste0(
          safe_name,
          "_SERPIN_TissueSimilarity_Heatmap.png"
        )
      ),
    plot = p,
    width = 9,
    height = 7,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}


# ------------------------------------------------------------------------------
# 25. CONSENSUS HEATMAP
#
# Displays median rho across the independent comparisons for each fly-human pair.
# This summarizes correlation evidence only; it does not combine raw expression.
# ------------------------------------------------------------------------------

consensus_plot_data <- consensus_summary %>%

  dplyr::mutate(

    Fly_Gene =
      factor(
        Fly_Gene,
        levels =
          rev(
            fly_serpins
          )
      ),

    Human_Gene =
      factor(
        Human_Gene,
        levels =
          human_serpins
      )
  )


consensus_heatmap <- ggplot2::ggplot(
  consensus_plot_data,
  ggplot2::aes(
    x = Human_Gene,
    y = Fly_Gene,
    fill = Median_Rho
  )
) +

  ggplot2::geom_tile(
    linewidth = 0.5
  ) +

  ggplot2::geom_text(
    ggplot2::aes(
      label =
        ifelse(
          is.na(
            Median_Rho
          ),
          "",
          sprintf(
            "%.2f",
            Median_Rho
          )
        )
    ),
    size = 3.5
  ) +

  ggplot2::scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0,
    name = "Median\nSpearman rho"
  ) +

  ggplot2::labs(
    title =
      "Consensus cross-species SERPIN tissue-expression similarity",
    subtitle =
      "Median Spearman rho across available FlyAtlas2/FCA × GTEx/HPA comparisons",
    x =
      "Human serpin",
    y =
      "Drosophila serpin"
  ) +

  ggplot2::theme_bw(
    base_size = 11
  ) +

  ggplot2::theme(

    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 10
      ),

    axis.text.y =
      ggplot2::element_text(
        size = 10
      ),

    axis.title =
      ggplot2::element_text(
        size = 11,
        face = "bold"
      ),

    plot.title =
      ggplot2::element_text(
        size = 14,
        face = "bold"
      ),

    plot.subtitle =
      ggplot2::element_text(
        size = 10
      ),

    legend.title =
      ggplot2::element_text(
        size = 10
      ),

    legend.text =
      ggplot2::element_text(
        size = 9
      )
  )


print(
  consensus_heatmap
)


ggplot2::ggsave(
  filename =
    file.path(
      figure_dir,
      "Consensus_FourDataset_SERPIN_TissueSimilarity_Heatmap.png"
    ),
  plot =
    consensus_heatmap,
  width =
    9,
  height =
    7,
  units =
    "in",
  dpi =
    600,
  bg =
    "white"
)


# ------------------------------------------------------------------------------
# 26. FINAL OUTPUT SUMMARY
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n",
  "FOUR-DATASET CROSS-SPECIES ANALYSIS COMPLETE\n",
  "============================================================\n"
)

cat(
  "\nInput files used:\n"
)

cat(
  "\nFCA:\n",
  fca_file,
  "\n"
)

cat(
  "\nFlyAtlas2:\n",
  flyatlas2_file,
  "\n"
)

cat(
  "\nGTEx:\n",
  gtex_file,
  "\n"
)

cat(
  "\nHPA:\n",
  hpa_file,
  "\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

cat(
  "\nFour independent comparisons:\n",
  "  1. FlyAtlas2 vs GTEx\n",
  "  2. FlyAtlas2 vs HPA\n",
  "  3. FCA vs GTEx\n",
  "  4. FCA vs HPA\n"
)

cat(
  "\nFCA testis excluded.\n"
)

cat(
  "\nFinal systems:\n"
)

print(
  system_order
)

cat(
  "\nTop consensus similarities:\n"
)

print(
  consensus_summary,
  n = Inf,
  width = Inf
)
