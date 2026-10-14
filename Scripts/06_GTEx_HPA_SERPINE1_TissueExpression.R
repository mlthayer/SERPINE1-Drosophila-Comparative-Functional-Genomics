# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 06_GTEx_HPA_SERPINE1_TissueExpression.R
# Run from the repository root.
# ============================================================
# 06_GTEx_HPA_Human_SERPIN_TissueExpression_REANALYSIS.R
#
# Human tissue-expression reanalysis for:
#   SERPINE1
#   SERPINI1
#   SERPINB3
#   SERPINB4
#   SERPINB11
#   SERPINB9
#
# GTEx v11:
#   Gene-level TPM from tissue-specific GCT files
#
# HPA:
#   Transcript-level TPM from transcript_rna_tissue.tsv.zip
#   Transcript TPM values are summed within each sample to
#   reconstruct gene-level TPM, then summarized by tissue.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(readr)
library(tibble)


# ============================================================
# 2. FILE PATHS
# ============================================================

# Set options(serpine1.choose_files = TRUE) in R/RStudio to select downloads.
# Missing inputs also offer file selection in interactive sessions.
choose_downloads <- function(label, url, multiple = FALSE) {
  message(label, " download page: ", url)
  if (!interactive()) {
    stop("Download the inputs and configure their input directory before running Rscript.")
  }
  if (.Platform$OS.type == "windows") {
    paths <- utils::choose.files(caption = paste("Select", label), multi = multiple)
  } else {
    paths <- character()
    repeat {
      paths <- c(paths, file.choose())
      if (!multiple || tolower(trimws(readline("Select another tissue file? [y/N]: "))) != "y") break
    }
  }
  if (!length(paths) || any(!file.exists(paths))) stop("No valid input files selected.")
  unique(paths)
}

gtex_dir <- Sys.getenv("SERPINE1_GTEX_INPUT", unset = "Expression/GTEx")

hpa_dir <- Sys.getenv("SERPINE1_HPA_INPUT", unset = "Expression/HPA")

output_dir <- file.path("Expression", "Human_SERPIN_reanalysis")

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
  recursive = TRUE,
  showWarnings = FALSE
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


# ============================================================
# 3. HUMAN SERPIN GENE SET
# ============================================================

human_serpins <- tribble(
  ~gene,       ~ensgid,
  "SERPINE1",  "ENSG00000106366",
  "SERPINI1",  "ENSG00000163536",
  "SERPINB3",  "ENSG00000057149",
  "SERPINB4",  "ENSG00000206073",
  "SERPINB11", "ENSG00000206072",
  "SERPINB9",  "ENSG00000170542"
)

cat("\nHuman SERPIN genes included:\n")
print(human_serpins)


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

safe_mean <- function(x) {

  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(
    x,
    na.rm = TRUE
  )
}


safe_median <- function(x) {

  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  median(
    x,
    na.rm = TRUE
  )
}


safe_sd <- function(x) {

  x <- as.numeric(x)

  if (sum(!is.na(x)) < 2) {
    return(NA_real_)
  }

  sd(
    x,
    na.rm = TRUE
  )
}


safe_iqr <- function(x) {

  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  IQR(
    x,
    na.rm = TRUE
  )
}


safe_min <- function(x) {

  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  min(
    x,
    na.rm = TRUE
  )
}


safe_max <- function(x) {

  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  max(
    x,
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# Clean GTEx filename into readable tissue name
# ------------------------------------------------------------

clean_gtex_tissue_name <- function(filename) {

  filename %>%
    basename() %>%
    str_remove("\\.gz$") %>%
    str_remove("\\.gct$") %>%
    str_remove("^gene_tpm_adult_gtex_v11_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}


# ============================================================
# PART A. GTEx MULTI-GENE ANALYSIS
# ============================================================


# ============================================================
# 5. IDENTIFY GTEx FILES
# ============================================================

gtex_files <- list.files(
  path = gtex_dir,
  pattern = "\\.gct(\\.gz)?$",
  full.names = TRUE
)

if (isTRUE(getOption("serpine1.choose_files", FALSE)) || length(gtex_files) == 0) {
  gtex_files <- choose_downloads(
    "GTEx v11 gene-level TPM GCT files (one per tissue)",
    "https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression",
    multiple = TRUE
  )
  if (any(!grepl("\\.gct(\\.gz)?$", gtex_files))) stop("Select .gct or .gct.gz files.")
}

cat(
  "\nNumber of GTEx files detected:",
  length(gtex_files),
  "\n\n"
)

print(
  basename(gtex_files)
)


# ============================================================
# 6. FUNCTION TO EXTRACT ALL HUMAN SERPINS FROM ONE GTEx FILE
# ============================================================

read_serpins_from_gtex <- function(file) {

  tissue_name <- clean_gtex_tissue_name(file)

  cat(
    "Reading GTEx tissue:",
    tissue_name,
    "\n"
  )

  dat <- read.delim(
    file,
    skip = 2,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(dat) < 3) {

    warning(
      "File ",
      basename(file),
      " has fewer than 3 columns and will be skipped."
    )

    return(NULL)
  }


  # ----------------------------------------------------------
  # Identify gene-symbol column
  # ----------------------------------------------------------

  possible_gene_columns <- c(
    "Description",
    "description",
    "gene_name",
    "Gene",
    "gene",
    "symbol",
    "Symbol"
  )

  gene_col <- intersect(
    possible_gene_columns,
    names(dat)
  )

  if (length(gene_col) == 0) {

    stop(
      "Could not identify a gene-symbol column in file: ",
      basename(file),
      "\nAvailable columns include: ",
      paste(
        head(
          names(dat),
          15
        ),
        collapse = ", "
      )
    )
  }

  gene_col <- gene_col[1]


  # ----------------------------------------------------------
  # Annotation columns that should NOT be treated as samples
  # ----------------------------------------------------------

  annotation_cols <- c(
    "Name",
    "name",
    "Description",
    "description",
    "gene_name",
    "Gene",
    "gene",
    "symbol",
    "Symbol"
  )

  sample_cols <- setdiff(
    names(dat),
    annotation_cols
  )


  # ----------------------------------------------------------
  # Extract each target gene
  # ----------------------------------------------------------

  result <- map_dfr(
    human_serpins$gene,
    function(target_gene) {

      gene_row <- dat %>%
        filter(
          .data[[gene_col]] == target_gene
        )

      if (nrow(gene_row) == 0) {

        warning(
          target_gene,
          " was not found in ",
          basename(file)
        )

        return(
          tibble(
            gene = target_gene,
            tissue = tissue_name,
            n_samples = 0,
            median_TPM = NA_real_,
            mean_TPM = NA_real_,
            sd_TPM = NA_real_,
            IQR_TPM = NA_real_,
            min_TPM = NA_real_,
            max_TPM = NA_real_
          )
        )
      }


      if (nrow(gene_row) > 1) {

        warning(
          "More than one row was found for ",
          target_gene,
          " in ",
          basename(file),
          ". Using the first matching gene row."
        )

        gene_row <- gene_row[
          1,
          ,
          drop = FALSE
        ]
      }


      expr_values <- suppressWarnings(
        as.numeric(
          unlist(
            gene_row[
              1,
              sample_cols
            ],
            use.names = FALSE
          )
        )
      )


      expr_values <- expr_values[
        !is.na(expr_values)
      ]


      if (length(expr_values) == 0) {

        return(
          tibble(
            gene = target_gene,
            tissue = tissue_name,
            n_samples = 0,
            median_TPM = NA_real_,
            mean_TPM = NA_real_,
            sd_TPM = NA_real_,
            IQR_TPM = NA_real_,
            min_TPM = NA_real_,
            max_TPM = NA_real_
          )
        )
      }


      if (
        any(
          expr_values < 0,
          na.rm = TRUE
        )
      ) {

        stop(
          "Negative TPM values detected for ",
          target_gene,
          " in ",
          basename(file)
        )
      }


      tibble(
        gene = target_gene,
        tissue = tissue_name,
        n_samples = length(expr_values),

        median_TPM =
          safe_median(expr_values),

        mean_TPM =
          safe_mean(expr_values),

        sd_TPM =
          safe_sd(expr_values),

        IQR_TPM =
          safe_iqr(expr_values),

        min_TPM =
          safe_min(expr_values),

        max_TPM =
          safe_max(expr_values)
      )
    }
  )

  result
}


# ============================================================
# 7. RUN GTEx EXTRACTION
# ============================================================

gtex_serpins <- map_dfr(
  gtex_files,
  read_serpins_from_gtex
) %>%
  mutate(

    log2_median_TPM =
      log2(
        median_TPM + 1
      ),

    log2_mean_TPM =
      log2(
        mean_TPM + 1
      )
  ) %>%
  arrange(
    gene,
    desc(median_TPM)
  )


cat(
  "\n========================================\n"
)

cat(
  "GTEx HUMAN SERPIN TISSUE SUMMARY\n"
)

cat(
  "========================================\n"
)

print(
  tibble::as_tibble(gtex_serpins),
  n = Inf
)


# ============================================================
# 8. SAVE GTEx TABLE
# ============================================================

write.csv(
  gtex_serpins,
  file.path(
    table_dir,
    "GTEx_v11_human_SERPIN_tissue_expression_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. GTEx FACETED EXPRESSION PLOT
# ============================================================

gtex_plot_data <- gtex_serpins %>%
  filter(
    !is.na(log2_median_TPM)
  )


p_gtex <- ggplot(
  gtex_plot_data,
  aes(
    x = log2_median_TPM,
    y = reorder(
      tissue,
      log2_median_TPM
    )
  )
) +
  geom_col(
    width = 0.7
  ) +
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 2
  ) +
  theme_bw(
    base_size = 10
  ) +
  labs(
    title =
      "GTEx v11 human SERPIN tissue expression",

    x =
      expression(
        log[2]~"(median TPM + 1)"
      ),

    y = NULL
  ) +
  theme(

    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),

    strip.text =
      element_text(
        face = "bold"
      ),

    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank()
  )


print(
  p_gtex
)


ggsave(
  filename = file.path(
    figure_dir,
    "GTEx_v11_human_SERPIN_tissue_expression.png"
  ),
  plot = p_gtex,
  width = 11,
  height = 13,
  dpi = 600,
  bg = "white"
)


# ============================================================
# PART B. HPA MULTI-GENE ANALYSIS
# ============================================================


# ============================================================
# 10. READ HPA TRANSCRIPT RNA DATA
# ============================================================

# CLI mode validates and exports GTEx independently of the HPA input.
if ("--gtex-only" %in% commandArgs(trailingOnly = TRUE)) {
  message("GTEx-only run completed for ", length(gtex_files),
          " available tissue files. HPA analysis was not run.")
  quit(save = "no", status = 0L)
}

hpa_zip <- file.path(
  hpa_dir,
  "transcript_rna_tissue.tsv.zip"
)


if (isTRUE(getOption("serpine1.choose_files", FALSE)) || !file.exists(hpa_zip)) {
  hpa_zip <- choose_downloads(
    "HPA transcript_rna_tissue.tsv.zip",
    "https://www.proteinatlas.org/about/download"
  )
}


zip_contents <- unzip(
  hpa_zip,
  list = TRUE
)


hpa_internal_file <- zip_contents$Name[
  grepl(
    "transcript_rna_tissue\\.tsv$",
    zip_contents$Name
  )
][1]


if (is.na(hpa_internal_file)) {

  stop(
    "Could not find transcript_rna_tissue.tsv inside the HPA ZIP archive."
  )
}


hpa_raw <- read.delim(
  unz(
    hpa_zip,
    hpa_internal_file
  ),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "\nHPA data dimensions:\n"
)

print(
  dim(hpa_raw)
)


cat(
  "\nFirst HPA columns:\n"
)

print(
  head(
    names(hpa_raw),
    20
  )
)


# ============================================================
# 11. CHECK REQUIRED HPA IDENTIFIER COLUMNS
# ============================================================

required_hpa_columns <- c(
  "ensgid",
  "enstid"
)


missing_hpa_columns <- setdiff(
  required_hpa_columns,
  names(hpa_raw)
)


if (
  length(
    missing_hpa_columns
  ) > 0
) {

  stop(
    "HPA file is missing required column(s): ",
    paste(
      missing_hpa_columns,
      collapse = ", "
    )
  )
}


# ============================================================
# 12. EXTRACT TARGET HUMAN SERPIN TRANSCRIPTS
# ============================================================

hpa_target <- hpa_raw %>%
  filter(
    ensgid %in%
      human_serpins$ensgid
  ) %>%
  left_join(
    human_serpins,
    by = "ensgid"
  )


cat(
  "\nNumber of HPA transcript rows recovered for each gene:\n"
)


hpa_transcript_counts <- hpa_target %>%
  count(
    gene,
    name = "n_transcript_rows"
  )


print(
  hpa_transcript_counts
)


# ------------------------------------------------------------
# Warn if one of the requested genes was not found
# ------------------------------------------------------------

missing_hpa_genes <- setdiff(
  human_serpins$gene,
  hpa_transcript_counts$gene
)


if (
  length(
    missing_hpa_genes
  ) > 0
) {

  warning(
    "The following genes were not found in the HPA transcript file: ",
    paste(
      missing_hpa_genes,
      collapse = ", "
    )
  )
}


# ============================================================
# 13. IDENTIFY TPM COLUMNS
# ============================================================

tpm_cols <- names(hpa_target)[
  str_detect(
    names(hpa_target),
    "^TPM\\."
  )
]


if (
  length(
    tpm_cols
  ) == 0
) {

  stop(
    "No TPM columns were detected in the HPA transcript file."
  )
}


cat(
  "\nNumber of HPA TPM sample columns:",
  length(tpm_cols),
  "\n"
)


# ============================================================
# 14. CONVERT HPA TRANSCRIPT DATA TO LONG FORMAT
# ============================================================

hpa_transcript_long <- hpa_target %>%
  select(
    gene,
    ensgid,
    enstid,
    all_of(
      tpm_cols
    )
  ) %>%
  pivot_longer(
    cols =
      all_of(
        tpm_cols
      ),

    names_to =
      "sample_column",

    values_to =
      "transcript_TPM"
  ) %>%
  mutate(

    transcript_TPM =
      as.numeric(
        transcript_TPM
      ),

    sample_info =
      str_remove(
        sample_column,
        "^TPM\\."
      ),

    sample_id =
      str_extract(
        sample_info,
        "[^.]+$"
      ),

    tissue =
      str_remove(
        sample_info,
        "\\.[^.]+$"
      )
  )


cat(
  "\nExample parsed HPA rows:\n"
)


print(
  head(
    hpa_transcript_long %>%
      select(
        gene,
        enstid,
        tissue,
        sample_id,
        transcript_TPM
      )
  )
)


# ============================================================
# 15. SUM TRANSCRIPT TPM TO GENE-LEVEL TPM WITHIN EACH SAMPLE
# ============================================================

hpa_sample_level <- hpa_transcript_long %>%
  group_by(
    gene,
    ensgid,
    tissue,
    sample_id
  ) %>%
  summarize(

    gene_TPM =
      sum(
        transcript_TPM,
        na.rm = TRUE
      ),

    n_transcripts =
      n_distinct(
        enstid
      ),

    .groups =
      "drop"
  )


cat(
  "\nHPA gene-level sample data:\n"
)


print(
  hpa_sample_level,
  n = 30
)


# ============================================================
# 16. SUMMARIZE HPA EXPRESSION BY GENE AND TISSUE
# ============================================================

hpa_serpins <- hpa_sample_level %>%
  group_by(
    gene,
    tissue
  ) %>%
  summarize(

    n_samples =
      n(),

    median_TPM =
      safe_median(
        gene_TPM
      ),

    mean_TPM =
      safe_mean(
        gene_TPM
      ),

    sd_TPM =
      safe_sd(
        gene_TPM
      ),

    IQR_TPM =
      safe_iqr(
        gene_TPM
      ),

    min_TPM =
      safe_min(
        gene_TPM
      ),

    max_TPM =
      safe_max(
        gene_TPM
      ),

    .groups =
      "drop"
  ) %>%
  mutate(

    log2_median_TPM =
      log2(
        median_TPM + 1
      ),

    log2_mean_TPM =
      log2(
        mean_TPM + 1
      )
  ) %>%
  arrange(
    gene,
    desc(median_TPM)
  )


cat(
  "\n========================================\n"
)

cat(
  "HPA HUMAN SERPIN TISSUE SUMMARY\n"
)

cat(
  "========================================\n"
)


print(
   tibble::as_tibble(hpa_serpins),
  n = Inf
)


# ============================================================
# 17. SAVE HPA TABLES
# ============================================================

write.csv(
  hpa_sample_level,
  file.path(
    table_dir,
    "HPA_human_SERPIN_sample_level_gene_TPM.csv"
  ),
  row.names = FALSE
)


write.csv(
  hpa_serpins,
  file.path(
    table_dir,
    "HPA_human_SERPIN_tissue_expression_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. HPA FACETED EXPRESSION PLOT
# ============================================================

p_hpa <- ggplot(
  hpa_serpins,
  aes(
    x = log2_median_TPM,
    y = reorder(
      tissue,
      log2_median_TPM
    )
  )
) +
  geom_col(
    width = 0.7
  ) +
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 2
  ) +
  theme_bw(
    base_size = 10
  ) +
  labs(

    title =
      "HPA human SERPIN tissue expression",

    x =
      expression(
        log[2]~"(median TPM + 1)"
      ),

    y = NULL
  ) +
  theme(

    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),

    strip.text =
      element_text(
        face = "bold"
      ),

    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank()
  )


print(
  p_hpa
)


ggsave(
  filename = file.path(
    figure_dir,
    "HPA_human_SERPIN_tissue_expression.png"
  ),
  plot = p_hpa,
  width = 11,
  height = 15,
  dpi = 600,
  bg = "white"
)


# ============================================================
# PART C. COMBINED LONG-FORM DATA
# ============================================================


# ============================================================
# 19. FORMAT GTEx FOR COMBINED DATASET
# ============================================================

gtex_combined <- gtex_serpins %>%
  transmute(

    database =
      "GTEx v11",

    gene,

    tissue,

    n_samples,

    median_TPM,

    mean_TPM,

    sd_TPM,

    IQR_TPM,

    min_TPM,

    max_TPM,

    log2_median_TPM,

    log2_mean_TPM
  )


# ============================================================
# 20. FORMAT HPA FOR COMBINED DATASET
# ============================================================

hpa_combined <- hpa_serpins %>%
  transmute(

    database =
      "HPA",

    gene,

    tissue,

    n_samples,

    median_TPM,

    mean_TPM,

    sd_TPM,

    IQR_TPM,

    min_TPM,

    max_TPM,

    log2_median_TPM,

    log2_mean_TPM
  )


# ============================================================
# 21. COMBINE WITHOUT AVERAGING DATABASES
# ============================================================

human_serpins_combined <- bind_rows(
  gtex_combined,
  hpa_combined
)


write.csv(
  human_serpins_combined,
  file.path(
    table_dir,
    "Human_SERPIN_GTEx_HPA_long_format.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART D. GENE × TISSUE MATRICES
# ============================================================


# ============================================================
# 22. GTEx GENE × TISSUE MATRIX
# ============================================================

gtex_matrix <- gtex_serpins %>%
  select(
    gene,
    tissue,
    log2_median_TPM
  ) %>%
  pivot_wider(
    names_from =
      tissue,

    values_from =
      log2_median_TPM
  )


write.csv(
  gtex_matrix,
  file.path(
    table_dir,
    "GTEx_human_SERPIN_gene_by_tissue_log2_matrix.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. HPA GENE × TISSUE MATRIX
# ============================================================

hpa_matrix <- hpa_serpins %>%
  select(
    gene,
    tissue,
    log2_median_TPM
  ) %>%
  pivot_wider(
    names_from =
      tissue,

    values_from =
      log2_median_TPM
  )


write.csv(
  hpa_matrix,
  file.path(
    table_dir,
    "HPA_human_SERPIN_gene_by_tissue_log2_matrix.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART E. WITHIN-GENE STANDARDIZATION
# ============================================================


# ============================================================
# 24. GTEx WITHIN-GENE Z-SCORES
# ============================================================

gtex_z <- gtex_serpins %>%
  group_by(
    gene
  ) %>%
  mutate(

    tissue_z =
      if (
        sum(
          !is.na(
            log2_median_TPM
          )
        ) > 1 &&
        sd(
          log2_median_TPM,
          na.rm = TRUE
        ) > 0
      ) {

        as.numeric(
          scale(
            log2_median_TPM
          )
        )

      } else {

        NA_real_
      }
  ) %>%
  ungroup()


# ============================================================
# 25. HPA WITHIN-GENE Z-SCORES
# ============================================================

hpa_z <- hpa_serpins %>%
  group_by(
    gene
  ) %>%
  mutate(

    tissue_z =
      if (
        sum(
          !is.na(
            log2_median_TPM
          )
        ) > 1 &&
        sd(
          log2_median_TPM,
          na.rm = TRUE
        ) > 0
      ) {

        as.numeric(
          scale(
            log2_median_TPM
          )
        )

      } else {

        NA_real_
      }
  ) %>%
  ungroup()


write.csv(
  gtex_z,
  file.path(
    table_dir,
    "GTEx_human_SERPIN_tissue_expression_zscores.csv"
  ),
  row.names = FALSE
)


write.csv(
  hpa_z,
  file.path(
    table_dir,
    "HPA_human_SERPIN_tissue_expression_zscores.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART F. TOP TISSUES FOR EACH GENE
# ============================================================


# ============================================================
# 26. TOP 10 GTEx TISSUES FOR EACH SERPIN
# ============================================================

gtex_top10 <- gtex_serpins %>%
  filter(
    !is.na(median_TPM)
  ) %>%
  group_by(
    gene
  ) %>%
  slice_max(
    order_by =
      median_TPM,

    n = 10,

    with_ties = FALSE
  ) %>%
  arrange(
    gene,
    desc(median_TPM)
  ) %>%
  ungroup()


cat(
  "\n========================================\n"
)

cat(
  "TOP GTEx TISSUES FOR EACH HUMAN SERPIN\n"
)

cat(
  "========================================\n"
)


print(
   tibble::as_tibble(gtex_top10 %>%
    select(
      gene,
      tissue,
      n_samples,
      median_TPM,
      mean_TPM,
      log2_median_TPM
    )),
  n = Inf
)


write.csv(
  gtex_top10,
  file.path(
    table_dir,
    "GTEx_human_SERPIN_top10_tissues.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 27. TOP 10 HPA TISSUES FOR EACH SERPIN
# ============================================================

hpa_top10 <- hpa_serpins %>%
  filter(
    !is.na(median_TPM)
  ) %>%
  group_by(
    gene
  ) %>%
  slice_max(
    order_by =
      median_TPM,

    n = 10,

    with_ties = FALSE
  ) %>%
  arrange(
    gene,
    desc(median_TPM)
  ) %>%
  ungroup()


cat(
  "\n========================================\n"
)

cat(
  "TOP HPA TISSUES FOR EACH HUMAN SERPIN\n"
)

cat(
  "========================================\n"
)


print(
   tibble::as_tibble(hpa_top10 %>%
    select(
      gene,
      tissue,
      n_samples,
      median_TPM,
      mean_TPM,
      log2_median_TPM
    )),
  n = Inf
)


write.csv(
  hpa_top10,
  file.path(
    table_dir,
    "HPA_human_SERPIN_top10_tissues.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART G. RECOVERY / QUALITY-CONTROL SUMMARY
# ============================================================


# ============================================================
# 28. GTEx GENE RECOVERY
# ============================================================

gtex_recovery <- gtex_serpins %>%
  group_by(
    gene
  ) %>%
  summarize(

    tissues_examined =
      n(),

    tissues_with_expression_data =
      sum(
        !is.na(
          median_TPM
        )
      ),

    .groups =
      "drop"
  )


cat(
  "\nGTEx gene recovery:\n"
)

print(
  gtex_recovery
)


# ============================================================
# 29. HPA GENE RECOVERY
# ============================================================

hpa_recovery <- hpa_serpins %>%
  group_by(
    gene
  ) %>%
  summarize(

    tissues_examined =
      n(),

    tissues_with_expression_data =
      sum(
        !is.na(
          median_TPM
        )
      ),

    .groups =
      "drop"
  )


cat(
  "\nHPA gene recovery:\n"
)

print(
  hpa_recovery
)


write.csv(
  gtex_recovery,
  file.path(
    table_dir,
    "GTEx_SERPIN_gene_recovery_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  hpa_recovery,
  file.path(
    table_dir,
    "HPA_SERPIN_gene_recovery_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 30. FINAL MESSAGE
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "HUMAN SERPIN REANALYSIS COMPLETE\n"
)

cat(
  "========================================\n"
)


cat(
  "\nTables saved to:\n",
  table_dir,
  "\n\nFigures saved to:\n",
  figure_dir,
  "\n"
)

# ============================================================
# PART H. GTEx-HPA MATCHED-TISSUE CONCORDANCE
# ============================================================


# ------------------------------------------------------------
# 31. Define conservative matched tissues
# ------------------------------------------------------------

tissue_mapping <- tibble::tibble(
  matched_tissue = c(
    "Adipose",
    "Adipose",
    "Heart",
    "Heart",
    "Liver",
    "Lung",
    "Ovary",
    "Pancreas",
    "Skeletal Muscle",
    "Small Intestine",
    "Testis"
  ),

  gtex_tissue = c(
    "Adipose Subcutaneous",
    "Adipose Visceral Omentum",
    "Heart Atrial Appendage",
    "Heart Left Ventricle",
    "Liver",
    "Lung",
    "Ovary",
    "Pancreas",
    "Muscle Skeletal",
    "Small Intestine Terminal Ileum",
    "Testis"
  ),

  hpa_tissue = c(
    "adipose tissue",
    "adipose tissue",
    "heart muscle",
    "heart muscle",
    "liver",
    "lung",
    "ovary",
    "pancreas",
    "skeletal muscle",
    "small intestine",
    "testis"
  )
)

cat("\nMatched tissue mapping:\n")

print(
  tibble::as_tibble(tissue_mapping),
  n = Inf
)


# ------------------------------------------------------------
# 32. Prepare GTEx matched tissues
# ------------------------------------------------------------

gtex_matched <- gtex_serpins %>%
  inner_join(
    tissue_mapping,
    by = c("tissue" = "gtex_tissue")
  ) %>%
  group_by(
    gene,
    matched_tissue
  ) %>%
  summarize(
    gtex_median_TPM = mean(
      median_TPM,
      na.rm = TRUE
    ),

    gtex_n_subtypes = n(),

    .groups = "drop"
  ) %>%
  mutate(
    gtex_log2_median_TPM =
      log2(gtex_median_TPM + 1)
  )


cat("\nGTEx matched tissue data:\n")

print(
  tibble::as_tibble(gtex_matched),
  n = Inf
)


# ------------------------------------------------------------
# 33. Prepare HPA matched tissues
# ------------------------------------------------------------

hpa_mapping <- tissue_mapping %>%
  select(
    matched_tissue,
    hpa_tissue
  ) %>%
  distinct()


hpa_matched <- hpa_serpins %>%
  inner_join(
    hpa_mapping,
    by = c("tissue" = "hpa_tissue")
  ) %>%
  group_by(
    gene,
    matched_tissue
  ) %>%
  summarize(
    hpa_median_TPM = mean(
      median_TPM,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    hpa_log2_median_TPM =
      log2(hpa_median_TPM + 1)
  )


cat("\nHPA matched tissue data:\n")

print(
  tibble::as_tibble(hpa_matched),
  n = Inf
)


# ------------------------------------------------------------
# 34. Join GTEx and HPA
# ------------------------------------------------------------

gtex_hpa_matched <- gtex_matched %>%
  inner_join(
    hpa_matched,
    by = c(
      "gene",
      "matched_tissue"
    )
  ) %>%
  arrange(
    gene,
    matched_tissue
  )


cat("\n========================================\n")
cat("GTEx-HPA MATCHED TISSUE TABLE\n")
cat("========================================\n")

print(
  tibble::as_tibble(gtex_hpa_matched),
  n = Inf
)


write.csv(
  gtex_hpa_matched,
  file.path(
    table_dir,
    "GTEx_HPA_SERPIN_matched_tissue_expression.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 35. Function for Pearson and Spearman correlations
# ------------------------------------------------------------

calculate_concordance <- function(df) {

  df <- df %>%
    filter(
      !is.na(gtex_log2_median_TPM),
      !is.na(hpa_log2_median_TPM)
    )

  n_matched <- nrow(df)

  if (n_matched < 3) {

    return(
      tibble::tibble(
        n_matched_tissues = n_matched,
        pearson_r = NA_real_,
        pearson_p = NA_real_,
        spearman_rho = NA_real_,
        spearman_p = NA_real_
      )
    )
  }

  pearson_test <- cor.test(
    df$gtex_log2_median_TPM,
    df$hpa_log2_median_TPM,
    method = "pearson"
  )

  spearman_test <- suppressWarnings(
    cor.test(
      df$gtex_log2_median_TPM,
      df$hpa_log2_median_TPM,
      method = "spearman",
      exact = FALSE
    )
  )

  tibble::tibble(
    n_matched_tissues = n_matched,

    pearson_r =
      unname(pearson_test$estimate),

    pearson_p =
      pearson_test$p.value,

    spearman_rho =
      unname(spearman_test$estimate),

    spearman_p =
      spearman_test$p.value
  )
}


# ------------------------------------------------------------
# 36. Calculate concordance for each SERPIN
# ------------------------------------------------------------

concordance_list <- split(
  gtex_hpa_matched,
  gtex_hpa_matched$gene
)


concordance_summary <- purrr::imap_dfr(
  concordance_list,

  function(df, gene_name) {

    result <- calculate_concordance(df)

    result$gene <- gene_name

    result
  }
) %>%
  select(
    gene,
    everything()
  ) %>%
  arrange(
    desc(pearson_r)
  )


cat("\n========================================\n")
cat("GTEx-HPA CONCORDANCE SUMMARY\n")
cat("========================================\n")

print(
  tibble::as_tibble(concordance_summary),
  n = Inf
)


write.csv(
  concordance_summary,
  file.path(
    table_dir,
    "GTEx_HPA_SERPIN_concordance_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 37. Calculate GTEx-HPA discrepancies
# ------------------------------------------------------------

gtex_hpa_discrepancy <- gtex_hpa_matched %>%
  mutate(
    log2_difference =
      gtex_log2_median_TPM -
      hpa_log2_median_TPM,

    absolute_log2_difference =
      abs(log2_difference)
  ) %>%
  arrange(
    gene,
    desc(absolute_log2_difference)
  )


write.csv(
  gtex_hpa_discrepancy,
  file.path(
    table_dir,
    "GTEx_HPA_SERPIN_tissue_discrepancies.csv"
  ),
  row.names = FALSE
)


cat("\nLargest GTEx-HPA discrepancies:\n")


largest_discrepancies <- gtex_hpa_discrepancy %>%
  group_by(gene) %>%
  slice_max(
    order_by = absolute_log2_difference,
    n = 3,
    with_ties = FALSE
  ) %>%
  ungroup()


largest_discrepancies_print <- largest_discrepancies %>%
  select(
    gene,
    matched_tissue,
    gtex_median_TPM,
    hpa_median_TPM,
    log2_difference
  )


print(
  tibble::as_tibble(largest_discrepancies_print),
  n = Inf
)


# ------------------------------------------------------------
# 38. Concordance scatter plots
# ------------------------------------------------------------

p_concordance <- ggplot(
  gtex_hpa_matched,
  aes(
    x = gtex_log2_median_TPM,
    y = hpa_log2_median_TPM
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  facet_wrap(
    ~ gene,
    ncol = 2
  ) +
  theme_bw(
    base_size = 10
  ) +
  labs(
    title =
      "GTEx-HPA concordance of human SERPIN tissue expression",

    x =
      "GTEx log2(median TPM + 1)",

    y =
      "HPA log2(median TPM + 1)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    strip.text = element_text(
      face = "bold"
    ),

    panel.grid.minor =
      element_blank()
  )


print(p_concordance)


ggsave(
  filename = file.path(
    figure_dir,
    "GTEx_HPA_SERPIN_concordance_scatterplots.png"
  ),

  plot = p_concordance,

  width = 10,
  height = 11,
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 39. Tissue-labelled scatter plots
# ------------------------------------------------------------

p_concordance_labeled <- ggplot(
  gtex_hpa_matched,
  aes(
    x = gtex_log2_median_TPM,
    y = hpa_log2_median_TPM,
    label = matched_tissue
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_text(
    size = 2.5,
    nudge_y = 0.12,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ gene,
    ncol = 2
  ) +
  theme_bw(
    base_size = 10
  ) +
  labs(
    title =
      "Matched GTEx and HPA tissue expression",

    x =
      "GTEx log2(median TPM + 1)",

    y =
      "HPA log2(median TPM + 1)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    strip.text = element_text(
      face = "bold"
    ),

    panel.grid.minor =
      element_blank()
  )


print(p_concordance_labeled)


ggsave(
  filename = file.path(
    figure_dir,
    "GTEx_HPA_SERPIN_concordance_labeled.png"
  ),

  plot = p_concordance_labeled,

  width = 11,
  height = 12,
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 40. Completion message
# ------------------------------------------------------------

cat("\n========================================\n")
cat("GTEx-HPA CONCORDANCE ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "\nCreated:\n",
  "GTEx_HPA_SERPIN_matched_tissue_expression.csv\n",
  "GTEx_HPA_SERPIN_concordance_summary.csv\n",
  "GTEx_HPA_SERPIN_tissue_discrepancies.csv\n",
  "GTEx_HPA_SERPIN_concordance_scatterplots.png\n",
  "GTEx_HPA_SERPIN_concordance_labeled.png\n"
)
