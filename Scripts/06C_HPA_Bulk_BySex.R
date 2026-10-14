# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Run from the repository root.
# ==============================================================================
# HPA BULK TISSUE SERPIN EXPRESSION BY SEX
#
# Purpose:
#   Build sex-specific HPA tissue-expression summaries for the six
#   human SERPIN genes from the previously generated sample-level table.
#
# Primary metric:
#   Median gene-level TPM within sex and tissue
#
# Input:
#   HPA_human_SERPIN_sample_level_gene_TPM.csv
#   rna_tissue_hpa_samples.tsv or .tsv.zip
#
# Outputs:
#   1. HPA sample-level SERPIN TPM with sex metadata
#   2. HPA tissue-expression summary by sex
#   3. HPA sample coverage by sex and tissue
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)


# ------------------------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------------------------

analysis_root <- "Expression/Human_SERPIN_reanalysis"

table_dir <- file.path(
  analysis_root,
  "Tables"
)

hpa_raw_dir <- Sys.getenv("SERPINE1_HPA_INPUT", unset = "Expression/HPA")


sample_level_file <- file.path(
  table_dir,
  "HPA_human_SERPIN_sample_level_gene_TPM.csv"
)


# ------------------------------------------------------------------------------
# 3. Locate HPA sample metadata
# ------------------------------------------------------------------------------

metadata_candidates <- list.files(
  hpa_raw_dir,
  pattern = "rna_tissue_hpa_samples.*\\.tsv(\\.zip)?$",
  full.names = TRUE,
  ignore.case = TRUE
)


if (
  length(metadata_candidates) == 0
) {

  stop(
    "Could not find HPA rna_tissue_hpa_samples TSV/ZIP file."
  )
}


metadata_file <- metadata_candidates[1]


cat(
  "\nUsing HPA metadata file:\n",
  metadata_file,
  "\n"
)


stopifnot(
  file.exists(sample_level_file)
)


# ------------------------------------------------------------------------------
# 4. Read sample-level expression
# ------------------------------------------------------------------------------

hpa_sample_level <- read_csv(
  sample_level_file,
  show_col_types = FALSE
)


cat("\nHPA sample-level columns:\n")

print(
  names(hpa_sample_level)
)


# ------------------------------------------------------------------------------
# 5. Read HPA sample metadata
# ------------------------------------------------------------------------------

if (
  str_detect(
    metadata_file,
    "\\.zip$"
  )
) {

  zip_contents <- unzip(
    metadata_file,
    list = TRUE
  )

  tsv_name <- zip_contents$Name[
    str_detect(
      zip_contents$Name,
      "\\.tsv$"
    )
  ][1]


  if (
    is.na(tsv_name)
  ) {

    stop(
      "No TSV file found inside HPA metadata ZIP."
    )
  }


  hpa_samples <- read.delim(
    unz(
      metadata_file,
      tsv_name
    ),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

} else {

  hpa_samples <- read.delim(
    metadata_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


cat("\nHPA sample metadata columns:\n")

print(
  names(hpa_samples)
)


# ------------------------------------------------------------------------------
# 6. Standardize metadata
# ------------------------------------------------------------------------------

hpa_samples_clean <- hpa_samples %>%

  transmute(

    sample_id =
      as.numeric(
        `HPA sample id`
      ),

    metadata_tissue =
      Tissue,

    Age =
      Age,

    Sex =
      Sex
  )


# ------------------------------------------------------------------------------
# 7. Standardize expression sample IDs
# ------------------------------------------------------------------------------

hpa_sample_level <- hpa_sample_level %>%

  mutate(
    sample_id =
      as.numeric(
        sample_id
      )
  )


# ------------------------------------------------------------------------------
# 8. Join sex metadata
# ------------------------------------------------------------------------------

hpa_expr_sex <- hpa_sample_level %>%

  left_join(
    hpa_samples_clean,
    by = "sample_id"
  )


# ------------------------------------------------------------------------------
# 9. Metadata join audit
# ------------------------------------------------------------------------------

join_audit <- hpa_expr_sex %>%

  summarise(

    Total_Rows =
      n(),

    Rows_With_Sex =
      sum(
        !is.na(Sex)
      ),

    Rows_Missing_Sex =
      sum(
        is.na(Sex)
      ),

    Percent_Matched =
      100 *
      Rows_With_Sex /
      Total_Rows
  )


cat("\nHPA metadata join audit:\n")

print(
  join_audit
)


# ------------------------------------------------------------------------------
# 10. Check tissue agreement
# ------------------------------------------------------------------------------

tissue_mismatch <- hpa_expr_sex %>%

  filter(
    !is.na(metadata_tissue),
    !is.na(tissue),
    str_to_lower(
      metadata_tissue
    ) !=
      str_to_lower(
        tissue
      )
  )


cat(
  "\nTissue-name mismatches:",
  nrow(
    tissue_mismatch
  ),
  "\n"
)


if (
  nrow(tissue_mismatch) > 0
) {

  print(
    tissue_mismatch %>%
      select(
        sample_id,
        tissue,
        metadata_tissue
      ) %>%
      distinct(),
    n = Inf
  )
}


# ------------------------------------------------------------------------------
# 11. Save joined sample-level table
# ------------------------------------------------------------------------------

write_csv(
  hpa_expr_sex,
  file.path(
    table_dir,
    "HPA_human_SERPIN_sample_level_TPM_by_sex.csv"
  )
)


# ==============================================================================
# 12. SEX-SPECIFIC TISSUE SUMMARY
# ==============================================================================

hpa_summary_by_sex <- hpa_expr_sex %>%

  filter(
    Sex %in% c(
      "Female",
      "Male"
    )
  ) %>%

  group_by(
    gene,
    tissue,
    Sex
  ) %>%

  summarise(

    N_Samples =
      n_distinct(
        sample_id
      ),

    Median_TPM =
      median(
        gene_TPM,
        na.rm = TRUE
      ),

    Mean_TPM =
      mean(
        gene_TPM,
        na.rm = TRUE
      ),

    SD_TPM =
      sd(
        gene_TPM,
        na.rm = TRUE
      ),

    IQR_TPM =
      IQR(
        gene_TPM,
        na.rm = TRUE
      ),

    Min_TPM =
      min(
        gene_TPM,
        na.rm = TRUE
      ),

    Max_TPM =
      max(
        gene_TPM,
        na.rm = TRUE
      ),

    Log2_Median_TPM =
      log2(
        Median_TPM + 1
      ),

    Log2_Mean_TPM =
      log2(
        Mean_TPM + 1
      ),

    .groups =
      "drop"
  ) %>%

  arrange(
    gene,
    Sex,
    desc(
      Median_TPM
    )
  )


# ------------------------------------------------------------------------------
# 13. Standardize column name
# ------------------------------------------------------------------------------

hpa_summary_by_sex <- hpa_summary_by_sex %>%

  rename(
    Gene =
      gene,

    Tissue =
      tissue
  )


write_csv(
  hpa_summary_by_sex,
  file.path(
    table_dir,
    "HPA_human_SERPIN_tissue_expression_summary_by_sex.csv"
  )
)


# ==============================================================================
# 14. COVERAGE TABLE
# ==============================================================================

hpa_coverage <- hpa_expr_sex %>%

  filter(
    Sex %in% c(
      "Female",
      "Male"
    )
  ) %>%

  group_by(
    tissue,
    Sex
  ) %>%

  summarise(

    N_Samples =
      n_distinct(
        sample_id
      ),

    .groups =
      "drop"
  ) %>%

  rename(
    Tissue =
      tissue
  ) %>%

  arrange(
    Tissue,
    Sex
  )


write_csv(
  hpa_coverage,
  file.path(
    table_dir,
    "HPA_tissue_sample_coverage_by_sex.csv"
  )
)


# ==============================================================================
# 15. SERPINE1 CHECK
# ==============================================================================

serpine1_check <- hpa_summary_by_sex %>%

  filter(
    Gene ==
      "SERPINE1"
  ) %>%

  select(
    Sex,
    Tissue,
    N_Samples,
    Median_TPM,
    Mean_TPM
  ) %>%

  arrange(
    Sex,
    desc(
      Median_TPM
    )
  )


cat("\n============================================\n")
cat("HPA SERPINE1 SEX-SPECIFIC SUMMARY\n")
cat("============================================\n")

print(
  serpine1_check,
  n = Inf
)


# ==============================================================================
# 16. FINAL CHECKS
# ==============================================================================

cat("\n============================================\n")
cat("HPA SEX-SPECIFIC BULK REANALYSIS COMPLETE\n")
cat("============================================\n")


cat(
  "\nNumber of tissues:",
  n_distinct(
    hpa_summary_by_sex$Tissue
  ),
  "\n"
)


cat(
  "Number of genes:",
  n_distinct(
    hpa_summary_by_sex$Gene
  ),
  "\n"
)


cat(
  "Sexes:",
  paste(
    unique(
      hpa_summary_by_sex$Sex
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "\nTables saved to:\n",
  table_dir,
  "\n"
)