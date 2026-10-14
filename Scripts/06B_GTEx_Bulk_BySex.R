# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Run from the repository root.
# ==============================================================================
# GTEx v11 BULK TISSUE SERPIN EXPRESSION BY SEX
#
# Purpose:
#   Rebuild GTEx v11 bulk-tissue expression summaries separately for
#   female and male donors.
#
# Genes:
#   SERPINE1
#   SERPINI1
#   SERPINB3
#   SERPINB4
#   SERPINB11
#   SERPINB9
#
# Primary tissue-level metric:
#   Median TPM within sex and tissue
#
# Outputs:
#   1. GTEx sample-level SERPIN TPM with subject sex
#   2. GTEx tissue-expression summary by sex
#   3. GTEx donor/sample coverage by sex and tissue
#   4. Gene-presence audit
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(purrr)
library(data.table)


# ------------------------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------------------------

gtex_dir <- Sys.getenv("SERPINE1_GTEX_INPUT", unset = "Expression/GTEx")


subject_metadata_file <- Sys.getenv("SERPINE1_GTEX_SUBJECT_INPUT", unset = file.path(
  gtex_dir, "Annotations",
  "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"
))


output_root <- "Expression/Human_SERPIN_reanalysis"


table_dir <- file.path(
  output_root,
  "Tables"
)


dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


stopifnot(
  dir.exists(gtex_dir),
  file.exists(subject_metadata_file)
)


# ------------------------------------------------------------------------------
# 3. Target genes
# ------------------------------------------------------------------------------

target_genes <- c(
  "SERPINE1",
  "SERPINI1",
  "SERPINB3",
  "SERPINB4",
  "SERPINB11",
  "SERPINB9"
)


# Optional Ensembl IDs for audit/reference
target_ensg <- c(
  SERPINE1   = "ENSG00000106366",
  SERPINI1   = "ENSG00000163536",
  SERPINB3   = "ENSG00000057149",
  SERPINB4   = "ENSG00000206073",
  SERPINB11  = "ENSG00000206072",
  SERPINB9   = "ENSG00000170542"
)


# ==============================================================================
# 4. READ GTEx SUBJECT METADATA
# ==============================================================================

gtex_subjects <- read.delim(
  subject_metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\nGTEx subject metadata columns:\n")
print(
  names(gtex_subjects)
)


# ------------------------------------------------------------------------------
# 5. Recode GTEx sex
#
# GTEx:
#   1 = Male
#   2 = Female
# ------------------------------------------------------------------------------

gtex_subjects <- gtex_subjects %>%

  mutate(

    Sex = case_when(
      SEX == 1 ~ "Male",
      SEX == 2 ~ "Female",
      TRUE ~ NA_character_
    )

  ) %>%

  select(
    SUBJID,
    Sex,
    AGE,
    DTHHRDY
  )


cat("\nGTEx subject counts by sex:\n")

gtex_subjects %>%
  count(
    Sex
  ) %>%
  print()


# ==============================================================================
# 6. LOCATE GTEx GCT FILES
# ==============================================================================

gct_files <- list.files(
  gtex_dir,
  pattern = "\\.gct(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)


cat(
  "\nNumber of GCT files found:",
  length(gct_files),
  "\n"
)


print(
  basename(
    gct_files
  )
)


if (
  length(gct_files) == 0
) {

  stop(
    "No GTEx GCT or GCT.GZ files were found in the GTEx directory."
  )
}


# ==============================================================================
# 7. TISSUE-NAME FUNCTION
#
# Uses filename text to assign the tissue names used in the existing
# sex-specific tissue-expression analysis.
# ==============================================================================

get_tissue_name <- function(filename) {

  x <- basename(filename) %>%
    str_to_lower()


  case_when(

    str_detect(
      x,
      "adipose.*subcutaneous"
    ) ~
      "Adipose Subcutaneous",


    str_detect(
      x,
      "adipose.*visceral.*omentum"
    ) ~
      "Adipose Visceral Omentum",


    str_detect(
      x,
      "artery.*aorta"
    ) ~
      "Artery Aorta",


    str_detect(
      x,
      "artery.*coronary"
    ) ~
      "Artery Coronary",


    str_detect(
      x,
      "artery.*tibial"
    ) ~
      "Artery Tibial",


    str_detect(
      x,
      "brain.*hypothalamus"
    ) ~
      "Brain Hypothalamus",


    str_detect(
      x,
      "nucleus.*accumbens"
    ) ~
      "Brain Nucleus Accumbens Basal Ganglia",


    str_detect(
      x,
      "spinal.*cord"
    ) ~
      "Brain Spinal Cord Cervical C-1",


    str_detect(
      x,
      "heart.*atrial"
    ) ~
      "Heart Atrial Appendage",


    str_detect(
      x,
      "heart.*left.*ventricle"
    ) ~
      "Heart Left Ventricle",


    str_detect(
      x,
      "liver"
    ) ~
      "Liver",


    str_detect(
      x,
      "lung"
    ) ~
      "Lung",


    str_detect(
      x,
      "muscle.*skeletal"
    ) ~
      "Muscle Skeletal",


    str_detect(
      x,
      "ovary"
    ) ~
      "Ovary",


    str_detect(
      x,
      "pancreas"
    ) ~
      "Pancreas",


    str_detect(
      x,
      "small.*intestine"
    ) ~
      "Small Intestine Terminal Ileum",


    str_detect(
      x,
      "testis"
    ) ~
      "Testis",


    str_detect(
      x,
      "whole.*blood"
    ) ~
      "Whole Blood",


    TRUE ~
      NA_character_
  )
}


# ==============================================================================
# 8. CHECK FILE-TO-TISSUE MAPPING
# ==============================================================================

file_map <- tibble(
  Source_File = basename(
    gct_files
  ),
  File_Path = gct_files,
  Tissue = map_chr(
    gct_files,
    get_tissue_name
  )
)


cat("\nGTEx file-to-tissue mapping:\n")

print(
  file_map,
  n = Inf
)


unmapped_files <- file_map %>%
  filter(
    is.na(Tissue)
  )


if (
  nrow(unmapped_files) > 0
) {

  cat(
    "\nWARNING: Some GCT files could not be mapped to a tissue:\n"
  )

  print(
    unmapped_files,
    n = Inf
  )

  stop(
    "Resolve unmapped GTEx filenames before continuing."
  )
}


write_csv(
  file_map,
  file.path(
    table_dir,
    "GTEx_v11_GCT_file_tissue_mapping.csv"
  )
)


# ==============================================================================
# 9. FUNCTION TO READ ONE GTEx TISSUE FILE
# ==============================================================================

read_gtex_tissue <- function(
    file_path,
    tissue_name
) {

  cat(
    "\nProcessing:",
    tissue_name,
    "\n"
  )


  # --------------------------------------------------------------------------
  # GTEx GCT files have two metadata/header lines before the column header.
  # --------------------------------------------------------------------------

  dat <- data.table::fread(
    file_path,
    skip = 2,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE
  )


  cat(
    "Rows:",
    nrow(dat),
    " Columns:",
    ncol(dat),
    "\n"
  )


  # --------------------------------------------------------------------------
  # Identify gene-symbol column
  # --------------------------------------------------------------------------

  possible_symbol_columns <- c(
    "Description",
    "description",
    "gene_name",
    "Gene",
    "gene",
    "symbol",
    "Symbol"
  )


  symbol_column <- intersect(
    possible_symbol_columns,
    names(dat)
  )


  if (
    length(symbol_column) == 0
  ) {

    stop(
      paste0(
        "No gene-symbol column detected in ",
        basename(file_path)
      )
    )
  }


  symbol_column <- symbol_column[1]


  # --------------------------------------------------------------------------
  # Identify annotation columns that should NOT be interpreted as samples
  # --------------------------------------------------------------------------

  annotation_columns <- unique(
    c(
      "Name",
      "Description",
      "description",
      "gene_name",
      "Gene",
      "gene",
      "symbol",
      "Symbol"
    )
  )


  annotation_columns <- intersect(
    annotation_columns,
    names(dat)
  )


  sample_columns <- setdiff(
    names(dat),
    annotation_columns
  )


  # --------------------------------------------------------------------------
  # Extract target genes
  # --------------------------------------------------------------------------

  gene_data <- dat %>%

    filter(
      .data[[symbol_column]] %in%
        target_genes
    ) %>%

    mutate(
      Gene =
        .data[[symbol_column]]
    )


  cat(
    "Target gene rows found:",
    nrow(gene_data),
    "\n"
  )


  # --------------------------------------------------------------------------
  # Convert to sample-level long format
  # --------------------------------------------------------------------------

  sample_long <- gene_data %>%

    select(
      Gene,
      all_of(
        sample_columns
      )
    ) %>%

    pivot_longer(
      cols =
        -Gene,
      names_to =
        "Sample_ID",
      values_to =
        "TPM"
    ) %>%

    mutate(

      TPM =
        as.numeric(
          TPM
        ),

      Tissue =
        tissue_name

    )


  # --------------------------------------------------------------------------
  # GTEx tissue sample:
  #
  # GTEX-XXXX-YYYY...
  #
  # Subject ID:
  #
  # GTEX-XXXX
  # --------------------------------------------------------------------------

  sample_long <- sample_long %>%

    mutate(

      SUBJID =
        str_extract(
          Sample_ID,
          "^GTEX-[^-]+"
        )

    )


  sample_long
}


# ==============================================================================
# 10. PROCESS ALL GTEx TISSUES
# ==============================================================================

gtex_sample_level <- map2_dfr(
  file_map$File_Path,
  file_map$Tissue,
  read_gtex_tissue
)


cat(
  "\nTotal SERPIN sample-level rows:",
  nrow(gtex_sample_level),
  "\n"
)


# ==============================================================================
# 11. JOIN SUBJECT SEX
# ==============================================================================

gtex_sample_level <- gtex_sample_level %>%

  left_join(
    gtex_subjects,
    by = "SUBJID"
  )


# ------------------------------------------------------------------------------
# Join audit
# ------------------------------------------------------------------------------

join_audit <- gtex_sample_level %>%

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


cat("\nMetadata join audit:\n")

print(
  join_audit
)


if (
  join_audit$Rows_Missing_Sex > 0
) {

  cat(
    "\nWARNING: Some GTEx sample rows did not match subject sex metadata.\n"
  )
}


# ==============================================================================
# 12. SAVE SAMPLE-LEVEL TABLE
# ==============================================================================

gtex_sample_level <- gtex_sample_level %>%

  arrange(
    Gene,
    Tissue,
    Sex,
    SUBJID,
    Sample_ID
  )


write_csv(
  gtex_sample_level,
  file.path(
    table_dir,
    "GTEx_v11_human_SERPIN_sample_level_TPM_by_sex.csv"
  )
)


# ==============================================================================
# 13. GENE-PRESENCE AUDIT
# ==============================================================================

gene_presence <- crossing(
  Tissue =
    unique(
      file_map$Tissue
    ),
  Gene =
    target_genes
) %>%

  left_join(

    gtex_sample_level %>%

      distinct(
        Tissue,
        Gene
      ) %>%

      mutate(
        Present =
          TRUE
      ),

    by = c(
      "Tissue",
      "Gene"
    )
  ) %>%

  mutate(
    Present =
      ifelse(
        is.na(Present),
        FALSE,
        Present
      )
  )


write_csv(
  gene_presence,
  file.path(
    table_dir,
    "GTEx_v11_SERPIN_gene_presence_by_tissue.csv"
  )
)


cat("\nGene-presence audit:\n")

gene_presence %>%
  count(
    Gene,
    Present
  ) %>%
  print(
    n = Inf
  )


# ==============================================================================
# 14. SEX-SPECIFIC TISSUE SUMMARY
#
# Primary summary for cross-species tissue analysis.
# ==============================================================================

gtex_summary_by_sex <- gtex_sample_level %>%

  filter(
    Sex %in% c(
      "Female",
      "Male"
    )
  ) %>%

  group_by(
    Gene,
    Tissue,
    Sex
  ) %>%

  summarise(

    N_Samples =
      n(),

    N_Donors =
      n_distinct(
        SUBJID
      ),

    Median_TPM =
      median(
        TPM,
        na.rm = TRUE
      ),

    Mean_TPM =
      mean(
        TPM,
        na.rm = TRUE
      ),

    SD_TPM =
      sd(
        TPM,
        na.rm = TRUE
      ),

    IQR_TPM =
      IQR(
        TPM,
        na.rm = TRUE
      ),

    Min_TPM =
      min(
        TPM,
        na.rm = TRUE
      ),

    Max_TPM =
      max(
        TPM,
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
    Gene,
    Sex,
    desc(
      Median_TPM
    )
  )


write_csv(
  gtex_summary_by_sex,
  file.path(
    table_dir,
    "GTEx_v11_human_SERPIN_tissue_expression_summary_by_sex.csv"
  )
)


# ==============================================================================
# 15. COVERAGE SUMMARY
# ==============================================================================

gtex_coverage <- gtex_sample_level %>%

  filter(
    Sex %in% c(
      "Female",
      "Male"
    )
  ) %>%

  group_by(
    Tissue,
    Sex
  ) %>%

  summarise(

    N_Samples =
      n_distinct(
        Sample_ID
      ),

    N_Donors =
      n_distinct(
        SUBJID
      ),

    .groups =
      "drop"
  ) %>%

  arrange(
    Tissue,
    Sex
  )


write_csv(
  gtex_coverage,
  file.path(
    table_dir,
    "GTEx_v11_tissue_donor_coverage_by_sex.csv"
  )
)


# ==============================================================================
# 16. FEMALE/Male SERPINE1 CHECK
# ==============================================================================

serpine1_check <- gtex_summary_by_sex %>%

  filter(
    Gene ==
      "SERPINE1"
  ) %>%

  select(
    Sex,
    Tissue,
    N_Samples,
    N_Donors,
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
cat("SERPINE1 SEX-SPECIFIC GTEx SUMMARY\n")
cat("============================================\n")

print(
  serpine1_check,
  n = Inf
)


# ==============================================================================
# 17. FINAL COVERAGE AUDIT
# ==============================================================================

cat("\n============================================\n")
cat("GTEx DONOR COVERAGE BY SEX AND TISSUE\n")
cat("============================================\n")

print(
  gtex_coverage,
  n = Inf
)


# ==============================================================================
# 18. FINAL CHECKS
# ==============================================================================

cat("\n============================================\n")
cat("GTEx SEX-SPECIFIC BULK REANALYSIS COMPLETE\n")
cat("============================================\n")


cat(
  "\nNumber of tissues:",
  n_distinct(
    gtex_summary_by_sex$Tissue
  ),
  "\n"
)


cat(
  "Number of genes:",
  n_distinct(
    gtex_summary_by_sex$Gene
  ),
  "\n"
)


cat(
  "Sexes:",
  paste(
    unique(
      gtex_summary_by_sex$Sex
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