# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Run from the repository root.
# ==============================================================================
# GTEx snRNA-seq SERPIN extraction
#
# Purpose:
#   Extract selected human SERPIN genes from the GTEx 8-tissue snRNA-seq atlas
#   and summarize expression by donor, sex, tissue, and functional system.
#
# Primary cross-species metric:
#   Percent of nuclei expressing each gene
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(rhdf5)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)


# ------------------------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------------------------

gtex_sc_dir <- Sys.getenv("SERPINE1_GTEX_SNRNA_DIR", unset = "Expression/GTEx")

gtex_sc_file <- Sys.getenv("SERPINE1_GTEX_SNRNA_INPUT", unset = file.path(
  gtex_sc_dir,
  "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
))

output_dir <- Sys.getenv("SERPINE1_GTEX_SNRNA_OUTPUT", unset = file.path(
  gtex_sc_dir,
  "GTEx_snRNA_SERPIN_Reanalysis"
))

tables_dir <- file.path(
  output_dir,
  "Tables"
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Select the downloaded atlas in R/RStudio; automated runs use the configured path.
if (isTRUE(getOption("serpine1.choose_files", FALSE)) || !file.exists(gtex_sc_file)) {
  message("Download GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad from:\n",
          "https://gtexportal.org/home/downloads/adult-gtex/single_cell")
  if (!interactive()) stop("Atlas missing. Set SERPINE1_GTEX_SNRNA_INPUT to the downloaded .h5ad file.")
  gtex_sc_file <- file.choose()
}
if (!file.exists(gtex_sc_file) || !grepl("\\.h5ad$", gtex_sc_file)) {
  stop("Select a valid downloaded .h5ad atlas.")
}


# ------------------------------------------------------------------------------
# 3. Human SERPIN genes
# ------------------------------------------------------------------------------

selected_genes <- c(
  "SERPINE1",
  "SERPINI1",
  "SERPINB3",
  "SERPINB4",
  "SERPINB11",
  "SERPINB9"
)


# ------------------------------------------------------------------------------
# 4. Decode AnnData categorical variable
# ------------------------------------------------------------------------------

decode_obs_category <- function(file, field) {

  codes <- h5read(
    file,
    paste0("/obs/", field)
  )

  categories <- h5read(
    file,
    paste0(
      "/obs/__categories/",
      field
    )
  )

  result <- rep(
    NA_character_,
    length(codes)
  )

  keep <- !is.na(codes) &
    codes >= 0

  result[keep] <-
    categories[
      codes[keep] + 1
    ]

  result
}


decode_var_category <- function(file, field) {

  codes <- h5read(
    file,
    paste0(
      "/raw/var/",
      field
    )
  )

  categories <- h5read(
    file,
    paste0(
      "/raw/var/__categories/",
      field
    )
  )

  result <- rep(
    NA_character_,
    length(codes)
  )

  keep <- !is.na(codes) &
    codes >= 0

  result[keep] <-
    categories[
      codes[keep] + 1
    ]

  result
}


# ------------------------------------------------------------------------------
# 5. Build GTEx cell metadata
# ------------------------------------------------------------------------------

cat("\nDecoding GTEx cell metadata...\n")

gtex_meta <- tibble(

  Cell_Index =
    seq_along(
      decode_obs_category(
        gtex_sc_file,
        "Sex"
      )
    ),

  Participant_ID =
    decode_obs_category(
      gtex_sc_file,
      "Participant ID"
    ),

  Sex =
    decode_obs_category(
      gtex_sc_file,
      "Sex"
    ),

  Age_bin =
    decode_obs_category(
      gtex_sc_file,
      "Age_bin"
    ),

  Tissue =
    decode_obs_category(
      gtex_sc_file,
      "Tissue"
    ),

  Tissue_Site_Detail =
    decode_obs_category(
      gtex_sc_file,
      "Tissue Site Detail"
    ),

  Broad_Cell_Type =
    decode_obs_category(
      gtex_sc_file,
      "Broad cell type"
    ),

  Cell_Type_L2 =
    decode_obs_category(
      gtex_sc_file,
      "Cell types level 2"
    ),

  Cell_Type_L3 =
    decode_obs_category(
      gtex_sc_file,
      "Cell types level 3"
    ),

  Granular_Cell_Type =
    decode_obs_category(
      gtex_sc_file,
      "Granular cell type"
    )
)


# ------------------------------------------------------------------------------
# 6. Cross-species functional-system mapping
#
# Breast and Skin intentionally excluded.
# ------------------------------------------------------------------------------

gtex_meta <- gtex_meta %>%
  mutate(

    System = case_when(

      Tissue == "Heart" ~
        "Heart",

      Tissue == "Skeletal muscle" ~
        "Skeletal muscle",

      Tissue %in% c(
        "Esophagus mucosa",
        "Esophagus muscularis"
      ) ~
        "Digestive",

      Tissue == "Lung" ~
        "Respiratory",

      Tissue == "Prostate" ~
        "Male reproductive",

      TRUE ~
        NA_character_
    )
  )


# ------------------------------------------------------------------------------
# 7. Save metadata audit table
# ------------------------------------------------------------------------------

write_csv(
  gtex_meta,
  file.path(
    tables_dir,
    "GTEx_snRNA_cell_metadata.csv"
  )
)


# ------------------------------------------------------------------------------
# 8. Inspect raw gene annotations
# ------------------------------------------------------------------------------

cat("\nDecoding gene symbols...\n")

gene_symbols <- decode_var_category(
  gtex_sc_file,
  "Approved symbol"
)

gene_check <- tibble(
  Gene = selected_genes,
  Present = selected_genes %in% gene_symbols,
  Gene_Index_R = match(
    selected_genes,
    gene_symbols
  ),
  Gene_Index_H5 = match(
    selected_genes,
    gene_symbols
  ) - 1
)

print(gene_check)

write_csv(
  gene_check,
  file.path(
    tables_dir,
    "GTEx_snRNA_SERPIN_gene_check.csv"
  )
)

if (
  any(
    !gene_check$Present
  )
) {

  stop(
    "One or more selected SERPIN genes were not found in /raw/var."
  )
}


# ------------------------------------------------------------------------------
# 9. Prepare sparse raw-count matrix information
#
# GTEx raw/X is CSR:
#   rows    = nuclei
#   columns = genes
# ------------------------------------------------------------------------------

indptr <- h5read(
  gtex_sc_file,
  "/raw/X/indptr"
)

n_cells <- length(indptr) - 1

if (
  n_cells != nrow(gtex_meta)
) {

  stop(
    "Number of raw matrix rows does not match GTEx metadata rows."
  )
}


target_gene_indices <- gene_check$Gene_Index_H5

names(
  target_gene_indices
) <- gene_check$Gene


# ------------------------------------------------------------------------------
# 10. Allocate selected-gene matrix
#
# Only 6 genes x ~209,000 nuclei, so this is memory-safe.
# ------------------------------------------------------------------------------

serpin_counts <- matrix(
  0,
  nrow = n_cells,
  ncol = length(selected_genes)
)

colnames(
  serpin_counts
) <- selected_genes


# ------------------------------------------------------------------------------
# 11. Extract selected genes from CSR matrix in chunks
# ------------------------------------------------------------------------------

chunk_size <- 5000

starts <- seq(
  1,
  n_cells,
  by = chunk_size
)

cat(
  "\nExtracting selected SERPIN genes from raw counts...\n"
)


for (
  chunk_start in starts
) {

  chunk_end <- min(
    chunk_start + chunk_size - 1,
    n_cells
  )

  cat(
    "Cells",
    chunk_start,
    "to",
    chunk_end,
    "of",
    n_cells,
    "\n"
  )

  # CSR offsets are zero-based
  entry_start <-
    indptr[chunk_start] + 1

  entry_end <-
    indptr[chunk_end + 1]


  if (
    entry_end < entry_start
  ) {

    next
  }


  indices_chunk <- h5read(
    gtex_sc_file,
    "/raw/X/indices",
    index = list(
      entry_start:entry_end
    )
  )

  data_chunk <- h5read(
    gtex_sc_file,
    "/raw/X/data",
    index = list(
      entry_start:entry_end
    )
  )


  row_lengths <- diff(
    indptr[
      chunk_start:(chunk_end + 1)
    ]
  )


  row_ids <- rep(
    chunk_start:chunk_end,
    times = row_lengths
  )


  for (
    gene_name in selected_genes
  ) {

    gene_index <-
      target_gene_indices[
        gene_name
      ]

    hits <-
      which(
        indices_chunk ==
          gene_index
      )


    if (
      length(hits) > 0
    ) {

      serpin_counts[
        row_ids[hits],
        gene_name
      ] <-
        data_chunk[hits]
    }
  }


  rm(
    indices_chunk,
    data_chunk,
    row_ids
  )

  gc()
}


# ------------------------------------------------------------------------------
# 12. Convert selected genes to long format
# ------------------------------------------------------------------------------

serpin_long <- as_tibble(
  serpin_counts
) %>%
  mutate(
    Cell_Index =
      seq_len(
        n()
      )
  ) %>%
  pivot_longer(
    cols =
      all_of(
        selected_genes
      ),
    names_to =
      "Gene",
    values_to =
      "Raw_Count"
  ) %>%
  left_join(
    gtex_meta,
    by =
      "Cell_Index"
  ) %>%
  mutate(
    Expressing =
      Raw_Count > 0
  )


# ------------------------------------------------------------------------------
# 13. Donor x tissue x gene summary
#
# Biological donor remains the replicate.
# ------------------------------------------------------------------------------

donor_tissue_summary <- serpin_long %>%
  filter(
    !is.na(Participant_ID),
    !is.na(Sex),
    !is.na(Tissue)
  ) %>%
  group_by(
    Participant_ID,
    Sex,
    Age_bin,
    Tissue,
    System,
    Gene
  ) %>%
  summarise(

    N_Nuclei =
      n(),

    N_Expressing =
      sum(
        Expressing,
        na.rm = TRUE
      ),

    Percent_Expressing =
      100 *
      N_Expressing /
      N_Nuclei,

    Mean_Raw_Count_All =
      mean(
        Raw_Count,
        na.rm = TRUE
      ),

    Mean_Raw_Count_Positive =
      ifelse(
        N_Expressing > 0,
        mean(
          Raw_Count[
            Expressing
          ],
          na.rm = TRUE
        ),
        0
      ),

    .groups =
      "drop"
  )


write_csv(
  donor_tissue_summary,
  file.path(
    tables_dir,
    "GTEx_snRNA_SERPIN_donor_tissue_summary.csv"
  )
)


# ------------------------------------------------------------------------------
# 14. Donor x functional system summary
#
# Esophagus mucosa + muscularis are pooled within donor as Digestive.
#
# Breast and Skin are excluded.
# ------------------------------------------------------------------------------

donor_system_summary <- serpin_long %>%
  filter(
    !is.na(Participant_ID),
    !is.na(Sex),
    !is.na(System)
  ) %>%
  group_by(
    Participant_ID,
    Sex,
    Age_bin,
    System,
    Gene
  ) %>%
  summarise(

    N_Nuclei =
      n(),

    N_Expressing =
      sum(
        Expressing,
        na.rm = TRUE
      ),

    Percent_Expressing =
      100 *
      N_Expressing /
      N_Nuclei,

    Mean_Raw_Count_All =
      mean(
        Raw_Count,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


write_csv(
  donor_system_summary,
  file.path(
    tables_dir,
    "GTEx_snRNA_SERPIN_donor_system_summary.csv"
  )
)


# ------------------------------------------------------------------------------
# 15. Sex-specific system summary
#
# Primary value for cross-species comparison:
# median donor % expressing
# ------------------------------------------------------------------------------

sex_system_summary <- donor_system_summary %>%
  group_by(
    Sex,
    System,
    Gene
  ) %>%
  summarise(

    N_Donors =
      n_distinct(
        Participant_ID
      ),

    Median_Percent_Expressing =
      median(
        Percent_Expressing,
        na.rm = TRUE
      ),

    Mean_Percent_Expressing =
      mean(
        Percent_Expressing,
        na.rm = TRUE
      ),

    SD_Percent_Expressing =
      sd(
        Percent_Expressing,
        na.rm = TRUE
      ),

    IQR_Percent_Expressing =
      IQR(
        Percent_Expressing,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


write_csv(
  sex_system_summary,
  file.path(
    tables_dir,
    "GTEx_snRNA_SERPIN_system_expression_by_sex.csv"
  )
)


# ------------------------------------------------------------------------------
# 16. Donor coverage table
# ------------------------------------------------------------------------------

donor_coverage <- donor_system_summary %>%
  distinct(
    Participant_ID,
    Sex,
    System
  ) %>%
  count(
    Sex,
    System,
    name =
      "N_Donors"
  ) %>%
  arrange(
    Sex,
    System
  )


print(
  donor_coverage
)

write_csv(
  donor_coverage,
  file.path(
    tables_dir,
    "GTEx_snRNA_donor_coverage_by_sex_system.csv"
  )
)


# ------------------------------------------------------------------------------
# 17. Final audit
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("GTEx snRNA-seq SERPIN extraction complete\n")
cat("============================================\n\n")

print(
  sex_system_summary,
  n = Inf
)

cat(
  "\nOutput folder:\n",
  output_dir,
  "\n"
)
