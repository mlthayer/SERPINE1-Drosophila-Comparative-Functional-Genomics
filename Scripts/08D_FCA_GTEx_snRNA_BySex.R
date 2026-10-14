# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Run from the repository root.
# ==============================================================================
# FCA vs GTEx snRNA-seq
# SEX-SPECIFIC CROSS-SPECIES SERPIN EXPRESSION SIMILARITY
#
# Primary metric:
#   Percent of cells/nuclei expressing gene
#
# Comparisons:
#   Female FCA vs Female GTEx snRNA-seq
#   Male   FCA vs Male   GTEx snRNA-seq
#
# Functional systems:
#   Heart
#   Skeletal muscle
#   Digestive
#   Respiratory
#   Male reproductive (male only)
#
# Breast and Skin are excluded.
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
library(ggplot2)


# ------------------------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------------------------

fca_root <- "Expression/FlyCellAtlas/Reanalysis_Output"

fca_file <- Sys.getenv("SERPINE1_FCA_SUMMARY_INPUT", unset = file.path(
  fca_root,
  "Tables",
  "FCA_SERPIN_SampleLevel_MaleFemale_Only.csv"
))


gtex_root <- "Expression/GTEx/GTEx_snRNA_SERPIN_Reanalysis"

gtex_file <- Sys.getenv("SERPINE1_GTEX_DONOR_INPUT", unset = file.path(
  gtex_root,
  "Tables",
  "GTEx_snRNA_SERPIN_donor_system_summary.csv"
))


output_root <- "Expression/Human_SERPIN_reanalysis/CrossSpecies_FCA_GTEx_snRNA_BySex"

table_dir <- file.path(
  output_root,
  "Tables"
)

figure_dir <- file.path(
  output_root,
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


stopifnot(
  file.exists(fca_file),
  file.exists(gtex_file)
)


# ------------------------------------------------------------------------------
# 3. Genes
# ------------------------------------------------------------------------------

fly_genes <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)

human_genes <- c(
  "SERPINE1",
  "SERPINI1",
  "SERPINB3",
  "SERPINB4",
  "SERPINB11",
  "SERPINB9"
)


# ------------------------------------------------------------------------------
# 4. Read data
# ------------------------------------------------------------------------------

fca <- read_csv(
  fca_file,
  show_col_types = FALSE
)

gtex <- read_csv(
  gtex_file,
  show_col_types = FALSE
)


cat("\nFCA columns:\n")
print(names(fca))

cat("\nGTEx columns:\n")
print(names(gtex))


# ------------------------------------------------------------------------------
# 5. Standardize FCA sex and tissue names
# ------------------------------------------------------------------------------

fca <- fca %>%
  mutate(

    Sex = case_when(
      str_to_lower(Sex) == "female" ~ "Female",
      str_to_lower(Sex) == "male"   ~ "Male",
      TRUE ~ NA_character_
    ),

    Tissue_clean = Tissue %>%
      str_to_lower() %>%
      str_replace_all("-", "_") %>%
      str_replace_all(" ", "_")
  )


# ------------------------------------------------------------------------------
# 6. Map FCA tissues to the systems represented in GTEx snRNA-seq
# ------------------------------------------------------------------------------

fca <- fca %>%
  mutate(

    System = case_when(

      Tissue_clean == "heart" ~
        "Heart",

      Tissue_clean %in% c(
        "body_wall",
        "leg",
        "wing",
        "haltere"
      ) ~
        "Skeletal muscle",

      Tissue_clean == "gut" ~
        "Digestive",

      Tissue_clean == "trachea" ~
        "Respiratory",

      Tissue_clean == "male_reproductive_glands" ~
        "Male reproductive",

      TRUE ~
        NA_character_
    )
  )


# ------------------------------------------------------------------------------
# 7. Create FCA biological/sample unit
#
# Prefer Sample_ID.
# If absent, use FCA_ID.
# If both absent, retain Source_File as the observational unit.
# ------------------------------------------------------------------------------

fca <- fca %>%
  mutate(

    Sample_ID_clean = ifelse(
      !is.na(Sample_ID) &
        Sample_ID != "",
      as.character(Sample_ID),
      NA_character_
    ),

    FCA_ID_clean = ifelse(
      !is.na(FCA_ID) &
        FCA_ID != "",
      as.character(FCA_ID),
      NA_character_
    ),

    Replicate_ID = case_when(

      !is.na(Sample_ID_clean) ~
        paste(
          Source_File,
          Sample_ID_clean,
          sep = "::"
        ),

      !is.na(FCA_ID_clean) ~
        paste(
          Source_File,
          FCA_ID_clean,
          sep = "::"
        ),

      TRUE ~
        as.character(Source_File)
    )
  )


# ------------------------------------------------------------------------------
# 8. Restrict FCA to relevant genes, sexes, and systems
# ------------------------------------------------------------------------------

fca_use <- fca %>%
  filter(
    Gene %in% fly_genes,
    Sex %in% c(
      "Female",
      "Male"
    ),
    !is.na(System)
  )


# ------------------------------------------------------------------------------
# 9. FCA sample x system x gene prevalence
#
# Sum cell counts across cell types BEFORE computing percentage.
# ------------------------------------------------------------------------------

fca_sample_system <- fca_use %>%
  group_by(
    Replicate_ID,
    Sex,
    System,
    Gene
  ) %>%
  summarise(

    N_Cells =
      sum(
        N_Cells,
        na.rm = TRUE
      ),

    N_Expressing =
      sum(
        N_Expressing,
        na.rm = TRUE
      ),

    Percent_Expressing =
      ifelse(
        N_Cells > 0,
        100 * N_Expressing / N_Cells,
        NA_real_
      ),

    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 10. FCA sex x system summary
#
# Median across available FCA observational/sample units.
# ------------------------------------------------------------------------------

fca_sex_system <- fca_sample_system %>%
  group_by(
    Sex,
    System,
    Gene
  ) %>%
  summarise(

    N_FCA_Units =
      n_distinct(
        Replicate_ID
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

    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 11. GTEx donor-level summaries
# ------------------------------------------------------------------------------

gtex_use <- gtex %>%
  filter(
    Gene %in% human_genes,
    Sex %in% c(
      "Female",
      "Male"
    ),
    System %in% c(
      "Heart",
      "Skeletal muscle",
      "Digestive",
      "Respiratory",
      "Male reproductive"
    )
  )


# ------------------------------------------------------------------------------
# 12. GTEx sex x system summary
#
# Median donor prevalence is the primary metric.
# ------------------------------------------------------------------------------

gtex_sex_system <- gtex_use %>%
  group_by(
    Sex,
    System,
    Gene
  ) %>%
  summarise(

    N_GTEx_Donors =
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

    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 13. Explicit sex-specific allowable systems
# ------------------------------------------------------------------------------

female_systems <- c(
  "Heart",
  "Skeletal muscle",
  "Digestive",
  "Respiratory"
)

male_systems <- c(
  "Heart",
  "Skeletal muscle",
  "Digestive",
  "Respiratory",
  "Male reproductive"
)


# ------------------------------------------------------------------------------
# 14. Save source summaries
# ------------------------------------------------------------------------------

write_csv(
  fca_sample_system,
  file.path(
    table_dir,
    "FCA_sample_system_expression.csv"
  )
)

write_csv(
  fca_sex_system,
  file.path(
    table_dir,
    "FCA_system_expression_by_sex.csv"
  )
)

write_csv(
  gtex_sex_system,
  file.path(
    table_dir,
    "GTEx_snRNA_system_expression_by_sex.csv"
  )
)


# ------------------------------------------------------------------------------
# 15. Create every fly-human gene pair
# ------------------------------------------------------------------------------

gene_pairs <- tidyr::crossing(
  Fly_Gene = fly_genes,
  Human_Gene = human_genes
)


# ------------------------------------------------------------------------------
# 16. Spearman correlation function
# ------------------------------------------------------------------------------

calculate_pair_correlation <- function(
    sex_value,
    fly_gene,
    human_gene
) {

  allowed_systems <- if (
    sex_value == "Female"
  ) {
    female_systems
  } else {
    male_systems
  }


  fly_profile <- fca_sex_system %>%
    filter(
      Sex == sex_value,
      Gene == fly_gene,
      System %in% allowed_systems
    ) %>%
    transmute(
      System,
      Fly_Value =
        Median_Percent_Expressing,
      N_FCA_Units
    )


  human_profile <- gtex_sex_system %>%
    filter(
      Sex == sex_value,
      Gene == human_gene,
      System %in% allowed_systems
    ) %>%
    transmute(
      System,
      Human_Value =
        Median_Percent_Expressing,
      N_GTEx_Donors
    )


  matched <- inner_join(
    fly_profile,
    human_profile,
    by = "System"
  ) %>%
    filter(
      !is.na(Fly_Value),
      !is.na(Human_Value)
    )


  n_systems <- nrow(
    matched
  )


  if (
    n_systems < 3
  ) {

    return(
      tibble(
        Sex = sex_value,
        Fly_Gene = fly_gene,
        Human_Gene = human_gene,
        N_Systems = n_systems,
        Rho = NA_real_,
        P_Value = NA_real_
      )
    )
  }


  test <- suppressWarnings(
    cor.test(
      matched$Fly_Value,
      matched$Human_Value,
      method = "spearman",
      exact = FALSE
    )
  )


  tibble(
    Sex = sex_value,
    Fly_Gene = fly_gene,
    Human_Gene = human_gene,
    N_Systems = n_systems,
    Rho = unname(
      test$estimate
    ),
    P_Value = test$p.value
  )
}


# ------------------------------------------------------------------------------
# 17. Calculate all 72 correlations
#
# 36 female + 36 male
# ------------------------------------------------------------------------------

correlation_results <- tidyr::crossing(
  Sex = c(
    "Female",
    "Male"
  ),
  gene_pairs
) %>%
  mutate(

    Result = pmap(
      list(
        Sex,
        Fly_Gene,
        Human_Gene
      ),
      calculate_pair_correlation
    )

  ) %>%
  select(
    Result
  ) %>%
  unnest(
    Result
  )


# ------------------------------------------------------------------------------
# 18. FDR correction within sex
# ------------------------------------------------------------------------------

correlation_results <- correlation_results %>%
  group_by(
    Sex
  ) %>%
  mutate(

    FDR =
      p.adjust(
        P_Value,
        method = "BH"
      )

  ) %>%
  ungroup() %>%
  mutate(

    Significance = case_when(
      !is.na(FDR) & FDR < 0.001 ~ "***",
      !is.na(FDR) & FDR < 0.01  ~ "**",
      !is.na(FDR) & FDR < 0.05  ~ "*",
      TRUE ~ ""
    ),

    Heatmap_Label = ifelse(
      is.na(Rho),
      "NA",
      sprintf(
        "%.2f%s",
        Rho,
        Significance
      )
    )
  )


write_csv(
  correlation_results,
  file.path(
    table_dir,
    "FCA_GTEx_snRNA_Spearman_by_sex.csv"
  )
)


# ------------------------------------------------------------------------------
# 19. Save detailed matched-system profiles
#
# IMPORTANT:
# Use argument names that do NOT duplicate dataframe column names.
# ------------------------------------------------------------------------------

matched_profiles <- tidyr::crossing(
  sex_value = c(
    "Female",
    "Male"
  ),
  fly_gene = fly_genes,
  human_gene = human_genes
) %>%
  purrr::pmap_dfr(

    function(
      sex_value,
      fly_gene,
      human_gene
    ) {

      allowed_systems <- if (
        sex_value == "Female"
      ) {
        female_systems
      } else {
        male_systems
      }


      fly_profile <- fca_sex_system %>%
        dplyr::filter(
          Sex == .env$sex_value,
          Gene == .env$fly_gene,
          System %in% .env$allowed_systems
        ) %>%
        dplyr::transmute(
          System,
          FCA_Percent_Expressing =
            Median_Percent_Expressing,
          N_FCA_Units
        )


      human_profile <- gtex_sex_system %>%
        dplyr::filter(
          Sex == .env$sex_value,
          Gene == .env$human_gene,
          System %in% .env$allowed_systems
        ) %>%
        dplyr::transmute(
          System,
          GTEx_Percent_Expressing =
            Median_Percent_Expressing,
          N_GTEx_Donors
        )


      # Each system should occur exactly once in each profile.
      stopifnot(
        !anyDuplicated(
          fly_profile$System
        ),
        !anyDuplicated(
          human_profile$System
        )
      )


      dplyr::inner_join(
        fly_profile,
        human_profile,
        by = "System",
        relationship = "one-to-one"
      ) %>%
        dplyr::mutate(
          Sex = sex_value,
          Fly_Gene = fly_gene,
          Human_Gene = human_gene,
          .before = 1
        )
    }
  )


readr::write_csv(
  matched_profiles,
  file.path(
    table_dir,
    "FCA_GTEx_snRNA_matched_system_profiles.csv"
  )
)


# ------------------------------------------------------------------------------
# 20. Check matched profile structure
# ------------------------------------------------------------------------------

matched_profile_check <- matched_profiles %>%
  dplyr::count(
    Sex,
    Fly_Gene,
    Human_Gene,
    System,
    name = "N_Rows"
  )


# There should NEVER be more than one row per gene-pair/system.
problem_rows <- matched_profile_check %>%
  dplyr::filter(
    N_Rows != 1
  )


cat(
  "\nNumber of duplicated matched-profile rows:",
  nrow(problem_rows),
  "\n"
)

if (
  nrow(problem_rows) > 0
) {

  print(
    problem_rows,
    n = Inf
  )

  stop(
    "Duplicate rows remain in matched_profiles."
  )
}


# ------------------------------------------------------------------------------
# 21. Factor ordering for figures
# ------------------------------------------------------------------------------

correlation_results <- correlation_results %>%
  dplyr::mutate(

    Fly_Gene = factor(
      Fly_Gene,
      levels = rev(
        fly_genes
      )
    ),

    Human_Gene = factor(
      Human_Gene,
      levels = human_genes
    )
  )


# ------------------------------------------------------------------------------
# 22. Heatmap function
# ------------------------------------------------------------------------------

make_heatmap <- function(
  data,
  sex_value
) {

  plot_data <- data %>%
    dplyr::filter(
      Sex == .env$sex_value
    )


  ggplot(
    plot_data,
    aes(
      x = Human_Gene,
      y = Fly_Gene,
      fill = Rho
    )
  ) +

    geom_tile(
      linewidth = 0.5
    ) +

    geom_text(
      aes(
        label = Heatmap_Label
      ),
      size = 3.2
    ) +

    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(
        -1,
        1
      ),
      name = "Spearman\nrho"
    ) +

    labs(
      title = paste0(
        sex_value,
        " FCA vs GTEx snRNA-seq"
      ),
      subtitle = paste0(
        "Cross-species similarity based on median % ",
        "cells/nuclei expressing each gene"
      ),
      x = "Human SERPIN",
      y = "Drosophila serpin"
    ) +

    theme_minimal(
      base_size = 10
    ) +

    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      axis.text.y = element_text(
        size = 9
      ),
      axis.title = element_text(
        size = 10
      ),
      plot.title = element_text(
        size = 12,
        face = "bold"
      ),
      plot.subtitle = element_text(
        size = 9
      ),
      panel.grid = element_blank()
    )
}


# ------------------------------------------------------------------------------
# 23. Female heatmap
# ------------------------------------------------------------------------------

female_heatmap <- make_heatmap(
  correlation_results,
  "Female"
)

ggsave(
  filename = file.path(
    figure_dir,
    "Female_FCA_GTEx_heatmap.png"
  ),
  plot = female_heatmap,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)


# ------------------------------------------------------------------------------
# 24. Male heatmap
# ------------------------------------------------------------------------------

male_heatmap <- make_heatmap(
  correlation_results,
  "Male"
)

ggsave(
  filename = file.path(
    figure_dir,
    "Male_FCA_GTEx_heatmap.png"
  ),
  plot = male_heatmap,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)


# ------------------------------------------------------------------------------
# 25. Print strongest associations
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("TOP FEMALE CORRELATIONS\n")
cat("============================================\n")

correlation_results %>%
  dplyr::filter(
    Sex == "Female"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Rho)
  ) %>%
  dplyr::select(
    Fly_Gene,
    Human_Gene,
    N_Systems,
    Rho,
    P_Value,
    FDR
  ) %>%
  print(
    n = 36
  )


cat("\n============================================\n")
cat("TOP MALE CORRELATIONS\n")
cat("============================================\n")

correlation_results %>%
  dplyr::filter(
    Sex == "Male"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Rho)
  ) %>%
  dplyr::select(
    Fly_Gene,
    Human_Gene,
    N_Systems,
    Rho,
    P_Value,
    FDR
  ) %>%
  print(
    n = 36
  )


# ------------------------------------------------------------------------------
# 26. Coverage audit
#
# Obtain coverage directly from the source summary tables rather than deriving
# it from all 36 gene-pair combinations.
# ------------------------------------------------------------------------------

fca_coverage <- fca_sex_system %>%
  dplyr::filter(
    System %in% c(
      female_systems,
      male_systems
    )
  ) %>%
  dplyr::group_by(
    Sex,
    System
  ) %>%
  dplyr::summarise(
    FCA_Units =
      max(
        N_FCA_Units,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


gtex_coverage <- gtex_sex_system %>%
  dplyr::filter(
    System %in% c(
      female_systems,
      male_systems
    )
  ) %>%
  dplyr::group_by(
    Sex,
    System
  ) %>%
  dplyr::summarise(
    GTEx_Donors =
      max(
        N_GTEx_Donors,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


coverage_audit <- dplyr::full_join(
  fca_coverage,
  gtex_coverage,
  by = c(
    "Sex",
    "System"
  ),
  relationship = "one-to-one"
) %>%
  dplyr::filter(
    !(Sex == "Female" &
        System == "Male reproductive")
  ) %>%
  dplyr::arrange(
    Sex,
    System
  )


readr::write_csv(
  coverage_audit,
  file.path(
    table_dir,
    "FCA_GTEx_system_coverage_audit.csv"
  )
)


cat("\n============================================\n")
cat("SYSTEM COVERAGE\n")
cat("============================================\n")

print(
  coverage_audit,
  n = Inf
)


# ------------------------------------------------------------------------------
# 27. Completion
# ------------------------------------------------------------------------------

cat("\n============================================\n")
cat("Corrected FCA vs GTEx snRNA-seq analysis complete\n")
cat("============================================\n")