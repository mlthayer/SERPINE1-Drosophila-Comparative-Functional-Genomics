# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 08B_SexSpecific_Bulk_TissueSimilarity.R
# Run from the repository root.
# ==============================================================================
# SEX-SPECIFIC CROSS-SPECIES BULK TISSUE EXPRESSION ANALYSIS
#
# Drosophila:
#   FlyAtlas2 adult tissue expression (FPKM)
#
# Human:
#   GTEx v11 bulk tissue expression (TPM)
#   HPA bulk tissue expression (TPM)
#
# Comparisons:
#   Female FlyAtlas2 vs Female GTEx
#   Male   FlyAtlas2 vs Male   GTEx
#   Female FlyAtlas2 vs Female HPA
#   Male   FlyAtlas2 vs Male   HPA
#
# Cross-species statistic:
#   Spearman rank correlation across matched functional systems
#
# IMPORTANT:
#   FPKM and TPM are NOT directly compared numerically.
#   Correlations evaluate relative tissue/system expression patterns.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)


# ==============================================================================
# 2. PATHS
# ==============================================================================

# FlyAtlas2 -------------------------------------------------------------

flyatlas_file <- file.path("Expression", "FlyAtlas2_Analysis", "Processed",
                           "05_FlyAtlas2_adult_gene_expression.csv")


# Human reanalysis ------------------------------------------------------

human_root <- "Expression/Human_SERPIN_reanalysis"

human_table_dir <- file.path(
  human_root,
  "Tables"
)


# GTEx sex-specific summary --------------------------------------------

gtex_file <- Sys.getenv("SERPINE1_GTEX_BYSEX_INPUT", unset = file.path(
  human_root, "Reference_Uploads",
  "GTEx_SERPIN_tissue_expression_by_sex.csv"
))


# HPA sex-specific summary ---------------------------------------------

hpa_file <- Sys.getenv("SERPINE1_HPA_BYSEX_INPUT", unset = file.path(
  human_root, "Reference_Uploads",
  "HPA_SERPIN_tissue_expression_by_sex.csv"
))


# Analysis output -------------------------------------------------------

output_root <- file.path(
  human_root,
  "CrossSpecies_Bulk_BySex"
)

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


# Verify files ----------------------------------------------------------

stopifnot(
  file.exists(flyatlas_file),
  file.exists(gtex_file),
  file.exists(hpa_file)
)


# ==============================================================================
# 3. TARGET GENES
# ==============================================================================

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


# ==============================================================================
# 4. FUNCTIONAL SYSTEM ORDER
# ==============================================================================

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


female_systems <- c(
  "Heart",
  "Skeletal muscle",
  "Liver-metabolic",
  "Digestive",
  "Renal/excretory",
  "Respiratory",
  "Ovary",
  "Sensory/neural"
)


male_systems <- c(
  "Heart",
  "Skeletal muscle",
  "Liver-metabolic",
  "Digestive",
  "Renal/excretory",
  "Respiratory",
  "Male reproductive",
  "Sensory/neural"
)


# ==============================================================================
# 5. READ INPUT DATA
# ==============================================================================

flyatlas <- read_csv(
  flyatlas_file,
  show_col_types = FALSE
)

gtex <- read_csv(
  gtex_file,
  show_col_types = FALSE
)

hpa <- read_csv(
  hpa_file,
  show_col_types = FALSE
)


# Uploaded prepared tissue summaries use median_TPM rather than Median_TPM.
# Normalize the column label only; keep expression values and mappings unchanged.
for (dataset_name in c("gtex", "hpa")) {
  dataset <- get(dataset_name)
  if (!"Median_TPM" %in% names(dataset) && "median_TPM" %in% names(dataset)) {
    names(dataset)[names(dataset) == "median_TPM"] <- "Median_TPM"
  }
  assign(dataset_name, dataset)
}

cat("\nFlyAtlas2 columns:\n")
print(names(flyatlas))

cat("\nGTEx columns:\n")
print(names(gtex))

cat("\nHPA columns:\n")
print(names(hpa))


# ==============================================================================
# 6. STANDARDIZE SEX
# ==============================================================================

standardize_sex <- function(x) {

  case_when(

    str_to_lower(
      as.character(x)
    ) %in% c(
      "female",
      "f"
    ) ~
      "Female",

    str_to_lower(
      as.character(x)
    ) %in% c(
      "male",
      "m"
    ) ~
      "Male",

    TRUE ~
      NA_character_
  )
}


flyatlas <- flyatlas %>%
  mutate(
    Sex =
      standardize_sex(
        Sex
      )
  )


gtex <- gtex %>%
  mutate(
    Sex =
      standardize_sex(
        Sex
      )
  )


hpa <- hpa %>%
  mutate(
    Sex =
      standardize_sex(
        Sex
      )
  )


# ==============================================================================
# 7. FLYATLAS2 SYSTEM MAPPING
#
# Conservative mappings retained from previous analysis.
# ==============================================================================

map_fly_system <- function(tissue) {

  x <- tissue %>%
    as.character() %>%
    str_to_lower()


  case_when(

    # Heart -------------------------------------------------------------

    str_detect(
      x,
      "heart|cardiac"
    ) ~
      "Heart",


    # Skeletal muscle ---------------------------------------------------

    str_detect(
      x,
      "carcass|body.?wall|skeletal.?muscle|thorax|leg|wing|haltere"
    ) ~
      "Skeletal muscle",


    # Liver-like metabolic tissues -------------------------------------

    str_detect(
      x,
      "fat.?body|oenocyte"
    ) ~
      "Liver-metabolic",


    # Digestive ---------------------------------------------------------

    str_detect(
      x,
      "midgut|hindgut|foregut|crop|gut|digest"
    ) ~
      "Digestive",


    # Renal / excretory -------------------------------------------------

    str_detect(
      x,
      "malpigh"
    ) ~
      "Renal/excretory",


    # Respiratory -------------------------------------------------------

    str_detect(
      x,
      "trache"
    ) ~
      "Respiratory",


    # Female reproductive ----------------------------------------------

    str_detect(
      x,
      "^ovary$|ovaries"
    ) ~
      "Ovary",


    # Male reproductive accessory tissue -------------------------------

    str_detect(
      x,
      "accessory.?gland|male.?reproductive|seminal.?vesicle"
    ) ~
      "Male reproductive",


    # Sensory / neural --------------------------------------------------

    str_detect(
      x,
      "antenna|brain|eye|head|sensory|proboscis|maxillary"
    ) ~
      "Sensory/neural",


    # All other tissues remain intentionally unmapped ------------------

    TRUE ~
      NA_character_
  )
}


# ==============================================================================
# 8. HUMAN SYSTEM MAPPING
#
# Same conservative classifier for GTEx and HPA.
#
# IMPORTANT exclusions:
#   Testis
#   Pancreas
#   Breast
#   Smooth muscle
#   Cerebral cortex
#   Choroid plexus
#   Salivary gland
#   Tongue
#   etc.
# ==============================================================================

map_human_system <- function(tissue) {

  x <- tissue %>%
    as.character() %>%
    str_to_lower()


  case_when(

    # Heart -------------------------------------------------------------

    str_detect(
      x,
      "heart|cardiac"
    ) ~
      "Heart",


    # Skeletal muscle only ---------------------------------------------

    str_detect(
      x,
      "skeletal.?muscle|muscle.?skeletal"
    ) ~
      "Skeletal muscle",


    # Liver / adipose ---------------------------------------------------

    str_detect(
      x,
      "liver|adipose"
    ) ~
      "Liver-metabolic",


    # Digestive ---------------------------------------------------------

    str_detect(
      x,
      paste0(
        "stomach|",
        "esophagus|oesophagus|",
        "small.?intestine|",
        "duodenum|jejunum|ileum|",
        "colon|large.?intestine|",
        "rectum|appendix|",
        "gallbladder|digestive"
      )
    ) ~
      "Digestive",


    # Renal -------------------------------------------------------------

    str_detect(
      x,
      "kidney|renal"
    ) ~
      "Renal/excretory",


    # Respiratory -------------------------------------------------------

    str_detect(
      x,
      "lung|respiratory"
    ) ~
      "Respiratory",


    # Ovary -------------------------------------------------------------

    str_detect(
      x,
      "^ovary$|ovaries"
    ) ~
      "Ovary",


    # Male reproductive accessory tissue -------------------------------

    str_detect(
      x,
      "prostate|seminal|epididym|vas.?deferens|male.?reproductive"
    ) ~
      "Male reproductive",


    # Sensory / neural --------------------------------------------------

    str_detect(
      x,
      "brain|nerve|neural|nervous|sensory"
    ) ~
      "Sensory/neural",


    # Everything else stays excluded -----------------------------------

    TRUE ~
      NA_character_
  )
}


# ==============================================================================
# 9. APPLY SYSTEM MAPPINGS
# ==============================================================================

flyatlas <- flyatlas %>%
  mutate(
    System =
      map_fly_system(
        Tissue
      )
  )


gtex <- gtex %>%
  mutate(
    System =
      map_human_system(
        Tissue
      )
  )


hpa <- hpa %>%
  mutate(
    System =
      map_human_system(
        Tissue
      )
  )


# ==============================================================================
# 10. MAPPING AUDITS
# ==============================================================================

fly_mapping_audit <- flyatlas %>%
  distinct(
    Tissue,
    System
  ) %>%
  arrange(
    System,
    Tissue
  )


gtex_mapping_audit <- gtex %>%
  distinct(
    Tissue,
    System
  ) %>%
  arrange(
    System,
    Tissue
  )


hpa_mapping_audit <- hpa %>%
  distinct(
    Tissue,
    System
  ) %>%
  arrange(
    System,
    Tissue
  )


write_csv(
  fly_mapping_audit,
  file.path(
    table_dir,
    "FlyAtlas2_tissue_system_mapping.csv"
  )
)

write_csv(
  gtex_mapping_audit,
  file.path(
    table_dir,
    "GTEx_tissue_system_mapping.csv"
  )
)

write_csv(
  hpa_mapping_audit,
  file.path(
    table_dir,
    "HPA_tissue_system_mapping.csv"
  )
)


cat("\n============================================\n")
cat("UNMAPPED FLYATLAS2 TISSUES\n")
cat("============================================\n")

fly_mapping_audit %>%
  filter(
    is.na(System)
  ) %>%
  print(
    n = Inf
  )


cat("\n============================================\n")
cat("UNMAPPED GTEx TISSUES\n")
cat("============================================\n")

gtex_mapping_audit %>%
  filter(
    is.na(System)
  ) %>%
  print(
    n = Inf
  )


cat("\n============================================\n")
cat("UNMAPPED HPA TISSUES\n")
cat("============================================\n")

hpa_mapping_audit %>%
  filter(
    is.na(System)
  ) %>%
  print(
    n = Inf
  )


# ==============================================================================
# 11. PREPARE FLYATLAS2 SYSTEM EXPRESSION
#
# Median FPKM across mapped FlyAtlas2 tissues within each system.
# ==============================================================================

fly_system <- flyatlas %>%

  filter(
    Gene %in% fly_genes,
    Sex %in% c(
      "Female",
      "Male"
    ),
    !is.na(System)
  ) %>%

  # Sex-appropriate systems only
  filter(

    (
      Sex == "Female" &
        System %in% female_systems
    ) |

    (
      Sex == "Male" &
        System %in% male_systems
    )

  ) %>%

  group_by(
    Sex,
    Gene,
    System
  ) %>%

  summarise(

    N_Source_Tissues =
      n_distinct(
        Tissue
      ),

    Median_FPKM =
      median(
        FPKM,
        na.rm = TRUE
      ),

    Mean_FPKM =
      mean(
        FPKM,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


# ==============================================================================
# 12. PREPARE GTEx SYSTEM EXPRESSION
#
# Median of sex-specific tissue medians within each system.
# ==============================================================================

gtex_system <- gtex %>%

  filter(
    Gene %in% human_genes,
    Sex %in% c(
      "Female",
      "Male"
    ),
    !is.na(System)
  ) %>%

  filter(

    (
      Sex == "Female" &
        System %in% female_systems
    ) |

    (
      Sex == "Male" &
        System %in% male_systems
    )

  ) %>%

  group_by(
    Sex,
    Gene,
    System
  ) %>%

  summarise(

    N_Source_Tissues =
      n_distinct(
        Tissue
      ),

    Median_TPM =
      median(
        Median_TPM,
        na.rm = TRUE
      ),

    Mean_Tissue_Median_TPM =
      mean(
        Median_TPM,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


# ==============================================================================
# 13. PREPARE HPA SYSTEM EXPRESSION
# ==============================================================================

hpa_system <- hpa %>%

  filter(
    Gene %in% human_genes,
    Sex %in% c(
      "Female",
      "Male"
    ),
    !is.na(System)
  ) %>%

  filter(

    (
      Sex == "Female" &
        System %in% female_systems
    ) |

    (
      Sex == "Male" &
        System %in% male_systems
    )

  ) %>%

  group_by(
    Sex,
    Gene,
    System
  ) %>%

  summarise(

    N_Source_Tissues =
      n_distinct(
        Tissue
      ),

    Median_TPM =
      median(
        Median_TPM,
        na.rm = TRUE
      ),

    Mean_Tissue_Median_TPM =
      mean(
        Median_TPM,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


# ==============================================================================
# 14. SAVE SYSTEM TABLES
# ==============================================================================

write_csv(
  fly_system,
  file.path(
    table_dir,
    "FlyAtlas2_system_expression_by_sex.csv"
  )
)

write_csv(
  gtex_system,
  file.path(
    table_dir,
    "GTEx_system_expression_by_sex.csv"
  )
)

write_csv(
  hpa_system,
  file.path(
    table_dir,
    "HPA_system_expression_by_sex.csv"
  )
)


# ==============================================================================
# 15. SYSTEM COVERAGE AUDIT
# ==============================================================================

coverage_audit <- bind_rows(

  fly_system %>%
    distinct(
      Sex,
      System
    ) %>%
    mutate(
      Dataset =
        "FlyAtlas2"
    ),

  gtex_system %>%
    distinct(
      Sex,
      System
    ) %>%
    mutate(
      Dataset =
        "GTEx"
    ),

  hpa_system %>%
    distinct(
      Sex,
      System
    ) %>%
    mutate(
      Dataset =
        "HPA"
    )

) %>%

  mutate(
    Present =
      TRUE
  ) %>%

  complete(
    Sex = c(
      "Female",
      "Male"
    ),
    System = system_order,
    Dataset = c(
      "FlyAtlas2",
      "GTEx",
      "HPA"
    ),
    fill = list(
      Present =
        FALSE
    )
  ) %>%

  arrange(
    Sex,
    System,
    Dataset
  )


write_csv(
  coverage_audit,
  file.path(
    table_dir,
    "Bulk_system_coverage_audit.csv"
  )
)


# ==============================================================================
# 16. CORRELATION FUNCTION
# ==============================================================================

calculate_bulk_correlation <- function(
    sex_value,
    fly_gene,
    human_gene,
    human_dataset
) {

  allowed_systems <- if (
    sex_value == "Female"
  ) {
    female_systems
  } else {
    male_systems
  }


  fly_profile <- fly_system %>%

    filter(
      Sex == .env$sex_value,
      Gene == .env$fly_gene,
      System %in% .env$allowed_systems
    ) %>%

    transmute(
      System,
      Fly_Value =
        Median_FPKM
    )


  if (
    human_dataset == "GTEx"
  ) {

    human_profile <- gtex_system %>%

      filter(
        Sex == .env$sex_value,
        Gene == .env$human_gene,
        System %in% .env$allowed_systems
      ) %>%

      transmute(
        System,
        Human_Value =
          Median_TPM
      )

  } else {

    human_profile <- hpa_system %>%

      filter(
        Sex == .env$sex_value,
        Gene == .env$human_gene,
        System %in% .env$allowed_systems
      ) %>%

      transmute(
        System,
        Human_Value =
          Median_TPM
      )
  }


  matched <- inner_join(
    fly_profile,
    human_profile,
    by = "System",
    relationship = "one-to-one"
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
        Dataset = human_dataset,
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
    Dataset = human_dataset,
    Sex = sex_value,
    Fly_Gene = fly_gene,
    Human_Gene = human_gene,
    N_Systems = n_systems,
    Rho =
      unname(
        test$estimate
      ),
    P_Value =
      test$p.value
  )
}


# ==============================================================================
# 17. CALCULATE ALL CORRELATIONS
#
# 6 fly genes x 6 human genes x 2 sexes x 2 datasets = 144 comparisons
# ==============================================================================

correlation_grid <- crossing(
  Dataset = c(
    "GTEx",
    "HPA"
  ),
  Sex = c(
    "Female",
    "Male"
  ),
  Fly_Gene =
    fly_genes,
  Human_Gene =
    human_genes
)


correlation_results <- correlation_grid %>%

  mutate(

    Result = pmap(
      list(
        Sex,
        Fly_Gene,
        Human_Gene,
        Dataset
      ),
      calculate_bulk_correlation
    )

  ) %>%

  select(
    Result
  ) %>%

  unnest(
    Result
  )


# ==============================================================================
# 18. FDR CORRECTION
#
# Correct separately within dataset and sex.
# ==============================================================================

correlation_results <- correlation_results %>%

  group_by(
    Dataset,
    Sex
  ) %>%

  mutate(

    FDR =
      p.adjust(
        P_Value,
        method = "BH"
      )

  ) %>%

  ungroup()


write_csv(
  correlation_results,
  file.path(
    table_dir,
    "Bulk_cross_species_Spearman_by_sex.csv"
  )
)


# ==============================================================================
# 19. SERPINE1-FOCUSED RESULTS
# ==============================================================================

serpine1_results <- correlation_results %>%

  filter(
    Human_Gene ==
      "SERPINE1"
  ) %>%

  arrange(
    Dataset,
    Sex,
    desc(
      Rho
    )
  )


write_csv(
  serpine1_results,
  file.path(
    table_dir,
    "SERPINE1_bulk_correlations_by_sex.csv"
  )
)


cat("\n============================================\n")
cat("SERPINE1 BULK CORRELATIONS\n")
cat("============================================\n")

print(
  serpine1_results,
  n = Inf
)


# ==============================================================================
# 20. SAVE MATCHED-SYSTEM DATA FOR EVERY COMPARISON
# ==============================================================================

make_matched_profile <- function(
    dataset_value,
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


  fly_profile <- fly_system %>%

    filter(
      Sex == .env$sex_value,
      Gene == .env$fly_gene,
      System %in% .env$allowed_systems
    ) %>%

    transmute(
      System,
      Fly_FPKM =
        Median_FPKM
    )


  human_source <- if (
    dataset_value == "GTEx"
  ) {
    gtex_system
  } else {
    hpa_system
  }


  human_profile <- human_source %>%

    filter(
      Sex == .env$sex_value,
      Gene == .env$human_gene,
      System %in% .env$allowed_systems
    ) %>%

    transmute(
      System,
      Human_TPM =
        Median_TPM
    )


  inner_join(
    fly_profile,
    human_profile,
    by = "System",
    relationship = "one-to-one"
  ) %>%

    mutate(
      Dataset =
        dataset_value,
      Sex =
        sex_value,
      Fly_Gene =
        fly_gene,
      Human_Gene =
        human_gene,
      .before = 1
    )
}


matched_profiles <- crossing(
  Dataset = c(
    "GTEx",
    "HPA"
  ),
  Sex = c(
    "Female",
    "Male"
  ),
  Fly_Gene =
    fly_genes,
  Human_Gene =
    human_genes
) %>%

  pmap_dfr(
    function(Dataset, Sex, Fly_Gene, Human_Gene) {
      make_matched_profile(Dataset, Sex, Fly_Gene, Human_Gene)
    }
  )


write_csv(
  matched_profiles,
  file.path(
    table_dir,
    "Bulk_cross_species_matched_system_profiles.csv"
  )
)


# ==============================================================================
# 21. SAFE Z-SCORE FUNCTION FOR PROFILE VISUALIZATION
# ==============================================================================

safe_z <- function(x) {

  result <- rep(
    NA_real_,
    length(x)
  )

  good <- !is.na(x)


  if (
    sum(good) < 2
  ) {

    return(
      result
    )
  }


  s <- sd(
    x[good],
    na.rm = TRUE
  )


  if (
    is.na(s) ||
    s == 0
  ) {

    result[good] <- 0

    return(
      result
    )
  }


  result[good] <-
    (
      x[good] -
        mean(
          x[good],
          na.rm = TRUE
        )
    ) /
    s


  result
}


# ==============================================================================
# 22. CREATE SERPINE1 PROFILE DATA
#
# Z-scores are calculated within each gene pair using only matched systems.
# ==============================================================================

serpine1_profile_data <- matched_profiles %>%

  filter(
    Human_Gene ==
      "SERPINE1"
  ) %>%

  group_by(
    Dataset,
    Sex,
    Fly_Gene
  ) %>%

  mutate(

    Fly_Z =
      safe_z(
        Fly_FPKM
      ),

    Human_Z =
      safe_z(
        Human_TPM
      )

  ) %>%

  ungroup()


# Convert to long plotting format

serpine1_profile_long <- serpine1_profile_data %>%

  select(
    Dataset,
    Sex,
    Fly_Gene,
    System,
    Fly_Z,
    Human_Z
  ) %>%

  pivot_longer(
    cols = c(
      Fly_Z,
      Human_Z
    ),
    names_to =
      "Profile",
    values_to =
      "Z"
  ) %>%

  mutate(

    Profile_Type =
      case_when(

        Profile == "Human_Z" ~
          "Human SERPINE1",

        Profile == "Fly_Z" ~
          "Drosophila serpin"
      )

  ) %>%

  left_join(

    serpine1_results %>%

      select(
        Dataset,
        Sex,
        Fly_Gene,
        N_Systems,
        Rho
      ),

    by = c(
      "Dataset",
      "Sex",
      "Fly_Gene"
    )
  ) %>%

  mutate(

    Facet_Label =
      paste0(
        Fly_Gene,
        "\n",
        "\u03C1 = ",
        sprintf(
          "%.2f",
          Rho
        ),
        "   n = ",
        N_Systems
      ),

    Fly_Gene =
      factor(
        Fly_Gene,
        levels =
          fly_genes
      ),

    Profile_Type =
      factor(
        Profile_Type,
        levels = c(
          "Human SERPINE1",
          "Drosophila serpin"
        )
      )
  )


write_csv(
  serpine1_profile_long,
  file.path(
    table_dir,
    "SERPINE1_bulk_profile_plot_data.csv"
  )
)


# ==============================================================================
# 23. PROFILE-PLOT FUNCTION
# ==============================================================================

make_bulk_profile_plot <- function(
    dat,
    dataset_value,
    sex_value
) {

  plot_data <- dat %>%

    filter(
      Dataset == .env$dataset_value,
      Sex == .env$sex_value
    )


  allowed_systems <- if (
    sex_value == "Female"
  ) {
    female_systems
  } else {
    male_systems
  }


  plot_data <- plot_data %>%

    mutate(
      System =
        factor(
          System,
          levels =
            allowed_systems
        )
    )


  facet_labels <- plot_data %>%

    distinct(
      Fly_Gene,
      Facet_Label
    )


  facet_label_vector <- setNames(
    facet_labels$Facet_Label,
    as.character(
      facet_labels$Fly_Gene
    )
  )


  ggplot(
    plot_data,
    aes(
      x = System,
      y = Z,
      group = Profile_Type,
      linetype = Profile_Type,
      shape = Profile_Type
    )
  ) +

    geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      linetype = "dotted"
    ) +

    geom_line(
      linewidth = 0.85,
      na.rm = TRUE
    ) +

    geom_point(
      size = 2.8,
      na.rm = TRUE
    ) +

    facet_wrap(
      ~ Fly_Gene,
      ncol = 2,
      labeller =
        labeller(
          Fly_Gene =
            facet_label_vector
        )
    ) +

    scale_linetype_manual(
      values = c(
        "Human SERPINE1" =
          "solid",
        "Drosophila serpin" =
          "dashed"
      )
    ) +

    scale_shape_manual(
      values = c(
        "Human SERPINE1" =
          16,
        "Drosophila serpin" =
          17
      )
    ) +

    labs(

      title =
        paste0(
          sex_value,
          " FlyAtlas2 vs ",
          dataset_value,
          " bulk tissue profiles"
        ),

      subtitle =
        "Cross-species SERPINE1 tissue-distribution similarity",

      x =
        "Functional system",

      y =
        paste0(
          "Relative tissue expression ",
          "(within-gene z-score)"
        ),

      linetype =
        NULL,

      shape =
        NULL,

      caption =
        paste0(
          "FlyAtlas2 expression is summarized as FPKM and human expression ",
          "as TPM. Values are standardized within gene across matched ",
          "functional systems before visualization. ",
          "\u03C1 is the Spearman correlation calculated from the unstandardized ",
          "system-level expression ranks; n is the number of matched systems."
        )
    ) +

    theme_classic(
      base_size = 10
    ) +

    theme(

      axis.text.x =
        element_text(
          angle = 35,
          hjust = 1,
          size = 9
        ),

      axis.text.y =
        element_text(
          size = 9
        ),

      axis.title =
        element_text(
          size = 10
        ),

      strip.text =
        element_text(
          size = 10,
          face = "bold"
        ),

      plot.title =
        element_text(
          size = 12,
          face = "bold"
        ),

      plot.subtitle =
        element_text(
          size = 10
        ),

      plot.caption =
        element_text(
          size = 8,
          hjust = 0
        ),

      legend.text =
        element_text(
          size = 9
        )
    )
}


# ==============================================================================
# 24. CREATE FOUR PRIMARY PROFILE FIGURES
# ==============================================================================

profile_conditions <- tribble(

  ~Dataset, ~Sex,      ~Filename,

  "GTEx",    "Female",
  "GTEx_Female_profiles.png",

  "GTEx",    "Male",
  "GTEx_Male_profiles.png",

  "HPA",     "Female",
  "HPA_Female_profiles.png",

  "HPA",     "Male",
  "HPA_Male_profiles.png"
)


for (
  i in seq_len(
    nrow(
      profile_conditions
    )
  )
) {

  dataset_value <-
    profile_conditions$Dataset[i]

  sex_value <-
    profile_conditions$Sex[i]

  filename <-
    profile_conditions$Filename[i]


  p <- make_bulk_profile_plot(
    dat =
      serpine1_profile_long,
    dataset_value =
      dataset_value,
    sex_value =
      sex_value
  )


  output_file <- file.path(
    figure_dir,
    filename
  )


  # Force genuine PNG output
  grDevices::png(
    filename =
      output_file,
    width =
      9.5,
    height =
      8.5,
    units =
      "in",
    res =
      600
  )

  print(
    p
  )

  grDevices::dev.off()


  cat(
    "\nSaved:",
    output_file,
    "\n"
  )
}


# ==============================================================================
# 25. CORRELATION HEATMAP FUNCTION
#
# Secondary comprehensive summary.
# No significance stars.
# ==============================================================================

make_correlation_heatmap <- function(
    dat,
    dataset_value,
    sex_value
) {

  plot_data <- dat %>%

    filter(
      Dataset == .env$dataset_value,
      Sex == .env$sex_value
    ) %>%

    mutate(

      Fly_Gene =
        factor(
          Fly_Gene,
          levels =
            rev(
              fly_genes
            )
        ),

      Human_Gene =
        factor(
          Human_Gene,
          levels =
            human_genes
        ),

      Label =
        ifelse(
          is.na(Rho),
          "NA",
          sprintf(
            "%.2f",
            Rho
          )
        )
    )


  ggplot(
    plot_data,
    aes(
      x =
        Human_Gene,
      y =
        Fly_Gene,
      fill =
        Rho
    )
  ) +

    geom_tile(
      linewidth =
        0.5
    ) +

    geom_text(
      aes(
        label =
          Label
      ),
      size =
        3.1
    ) +

    scale_fill_gradient2(
      low =
        "#3B4CC0",
      mid =
        "white",
      high =
        "#B40426",
      midpoint =
        0,
      limits =
        c(
          -1,
          1
        ),
      na.value =
        "grey90",
      name =
        "Spearman\nrho"
    ) +

    labs(

      title =
        paste0(
          sex_value,
          " FlyAtlas2 vs ",
          dataset_value
        ),

      subtitle =
        "Bulk tissue-expression pattern similarity",

      x =
        "Human SERPIN",

      y =
        "Drosophila serpin"
    ) +

    theme_minimal(
      base_size =
        10
    ) +

    theme(

      panel.grid =
        element_blank(),

      axis.text.x =
        element_text(
          angle =
            45,
          hjust =
            1,
          size =
            9
        ),

      axis.text.y =
        element_text(
          size =
            9
        ),

      plot.title =
        element_text(
          size =
            12,
          face =
            "bold"
        ),

      plot.subtitle =
        element_text(
          size =
            9
        )
    )
}


# ==============================================================================
# 26. SAVE FOUR HEATMAPS
# ==============================================================================

heatmap_conditions <- tribble(

  ~Dataset, ~Sex,      ~Filename,

  "GTEx",    "Female",
  "GTEx_Female_heatmap.png",

  "GTEx",    "Male",
  "GTEx_Male_heatmap.png",

  "HPA",     "Female",
  "HPA_Female_heatmap.png",

  "HPA",     "Male",
  "HPA_Male_heatmap.png"
)


for (
  i in seq_len(
    nrow(
      heatmap_conditions
    )
  )
) {

  dataset_value <-
    heatmap_conditions$Dataset[i]

  sex_value <-
    heatmap_conditions$Sex[i]

  filename <-
    heatmap_conditions$Filename[i]


  p <- make_correlation_heatmap(
    dat =
      correlation_results,
    dataset_value =
      dataset_value,
    sex_value =
      sex_value
  )


  output_file <- file.path(
    figure_dir,
    filename
  )


  grDevices::png(
    filename =
      output_file,
    width =
      8,
    height =
      6,
    units =
      "in",
    res =
      600
  )

  print(
    p
  )

  grDevices::dev.off()
}


# ==============================================================================
# 27. OUTPUT VERIFICATION
# ==============================================================================

png_files <- list.files(
  figure_dir,
  pattern =
    "\\.png$",
  full.names =
    TRUE,
  ignore.case =
    TRUE
)


cat("\n============================================\n")
cat("PNG FIGURES CREATED\n")
cat("============================================\n")

print(
  png_files
)


cat(
  "\nNumber of PNG files:",
  length(
    png_files
  ),
  "\n"
)


# ==============================================================================
# 28. TOP SERPINE1 RESULTS BY DATASET AND SEX
# ==============================================================================

cat("\n============================================\n")
cat("TOP SERPINE1 MATCHES\n")
cat("============================================\n")


serpine1_results %>%

  group_by(
    Dataset,
    Sex
  ) %>%

  arrange(
    desc(
      Rho
    ),
    .by_group =
      TRUE
  ) %>%

  select(
    Dataset,
    Sex,
    Fly_Gene,
    N_Systems,
    Rho,
    P_Value,
    FDR
  ) %>%

  print(
    n = Inf
  )


# ==============================================================================
# 29. FINAL MESSAGE
# ==============================================================================

cat("\n============================================\n")
cat("SEX-SPECIFIC BULK TISSUE ANALYSIS COMPLETE\n")
cat("============================================\n")

cat(
  "\nTables saved to:\n",
  table_dir,
  "\n"
)

cat(
  "\nFigures saved to:\n",
  figure_dir,
  "\n"
)
