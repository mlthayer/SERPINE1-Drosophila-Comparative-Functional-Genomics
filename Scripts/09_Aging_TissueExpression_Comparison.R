# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 09_Aging_TissueExpression_Comparison.R
# Run from the repository root.
# ============================================================
# 09_Aging_TissueExpression_Comparison_REANALYSIS_v2.R
#
# Cross-species SERPIN aging trajectory comparison
# AFCA Drosophila aging data + GTEx v11 human aging data
#
# REVISION v2
#
# Major changes:
#
# 1. Seven matched biological systems retained:
#       Metabolic
#       Digestive
#       Cardiovascular
#       Respiratory
#       Excretory
#       Secretory
#       Immune defense
#
# 2. AMP-producing cells retained and mapped to
#    Immune defense.
#
# 3. Human Spleen + Whole Blood mapped to
#    Immune defense.
#
# 4. Species ages mapped to a standardized within-species
#    adult-aging axis (0 to 1).
#
# 5. Signed aging changes retained.
#
# 6. Two different z-score frameworks are used:
#
#    OVERALL Z-SCORE:
#      calculated across all system x age points for each
#      fly-human gene pair.
#      Used only for overall multisystem similarity.
#
#    SYSTEM Z-SCORE:
#      calculated separately within each biological system.
#      Used only for system-specific similarity.
#
# 7. Spearman correlations use original signed values.
#
# 8. No inferential P-values are calculated from interpolated
#    trajectory points because these are not independent
#    biological observations.
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(purrr)
library(readr)


# ============================================================
# 2. PROJECT PATHS
# ============================================================

analysis_dir <- "."


aging_root <- file.path(
  analysis_dir,
  "Aging"
)


# ------------------------------------------------------------
# AFCA input
# ------------------------------------------------------------

afca_file <- file.path(
  aging_root,
  "Fly_Aging_tissue_expression",
  "Data",
  "Processed",
  "AFCA_serpin_aging_change_summary_table.csv"
)


# ------------------------------------------------------------
# GTEx input
# ------------------------------------------------------------

gtex_file <- file.path(
  aging_root,
  "Human_Aging",
  "Processed",
  "GTEx_v11_SERPIN_change_from_20s.csv"
)


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

output_dir <- file.path(
  aging_root,
  "CrossSpecies_Aging_Reanalysis_v2"
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


cat(
  "\nAFCA input:\n",
  afca_file,
  "\n"
)

cat(
  "\nGTEx input:\n",
  gtex_file,
  "\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)


if (!file.exists(afca_file)) {

  stop(
    "\nAFCA input file not found:\n",
    afca_file
  )
}


if (!file.exists(gtex_file)) {

  stop(
    "\nGTEx input file not found:\n",
    gtex_file
  )
}


# ============================================================
# 3. ANALYSIS SETTINGS
# ============================================================

fly_genes <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)


# These are the human SERPIN genes present in the
# GTEx aging dataset used for this analysis.

human_genes <- c(
  "SERPINB1",
  "SERPINB3",
  "SERPINB9",
  "SERPINE1",
  "SERPINI1"
)


system_order <- c(
  "Metabolic",
  "Digestive",
  "Cardiovascular",
  "Respiratory",
  "Excretory",
  "Secretory",
  "Immune defense"
)


# ------------------------------------------------------------
# Alary muscle sensitivity setting
#
# FALSE = primary analysis excludes alary muscle
# TRUE  = include alary muscle as Cardiovascular
# ------------------------------------------------------------

include_alary_muscle <- FALSE


# ------------------------------------------------------------
# Shared relative-age grid
# ------------------------------------------------------------

common_age_grid <- seq(
  0,
  1,
  by = 0.25
)


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================


mean_or_na <- function(x) {

  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(
    x,
    na.rm = TRUE
  )
}


median_or_na <- function(x) {

  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  median(
    x,
    na.rm = TRUE
  )
}


weighted_mean_or_na <- function(x, w) {

  x <- as.numeric(x)
  w <- as.numeric(w)

  keep <- (
    !is.na(x) &
      !is.na(w) &
      w > 0
  )

  if (!any(keep)) {
    return(NA_real_)
  }

  weighted.mean(
    x[keep],
    w[keep],
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# Safe z-score
# ------------------------------------------------------------

zscore_safe <- function(x) {

  x <- as.numeric(x)

  if (all(is.na(x))) {

    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  x_sd <- sd(
    x,
    na.rm = TRUE
  )

  if (
    is.na(x_sd) ||
      x_sd == 0
  ) {

    return(
      rep(
        0,
        length(x)
      )
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) / x_sd
}


# ------------------------------------------------------------
# Correlation without inferential P-value
# ------------------------------------------------------------

safe_cor <- function(
  x,
  y,
  method = "spearman"
) {

  keep <- complete.cases(
    x,
    y
  )

  x <- as.numeric(
    x[keep]
  )

  y <- as.numeric(
    y[keep]
  )

  if (length(x) < 3) {
    return(NA_real_)
  }

  if (
    sd(x) == 0 ||
      sd(y) == 0
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x,
      y,
      method = method
    )
  )
}


# ------------------------------------------------------------
# Cosine similarity
# ------------------------------------------------------------

cosine_similarity <- function(
  x,
  y
) {

  keep <- complete.cases(
    x,
    y
  )

  x <- as.numeric(
    x[keep]
  )

  y <- as.numeric(
    y[keep]
  )

  if (length(x) < 2) {
    return(NA_real_)
  }

  denominator <- sqrt(
    sum(x^2) *
      sum(y^2)
  )

  if (denominator == 0) {
    return(NA_real_)
  }

  sum(
    x * y
  ) / denominator
}


# ------------------------------------------------------------
# RMSE
# ------------------------------------------------------------

rmse <- function(
  x,
  y
) {

  keep <- complete.cases(
    x,
    y
  )

  if (!any(keep)) {
    return(NA_real_)
  }

  sqrt(
    mean(
      (
        x[keep] -
          y[keep]
      )^2
    )
  )
}


# ------------------------------------------------------------
# MAE
# ------------------------------------------------------------

mae <- function(
  x,
  y
) {

  keep <- complete.cases(
    x,
    y
  )

  if (!any(keep)) {
    return(NA_real_)
  }

  mean(
    abs(
      x[keep] -
        y[keep]
    )
  )
}


# ------------------------------------------------------------
# Direction concordance
#
# Fraction of informative points where fly and human
# aging changes have the same sign.
# ------------------------------------------------------------

direction_concordance <- function(
  x,
  y,
  tolerance = 1e-08
) {

  keep <- complete.cases(
    x,
    y
  )

  x <- x[keep]
  y <- y[keep]


  informative <- (
    abs(x) > tolerance |
      abs(y) > tolerance
  )


  x <- x[informative]
  y <- y[informative]


  if (length(x) == 0) {
    return(NA_real_)
  }


  mean(
    sign(x) ==
      sign(y)
  )
}


# ------------------------------------------------------------
# Linear interpolation
# ------------------------------------------------------------

interpolate_profile <- function(
  age,
  value,
  grid = common_age_grid
) {

  keep <- complete.cases(
    age,
    value
  )

  age <- age[keep]
  value <- value[keep]


  if (length(age) < 2) {

    return(
      tibble(
        Relative_Age = grid,
        Interpolated_Change = NA_real_
      )
    )
  }


  tmp <- tibble(
    age = age,
    value = value
  ) %>%

    group_by(
      age
    ) %>%

    summarize(
      value =
        mean_or_na(
          value
        ),
      .groups =
        "drop"
    ) %>%

    arrange(
      age
    )


  interpolated <- approx(
    x =
      tmp$age,
    y =
      tmp$value,
    xout =
      grid,
    method =
      "linear",
    rule =
      1,
    ties =
      mean
  )


  tibble(
    Relative_Age =
      interpolated$x,
    Interpolated_Change =
      interpolated$y
  )
}


# ============================================================
# 5. READ INPUT DATA
# ============================================================

afca <- read.csv(
  afca_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


gtex <- read.csv(
  gtex_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "\nAFCA rows:",
  nrow(afca),
  "\n"
)

cat(
  "GTEx rows:",
  nrow(gtex),
  "\n"
)


# ============================================================
# 6. INPUT VALIDATION
# ============================================================

required_afca <- c(
  "Gene",
  "cell_type_clean",
  "30d_vs_5d",
  "50d_vs_5d",
  "70d_vs_5d"
)


required_gtex <- c(
  "Gene",
  "SMTSD",
  "AGE",
  "Total_N",
  "change_from_20s"
)


missing_afca <- setdiff(
  required_afca,
  names(afca)
)


missing_gtex <- setdiff(
  required_gtex,
  names(gtex)
)


if (length(missing_afca) > 0) {

  stop(
    "AFCA missing required columns: ",
    paste(
      missing_afca,
      collapse = ", "
    )
  )
}


if (length(missing_gtex) > 0) {

  stop(
    "GTEx missing required columns: ",
    paste(
      missing_gtex,
      collapse = ", "
    )
  )
}


# ============================================================
# 7. GENE COVERAGE AUDIT
# ============================================================

cat(
  "\nFly genes found:\n"
)

print(
  sort(
    unique(
      afca$Gene
    )
  )
)


cat(
  "\nHuman genes found:\n"
)

print(
  sort(
    unique(
      gtex$Gene
    )
  )
)


missing_fly_genes <- setdiff(
  fly_genes,
  unique(
    afca$Gene
  )
)


missing_human_genes <- setdiff(
  human_genes,
  unique(
    gtex$Gene
  )
)


if (length(missing_fly_genes) > 0) {

  warning(
    "Missing fly genes: ",
    paste(
      missing_fly_genes,
      collapse = ", "
    )
  )
}


if (length(missing_human_genes) > 0) {

  warning(
    "Missing human genes: ",
    paste(
      missing_human_genes,
      collapse = ", "
    )
  )
}


# ============================================================
# 8. FUNCTIONAL-SYSTEM MAPPING
# ============================================================

fly_system_map <- tribble(

  ~cell_type_clean,
  ~Functional_System,
  ~Primary_Include,
  ~Mapping_Type,

  "Fat body body",
  "Metabolic",
  TRUE,
  "Functional",

  "Fat body head",
  "Metabolic",
  TRUE,
  "Functional",

  "Oenocyte",
  "Metabolic",
  TRUE,
  "Functional",

  "Midgut enterocyte",
  "Digestive",
  TRUE,
  "Functional",

  "Differentiating enterocyte",
  "Digestive",
  TRUE,
  "Functional",

  "Hindgut",
  "Digestive",
  TRUE,
  "Functional",

  "Heart muscle",
  "Cardiovascular",
  TRUE,
  "Anatomical/functional",

  "Alary muscle",
  "Cardiovascular",
  include_alary_muscle,
  "Functional proxy / sensitivity",

  "Tracheal cell",
  "Respiratory",
  TRUE,
  "Functional",

  "Malpighian tubule principal",
  "Excretory",
  TRUE,
  "Functional",

  "Malpighian tubule stellate",
  "Excretory",
  TRUE,
  "Functional",

  "Salivary gland",
  "Secretory",
  TRUE,
  "Functional",

  "Amp-producing cell",
  "Immune defense",
  TRUE,
  "Functional immune-defense proxy"
)


human_system_map <- tribble(

  ~SMTSD,
  ~Functional_System,
  ~Mapping_Type,

  "Adipose - Subcutaneous",
  "Metabolic",
  "Functional",

  "Adipose - Visceral (Omentum)",
  "Metabolic",
  "Functional",

  "Liver",
  "Metabolic",
  "Functional",

  "Small Intestine - Terminal Ileum",
  "Digestive",
  "Functional",

  "Colon - Sigmoid",
  "Digestive",
  "Functional",

  "Colon - Transverse",
  "Digestive",
  "Functional",

  "Heart - Left Ventricle",
  "Cardiovascular",
  "Anatomical/functional",

  "Heart - Atrial Appendage",
  "Cardiovascular",
  "Anatomical/functional",

  "Lung",
  "Respiratory",
  "Functional",

  "Kidney - Cortex",
  "Excretory",
  "Functional",

  "Kidney - Medulla",
  "Excretory",
  "Functional",

  "Minor Salivary Gland",
  "Secretory",
  "Functional",

  "Spleen",
  "Immune defense",
  "Systemic immune proxy",

  "Whole Blood",
  "Immune defense",
  "Systemic immune proxy"
)


write.csv(
  fly_system_map,
  file.path(
    table_dir,
    "Fly_Aging_Functional_System_Mapping.csv"
  ),
  row.names = FALSE
)


write.csv(
  human_system_map,
  file.path(
    table_dir,
    "Human_Aging_Functional_System_Mapping.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. MAPPING AUDIT
# ============================================================

unmapped_fly <- setdiff(
  unique(
    afca$cell_type_clean
  ),
  fly_system_map$cell_type_clean
)


unmapped_human <- setdiff(
  unique(
    gtex$SMTSD
  ),
  human_system_map$SMTSD
)


if (length(unmapped_fly) > 0) {

  cat(
    "\nAFCA cell types not used in mapping:\n"
  )

  print(
    unmapped_fly
  )
}


if (length(unmapped_human) > 0) {

  cat(
    "\nGTEx tissues not used in matched-system analysis:\n"
  )

  print(
    unmapped_human
  )
}


# ============================================================
# 10. PROCESS AFCA AGING DATA
# ============================================================

afca_long <- afca %>%

  filter(
    Gene %in%
      fly_genes
  ) %>%

  inner_join(
    fly_system_map,
    by =
      "cell_type_clean"
  ) %>%

  filter(
    Primary_Include
  ) %>%

  pivot_longer(

    cols =
      c(
        `30d_vs_5d`,
        `50d_vs_5d`,
        `70d_vs_5d`
      ),

    names_to =
      "Age_Contrast",

    values_to =
      "Aging_Change"

  ) %>%

  mutate(

    Fly_Age_Day = case_when(

      Age_Contrast ==
        "30d_vs_5d" ~ 30,

      Age_Contrast ==
        "50d_vs_5d" ~ 50,

      Age_Contrast ==
        "70d_vs_5d" ~ 70,

      TRUE ~
        NA_real_

    ),

    Aging_Change =
      as.numeric(
        Aging_Change
      ),

    Relative_Age =
      (
        Fly_Age_Day -
          5
      ) /
      (
        70 -
          5
      )

  )


# ============================================================
# 11. ADD FLY YOUNG-ADULT BASELINE
# ============================================================

afca_baseline <- afca_long %>%

  distinct(
    Gene,
    cell_type_clean,
    Functional_System,
    Mapping_Type
  ) %>%

  mutate(

    Age_Contrast =
      "5d_reference",

    Fly_Age_Day =
      5,

    Relative_Age =
      0,

    Aging_Change =
      0

  )


afca_long <- bind_rows(
  afca_baseline,
  afca_long
) %>%

  arrange(
    Gene,
    Functional_System,
    cell_type_clean,
    Relative_Age
  )


write.csv(
  afca_long,
  file.path(
    table_dir,
    "AFCA_Aging_CellType_Long.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. AFCA SYSTEM-LEVEL TRAJECTORIES
# ============================================================

afca_system <- afca_long %>%

  group_by(
    Gene,
    Functional_System,
    Relative_Age,
    Fly_Age_Day
  ) %>%

  summarize(

    Fly_Median_Change =
      median_or_na(
        Aging_Change
      ),

    Fly_Mean_Change =
      mean_or_na(
        Aging_Change
      ),

    N_Fly_CellTypes =
      sum(
        !is.na(
          Aging_Change
        )
      ),

    .groups =
      "drop"
  ) %>%

  mutate(

    Functional_System =
      factor(
        Functional_System,
        levels =
          system_order
      )

  ) %>%

  arrange(
    Gene,
    Functional_System,
    Relative_Age
  )


write.csv(
  afca_system,
  file.path(
    table_dir,
    "AFCA_Aging_System_Trajectories.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. PROCESS GTEx AGING DATA
# ============================================================

gtex_long <- gtex %>%

  filter(
    Gene %in%
      human_genes
  ) %>%

  inner_join(
    human_system_map,
    by =
      "SMTSD"
  ) %>%

  mutate(

    Human_Age_Midpoint = case_when(

      AGE ==
        "30-39" ~ 35,

      AGE ==
        "40-49" ~ 45,

      AGE ==
        "50-59" ~ 55,

      AGE ==
        "60-69" ~ 65,

      AGE ==
        "70-79" ~ 75,

      TRUE ~
        NA_real_

    ),

    change_from_20s =
      as.numeric(
        change_from_20s
      ),

    Total_N =
      as.numeric(
        Total_N
      ),

    Relative_Age =
      (
        Human_Age_Midpoint -
          25
      ) /
      (
        75 -
          25
      )

  ) %>%

  filter(
    !is.na(
      Relative_Age
    )
  )


# ============================================================
# 14. ADD HUMAN YOUNG-ADULT BASELINE
# ============================================================

gtex_baseline <- gtex_long %>%

  distinct(
    Gene,
    SMTSD,
    Functional_System,
    Mapping_Type
  ) %>%

  mutate(

    AGE =
      "20s_reference",

    Human_Age_Midpoint =
      25,

    Relative_Age =
      0,

    change_from_20s =
      0,

    Total_N =
      NA_real_

  )


gtex_long_with_baseline <- bind_rows(
  gtex_baseline,
  gtex_long
) %>%

  arrange(
    Gene,
    Functional_System,
    SMTSD,
    Relative_Age
  )


write.csv(
  gtex_long_with_baseline,
  file.path(
    table_dir,
    "GTEx_Aging_Tissue_Long.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. GTEx SYSTEM-LEVEL TRAJECTORIES
# ============================================================

gtex_system_observed <- gtex_long %>%

  group_by(
    Gene,
    Functional_System,
    Relative_Age,
    Human_Age_Midpoint,
    AGE
  ) %>%

  summarize(

    Human_Weighted_Mean_Change =
      weighted_mean_or_na(
        change_from_20s,
        Total_N
      ),

    Human_Median_Change =
      median_or_na(
        change_from_20s
      ),

    Human_Mean_Change =
      mean_or_na(
        change_from_20s
      ),

    N_GTEx_Tissues =
      n_distinct(
        SMTSD
      ),

    Total_GTEx_N =
      sum(
        Total_N,
        na.rm =
          TRUE
      ),

    .groups =
      "drop"
  )


gtex_system_baseline <- gtex_system_observed %>%

  distinct(
    Gene,
    Functional_System
  ) %>%

  mutate(

    Relative_Age =
      0,

    Human_Age_Midpoint =
      25,

    AGE =
      "20s_reference",

    Human_Weighted_Mean_Change =
      0,

    Human_Median_Change =
      0,

    Human_Mean_Change =
      0,

    N_GTEx_Tissues =
      NA_integer_,

    Total_GTEx_N =
      NA_real_

  )


gtex_system <- bind_rows(
  gtex_system_baseline,
  gtex_system_observed
) %>%

  mutate(

    Functional_System =
      factor(
        Functional_System,
        levels =
          system_order
      )

  ) %>%

  arrange(
    Gene,
    Functional_System,
    Relative_Age
  )


write.csv(
  gtex_system,
  file.path(
    table_dir,
    "GTEx_Aging_System_Trajectories.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. INTERPOLATE FLY TRAJECTORIES
# ============================================================

fly_interpolated <- afca_system %>%

  group_by(
    Gene,
    Functional_System
  ) %>%

  group_modify(
    ~ interpolate_profile(

      age =
        .x$Relative_Age,

      value =
        .x$Fly_Median_Change,

      grid =
        common_age_grid

    )
  ) %>%

  ungroup() %>%

  rename(
    Fly_Change =
      Interpolated_Change
  )


write.csv(
  fly_interpolated,
  file.path(
    table_dir,
    "AFCA_Aging_Interpolated_Trajectories.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. INTERPOLATE HUMAN TRAJECTORIES
# ============================================================

human_interpolated <- gtex_system %>%

  group_by(
    Gene,
    Functional_System
  ) %>%

  group_modify(
    ~ interpolate_profile(

      age =
        .x$Relative_Age,

      value =
        .x$Human_Weighted_Mean_Change,

      grid =
        common_age_grid

    )
  ) %>%

  ungroup() %>%

  rename(
    Human_Change =
      Interpolated_Change
  )


write.csv(
  human_interpolated,
  file.path(
    table_dir,
    "GTEx_Aging_Interpolated_Trajectories.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. BUILD ALL FLY × HUMAN GENE PAIRS
# ============================================================

gene_pairs <- expand_grid(

  Fly_Gene =
    fly_genes,

  Human_Gene =
    human_genes

)


# ============================================================
# 19. MATCH FLY AND HUMAN TRAJECTORIES
# ============================================================

matched_profiles <- gene_pairs %>%

  left_join(

    fly_interpolated %>%

      rename(
        Fly_Gene =
          Gene
      ),

    by =
      "Fly_Gene"

  ) %>%

  left_join(

    human_interpolated %>%

      rename(
        Human_Gene =
          Gene
      ),

    by = c(
      "Human_Gene",
      "Functional_System",
      "Relative_Age"
    )

  ) %>%

  filter(
    !is.na(
      Fly_Change
    ),
    !is.na(
      Human_Change
    )
  ) %>%

  mutate(

    Functional_System =
      factor(
        Functional_System,
        levels =
          system_order
      )

  ) %>%

  arrange(
    Fly_Gene,
    Human_Gene,
    Functional_System,
    Relative_Age
  )


# ============================================================
# 20. OVERALL GENE-PAIR Z-SCORES
#
# These are calculated across ALL matched system x age
# observations for each fly-human gene pair.
#
# Used only for overall multisystem similarity.
# ============================================================

matched_profiles <- matched_profiles %>%

  group_by(
    Fly_Gene,
    Human_Gene
  ) %>%

  mutate(

    Fly_Z_Overall =
      zscore_safe(
        Fly_Change
      ),

    Human_Z_Overall =
      zscore_safe(
        Human_Change
      )

  ) %>%

  ungroup()


# ============================================================
# 21. WITHIN-SYSTEM Z-SCORES
#
# These are calculated independently inside each:
#
# Fly_Gene × Human_Gene × Functional_System
#
# Used only for system-specific trajectory comparison.
# ============================================================

matched_profiles <- matched_profiles %>%

  group_by(
    Fly_Gene,
    Human_Gene,
    Functional_System
  ) %>%

  mutate(

    Fly_Z_System =
      zscore_safe(
        Fly_Change
      ),

    Human_Z_System =
      zscore_safe(
        Human_Change
      )

  ) %>%

  ungroup()


# ------------------------------------------------------------
# Raw pointwise diagnostics
# ------------------------------------------------------------

matched_profiles <- matched_profiles %>%

  mutate(

    Signed_Difference =
      Fly_Change -
        Human_Change,

    Absolute_Difference =
      abs(
        Fly_Change -
          Human_Change
      ),

    System_Z_Difference =
      Fly_Z_System -
        Human_Z_System,

    Overall_Z_Difference =
      Fly_Z_Overall -
        Human_Z_Overall

  )


write.csv(
  matched_profiles,
  file.path(
    table_dir,
    "CrossSpecies_Aging_Matched_Trajectories_Long_v2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. SYSTEM-SPECIFIC SIMILARITY
#
# IMPORTANT:
#
# Spearman:
#   raw signed trajectory values
#
# Pearson / cosine / RMSE / MAE:
#   within-system z-scores
#
# ============================================================

system_similarity <- matched_profiles %>%

  group_by(
    Fly_Gene,
    Human_Gene,
    Functional_System
  ) %>%

  summarize(

    N_Interpolated_Age_Points =
      sum(
        complete.cases(
          Fly_Change,
          Human_Change
        )
      ),

    # ------------------------------------------
    # Raw trajectory rank agreement
    # ------------------------------------------

    Spearman_Rho =
      safe_cor(
        Fly_Change,
        Human_Change,
        method =
          "spearman"
      ),

    # ------------------------------------------
    # Within-system standardized shape
    # ------------------------------------------

    Pearson_R_System =
      safe_cor(
        Fly_Z_System,
        Human_Z_System,
        method =
          "pearson"
      ),

    Cosine_System =
      cosine_similarity(
        Fly_Z_System,
        Human_Z_System
      ),

    Z_RMSE_System =
      rmse(
        Fly_Z_System,
        Human_Z_System
      ),

    Z_MAE_System =
      mae(
        Fly_Z_System,
        Human_Z_System
      ),

    # ------------------------------------------
    # Original signed direction
    # ------------------------------------------

    Direction_Concordance =
      direction_concordance(
        Fly_Change,
        Human_Change
      ),

    # ------------------------------------------
    # Descriptive raw-change differences
    # ------------------------------------------

    Mean_Signed_Difference =
      mean_or_na(
        Fly_Change -
          Human_Change
      ),

    Mean_Absolute_Difference =
      mean_or_na(
        abs(
          Fly_Change -
            Human_Change
        )
      ),

    .groups =
      "drop"
  ) %>%

  arrange(
    Human_Gene,
    Fly_Gene,
    Functional_System
  )


write.csv(
  system_similarity,
  file.path(
    table_dir,
    "CrossSpecies_Aging_System_Similarity_v2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. OVERALL MULTISYSTEM SIMILARITY
#
# Uses overall gene-pair z-scores.
# ============================================================

overall_similarity <- matched_profiles %>%

  group_by(
    Fly_Gene,
    Human_Gene
  ) %>%

  summarize(

    N_Systems =
      n_distinct(
        Functional_System
      ),

    N_Interpolated_Profile_Points =
      sum(
        complete.cases(
          Fly_Change,
          Human_Change
        )
      ),

    # ------------------------------------------
    # Raw rank-profile agreement
    # ------------------------------------------

    Spearman_Rho =
      safe_cor(
        Fly_Change,
        Human_Change,
        method =
          "spearman"
      ),

    # ------------------------------------------
    # Overall standardized profile agreement
    # ------------------------------------------

    Pearson_R_Overall =
      safe_cor(
        Fly_Z_Overall,
        Human_Z_Overall,
        method =
          "pearson"
      ),

    Cosine_Overall =
      cosine_similarity(
        Fly_Z_Overall,
        Human_Z_Overall
      ),

    Z_RMSE_Overall =
      rmse(
        Fly_Z_Overall,
        Human_Z_Overall
      ),

    Z_MAE_Overall =
      mae(
        Fly_Z_Overall,
        Human_Z_Overall
      ),

    # ------------------------------------------
    # Direction
    # ------------------------------------------

    Direction_Concordance =
      direction_concordance(
        Fly_Change,
        Human_Change
      ),

    .groups =
      "drop"
  )


# ============================================================
# 24. FIXED DESCRIPTIVE INTEGRATED SCORE
#
# No candidate-dependent min-max scaling.
#
# Spearman:
#       -1 = 0
#        0 = 0.5
#       +1 = 1
#
# Cosine:
#       -1 = 0
#       +1 = 1
#
# RMSE:
#       exp(-RMSE)
#
# Direction concordance:
#       already 0 to 1
#
# This is DESCRIPTIVE.
# Individual metrics should remain available.
# ============================================================

overall_similarity <- overall_similarity %>%

  mutate(

    Spearman_Similarity =
      (
        Spearman_Rho +
          1
      ) / 2,

    Cosine_Similarity_01 =
      (
        Cosine_Overall +
          1
      ) / 2,

    RMSE_Similarity =
      exp(
        -Z_RMSE_Overall
      ),

    Integrated_Aging_Similarity =
      rowMeans(
        cbind(
          Spearman_Similarity,
          Cosine_Similarity_01,
          RMSE_Similarity,
          Direction_Concordance
        ),
        na.rm =
          TRUE
      )

  ) %>%

  arrange(
    Human_Gene,
    desc(
      Integrated_Aging_Similarity
    )
  )


write.csv(
  overall_similarity,
  file.path(
    table_dir,
    "CrossSpecies_Aging_Overall_Similarity_v2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 25. RANK FLY SERPINS FOR EACH HUMAN SERPIN
# ============================================================

aging_rankings <- overall_similarity %>%

  group_by(
    Human_Gene
  ) %>%

  arrange(
    desc(
      Integrated_Aging_Similarity
    ),
    .by_group =
      TRUE
  ) %>%

  mutate(

    Aging_Rank =
      row_number()

  ) %>%

  ungroup()


write.csv(
  aging_rankings,
  file.path(
    table_dir,
    "CrossSpecies_Aging_Rankings_v2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 26. SERPINE1-SPECIFIC RANKING
# ============================================================

serpine1_results <- aging_rankings %>%

  filter(
    Human_Gene ==
      "SERPINE1"
  )


write.csv(
  serpine1_results,
  file.path(
    table_dir,
    "SERPINE1_Aging_Similarity_Ranking_v2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 27. SERPINE1 SYSTEM-SPECIFIC RESULTS
# ============================================================

serpine1_system_results <- system_similarity %>%

  filter(
    Human_Gene ==
      "SERPINE1"
  ) %>%

  arrange(
    Fly_Gene,
    Functional_System
  )


write.csv(
  serpine1_system_results,
  file.path(
    table_dir,
    "SERPINE1_Aging_System_Similarity_v2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 28. COVERAGE AUDIT — FLY
# ============================================================

fly_coverage <- afca_system %>%

  group_by(
    Gene,
    Functional_System
  ) %>%

  summarize(

    N_Observed_Ages =
      n_distinct(
        Fly_Age_Day[
          !is.na(
            Fly_Median_Change
          )
        ]
      ),

    Min_CellTypes =
      min(
        N_Fly_CellTypes,
        na.rm =
          TRUE
      ),

    Max_CellTypes =
      max(
        N_Fly_CellTypes,
        na.rm =
          TRUE
      ),

    .groups =
      "drop"
  )


write.csv(
  fly_coverage,
  file.path(
    table_dir,
    "AFCA_Aging_Coverage_Audit_v2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 29. COVERAGE AUDIT — HUMAN
# ============================================================

human_coverage <- gtex_system %>%

  filter(
    Relative_Age >
      0
  ) %>%

  group_by(
    Gene,
    Functional_System
  ) %>%

  summarize(

    N_Observed_Age_Bins =
      n_distinct(
        AGE[
          !is.na(
            Human_Weighted_Mean_Change
          )
        ]
      ),

    Total_Samples_Across_Bins =
      sum(
        Total_GTEx_N,
        na.rm =
          TRUE
      ),

    .groups =
      "drop"
  )


write.csv(
  human_coverage,
  file.path(
    table_dir,
    "GTEx_Aging_Coverage_Audit_v2.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 30. FIGURE:
#     SERPINE1 STANDARDIZED TRAJECTORIES
#
# Uses WITHIN-SYSTEM z-scores because each panel represents
# a distinct biological system.
# ============================================================

serpine1_profile <- matched_profiles %>%

  filter(
    Human_Gene ==
      "SERPINE1"
  )


p_serpine1_system <- ggplot(
  serpine1_profile,
  aes(
    x =
      Relative_Age
  )
) +

  geom_line(
    aes(
      y =
        Fly_Z_System,
      linetype =
        "Drosophila serpin"
    ),
    linewidth =
      0.7
  ) +

  geom_point(
    aes(
      y =
        Fly_Z_System
    ),
    size =
      1.5
  ) +

  geom_line(
    aes(
      y =
        Human_Z_System,
      linetype =
        "Human SERPINE1"
    ),
    linewidth =
      0.7
  ) +

  geom_point(
    aes(
      y =
        Human_Z_System
    ),
    shape =
      17,
    size =
      1.6
  ) +

  facet_grid(
    Fly_Gene ~ Functional_System
  ) +

  scale_linetype_manual(

    name =
      NULL,

    values = c(

      "Drosophila serpin" =
        "solid",

      "Human SERPINE1" =
        "dashed"

    )
  ) +

  geom_hline(

    yintercept =
      0,

    linewidth =
      0.25

  ) +

  theme_bw(
    base_size =
      9
  ) +

  labs(

    title =
      "Cross-species aging trajectories: Drosophila serpins vs SERPINE1",

    subtitle =
      "Aging trajectories standardized independently within each functional system",

    x =
      "Relative adult-aging position",

    y =
      "Within-system aging change z-score"

  ) +

  theme(

    legend.position =
      "top",

    axis.text =
      element_text(
        size =
          8
      ),

    axis.title =
      element_text(
        size =
          9
      ),

    strip.text =
      element_text(
        size =
          8
      )

  )


ggsave(
  file.path(
    figure_dir,
    "SERPINE1_Aging_Trajectories_WithinSystem_Z_v2.png"
  ),
  p_serpine1_system,
  width =
    13,
  height =
    9,
  dpi =
    600,
  bg =
    "white"
)


# ============================================================
# 31. FIGURE:
#     OVERALL SERPINE1 PROFILES
#
# Uses overall gene-pair z-score.
# ============================================================

p_serpine1_overall <- ggplot(
  serpine1_profile,
  aes(
    x =
      Relative_Age
  )
) +

  geom_line(
    aes(
      y =
        Fly_Z_Overall,
      group =
        interaction(
          Fly_Gene,
          Functional_System
        ),
      linetype =
        "Drosophila serpin"
    ),
    linewidth =
      0.7
  ) +

  geom_point(
    aes(
      y =
        Fly_Z_Overall
    ),
    size =
      1.4
  ) +

  geom_line(
    aes(
      y =
        Human_Z_Overall,
      group =
        interaction(
          Fly_Gene,
          Functional_System
        ),
      linetype =
        "Human SERPINE1"
    ),
    linewidth =
      0.7
  ) +

  geom_point(
    aes(
      y =
        Human_Z_Overall
    ),
    shape =
      17,
    size =
      1.5
  ) +

  facet_grid(
    Fly_Gene ~ Functional_System
  ) +

  scale_linetype_manual(

    name =
      NULL,

    values = c(

      "Drosophila serpin" =
        "solid",

      "Human SERPINE1" =
        "dashed"

    )
  ) +

  geom_hline(
    yintercept =
      0,
    linewidth =
      0.25
  ) +

  theme_bw(
    base_size =
      9
  ) +

  labs(

    title =
      "Overall cross-species aging profiles: Drosophila serpins vs SERPINE1",

    subtitle =
      "Standardization across all matched system-by-age observations within each gene pair",

    x =
      "Relative adult-aging position",

    y =
      "Overall aging-profile z-score"

  ) +

  theme(

    legend.position =
      "top",

    axis.text =
      element_text(
        size =
          8
      ),

    axis.title =
      element_text(
        size =
          9
      ),

    strip.text =
      element_text(
        size =
          8
      )

  )


ggsave(
  file.path(
    figure_dir,
    "SERPINE1_Aging_Trajectories_Overall_Z_v2.png"
  ),
  p_serpine1_overall,
  width =
    13,
  height =
    9,
  dpi =
    600,
  bg =
    "white"
)


# ============================================================
# 32. FIGURE:
#     OVERALL AGING SIMILARITY RANKINGS
# ============================================================

rank_plot_data <- aging_rankings %>%

  mutate(

    Fly_Gene =
      factor(
        Fly_Gene,
        levels =
          fly_genes
      ),

    Human_Gene =
      factor(
        Human_Gene,
        levels =
          human_genes
      )

  )


p_rank <- ggplot(
  rank_plot_data,
  aes(
    x =
      Fly_Gene,
    y =
      Integrated_Aging_Similarity
  )
) +

  geom_col(
    width =
      0.7
  ) +

  facet_wrap(
    ~Human_Gene,
    ncol =
      3
  ) +

  coord_cartesian(
    ylim =
      c(
        0,
        1
      )
  ) +

  theme_bw(
    base_size =
      10
  ) +

  labs(

    title =
      "Cross-species SERPIN aging-profile similarity",

    subtitle =
      "Descriptive integration of aging trajectory shape and direction",

    x =
      "Drosophila serpin",

    y =
      "Integrated aging similarity"

  ) +

  theme(

    axis.text.x =
      element_text(
        angle =
          45,
        hjust =
          1,
        size =
          8
      ),

    axis.text.y =
      element_text(
        size =
          8
      ),

    axis.title =
      element_text(
        size =
          9
      ),

    strip.text =
      element_text(
        face =
          "bold",
        size =
          9
      )

  )


ggsave(
  file.path(
    figure_dir,
    "CrossSpecies_SERPIN_Aging_Similarity_v2.png"
  ),
  p_rank,
  width =
    10,
  height =
    7,
  dpi =
    600,
  bg =
    "white"
)


# ============================================================
# 33. SESSION INFORMATION
# ============================================================

capture.output(
  sessionInfo(),
  file =
    file.path(
      output_dir,
      "R_sessionInfo_v2.txt"
    )
)


# ============================================================
# 34. FINAL CONSOLE SUMMARY
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "AGING REANALYSIS v2 COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

cat(
  "\nPrimary tables:\n"
)

cat(
  "\n1. CrossSpecies_Aging_Matched_Trajectories_Long_v2.csv"
)

cat(
  "\n2. CrossSpecies_Aging_System_Similarity_v2.csv"
)

cat(
  "\n3. CrossSpecies_Aging_Overall_Similarity_v2.csv"
)

cat(
  "\n4. CrossSpecies_Aging_Rankings_v2.csv"
)

cat(
  "\n5. SERPINE1_Aging_Similarity_Ranking_v2.csv"
)

cat(
  "\n6. SERPINE1_Aging_System_Similarity_v2.csv\n"
)

cat(
  "\nNo inferential P-values were calculated from interpolated trajectory points.\n"
)

cat(
  "\nDone.\n"
)
