# SERPINE1-Drosophila Comparative Functional Genomics Workflow
# 09_Aging_TissueExpression_Comparison.R
# Run from the repository root.
# ============================================================
# Sex-stratified aging expression analysis
#
# This script combines:
#   A. Drosophila AFCA aging preprocessing
#   B. Human GTEx aging preprocessing
#   C. Same-sex cross-species trajectory comparison
#
# Female fly profiles are compared only with female human profiles,
# and male fly profiles are compared only with male human profiles.
#
# Inputs are supplied with environment variables or selected
# interactively. No user-specific paths are stored in this script.
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(purrr)
library(readr)

# ------------------------------------------------------------
# 1. Input/output configuration
# ------------------------------------------------------------

resolve_input <- function(label, variable) {
  chosen <- Sys.getenv(variable, unset = "")

  if (!nzchar(chosen)) {
    if (!interactive()) {
      stop(
        "Missing input configuration: set ", variable,
        " to the path for ", label, ".",
        call. = FALSE
      )
    }
    message("Select file for: ", label)
    chosen <- file.choose()
  }

  if (!file.exists(chosen) || dir.exists(chosen)) {
    stop("Input file not found for ", variable, ": ", chosen, call. = FALSE)
  }

  normalizePath(chosen, winslash = "/", mustWork = TRUE)
}

afca_input <- resolve_input(
  "AFCA serpin age-by-sex cell-type summary",
  "SERPINE1_AFCA_AGING_INPUT"
)

gtex_input <- resolve_input(
  "GTEx SERPIN age-by-sex tissue summary",
  "SERPINE1_GTEX_AGING_INPUT"
)

default_output_dir <- file.path(
  ".",
  "Aging",
  "Sex_Stratified_Aging_Analysis"
)

output_dir <- Sys.getenv(
  "SERPINE1_AGING_OUTPUT",
  unset = default_output_dir
)

table_dir <- file.path(output_dir, "Tables")
figure_dir <- file.path(output_dir, "Figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("AFCA input: ", afca_input)
message("GTEx input: ", gtex_input)
message("Output directory: ", output_dir)

# ------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------

fly_genes <- c(
  "Spn42Da", "Spn42Db", "Spn42Dc",
  "Spn42Dd", "Spn28B", "Spn28F"
)

human_genes <- c(
  "SERPINB1", "SERPINB3", "SERPINB9",
  "SERPINE1", "SERPINI1"
)

sex_order <- c("Female", "Male")

system_order <- c(
  "Metabolic",
  "Digestive",
  "Cardiovascular",
  "Respiratory",
  "Excretory",
  "Secretory",
  "Immune defense"
)

common_age_grid <- seq(0, 1, by = 0.25)

age_midpoints <- c(
  "20-29" = 25,
  "30-39" = 35,
  "40-49" = 45,
  "50-59" = 55,
  "60-69" = 65,
  "70-79" = 75
)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

mean_or_na <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

median_or_na <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

weighted_mean_or_na <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  weighted.mean(x[keep], w[keep], na.rm = TRUE)
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- complete.cases(x, y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])

  if (length(x) < 3) return(NA_real_)
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)

  suppressWarnings(cor(x, y, method = method))
}

cosine_similarity <- function(x, y) {
  keep <- complete.cases(x, y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])

  if (length(x) < 2) return(NA_real_)

  d <- sqrt(sum(x^2) * sum(y^2))
  if (d == 0) return(NA_real_)

  sum(x * y) / d
}

rmse <- function(x, y) {
  keep <- complete.cases(x, y)
  if (!any(keep)) return(NA_real_)
  sqrt(mean((x[keep] - y[keep])^2))
}

mae <- function(x, y) {
  keep <- complete.cases(x, y)
  if (!any(keep)) return(NA_real_)
  mean(abs(x[keep] - y[keep]))
}

direction_concordance <- function(x, y, tolerance = 1e-8) {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  informative <- abs(x) > tolerance | abs(y) > tolerance
  x <- x[informative]
  y <- y[informative]

  if (length(x) == 0) return(NA_real_)

  mean(sign(x) == sign(y))
}

interpolate_profile <- function(age, value, grid = common_age_grid) {
  keep <- complete.cases(age, value)
  age <- as.numeric(age[keep])
  value <- as.numeric(value[keep])

  if (length(age) < 2) {
    return(
      tibble(
        Relative_Age = grid,
        Interpolated_Change = NA_real_
      )
    )
  }

  tmp <- tibble(age = age, value = value) %>%
    group_by(age) %>%
    summarise(
      value = mean_or_na(value),
      .groups = "drop"
    ) %>%
    arrange(age)

  interpolated <- approx(
    x = tmp$age,
    y = tmp$value,
    xout = grid,
    method = "linear",
    rule = 1,
    ties = mean
  )

  tibble(
    Relative_Age = interpolated$x,
    Interpolated_Change = interpolated$y
  )
}

check_required_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      "\nAvailable columns:\n",
      paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
}

# ------------------------------------------------------------
# 4. Functional-system mappings
# ------------------------------------------------------------

fly_system_map <- tribble(
  ~cell_type, ~cell_type_clean, ~Functional_System,
  "adult fat body_body", "Fat body body", "Metabolic",
  "adult fat body_head", "Fat body head", "Metabolic",
  "adult oenocyte", "Oenocyte", "Metabolic",
  "adult midgut enterocyte", "Midgut enterocyte", "Digestive",
  "adult differentiating enterocyte", "Differentiating enterocyte", "Digestive",
  "adult hindgut", "Hindgut", "Digestive",
  "cardiomyocyte, working adult heart muscle (non-ostia)", "Heart muscle", "Cardiovascular",
  "adult heart ventral longitudinal muscle", "Heart muscle", "Cardiovascular",
  "adult tracheal cell", "Tracheal cell", "Respiratory",
  "adult Malpighian tubule principal cell", "Malpighian tubule principal", "Excretory",
  "adult Malpighian tubule principal cell of initial segment", "Malpighian tubule principal", "Excretory",
  "adult Malpighian tubule principal cell of lower segment", "Malpighian tubule principal", "Excretory",
  "adult Malpighian tubule principal cell of lower ureter", "Malpighian tubule principal", "Excretory",
  "adult Malpighian tubule stellate cell of main segment", "Malpighian tubule stellate", "Excretory",
  "adult salivary gland", "Salivary gland", "Secretory",
  "antimicrobial peptide-producing cell", "AMP-producing cell", "Immune defense"
)

human_system_map <- tribble(
  ~SMTSD, ~Functional_System,
  "Adipose - Subcutaneous", "Metabolic",
  "Adipose - Visceral (Omentum)", "Metabolic",
  "Liver", "Metabolic",
  "Small Intestine - Terminal Ileum", "Digestive",
  "Colon - Sigmoid", "Digestive",
  "Colon - Transverse", "Digestive",
  "Heart - Left Ventricle", "Cardiovascular",
  "Heart - Atrial Appendage", "Cardiovascular",
  "Lung", "Respiratory",
  "Kidney - Cortex", "Excretory",
  "Kidney - Medulla", "Excretory",
  "Minor Salivary Gland", "Secretory",
  "Spleen", "Immune defense",
  "Whole Blood", "Immune defense"
)

write_csv(
  fly_system_map,
  file.path(table_dir, "Fly_Aging_Functional_System_Mapping.csv")
)

write_csv(
  human_system_map,
  file.path(table_dir, "Human_Aging_Functional_System_Mapping.csv")
)

# ============================================================
# A. DROSOPHILA AFCA PREPROCESSING
# ============================================================

afca_raw <- read_csv(afca_input, show_col_types = FALSE)

check_required_columns(
  afca_raw,
  c(
    "Gene", "age", "sex", "tissue", "cell_type_broad", "cell_type",
    "mean_expression", "median_expression", "percent_expressing", "n_cells"
  ),
  "AFCA aging input"
)

afca <- afca_raw %>%
  mutate(
    Sex = case_when(
      str_to_lower(as.character(sex)) == "female" ~ "Female",
      str_to_lower(as.character(sex)) == "male" ~ "Male",
      TRUE ~ "Other"
    ),
    Age_Day = as.numeric(age)
  ) %>%
  filter(
    Gene %in% fly_genes,
    Sex %in% sex_order,
    Age_Day %in% c(5, 30, 50, 70)
  )

afca_mapping_audit <- afca %>%
  distinct(cell_type) %>%
  anti_join(fly_system_map, by = "cell_type") %>%
  arrange(cell_type)

write_csv(
  afca_mapping_audit,
  file.path(table_dir, "AFCA_Unmapped_Cell_Types.csv")
)

afca_mapped <- afca %>%
  inner_join(fly_system_map, by = "cell_type")

afca_celltype_age <- afca_mapped %>%
  group_by(
    Gene, Sex, Age_Day,
    Functional_System, cell_type_clean
  ) %>%
  summarise(
    Mean_Expression = weighted_mean_or_na(mean_expression, n_cells),
    Median_Expression = weighted_mean_or_na(median_expression, n_cells),
    Percent_Expressing = weighted_mean_or_na(percent_expressing, n_cells),
    Total_Cells = sum(n_cells, na.rm = TRUE),
    N_Source_Rows = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Log2_Mean_Expression = log2(Mean_Expression + 1)
  )

afca_baseline <- afca_celltype_age %>%
  filter(Age_Day == 5) %>%
  select(
    Gene,
    Sex,
    Functional_System,
    cell_type_clean,
    Baseline_5d_Log2_Mean_Expression = Log2_Mean_Expression,
    Baseline_5d_Mean_Expression = Mean_Expression
  )

afca_change <- afca_celltype_age %>%
  left_join(
    afca_baseline,
    by = c(
      "Gene", "Sex",
      "Functional_System", "cell_type_clean"
    )
  ) %>%
  mutate(
    Aging_Change =
      Log2_Mean_Expression -
      Baseline_5d_Log2_Mean_Expression,
    Age_Contrast = if_else(
      Age_Day == 5,
      "5d_reference",
      paste0(Age_Day, "d_vs_5d")
    ),
    Relative_Age = (Age_Day - 5) / (70 - 5)
  ) %>%
  arrange(
    Gene, Sex, Functional_System,
    cell_type_clean, Age_Day
  )

afca_coverage <- afca_change %>%
  group_by(
    Gene, Sex, Functional_System, cell_type_clean
  ) %>%
  summarise(
    N_Ages = n_distinct(Age_Day[!is.na(Aging_Change)]),
    .groups = "drop"
  )

write_csv(
  afca_celltype_age,
  file.path(table_dir, "AFCA_Expression_BySex_Age_CellType.csv")
)

write_csv(
  afca_change,
  file.path(table_dir, "AFCA_Aging_Change_BySex.csv")
)

write_csv(
  afca_coverage,
  file.path(table_dir, "AFCA_Aging_Coverage_CellType_BySex.csv")
)

# ============================================================
# B. HUMAN GTEx PREPROCESSING
# ============================================================

gtex_raw <- read_csv(gtex_input, show_col_types = FALSE)

check_required_columns(
  gtex_raw,
  c(
    "Gene", "SMTSD", "AGE", "SEX",
    "Mean_TPM", "Median_TPM",
    "SD", "SE", "Mean_log2TPM", "N"
  ),
  "GTEx aging input"
)

gtex <- gtex_raw %>%
  mutate(
    Sex = case_when(
      str_to_lower(as.character(SEX)) == "female" ~ "Female",
      str_to_lower(as.character(SEX)) == "male" ~ "Male",
      TRUE ~ as.character(SEX)
    ),
    Human_Age_Midpoint = unname(age_midpoints[AGE])
  ) %>%
  filter(
    Gene %in% human_genes,
    Sex %in% sex_order,
    !is.na(Human_Age_Midpoint)
  )

gtex_mapping_audit <- gtex %>%
  distinct(SMTSD) %>%
  anti_join(human_system_map, by = "SMTSD") %>%
  arrange(SMTSD)

write_csv(
  gtex_mapping_audit,
  file.path(table_dir, "GTEx_Unmapped_Tissues.csv")
)

gtex_mapped <- gtex %>%
  inner_join(human_system_map, by = "SMTSD")

gtex_baseline <- gtex_mapped %>%
  filter(AGE == "20-29") %>%
  select(
    Gene,
    SMTSD,
    Sex,
    Baseline_20s_Mean_log2TPM = Mean_log2TPM,
    Baseline_20s_Mean_TPM = Mean_TPM,
    Baseline_20s_N = N
  )

gtex_change <- gtex_mapped %>%
  left_join(
    gtex_baseline,
    by = c("Gene", "SMTSD", "Sex")
  ) %>%
  mutate(
    change_from_20s =
      as.numeric(Mean_log2TPM) -
      as.numeric(Baseline_20s_Mean_log2TPM),
    Relative_Age =
      (Human_Age_Midpoint - 25) / (75 - 25),
    Total_N = as.numeric(N)
  ) %>%
  arrange(
    Gene, Sex, Functional_System,
    SMTSD, Human_Age_Midpoint
  )

gtex_coverage <- gtex_change %>%
  group_by(
    Gene, Sex, Functional_System, SMTSD
  ) %>%
  summarise(
    N_Age_Bins =
      n_distinct(AGE[!is.na(change_from_20s)]),
    Total_Samples =
      sum(Total_N, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  gtex_change,
  file.path(table_dir, "GTEx_Aging_Change_BySex.csv")
)

write_csv(
  gtex_coverage,
  file.path(table_dir, "GTEx_Aging_Coverage_Tissue_BySex.csv")
)

# ============================================================
# C. SAME-SEX CROSS-SPECIES AGING COMPARISON
# ============================================================

fly_system <- afca_change %>%
  group_by(
    Gene, Sex, Functional_System,
    Relative_Age, Age_Day
  ) %>%
  summarise(
    Fly_Median_Change = median_or_na(Aging_Change),
    Fly_Mean_Change = mean_or_na(Aging_Change),
    N_Fly_CellTypes =
      n_distinct(cell_type_clean[!is.na(Aging_Change)]),
    .groups = "drop"
  ) %>%
  mutate(
    Functional_System =
      factor(Functional_System, levels = system_order),
    Sex = factor(Sex, levels = sex_order)
  ) %>%
  arrange(Gene, Sex, Functional_System, Relative_Age)

human_system <- gtex_change %>%
  group_by(
    Gene, Sex, Functional_System,
    Relative_Age, Human_Age_Midpoint, AGE
  ) %>%
  summarise(
    Human_Weighted_Mean_Change =
      weighted_mean_or_na(change_from_20s, Total_N),
    Human_Median_Change =
      median_or_na(change_from_20s),
    Human_Mean_Change =
      mean_or_na(change_from_20s),
    N_GTEx_Tissues =
      n_distinct(SMTSD[!is.na(change_from_20s)]),
    Total_GTEx_N =
      sum(Total_N, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Functional_System =
      factor(Functional_System, levels = system_order),
    Sex = factor(Sex, levels = sex_order)
  ) %>%
  arrange(Gene, Sex, Functional_System, Relative_Age)

write_csv(
  fly_system,
  file.path(table_dir, "AFCA_Aging_System_Trajectories_BySex.csv")
)

write_csv(
  human_system,
  file.path(table_dir, "GTEx_Aging_System_Trajectories_BySex.csv")
)

fly_interpolated <- fly_system %>%
  group_by(Gene, Sex, Functional_System) %>%
  group_modify(
    ~interpolate_profile(
      .x$Relative_Age,
      .x$Fly_Median_Change,
      common_age_grid
    )
  ) %>%
  ungroup() %>%
  rename(Fly_Change = Interpolated_Change)

human_interpolated <- human_system %>%
  group_by(Gene, Sex, Functional_System) %>%
  group_modify(
    ~interpolate_profile(
      .x$Relative_Age,
      .x$Human_Weighted_Mean_Change,
      common_age_grid
    )
  ) %>%
  ungroup() %>%
  rename(Human_Change = Interpolated_Change)

write_csv(
  fly_interpolated,
  file.path(table_dir, "AFCA_Aging_Interpolated_BySex.csv")
)

write_csv(
  human_interpolated,
  file.path(table_dir, "GTEx_Aging_Interpolated_BySex.csv")
)

gene_pairs <- expand_grid(
  Fly_Gene = fly_genes,
  Human_Gene = human_genes,
  Sex = sex_order
)

matched_profiles <- gene_pairs %>%
  left_join(
    fly_interpolated %>% rename(Fly_Gene = Gene),
    by = c("Fly_Gene", "Sex")
  ) %>%
  left_join(
    human_interpolated %>% rename(Human_Gene = Gene),
    by = c(
      "Human_Gene", "Sex",
      "Functional_System", "Relative_Age"
    )
  ) %>%
  filter(
    !is.na(Fly_Change),
    !is.na(Human_Change)
  ) %>%
  mutate(
    Functional_System =
      factor(Functional_System, levels = system_order),
    Sex = factor(Sex, levels = sex_order)
  ) %>%
  arrange(
    Sex, Fly_Gene, Human_Gene,
    Functional_System, Relative_Age
  )

matched_profiles <- matched_profiles %>%
  group_by(Fly_Gene, Human_Gene, Sex) %>%
  mutate(
    Fly_Z_Overall = zscore_safe(Fly_Change),
    Human_Z_Overall = zscore_safe(Human_Change)
  ) %>%
  ungroup() %>%
  group_by(
    Fly_Gene, Human_Gene,
    Sex, Functional_System
  ) %>%
  mutate(
    Fly_Z_System = zscore_safe(Fly_Change),
    Human_Z_System = zscore_safe(Human_Change)
  ) %>%
  ungroup() %>%
  mutate(
    Signed_Difference = Fly_Change - Human_Change,
    Absolute_Difference = abs(Fly_Change - Human_Change),
    System_Z_Difference = Fly_Z_System - Human_Z_System,
    Overall_Z_Difference = Fly_Z_Overall - Human_Z_Overall
  )

write_csv(
  matched_profiles,
  file.path(
    table_dir,
    "CrossSpecies_Aging_Matched_Trajectories_BySex.csv"
  )
)

system_similarity <- matched_profiles %>%
  group_by(
    Fly_Gene, Human_Gene,
    Sex, Functional_System
  ) %>%
  summarise(
    N_Interpolated_Age_Points =
      sum(complete.cases(Fly_Change, Human_Change)),
    Spearman_Rho =
      safe_cor(Fly_Change, Human_Change, "spearman"),
    Pearson_R_System =
      safe_cor(Fly_Z_System, Human_Z_System, "pearson"),
    Cosine_System =
      cosine_similarity(Fly_Z_System, Human_Z_System),
    Z_RMSE_System =
      rmse(Fly_Z_System, Human_Z_System),
    Z_MAE_System =
      mae(Fly_Z_System, Human_Z_System),
    Direction_Concordance =
      direction_concordance(Fly_Change, Human_Change),
    Mean_Signed_Difference =
      mean_or_na(Fly_Change - Human_Change),
    Mean_Absolute_Difference =
      mean_or_na(abs(Fly_Change - Human_Change)),
    .groups = "drop"
  ) %>%
  arrange(
    Human_Gene, Sex,
    Fly_Gene, Functional_System
  )

overall_similarity <- matched_profiles %>%
  group_by(Fly_Gene, Human_Gene, Sex) %>%
  summarise(
    N_Systems = n_distinct(Functional_System),
    N_Interpolated_Profile_Points =
      sum(complete.cases(Fly_Change, Human_Change)),
    Spearman_Rho =
      safe_cor(Fly_Change, Human_Change, "spearman"),
    Pearson_R_Overall =
      safe_cor(Fly_Z_Overall, Human_Z_Overall, "pearson"),
    Cosine_Overall =
      cosine_similarity(Fly_Z_Overall, Human_Z_Overall),
    Z_RMSE_Overall =
      rmse(Fly_Z_Overall, Human_Z_Overall),
    Z_MAE_Overall =
      mae(Fly_Z_Overall, Human_Z_Overall),
    Direction_Concordance =
      direction_concordance(Fly_Change, Human_Change),
    .groups = "drop"
  ) %>%
  mutate(
    Spearman_Similarity = (Spearman_Rho + 1) / 2,
    Cosine_Similarity_01 = (Cosine_Overall + 1) / 2,
    RMSE_Similarity = exp(-Z_RMSE_Overall),
    Integrated_Aging_Similarity = apply(
      cbind(
        Spearman_Similarity,
        Cosine_Similarity_01,
        RMSE_Similarity,
        Direction_Concordance
      ),
      1,
      mean_or_na
    )
  ) %>%
  group_by(Human_Gene, Sex) %>%
  mutate(
    Aging_Rank = rank(
      -Integrated_Aging_Similarity,
      ties.method = "min",
      na.last = "keep"
    )
  ) %>%
  ungroup() %>%
  arrange(Human_Gene, Sex, Aging_Rank, Fly_Gene)

write_csv(
  system_similarity,
  file.path(
    table_dir,
    "CrossSpecies_Aging_System_Similarity_BySex.csv"
  )
)

write_csv(
  overall_similarity,
  file.path(
    table_dir,
    "CrossSpecies_Aging_Overall_Similarity_BySex.csv"
  )
)

# ------------------------------------------------------------
# SERPINE1-focused exports
# ------------------------------------------------------------

serpine1_results <- overall_similarity %>%
  filter(Human_Gene == "SERPINE1") %>%
  arrange(Sex, Aging_Rank, Fly_Gene)

serpine1_system_results <- system_similarity %>%
  filter(Human_Gene == "SERPINE1") %>%
  arrange(Sex, Fly_Gene, Functional_System)

write_csv(
  serpine1_results,
  file.path(
    table_dir,
    "SERPINE1_Aging_Similarity_Ranking_BySex.csv"
  )
)

write_csv(
  serpine1_system_results,
  file.path(
    table_dir,
    "SERPINE1_Aging_System_Similarity_BySex.csv"
  )
)

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------

p_rank <- serpine1_results %>%
  mutate(
    Sex = factor(Sex, levels = sex_order)
  ) %>%
  ggplot(
    aes(
      x = Integrated_Aging_Similarity,
      y = reorder(Fly_Gene, Integrated_Aging_Similarity)
    )
  ) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%.3f", Integrated_Aging_Similarity)),
    hjust = -0.15,
    size = 3
  ) +
  facet_wrap(~Sex, nrow = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_bw(base_size = 10) +
  labs(
    title = "Integrated cross-species aging similarity to human SERPINE1",
    subtitle = "Female and male comparisons are calculated separately",
    x = "Integrated aging similarity score",
    y = NULL
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

ggsave(
  file.path(
    figure_dir,
    "SERPINE1_Aging_Similarity_BySex.png"
  ),
  p_rank,
  width = 9,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

serpine1_profile <- matched_profiles %>%
  filter(Human_Gene == "SERPINE1")

serpine1_human <- serpine1_profile %>%
  distinct(
    Sex,
    Functional_System,
    Relative_Age,
    Human_Z_System
  )

p_trajectory <- ggplot() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    color = "grey60"
  ) +
  geom_line(
    data = serpine1_profile,
    aes(
      x = Relative_Age,
      y = Fly_Z_System,
      color = Fly_Gene,
      group = Fly_Gene
    ),
    linewidth = 0.7
  ) +
  geom_point(
    data = serpine1_profile,
    aes(
      x = Relative_Age,
      y = Fly_Z_System,
      color = Fly_Gene
    ),
    size = 1.2
  ) +
  geom_line(
    data = serpine1_human,
    aes(
      x = Relative_Age,
      y = Human_Z_System
    ),
    color = "black",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  facet_grid(Sex ~ Functional_System) +
  theme_bw(base_size = 8) +
  labs(
    title = "Cross-species aging trajectories compared with human SERPINE1",
    subtitle = "Colored lines: Drosophila serpins; dashed black line: human SERPINE1",
    x = "Relative adult-aging position",
    y = "Within-system standardized aging change",
    color = "Drosophila serpin"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    figure_dir,
    "SERPINE1_Aging_Trajectories_BySex.png"
  ),
  p_trajectory,
  width = 14,
  height = 7.5,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# Coverage audit and reproducibility record
# ------------------------------------------------------------

fly_system_coverage <- fly_system %>%
  group_by(Gene, Sex, Functional_System) %>%
  summarise(
    N_Observed_Ages =
      n_distinct(Age_Day[!is.na(Fly_Median_Change)]),
    Min_CellTypes =
      min(N_Fly_CellTypes, na.rm = TRUE),
    Max_CellTypes =
      max(N_Fly_CellTypes, na.rm = TRUE),
    .groups = "drop"
  )

human_system_coverage <- human_system %>%
  group_by(Gene, Sex, Functional_System) %>%
  summarise(
    N_Observed_Age_Bins =
      n_distinct(AGE[!is.na(Human_Weighted_Mean_Change)]),
    Total_Samples_Across_Bins =
      sum(Total_GTEx_N, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  fly_system_coverage,
  file.path(table_dir, "AFCA_Aging_Coverage_BySex.csv")
)

write_csv(
  human_system_coverage,
  file.path(table_dir, "GTEx_Aging_Coverage_BySex.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("SEX-STRATIFIED CROSS-SPECIES AGING ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nPrimary downstream input:\n")
cat(
  file.path(
    table_dir,
    "CrossSpecies_Aging_Overall_Similarity_BySex.csv"
  ),
  "\n"
)
cat(
  "\nNo inferential P-values are calculated from interpolated trajectory points.\n"
)
