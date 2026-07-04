# ============================================================
# 09_Aging_TissueExpression_Comparison.R
# Matched tissue aging comparison: AFCA + GTEx SERPINE1
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)

# Prefer an environment variable for reproducibility across machines.
# If it is not set, fall back to the original local project path.
project_dir <- Sys.getenv(
  "xxxxxx_xxxxxxxx_PROJECT_DIR",
  unset = "D:/xxxxxx-xxxxxxxx xxx/xxxxxx-xxxxxxxx xxx"
)

chapter3_dir <- file.path(project_dir, "Bioinformatics chapter three")

expression_dir <- file.path(
  project_dir,
  "Bioinformatics_tissue expression collection",
  "Bioinformatics_tissue_expression"
)

aging_dir <- file.path(expression_dir, "Aging_tissue_expression")

human_aging_dir <- file.path(
  expression_dir,
  "Data",
  "Human_Aging"
)

integrated_dir <- file.path(chapter3_dir, "Integrated_Evidence")

integrated_processed_dir <- file.path(
  integrated_dir,
  "Data",
  "Processed"
)

fig_dir <- file.path(integrated_dir, "Figures", "Aging")
tables_dir <- file.path(integrated_dir, "Tables", "Aging")

dir.create(integrated_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

selected_serpins <- c(
  "Spn42Da", "Spn42Db", "Spn42Dc",
  "Spn42Dd", "Spn28B", "Spn28F"
)

if (length(selected_serpins) == 0) {
  stop("selected_serpins is empty. Add at least one candidate gene before running aging comparison.")
}

afca_aging_file <- file.path(
  aging_dir,
  "Data",
  "Processed",
  "AFCA_serpin_aging_change_summary_table.csv"
)

gtex_aging_file <- file.path(
  human_aging_dir,
  "Processed",
  "GTEx_v11_SERPIN_change_from_20s.csv"
)

required_files <- c(
  AFCA = afca_aging_file,
  GTEx = gtex_aging_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required aging file(s): ",
    paste(names(missing_files), basename(missing_files), sep = " = ", collapse = "; ")
  )
}

afca_aging <- read.csv(afca_aging_file, stringsAsFactors = FALSE)
gtex_aging <- read.csv(gtex_aging_file, stringsAsFactors = FALSE)

required_afca_columns <- c("Gene", "cell_type_clean")
required_gtex_columns <- c("Gene", "SMTSD", "AGE", "change_from_20s")

missing_afca_columns <- setdiff(required_afca_columns, names(afca_aging))
missing_gtex_columns <- setdiff(required_gtex_columns, names(gtex_aging))

if (length(missing_afca_columns) > 0) {
  stop(
    "AFCA aging table is missing required column(s): ",
    paste(missing_afca_columns, collapse = ", ")
  )
}

if (length(missing_gtex_columns) > 0) {
  stop(
    "GTEx aging table is missing required column(s): ",
    paste(missing_gtex_columns, collapse = ", ")
  )
}

fly_system_map <- data.frame(
  cell_type_clean = c(
    "Fat body body",
    "Fat body head",
    "Midgut enterocyte",
    "Differentiating enterocyte",
    "Hindgut",
    "Heart muscle",
    "Tracheal cell",
    "Malpighian tubule principal cell",
    "Malpighian tubule stellate cell",
    "Salivary gland",
    "Oenocyte"
  ),
  Functional_System = c(
    "Metabolic",
    "Metabolic",
    "Digestive",
    "Digestive",
    "Digestive",
    "Cardiovascular",
    "Respiratory",
    "Excretory",
    "Excretory",
    "Secretory",
    "Metabolic"
  ),
  stringsAsFactors = FALSE
)

human_system_map <- data.frame(
  SMTSD = c(
    "Adipose - Subcutaneous",
    "Adipose - Visceral (Omentum)",
    "Liver",
    "Small Intestine - Terminal Ileum",
    "Colon - Transverse",
    "Heart - Left Ventricle",
    "Heart - Atrial Appendage",
    "Lung",
    "Kidney - Cortex",
    "Kidney - Medulla",
    "Minor Salivary Gland"
  ),
  Functional_System = c(
    "Metabolic",
    "Metabolic",
    "Metabolic",
    "Digestive",
    "Digestive",
    "Cardiovascular",
    "Cardiovascular",
    "Respiratory",
    "Excretory",
    "Excretory",
    "Secretory"
  ),
  stringsAsFactors = FALSE
)

system_order <- c(
  "Metabolic",
  "Digestive",
  "Cardiovascular",
  "Respiratory",
  "Excretory",
  "Secretory"
)

age_group_order <- c(
  "Middle_age",
  "Older_age",
  "Oldest_age"
)

if (anyDuplicated(fly_system_map$cell_type_clean) > 0) {
  stop("fly_system_map contains duplicate cell_type_clean values.")
}

if (anyDuplicated(human_system_map$SMTSD) > 0) {
  stop("human_system_map contains duplicate SMTSD values.")
}

mapped_systems <- unique(c(fly_system_map$Functional_System, human_system_map$Functional_System))
missing_systems <- setdiff(mapped_systems, system_order)

if (length(missing_systems) > 0) {
  stop(
    "system_order is missing mapped functional system(s): ",
    paste(missing_systems, collapse = ", ")
  )
}

mean_or_na <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

median_or_na <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  median(x, na.rm = TRUE)
}

scale_01 <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  if (x_max == x_min) {
    return(rep(0.5, length(x)))
  }

  (x - x_min) / (x_max - x_min)
}

safe_cor <- function(x, y, method = "pearson") {
  keep <- complete.cases(x, y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])

  if (length(x) < 3) {
    return(NA_real_)
  }

  if (sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }

  suppressWarnings(cor(x, y, method = method))
}

cosine_similarity <- function(x, y) {
  keep <- complete.cases(x, y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])

  if (length(x) < 2) {
    return(NA_real_)
  }

  denominator <- sqrt(sum(x^2) * sum(y^2))

  if (denominator == 0) {
    return(NA_real_)
  }

  sum(x * y) / denominator
}

euclidean_distance <- function(x, y) {
  keep <- complete.cases(x, y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])

  if (length(x) < 2) {
    return(NA_real_)
  }

  sqrt(sum((x - y)^2))
}

# ------------------------------------------------------------
# Process AFCA aging changes
# ------------------------------------------------------------

afca_change_cols <- grep(
  "vs_5d|X[0-9]+d_vs_5d",
  names(afca_aging),
  value = TRUE
)

if (length(afca_change_cols) == 0) {
  stop("No AFCA age-change columns were detected using pattern: vs_5d|X[0-9]+d_vs_5d")
}

message("AFCA age-change columns: ", paste(afca_change_cols, collapse = ", "))

afca_candidates <- afca_aging %>%
  filter(Gene %in% selected_serpins)

missing_serpins <- setdiff(selected_serpins, unique(afca_candidates$Gene))

if (length(missing_serpins) > 0) {
  warning(
    "The following selected serpin(s) were not found in AFCA aging data: ",
    paste(missing_serpins, collapse = ", ")
  )
}

unmapped_fly_cell_types <- setdiff(unique(afca_candidates$cell_type_clean), fly_system_map$cell_type_clean)

if (length(unmapped_fly_cell_types) > 0) {
  warning(
    "Selected AFCA cell type(s) not represented in fly_system_map will be excluded: ",
    paste(unmapped_fly_cell_types, collapse = ", ")
  )
}

afca_system_long <- afca_candidates %>%
  inner_join(
    fly_system_map,
    by = "cell_type_clean"
  ) %>%
  pivot_longer(
    cols = all_of(afca_change_cols),
    names_to = "Age_Comparison",
    values_to = "Fly_Aging_Change"
  ) %>%
  mutate(
    Fly_Aging_Change = as.numeric(Fly_Aging_Change),
    Age_Comparison = str_replace(Age_Comparison, "^X", ""),
    Age_Comparison = str_replace(Age_Comparison, "_vs_5d", ""),
    Age_Comparison = str_replace(Age_Comparison, "d", ""),
    Fly_Age_Day = as.numeric(Age_Comparison),
    Age_Group = case_when(
      Fly_Age_Day == 30 ~ "Middle_age",
      Fly_Age_Day == 50 ~ "Older_age",
      Fly_Age_Day == 70 ~ "Oldest_age",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(afca_system_long$Fly_Aging_Change))) {
  warning(
    "One or more AFCA aging-change values are NA after numeric conversion. ",
    "Missing values will be ignored in summary statistics."
  )
}

unknown_afca_ages <- afca_system_long %>%
  filter(is.na(Age_Group)) %>%
  distinct(Age_Comparison) %>%
  pull(Age_Comparison)

if (length(unknown_afca_ages) > 0) {
  warning(
    "The following AFCA age comparison(s) were not mapped to an Age_Group and will be excluded: ",
    paste(unknown_afca_ages, collapse = ", ")
  )
}

afca_system_long <- afca_system_long %>%
  filter(!is.na(Age_Group))

afca_system_profile <- afca_system_long %>%
  group_by(
    Gene,
    Functional_System,
    Age_Group
  ) %>%
  summarize(
    Fly_Mean_Aging_Change = mean_or_na(Fly_Aging_Change),
    Fly_Median_Aging_Change = median_or_na(Fly_Aging_Change),
    Fly_Mean_Abs_Aging_Change = mean_or_na(abs(Fly_Aging_Change)),
    N_Fly_Celltypes = n(),
    .groups = "drop"
  ) %>%
  arrange(Gene, Functional_System, Age_Group)

write.csv(
  afca_system_profile,
  file.path(integrated_processed_dir, "AFCA_Matched_Functional_System_Aging_Profile.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Process GTEx SERPINE1 aging changes
# ------------------------------------------------------------

gtex_serpine1 <- gtex_aging %>%
  filter(Gene == "SERPINE1") %>%
  mutate(change_from_20s = as.numeric(change_from_20s))

if (nrow(gtex_serpine1) == 0) {
  warning("No GTEx aging rows were found for Gene == 'SERPINE1'.")
}

if (any(is.na(gtex_serpine1$change_from_20s))) {
  warning(
    "One or more GTEx SERPINE1 change_from_20s values are NA after numeric conversion. ",
    "Missing values will be ignored in summary statistics."
  )
}

unmapped_human_tissues <- setdiff(unique(gtex_serpine1$SMTSD), human_system_map$SMTSD)

if (length(unmapped_human_tissues) > 0) {
  warning(
    "SERPINE1 GTEx tissue(s) not represented in human_system_map will be excluded: ",
    paste(unmapped_human_tissues, collapse = ", ")
  )
}

gtex_system_long <- gtex_serpine1 %>%
  inner_join(
    human_system_map,
    by = "SMTSD"
  ) %>%
  mutate(
    Age_Group = case_when(
      AGE %in% c("40-49", "50-59") ~ "Middle_age",
      AGE == "60-69" ~ "Older_age",
      AGE == "70-79" ~ "Oldest_age",
      TRUE ~ NA_character_
    )
  )

unknown_gtex_ages <- gtex_system_long %>%
  filter(is.na(Age_Group)) %>%
  distinct(AGE) %>%
  pull(AGE)

if (length(unknown_gtex_ages) > 0) {
  warning(
    "The following GTEx AGE value(s) were not mapped to an Age_Group and will be excluded: ",
    paste(unknown_gtex_ages, collapse = ", ")
  )
}

gtex_system_profile <- gtex_system_long %>%
  filter(!is.na(Age_Group)) %>%
  group_by(
    Functional_System,
    Age_Group
  ) %>%
  summarize(
    SERPINE1_Mean_Aging_Change = mean_or_na(change_from_20s),
    SERPINE1_Median_Aging_Change = median_or_na(change_from_20s),
    SERPINE1_Mean_Abs_Aging_Change = mean_or_na(abs(change_from_20s)),
    N_GTEx_Tissues = n(),
    .groups = "drop"
  ) %>%
  arrange(Functional_System, Age_Group)

write.csv(
  gtex_system_profile,
  file.path(integrated_processed_dir, "GTEx_SERPINE1_Matched_Functional_System_Aging_Profile.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Merge fly and human aging profiles
# ------------------------------------------------------------

duplicate_afca_keys <- afca_system_profile %>%
  count(Gene, Functional_System, Age_Group) %>%
  filter(n > 1)

if (nrow(duplicate_afca_keys) > 0) {
  stop("afca_system_profile contains duplicate Gene/Functional_System/Age_Group rows.")
}

duplicate_gtex_keys <- gtex_system_profile %>%
  count(Functional_System, Age_Group) %>%
  filter(n > 1)

if (nrow(duplicate_gtex_keys) > 0) {
  stop("gtex_system_profile contains duplicate Functional_System/Age_Group rows.")
}

joined_aging_comparison <- afca_system_profile %>%
  left_join(
    gtex_system_profile,
    by = c(
      "Functional_System",
      "Age_Group"
    )
  )

unmatched_fly_rows <- joined_aging_comparison %>%
  filter(is.na(SERPINE1_Mean_Aging_Change)) %>%
  distinct(Functional_System, Age_Group)

if (nrow(unmatched_fly_rows) > 0) {
  warning(
    "No matching GTEx SERPINE1 aging profile was joined for ",
    nrow(unmatched_fly_rows),
    " Functional_System/Age_Group combination(s)."
  )
}

matched_aging_comparison <- joined_aging_comparison %>%
  filter(
    !is.na(Fly_Mean_Aging_Change),
    !is.na(SERPINE1_Mean_Aging_Change)
  ) %>%
  mutate(
    Functional_System = factor(
      Functional_System,
      levels = system_order
    ),
    Age_Group = factor(
      Age_Group,
      levels = age_group_order
    )
  ) %>%
  arrange(Gene, Functional_System, Age_Group)

if (nrow(matched_aging_comparison) == 0) {
  warning("No matched fly-human aging comparison rows remained after filtering missing mean aging changes.")
}

write.csv(
  matched_aging_comparison,
  file.path(integrated_processed_dir, "Matched_Tissue_CrossSpecies_Aging_Profile_Long.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Calculate aging similarity
# ------------------------------------------------------------

matched_aging_input <- matched_aging_comparison %>%
  mutate(
    Fly_Mean_Aging_Change = as.numeric(Fly_Mean_Aging_Change),
    SERPINE1_Mean_Aging_Change = as.numeric(SERPINE1_Mean_Aging_Change),
    Fly_Mean_Abs_Aging_Change = as.numeric(Fly_Mean_Abs_Aging_Change),
    SERPINE1_Mean_Abs_Aging_Change = as.numeric(SERPINE1_Mean_Abs_Aging_Change)
  ) %>%
  group_by(Gene, Functional_System, Age_Group) %>%
  summarize(
    Fly_Mean_Aging_Change = mean_or_na(Fly_Mean_Aging_Change),
    SERPINE1_Mean_Aging_Change = mean_or_na(SERPINE1_Mean_Aging_Change),
    Fly_Mean_Abs_Aging_Change = mean_or_na(Fly_Mean_Abs_Aging_Change),
    SERPINE1_Mean_Abs_Aging_Change = mean_or_na(SERPINE1_Mean_Abs_Aging_Change),
    .groups = "drop"
  )

matched_aging_similarity <- matched_aging_input %>%
  group_by(Gene) %>%
  summarize(
    N_Comparisons = sum(complete.cases(Fly_Mean_Aging_Change, SERPINE1_Mean_Aging_Change)),
    Aging_Tissue_Pearson = safe_cor(
      Fly_Mean_Aging_Change,
      SERPINE1_Mean_Aging_Change,
      "pearson"
    ),
    Aging_Tissue_Spearman = safe_cor(
      Fly_Mean_Aging_Change,
      SERPINE1_Mean_Aging_Change,
      "spearman"
    ),
    Aging_Tissue_Cosine = cosine_similarity(
      Fly_Mean_Aging_Change,
      SERPINE1_Mean_Aging_Change
    ),
    Aging_Tissue_Euclidean = euclidean_distance(
      Fly_Mean_Aging_Change,
      SERPINE1_Mean_Aging_Change
    ),
    Aging_Tissue_MeanAbs_Difference = mean_or_na(
      abs(
        Fly_Mean_Abs_Aging_Change -
          SERPINE1_Mean_Abs_Aging_Change
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Aging_Tissue_Pearson_scaled = scale_01(Aging_Tissue_Pearson),
    Aging_Tissue_Spearman_scaled = scale_01(Aging_Tissue_Spearman),
    Aging_Tissue_Cosine_scaled = scale_01(Aging_Tissue_Cosine),
    Aging_Tissue_Euclidean_scaled = 1 - scale_01(Aging_Tissue_Euclidean),
    Aging_Tissue_MeanAbs_Difference_scaled =
      1 - scale_01(Aging_Tissue_MeanAbs_Difference)
  ) %>%
  rowwise() %>%
  mutate(
    Matched_Tissue_Aging_Similarity = mean_or_na(
      c(
        Aging_Tissue_Pearson_scaled,
        Aging_Tissue_Spearman_scaled,
        Aging_Tissue_Cosine_scaled,
        Aging_Tissue_Euclidean_scaled,
        Aging_Tissue_MeanAbs_Difference_scaled
      )
    )
  ) %>%
  ungroup() %>%
  arrange(desc(Matched_Tissue_Aging_Similarity), Gene)

low_comparison_genes <- matched_aging_similarity %>%
  filter(N_Comparisons < 3) %>%
  pull(Gene)

if (length(low_comparison_genes) > 0) {
  warning(
    "The following gene(s) have fewer than 3 complete matched comparisons, so correlation metrics may be NA: ",
    paste(low_comparison_genes, collapse = ", ")
  )
}

write.csv(
  matched_aging_similarity,
  file.path(integrated_processed_dir, "Matched_Tissue_CrossSpecies_Aging_Similarity_to_SERPINE1.csv"),
  row.names = FALSE
)

write.csv(
  matched_aging_similarity,
  file.path(tables_dir, "Matched_Tissue_CrossSpecies_Aging_Similarity_to_SERPINE1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Aging similarity ranking plot
# ------------------------------------------------------------

plot_data <- matched_aging_similarity %>%
  mutate(
    Matched_Tissue_Aging_Similarity = as.numeric(Matched_Tissue_Aging_Similarity)
  ) %>%
  filter(!is.na(Matched_Tissue_Aging_Similarity))

if (nrow(plot_data) == 0) {
  warning("Aging similarity ranking plot was not generated because no non-missing scores are available.")
} else {
  x_min <- min(0, min(plot_data$Matched_Tissue_Aging_Similarity, na.rm = TRUE))
  x_max <- max(plot_data$Matched_Tissue_Aging_Similarity, na.rm = TRUE)

  if (x_max == x_min) {
    x_limits <- c(x_min - 0.05, x_max + 0.05)
  } else {
    x_padding <- (x_max - x_min) * 0.15
    x_limits <- c(x_min, x_max + x_padding)
  }

  p_aging_rank <- ggplot(
    plot_data,
    aes(
      x = Matched_Tissue_Aging_Similarity,
      y = reorder(Gene, Matched_Tissue_Aging_Similarity)
    )
  ) +
    geom_col(width = 0.7, fill = "#B2182B") +
    geom_text(
      aes(label = round(Matched_Tissue_Aging_Similarity, 3)),
      hjust = -0.15,
      size = 3.5
    ) +
    coord_cartesian(
      xlim = x_limits
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Matched tissue aging similarity to human SERPINE1",
      subtitle = "Aging changes compared across matched biological systems",
      x = "Matched tissue aging similarity score",
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  print(p_aging_rank)

  ggsave(
    file.path(fig_dir, "Matched_Tissue_Aging_Similarity_to_SERPINE1.png"),
    p_aging_rank,
    width = 7,
    height = 4.5,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(fig_dir, "Matched_Tissue_Aging_Similarity_to_SERPINE1.pdf"),
    p_aging_rank,
    width = 7,
    height = 4.5
  )
}
