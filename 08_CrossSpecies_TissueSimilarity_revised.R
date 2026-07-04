# ============================================================
# 08_CrossSpecies_TissueSimilarity.R
# Cross-species tissue similarity: FlyAtlas2 vs GTEx SERPINE1
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# Prefer an environment variable for reproducibility across machines.
# If it is not set, fall back to the original local project path.
project_dir <- Sys.getenv(
  "xxxxxx_xxxxxxxx_PROJECT_DIR",
  unset = "D:/xxxxxx-xxxxxxxx xxx/xxxxxx-xxxxxxxx xxx"
)

expression_dir <- file.path(
  project_dir,
  "Bioinformatics_tissue expression collection",
  "Bioinformatics_tissue_expression"
)

processed_dir <- file.path(expression_dir, "Data", "Processed")
fig_dir <- file.path(expression_dir, "Figures", "TissueExpression")
tables_dir <- file.path(expression_dir, "Tables", "TissueExpression")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

flyatlas_file <- file.path(processed_dir, "FlyAtlas2_serpin_tissue_summary.csv")
gtex_file <- file.path(processed_dir, "human_SERPINE1_GTEx_tissue_expression_summary.csv")

required_files <- c(
  FlyAtlas2 = flyatlas_file,
  GTEx = gtex_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required tissue-expression file(s): ",
    paste(names(missing_files), basename(missing_files), sep = " = ", collapse = "; ")
  )
}

flyatlas <- read.csv(flyatlas_file, stringsAsFactors = FALSE)
gtex <- read.csv(gtex_file, stringsAsFactors = FALSE)

required_flyatlas_columns <- c("Gene", "Tissue", "MedianFPKM")
required_gtex_columns <- c("gene", "tissue", "median_TPM")

missing_flyatlas_columns <- setdiff(required_flyatlas_columns, names(flyatlas))
missing_gtex_columns <- setdiff(required_gtex_columns, names(gtex))

if (length(missing_flyatlas_columns) > 0) {
  stop(
    "FlyAtlas2 input table is missing required column(s): ",
    paste(missing_flyatlas_columns, collapse = ", ")
  )
}

if (length(missing_gtex_columns) > 0) {
  stop(
    "GTEx input table is missing required column(s): ",
    paste(missing_gtex_columns, collapse = ", ")
  )
}

selected_serpins <- c(
  "Spn42Da", "Spn42Db", "Spn42Dc",
  "Spn42Dd", "Spn28B", "Spn28F"
)

if (length(selected_serpins) == 0) {
  stop("selected_serpins is empty. Add at least one candidate gene before calculating tissue similarity.")
}

fly_system_map <- data.frame(
  Tissue = c(
    "Fat body",
    "Midgut",
    "Hindgut",
    "Heart",
    "Trachea",
    "Brain",
    "CNS",
    "Malpighian Tubules",
    "Ovary",
    "Testis",
    "Salivary gland",
    "Whole body"
  ),
  Functional_System = c(
    "Metabolic",
    "Digestive",
    "Digestive",
    "Cardiovascular",
    "Respiratory",
    "Nervous",
    "Nervous",
    "Excretory",
    "Reproductive",
    "Reproductive",
    "Secretory",
    "Immune/systemic"
  ),
  stringsAsFactors = FALSE
)

human_system_map <- data.frame(
  tissue = c(
    "adipose subcutaneous",
    "adipose visceral omentum",
    "liver",
    "small intestine terminal ileum",
    "heart left ventricle",
    "heart atrial appendage",
    "lung",
    "brain hypothalamus",
    "brain nucleus accumbens basal ganglia",
    "brain spinal cord cervical c-1",
    "whole blood",
    "ovary",
    "testis",
    "pancreas"
  ),
  Functional_System = c(
    "Metabolic",
    "Metabolic",
    "Metabolic",
    "Digestive",
    "Cardiovascular",
    "Cardiovascular",
    "Respiratory",
    "Nervous",
    "Nervous",
    "Nervous",
    "Immune/systemic",
    "Reproductive",
    "Reproductive",
    "Secretory"
  ),
  stringsAsFactors = FALSE
)

system_order <- c(
  "Metabolic",
  "Digestive",
  "Cardiovascular",
  "Respiratory",
  "Nervous",
  "Excretory",
  "Reproductive",
  "Secretory",
  "Immune/systemic"
)

gene_order <- c(
  "Human_SERPINE1",
  "Spn42Da",
  "Spn42Dd",
  "Spn42Db",
  "Spn42Dc",
  "Spn28F",
  "Spn28B"
)

if (anyDuplicated(fly_system_map$Tissue) > 0) {
  stop("fly_system_map contains duplicate Tissue values.")
}

if (anyDuplicated(human_system_map$tissue) > 0) {
  stop("human_system_map contains duplicate tissue values.")
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

safe_z <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  x_sd <- sd(x, na.rm = TRUE)

  if (is.na(x_sd) || x_sd == 0) {
    return(rep(0, length(x)))
  }

  as.numeric(scale(x))
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

safe_euclidean <- function(x, y) {
  keep <- complete.cases(x, y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])

  if (length(x) < 2) {
    return(NA_real_)
  }

  sqrt(sum((x - y)^2))
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

flyatlas_clean <- flyatlas %>%
  mutate(MedianFPKM = as.numeric(MedianFPKM))

gtex_clean <- gtex %>%
  mutate(median_TPM = as.numeric(median_TPM))

if (any(is.na(flyatlas_clean$MedianFPKM))) {
  warning(
    "One or more FlyAtlas2 MedianFPKM values are NA after numeric conversion. ",
    "Missing values will be ignored when calculating medians."
  )
}

if (any(is.na(gtex_clean$median_TPM))) {
  warning(
    "One or more GTEx median_TPM values are NA after numeric conversion. ",
    "Missing values will be ignored when calculating medians."
  )
}

if (any(flyatlas_clean$MedianFPKM < 0, na.rm = TRUE)) {
  stop("FlyAtlas2 MedianFPKM contains negative values, so log2(Expression + 1) cannot be safely computed.")
}

if (any(gtex_clean$median_TPM < 0, na.rm = TRUE)) {
  stop("GTEx median_TPM contains negative values, so log2(Expression + 1) cannot be safely computed.")
}

gtex_serpine1 <- gtex_clean %>%
  filter(gene == "SERPINE1")

if (nrow(gtex_serpine1) == 0) {
  warning("No GTEx rows were found for gene == 'SERPINE1'.")
}

unmapped_human_tissues <- setdiff(unique(gtex_serpine1$tissue), human_system_map$tissue)

if (length(unmapped_human_tissues) > 0) {
  warning(
    "SERPINE1 GTEx tissue(s) not represented in human_system_map will be excluded: ",
    paste(unmapped_human_tissues, collapse = ", ")
  )
}

human_system_profile <- gtex_serpine1 %>%
  inner_join(human_system_map, by = "tissue") %>%
  group_by(Functional_System) %>%
  summarize(
    Expression = median_or_na(median_TPM),
    .groups = "drop"
  ) %>%
  mutate(
    Gene = "Human_SERPINE1",
    Species = "Human",
    Unit = "GTEx median TPM"
  )

flyatlas_candidates <- flyatlas_clean %>%
  filter(Gene %in% selected_serpins)

missing_serpins <- setdiff(selected_serpins, unique(flyatlas_candidates$Gene))

if (length(missing_serpins) > 0) {
  warning(
    "The following selected serpin(s) were not found in FlyAtlas2: ",
    paste(missing_serpins, collapse = ", ")
  )
}

unmapped_fly_tissues <- setdiff(unique(flyatlas_candidates$Tissue), fly_system_map$Tissue)

if (length(unmapped_fly_tissues) > 0) {
  warning(
    "Selected FlyAtlas2 tissue(s) not represented in fly_system_map will be excluded: ",
    paste(unmapped_fly_tissues, collapse = ", ")
  )
}

fly_system_profile <- flyatlas_candidates %>%
  inner_join(fly_system_map, by = "Tissue") %>%
  group_by(Gene, Functional_System) %>%
  summarize(
    Expression = median_or_na(MedianFPKM),
    .groups = "drop"
  ) %>%
  mutate(
    Species = "Drosophila",
    Unit = "FlyAtlas2 median FPKM"
  )

system_profile <- bind_rows(
  human_system_profile,
  fly_system_profile
) %>%
  select(Gene, Species, Unit, Functional_System, Expression) %>%
  complete(
    Gene,
    Functional_System = system_order,
    fill = list(Expression = 0)
  ) %>%
  mutate(
    Species = if_else(Gene == "Human_SERPINE1", "Human", "Drosophila"),
    Unit = if_else(Gene == "Human_SERPINE1", "GTEx median TPM", "FlyAtlas2 median FPKM"),
    log2_expression = log2(Expression + 1),
    Gene = factor(Gene, levels = gene_order),
    Functional_System = factor(Functional_System, levels = system_order)
  ) %>%
  arrange(Gene, Functional_System)

if (nrow(system_profile) == 0) {
  warning("No cross-species functional-system profile rows were generated.")
}

write.csv(
  system_profile,
  file.path(processed_dir, "CrossSpecies_TissueExpression_FunctionalSystem_Profile.csv"),
  row.names = FALSE
)

write.csv(
  system_profile,
  file.path(tables_dir, "CrossSpecies_TissueExpression_FunctionalSystem_Profile.csv"),
  row.names = FALSE
)

heatmap_data <- system_profile %>%
  filter(!is.na(log2_expression))

if (nrow(heatmap_data) == 0) {
  warning("Cross-species heatmap was not generated because no non-missing log2 expression values are available.")
} else {
  p_cross_species <- ggplot(
    heatmap_data,
    aes(
      x = Functional_System,
      y = Gene,
      fill = log2_expression
    )
  ) +
    geom_tile(color = "grey85", linewidth = 0.5) +
    geom_text(aes(label = round(log2_expression, 1)), size = 3, na.rm = TRUE) +
    facet_grid(
      Species ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_fill_gradientn(
      colors = c(
        "#F7FBFF",
        "#D0E1F2",
        "#92C5DE",
        "#4393C3",
        "#2166AC",
        "#053061"
      ),
      name = expression(log[2]~"(expression + 1)"),
      na.value = "grey95"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Cross-species functional-system tissue expression",
      subtitle = "Human: GTEx median TPM; Drosophila: FlyAtlas2 median FPKM",
      x = NULL,
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text.y = element_text(face = "bold")
    )

  print(p_cross_species)

  ggsave(
    file.path(fig_dir, "CrossSpecies_TissueExpression_FunctionalSystem_Heatmap.png"),
    p_cross_species,
    width = 10,
    height = 6,
    dpi = 600,
    bg = "white"
  )
}

# ------------------------------------------------------------
# Similarity metrics
# ------------------------------------------------------------

system_profile_for_similarity <- system_profile %>%
  mutate(
    Gene = as.character(Gene),
    Functional_System = as.character(Functional_System),
    log2_expression = as.numeric(log2_expression)
  )

if (!"Human_SERPINE1" %in% system_profile_for_similarity$Gene) {
  stop("system_profile does not contain Human_SERPINE1, which is required as the similarity reference.")
}

duplicate_profile_rows <- system_profile_for_similarity %>%
  count(Gene, Functional_System) %>%
  filter(n > 1)

if (nrow(duplicate_profile_rows) > 0) {
  warning(
    "Duplicate Gene/Functional_System rows were found in system_profile. ",
    "Mean log2_expression will be used before calculating similarity."
  )
}

system_profile_z <- system_profile_for_similarity %>%
  group_by(Gene, Functional_System) %>%
  summarize(
    log2_expression = mean_or_na(log2_expression),
    .groups = "drop"
  ) %>%
  group_by(Gene) %>%
  mutate(Expression_z = safe_z(log2_expression)) %>%
  ungroup()

wide_z <- system_profile_z %>%
  select(Gene, Functional_System, Expression_z) %>%
  complete(
    Gene,
    Functional_System = system_order
  ) %>%
  pivot_wider(
    names_from = Functional_System,
    values_from = Expression_z
  )

missing_wide_columns <- setdiff(system_order, names(wide_z))

if (length(missing_wide_columns) > 0) {
  stop(
    "The wide tissue-similarity table is missing functional-system column(s): ",
    paste(missing_wide_columns, collapse = ", ")
  )
}

serpine1_rows <- wide_z %>%
  filter(Gene == "Human_SERPINE1")

if (nrow(serpine1_rows) != 1) {
  stop("Expected exactly one Human_SERPINE1 row in the wide tissue-similarity table.")
}

serpine1_vec <- as.numeric(serpine1_rows[1, system_order])

if (all(is.na(serpine1_vec))) {
  stop("Human_SERPINE1 has no usable functional-system expression values for similarity calculation.")
}

tissue_similarity <- wide_z %>%
  filter(Gene != "Human_SERPINE1") %>%
  rowwise() %>%
  mutate(
    Pearson = safe_cor(c_across(all_of(system_order)), serpine1_vec, "pearson"),
    Spearman = safe_cor(c_across(all_of(system_order)), serpine1_vec, "spearman"),
    Cosine = cosine_similarity(c_across(all_of(system_order)), serpine1_vec),
    Euclidean = safe_euclidean(c_across(all_of(system_order)), serpine1_vec)
  ) %>%
  ungroup() %>%
  mutate(
    Pearson_scaled = scale_01(Pearson),
    Spearman_scaled = scale_01(Spearman),
    Cosine_scaled = scale_01(Cosine),
    Euclidean_scaled = 1 - scale_01(Euclidean)
  ) %>%
  rowwise() %>%
  mutate(
    Integrated_Tissue_Expression_Similarity = mean_or_na(
      c(
        Pearson_scaled,
        Spearman_scaled,
        Cosine_scaled,
        Euclidean_scaled
      )
    )
  ) %>%
  ungroup() %>%
  arrange(desc(Integrated_Tissue_Expression_Similarity), Gene)

if (nrow(tissue_similarity) == 0) {
  warning("No tissue-similarity rows were generated.")
}

print(tissue_similarity)

write.csv(
  tissue_similarity,
  file.path(processed_dir, "TissueExpression_Similarity_To_SERPINE1.csv"),
  row.names = FALSE
)

write.csv(
  tissue_similarity,
  file.path(tables_dir, "TissueExpression_Similarity_To_SERPINE1.csv"),
  row.names = FALSE
)

similarity_plot_data <- tissue_similarity %>%
  mutate(
    Integrated_Tissue_Expression_Similarity =
      as.numeric(Integrated_Tissue_Expression_Similarity)
  ) %>%
  filter(!is.na(Integrated_Tissue_Expression_Similarity))

if (nrow(similarity_plot_data) == 0) {
  warning("Tissue-expression similarity ranking plot was not generated because no non-missing scores are available.")
} else {
  x_min <- min(0, min(similarity_plot_data$Integrated_Tissue_Expression_Similarity, na.rm = TRUE))
  x_max <- max(similarity_plot_data$Integrated_Tissue_Expression_Similarity, na.rm = TRUE)

  if (x_max == x_min) {
    x_limits <- c(x_min - 0.05, x_max + 0.05)
  } else {
    x_padding <- (x_max - x_min) * 0.15
    x_limits <- c(x_min, x_max + x_padding)
  }

  p_similarity <- ggplot(
    similarity_plot_data,
    aes(
      x = Integrated_Tissue_Expression_Similarity,
      y = reorder(Gene, Integrated_Tissue_Expression_Similarity)
    )
  ) +
    geom_col(width = 0.7, fill = "#2166AC") +
    geom_text(
      aes(label = round(Integrated_Tissue_Expression_Similarity, 3)),
      hjust = -0.15,
      size = 3.5
    ) +
    coord_cartesian(
      xlim = x_limits
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Functional-system tissue-expression similarity to human SERPINE1",
      x = "Integrated tissue-expression similarity score",
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  print(p_similarity)

  ggsave(
    file.path(fig_dir, "TissueExpression_Similarity_To_SERPINE1_Ranking.png"),
    p_similarity,
    width = 7,
    height = 4.5,
    dpi = 600,
    bg = "white"
  )
}
