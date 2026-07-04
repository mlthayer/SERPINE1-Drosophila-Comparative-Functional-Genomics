# ============================================================
# 05_FlyAtlas2_TissueExpression.R
# FlyAtlas2 tissue expression analysis for candidate serpins
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Prefer an environment variable for reproducibility across machines.
# If it is not set, fall back to the original local project path.
project_dir <- Sys.getenv(
  "xxxxx_xxxxxxxx_PROJECT_DIR",
  unset = "D:/xxxxx-xxxxxxxx xxx/xxxxx-xxxxxxxx xxx"
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

if (!file.exists(flyatlas_file)) {
  stop("FlyAtlas2 input file not found: ", flyatlas_file)
}

flyatlas <- read.csv(
  flyatlas_file,
  stringsAsFactors = FALSE
)

if (nrow(flyatlas) == 0) {
  warning("FlyAtlas2 input table has zero rows. Output tables will be empty.")
}

required_columns <- c("Gene", "Tissue", "MeanFPKM", "MedianFPKM")
missing_columns <- setdiff(required_columns, names(flyatlas))

if (length(missing_columns) > 0) {
  stop(
    "FlyAtlas2 input table is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

selected_serpins <- c(
  "Spn42Da", "Spn42Db", "Spn42Dc",
  "Spn42Dd", "Spn28B", "Spn28F"
)

if (length(selected_serpins) == 0) {
  stop("selected_serpins is empty. Add at least one candidate gene before running FlyAtlas2 analysis.")
}

relevant_tissues <- c(
  "Fat body",
  "Midgut",
  "Hindgut",
  "Heart",
  "Carcass",
  "Trachea",
  "Malpighian Tubules",
  "Salivary gland",
  "Brain",
  "Ovary",
  "Testis",
  "Whole body"
)

expression_threshold_fpkm <- 1

safe_mean <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  median(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  max(x, na.rm = TRUE)
}

flyatlas_clean <- flyatlas %>%
  mutate(
    MeanFPKM = as.numeric(MeanFPKM),
    MedianFPKM = as.numeric(MedianFPKM)
  )

if (any(is.na(flyatlas_clean$MeanFPKM)) || any(is.na(flyatlas_clean$MedianFPKM))) {
  warning(
    "One or more FlyAtlas2 FPKM values are NA after numeric conversion. ",
    "Missing values will be ignored in summary statistics and shown as blank heatmap cells."
  )
}

if (any(flyatlas_clean$MedianFPKM < 0 | flyatlas_clean$MeanFPKM < 0, na.rm = TRUE)) {
  stop("FlyAtlas2 FPKM values contain negatives, so log2(MedianFPKM + 1) cannot be safely computed.")
}

missing_serpins <- setdiff(selected_serpins, unique(flyatlas_clean$Gene))

if (length(missing_serpins) > 0) {
  warning(
    "The following selected serpin(s) were not found in FlyAtlas2: ",
    paste(missing_serpins, collapse = ", ")
  )
}

missing_relevant_tissues <- setdiff(relevant_tissues, unique(flyatlas_clean$Tissue))

if (length(missing_relevant_tissues) > 0) {
  warning(
    "The following relevant tissue(s) were not found in FlyAtlas2: ",
    paste(missing_relevant_tissues, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Summary table
# ------------------------------------------------------------

flyatlas_summary <- flyatlas_clean %>%
  filter(Gene %in% selected_serpins) %>%
  group_by(Gene) %>%
  summarize(
    Tissues_Expressed = sum(MedianFPKM > expression_threshold_fpkm, na.rm = TRUE),
    Total_Tissues = n(),
    Percent_Tissues = if_else(
      Total_Tissues > 0,
      100 * Tissues_Expressed / Total_Tissues,
      NA_real_
    ),
    Mean_FPKM = safe_mean(MeanFPKM),
    Median_FPKM = safe_median(MedianFPKM),
    Max_FPKM = safe_max(MedianFPKM),
    .groups = "drop"
  ) %>%
  arrange(desc(Percent_Tissues), desc(Median_FPKM), Gene)

if (nrow(flyatlas_summary) == 0) {
  warning("No FlyAtlas2 summary rows were generated for selected serpins.")
}

print(flyatlas_summary)

write.csv(
  flyatlas_summary,
  file.path(processed_dir, "FlyAtlas2_tissue_expression_summary_for_candidates.csv"),
  row.names = FALSE
)

write.csv(
  flyatlas_summary,
  file.path(tables_dir, "FlyAtlas2_tissue_expression_summary_for_candidates.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Heatmap data
# ------------------------------------------------------------

flyatlas_heatmap_data <- flyatlas_clean %>%
  filter(
    Gene %in% selected_serpins,
    Tissue %in% relevant_tissues
  ) %>%
  mutate(
    log2_MedianFPKM = log2(MedianFPKM + 1),
    Gene = factor(Gene, levels = selected_serpins),
    Tissue = factor(Tissue, levels = rev(relevant_tissues))
  )

if (nrow(flyatlas_heatmap_data) == 0) {
  warning("No FlyAtlas2 heatmap rows matched selected serpins and relevant tissues.")
}

write.csv(
  flyatlas_heatmap_data,
  file.path(processed_dir, "FlyAtlas2_candidate_serpin_tissue_heatmap_data.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------

if (nrow(flyatlas_heatmap_data) == 0) {
  warning("FlyAtlas2 heatmap was not generated because there are no matching rows.")
} else {
  p_flyatlas <- ggplot(
    flyatlas_heatmap_data,
    aes(
      x = Gene,
      y = Tissue,
      fill = log2_MedianFPKM
    )
  ) +
    geom_tile(color = "grey85", linewidth = 0.5) +
    geom_text(
      aes(label = round(log2_MedianFPKM, 1)),
      size = 3,
      na.rm = TRUE
    ) +
    scale_fill_gradientn(
      colors = c(
        "#FFFFCC",
        "#C2E699",
        "#78C679",
        "#31A354",
        "#006837"
      ),
      name = expression(log[2]~"(median FPKM + 1)"),
      na.value = "grey95"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "FlyAtlas2 tissue expression of candidate Drosophila serpins",
      x = NULL,
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  print(p_flyatlas)

  ggsave(
    file.path(fig_dir, "FlyAtlas2_candidate_serpin_tissue_heatmap.png"),
    p_flyatlas,
    width = 8,
    height = 6,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(fig_dir, "FlyAtlas2_candidate_serpin_tissue_heatmap.pdf"),
    p_flyatlas,
    width = 8,
    height = 6
  )
}
