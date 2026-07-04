# ============================================================
# 06_GTEx_HPA_SERPINE1_TissueExpression.R
# Human SERPINE1 tissue expression using GTEx and HPA
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
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

gtex_file <- file.path(processed_dir, "human_SERPINE1_GTEx_tissue_expression_summary.csv")
hpa_file <- file.path(processed_dir, "HPA_SERPINE1_tissue_summary.csv")

required_files <- c(
  GTEx = gtex_file,
  HPA = hpa_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required tissue-expression file(s): ",
    paste(names(missing_files), basename(missing_files), sep = " = ", collapse = "; ")
  )
}

gtex <- read.csv(gtex_file, stringsAsFactors = FALSE)
hpa <- read.csv(hpa_file, stringsAsFactors = FALSE)

if (nrow(gtex) == 0) {
  warning("GTEx input table has zero rows.")
}

if (nrow(hpa) == 0) {
  warning("HPA input table has zero rows.")
}

required_gtex_columns <- c("gene", "tissue", "median_TPM", "mean_TPM")
missing_gtex_columns <- setdiff(required_gtex_columns, names(gtex))

if (length(missing_gtex_columns) > 0) {
  stop(
    "GTEx input table is missing required column(s): ",
    paste(missing_gtex_columns, collapse = ", ")
  )
}

if (ncol(hpa) == 0) {
  stop("HPA input table has zero columns and cannot be cleaned.")
}

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

# ------------------------------------------------------------
# GTEx SERPINE1 summary
# ------------------------------------------------------------

gtex_clean <- gtex %>%
  mutate(
    median_TPM = as.numeric(median_TPM),
    mean_TPM = as.numeric(mean_TPM)
  )

if (any(is.na(gtex_clean$median_TPM)) || any(is.na(gtex_clean$mean_TPM))) {
  warning(
    "One or more GTEx TPM values are NA after numeric conversion. ",
    "Missing values will be retained in output tables and ignored by plots where needed."
  )
}

if (any(gtex_clean$median_TPM < 0 | gtex_clean$mean_TPM < 0, na.rm = TRUE)) {
  stop("GTEx TPM values contain negatives, so log2(TPM + 1) cannot be safely computed.")
}

gtex_serpine1 <- gtex_clean %>%
  filter(gene == "SERPINE1") %>%
  mutate(
    log2_median_TPM = log2(median_TPM + 1),
    log2_mean_TPM = log2(mean_TPM + 1)
  ) %>%
  arrange(desc(median_TPM), tissue)

if (nrow(gtex_serpine1) == 0) {
  warning("No GTEx rows were found for gene == 'SERPINE1'.")
}

print(gtex_serpine1)

write.csv(
  gtex_serpine1,
  file.path(processed_dir, "GTEx_SERPINE1_tissue_expression_summary.csv"),
  row.names = FALSE
)

write.csv(
  gtex_serpine1,
  file.path(tables_dir, "GTEx_SERPINE1_tissue_expression_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# GTEx SERPINE1 plot
# ------------------------------------------------------------

gtex_plot_data <- gtex_serpine1 %>%
  filter(!is.na(log2_median_TPM))

if (nrow(gtex_plot_data) == 0) {
  warning("GTEx SERPINE1 plot was not generated because no non-missing log2 median TPM values are available.")
} else {
  p_gtex <- ggplot(
    gtex_plot_data,
    aes(
      x = log2_median_TPM,
      y = reorder(tissue, log2_median_TPM)
    )
  ) +
    geom_col(width = 0.7, fill = "#2166AC") +
    theme_bw(base_size = 11) +
    labs(
      title = "GTEx SERPINE1 tissue expression",
      x = expression(log[2]~"(median TPM + 1)"),
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  print(p_gtex)

  ggsave(
    file.path(fig_dir, "GTEx_SERPINE1_tissue_expression.png"),
    p_gtex,
    width = 8,
    height = 7,
    dpi = 600,
    bg = "white"
  )
}

# ------------------------------------------------------------
# HPA clean table
# ------------------------------------------------------------

hpa_clean <- hpa %>%
  rename_with(~ str_replace_all(.x, "\\.", "_")) %>%
  rename_with(~ str_replace_all(.x, " ", "_")) %>%
  rename_with(~ str_replace_all(.x, "_+", "_")) %>%
  rename_with(~ str_replace_all(.x, "^_|_$", ""))

duplicated_columns <- names(hpa_clean)[duplicated(names(hpa_clean))]

if (length(duplicated_columns) > 0) {
  stop(
    "HPA column-name cleaning created duplicate column name(s): ",
    paste(unique(duplicated_columns), collapse = ", ")
  )
}

write.csv(
  hpa_clean,
  file.path(processed_dir, "HPA_SERPINE1_tissue_expression_validation_clean.csv"),
  row.names = FALSE
)

write.csv(
  hpa_clean,
  file.path(tables_dir, "HPA_SERPINE1_tissue_expression_validation_clean.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Optional HPA numeric summary
# ------------------------------------------------------------

hpa_numeric_cols <- names(hpa_clean)[vapply(hpa_clean, is.numeric, logical(1))]

if (length(hpa_numeric_cols) > 0) {
  hpa_value_col <- hpa_numeric_cols[1]

  hpa_summary <- hpa_clean %>%
    summarize(
      HPA_Value_Column = hpa_value_col,
      HPA_N_Tissues = n(),
      HPA_N_Nonmissing = sum(!is.na(.data[[hpa_value_col]])),
      HPA_Mean_Value = safe_mean(.data[[hpa_value_col]]),
      HPA_Median_Value = safe_median(.data[[hpa_value_col]]),
      HPA_Max_Value = safe_max(.data[[hpa_value_col]])
    )

  print(hpa_summary)

  write.csv(
    hpa_summary,
    file.path(processed_dir, "HPA_SERPINE1_expression_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    hpa_summary,
    file.path(tables_dir, "HPA_SERPINE1_expression_summary.csv"),
    row.names = FALSE
  )
} else {
  warning("No numeric HPA columns were detected, so HPA numeric summary was not generated.")
}
