# ============================================================
# 07_FlyCellAtlas_CellExpression.R
# Fly Cell Atlas functional cell-class expression
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

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
fig_dir <- file.path(expression_dir, "Figures", "CellExpression")
tables_dir <- file.path(expression_dir, "Tables", "CellExpression")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

flycell_file <- file.path(processed_dir, "FlyCellAtlas_serpin_celltype_expression_summary.csv")

if (!file.exists(flycell_file)) {
  stop("Fly Cell Atlas input file not found: ", flycell_file)
}

flycell <- read.csv(
  flycell_file,
  stringsAsFactors = FALSE
)

if (nrow(flycell) == 0) {
  warning("Fly Cell Atlas input table has zero rows. Output tables will be empty.")
}

if (ncol(flycell) == 0) {
  stop("Fly Cell Atlas input table has zero columns.")
}

selected_serpins <- c(
  "Spn42Da", "Spn42Db", "Spn42Dc",
  "Spn42Dd", "Spn28B", "Spn28F"
)

if (length(selected_serpins) == 0) {
  stop("selected_serpins is empty. Add at least one candidate gene before running Fly Cell Atlas analysis.")
}

flycell_clean <- flycell %>%
  rename_with(~ str_replace_all(.x, "\\.", "_")) %>%
  rename_with(~ str_replace_all(.x, " ", "_")) %>%
  rename_with(~ str_replace_all(.x, "_+", "_")) %>%
  rename_with(~ str_replace_all(.x, "^_|_$", ""))

duplicated_columns <- names(flycell_clean)[duplicated(names(flycell_clean))]

if (length(duplicated_columns) > 0) {
  stop(
    "Fly Cell Atlas column-name cleaning created duplicate column name(s): ",
    paste(unique(duplicated_columns), collapse = ", ")
  )
}

gene_col <- "Gene"

if (!gene_col %in% names(flycell_clean)) {
  stop("Fly Cell Atlas table is missing required column: Gene")
}

detect_single_column <- function(data, pattern, label) {
  matches <- names(data)[
    str_detect(names(data), regex(pattern, ignore_case = TRUE))
  ]

  if (length(matches) == 0) {
    stop("Could not detect a ", label, " column using pattern: ", pattern)
  }

  if (length(matches) > 1) {
    warning(
      "Multiple possible ",
      label,
      " columns detected: ",
      paste(matches, collapse = ", "),
      ". Using: ",
      matches[1]
    )
  }

  matches[1]
}

celltype_col <- detect_single_column(
  flycell_clean,
  "cell|annotation|type",
  "cell type"
)

mean_expr_col <- detect_single_column(
  flycell_clean,
  "mean.*expression|mean_log|expression",
  "mean expression"
)

pct_col <- detect_single_column(
  flycell_clean,
  "percent|pct",
  "percent expressing"
)

classify_cell_class <- function(x) {
  x <- as.character(x)

  case_when(
    is.na(x) ~ NA_character_,
    str_detect(x, regex("fat body", ignore_case = TRUE)) ~ "Fat body",
    str_detect(x, regex("enterocyte|midgut|hindgut|intestinal|gut", ignore_case = TRUE)) ~ "Gut",
    str_detect(x, regex("muscle|alary|heart|cardio", ignore_case = TRUE)) ~ "Muscle/heart",
    str_detect(x, regex("trache", ignore_case = TRUE)) ~ "Tracheal",
    str_detect(x, regex("malpighian|tubule", ignore_case = TRUE)) ~ "Malpighian",
    str_detect(x, regex("neuron|kenyon|brain|CNS|optic|lamina", ignore_case = TRUE)) ~ "Neuron",
    str_detect(x, regex("glia|glial", ignore_case = TRUE)) ~ "Glia",
    str_detect(x, regex("oenocyte", ignore_case = TRUE)) ~ "Oenocyte",
    str_detect(x, regex("salivary", ignore_case = TRUE)) ~ "Salivary",
    str_detect(x, regex("hemocyte|immune|antimicrobial|AMP", ignore_case = TRUE)) ~ "Immune",
    str_detect(x, regex("germline|ovary|testis|sperm|follicle|egg|oocyte", ignore_case = TRUE)) ~ "Reproductive",
    TRUE ~ "Other"
  )
}

cell_class_order <- c(
  "Fat body",
  "Gut",
  "Muscle/heart",
  "Tracheal",
  "Malpighian",
  "Neuron",
  "Glia",
  "Oenocyte",
  "Salivary",
  "Immune",
  "Reproductive"
)

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

mean_or_na <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
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

tau_specificity <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]

  if (length(x) <= 1) {
    return(NA_real_)
  }

  if (any(x < 0, na.rm = TRUE)) {
    warning("Negative expression values were encountered in tau_specificity().")
  }

  x_max <- max(x, na.rm = TRUE)

  if (is.na(x_max) || x_max == 0) {
    return(0)
  }

  x_norm <- x / x_max
  sum(1 - x_norm) / (length(x) - 1)
}

scale_centered_count <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  x_median <- median(x, na.rm = TRUE)
  distance <- abs(x - x_median)
  max_distance <- max(distance, na.rm = TRUE)

  if (is.na(max_distance) || max_distance == 0) {
    return(rep(0.5, length(x)))
  }

  1 - distance / max_distance
}

# ------------------------------------------------------------
# Collapse to functional cell classes
# ------------------------------------------------------------

flycell_candidates <- flycell_clean %>%
  filter(.data[[gene_col]] %in% selected_serpins) %>%
  mutate(
    Celltype_clean = as.character(.data[[celltype_col]]),
    Functional_Cell_Class = classify_cell_class(Celltype_clean),
    Mean_Cell_Expression = as.numeric(.data[[mean_expr_col]]),
    Percent_Expressing = as.numeric(.data[[pct_col]])
  )

missing_serpins <- setdiff(selected_serpins, unique(flycell_candidates[[gene_col]]))

if (length(missing_serpins) > 0) {
  warning(
    "The following selected serpin(s) were not found in Fly Cell Atlas: ",
    paste(missing_serpins, collapse = ", ")
  )
}

if (any(is.na(flycell_candidates$Mean_Cell_Expression))) {
  warning(
    "One or more mean expression values are NA after numeric conversion. ",
    "Missing values will be ignored in summaries."
  )
}

if (any(is.na(flycell_candidates$Percent_Expressing))) {
  warning(
    "One or more percent expressing values are NA after numeric conversion. ",
    "Missing values will be ignored in summaries."
  )
}

if (any(flycell_candidates$Percent_Expressing < 0 | flycell_candidates$Percent_Expressing > 100, na.rm = TRUE)) {
  warning("One or more Percent_Expressing values are outside the expected 0-100 range.")
}

if (any(flycell_candidates$Mean_Cell_Expression < 0, na.rm = TRUE)) {
  stop("Mean cell expression contains negative values, which are not appropriate for the pseudo-log expression color scale.")
}

unclassified_count <- sum(
  is.na(flycell_candidates$Functional_Cell_Class) |
    flycell_candidates$Functional_Cell_Class == "Other"
)

if (unclassified_count > 0) {
  warning(
    unclassified_count,
    " selected Fly Cell Atlas row(s) were unclassified or classified as Other and will be excluded."
  )
}

flycell_class <- flycell_candidates %>%
  filter(!is.na(Functional_Cell_Class), Functional_Cell_Class != "Other") %>%
  group_by(Gene, Functional_Cell_Class) %>%
  summarize(
    Cell_Mean_Expression = safe_median(Mean_Cell_Expression),
    Cell_Percent_Expressing = safe_median(Percent_Expressing),
    Cell_Max_Expression = safe_max(Mean_Cell_Expression),
    Cell_Max_Percent = safe_max(Percent_Expressing),
    N_Celltypes = n(),
    .groups = "drop"
  ) %>%
  arrange(Gene, Functional_Cell_Class)

if (nrow(flycell_class) == 0) {
  warning("No Fly Cell Atlas rows remained after filtering to selected serpins and functional cell classes.")
}

print(flycell_class)

write.csv(
  flycell_class,
  file.path(processed_dir, "FlyCellAtlas_Functional_Cell_Class_Summary.csv"),
  row.names = FALSE
)

write.csv(
  flycell_class,
  file.path(tables_dir, "FlyCellAtlas_Functional_Cell_Class_Summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Bubble plot
# ------------------------------------------------------------

flycell_plot_data <- flycell_class %>%
  mutate(
    Cell_Max_Percent = as.numeric(Cell_Max_Percent),
    Cell_Max_Expression = as.numeric(Cell_Max_Expression),
    Gene = factor(Gene, levels = selected_serpins),
    Functional_Cell_Class = factor(
      Functional_Cell_Class,
      levels = rev(cell_class_order)
    )
  )

unknown_cell_classes <- setdiff(
  unique(as.character(flycell_class$Functional_Cell_Class)),
  cell_class_order
)

if (length(unknown_cell_classes) > 0) {
  warning(
    "The following Functional_Cell_Class value(s) are not in cell_class_order and may appear as NA in plots/scores: ",
    paste(unknown_cell_classes, collapse = ", ")
  )
}

if (nrow(flycell_plot_data) == 0) {
  warning("Fly Cell Atlas bubble plot was not generated because there are no plotting rows.")
} else {
  p_flycell <- ggplot(
    flycell_plot_data,
    aes(
      x = Gene,
      y = Functional_Cell_Class
    )
  ) +
    geom_point(
      aes(
        size = Cell_Max_Percent,
        fill = Cell_Max_Expression
      ),
      shape = 21,
      color = "grey20",
      stroke = 0.35
    ) +
    scale_size_continuous(
      range = c(2, 10),
      limits = c(0, 100),
      name = "% expressing"
    ) +
    scale_fill_gradientn(
      colors = c(
        "#F7FCF5",
        "#E5F5E0",
        "#C7E9C0",
        "#A1D99B",
        "#74C476",
        "#31A354",
        "#006D2C"
      ),
      trans = scales::pseudo_log_trans(base = 2),
      name = "Mean\nexpression",
      na.value = "grey90"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Fly Cell Atlas functional cell-class expression",
      subtitle = "Dot size indicates percent expressing; color indicates expression",
      x = NULL,
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  print(p_flycell)

  ggsave(
    file.path(fig_dir, "FlyCellAtlas_Functional_Cell_Class_BubblePlot.png"),
    p_flycell,
    width = 9,
    height = 7,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(fig_dir, "FlyCellAtlas_Functional_Cell_Class_BubblePlot.pdf"),
    p_flycell,
    width = 9,
    height = 7
  )
}

# ------------------------------------------------------------
# Cell specificity score
# ------------------------------------------------------------

cell_specificity <- flycell_class %>%
  filter(Functional_Cell_Class %in% cell_class_order) %>%
  select(
    Gene,
    Functional_Cell_Class,
    Cell_Mean_Expression,
    Cell_Max_Expression,
    Cell_Max_Percent
  ) %>%
  group_by(Gene, Functional_Cell_Class) %>%
  summarize(
    Cell_Mean_Expression = mean_or_na(Cell_Mean_Expression),
    Cell_Max_Expression = safe_max(Cell_Max_Expression),
    Cell_Max_Percent = safe_max(Cell_Max_Percent),
    .groups = "drop"
  ) %>%
  complete(
    Gene = selected_serpins,
    Functional_Cell_Class = cell_class_order,
    fill = list(
      Cell_Mean_Expression = 0,
      Cell_Max_Expression = 0,
      Cell_Max_Percent = 0
    )
  ) %>%
  group_by(Gene) %>%
  summarize(
    Cell_Tau_Specificity = tau_specificity(Cell_Mean_Expression),
    Cell_Max_Class_Expression = safe_max(Cell_Max_Expression),
    Cell_Strong_Class_Count = sum(Cell_Max_Percent > 5, na.rm = TRUE),
    Cell_Detected_Class_Count = sum(Cell_Max_Percent > 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Cell_Tau_scaled = scale_01(Cell_Tau_Specificity),
    Cell_Max_Expression_scaled = scale_01(Cell_Max_Class_Expression),
    Cell_Strong_Class_Count_scaled = scale_centered_count(Cell_Strong_Class_Count)
  ) %>%
  rowwise() %>%
  mutate(
    Cell_Expression_Score = mean_or_na(
      c(
        Cell_Tau_scaled,
        Cell_Max_Expression_scaled,
        Cell_Strong_Class_Count_scaled
      )
    )
  ) %>%
  ungroup() %>%
  arrange(desc(Cell_Expression_Score), Gene)

if (nrow(cell_specificity) == 0) {
  warning("No cell specificity score rows were generated.")
}

print(cell_specificity)

write.csv(
  cell_specificity,
  file.path(processed_dir, "FlyCellAtlas_Cell_Specificity_Score.csv"),
  row.names = FALSE
)

write.csv(
  cell_specificity,
  file.path(tables_dir, "FlyCellAtlas_Cell_Specificity_Score.csv"),
  row.names = FALSE
)
