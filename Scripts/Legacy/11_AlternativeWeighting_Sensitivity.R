# ============================================================
# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 11_SensitivityAnalysis.R
# Alternative weighting models for integrated ranking
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Run from the repository root; all analysis paths are relative to it.
project_dir <- "."

analysis_dir <- project_dir
integrated_dir <- file.path(analysis_dir, "Integrated_Evidence")

integrated_processed_dir <- file.path(integrated_dir, "Data", "Processed")
fig_dir <- file.path(integrated_dir, "Figures")
tables_dir <- file.path(integrated_dir, "Tables")

dir.create(integrated_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

weighted_integrated_file <- file.path(
  integrated_processed_dir,
  "Final_Weighted_Integrated_Evidence_Table.csv"
)

if (!file.exists(weighted_integrated_file)) {
  stop("Weighted integrated evidence file not found: ", weighted_integrated_file)
}

weighted_integrated <- read.csv(
  weighted_integrated_file,
  stringsAsFactors = FALSE
)

if (nrow(weighted_integrated) == 0) {
  warning("Weighted integrated evidence table has zero rows.")
}

required_columns <- c(
  "Gene",
  "Orthology_Score",
  "Structure_Score",
  "Tissue_Expression_Score",
  "Cell_Expression_Score_Final",
  "Aging_Score"
)

missing_columns <- setdiff(required_columns, names(weighted_integrated))

if (length(missing_columns) > 0) {
  stop(
    "Weighted integrated evidence table is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(weighted_integrated$Gene) > 0) {
  duplicated_genes <- unique(weighted_integrated$Gene[duplicated(weighted_integrated$Gene)])
  stop(
    "Weighted integrated evidence table contains duplicate Gene values: ",
    paste(duplicated_genes, collapse = ", ")
  )
}

sum_or_na <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  sum(x, na.rm = TRUE)
}

# ------------------------------------------------------------
# Prepare core score table
# ------------------------------------------------------------

score_table <- weighted_integrated %>%
  transmute(
    Gene,
    Orthology_Score = as.numeric(Orthology_Score),
    Structure_Score = as.numeric(Structure_Score),
    Tissue_Expression_Score = as.numeric(Tissue_Expression_Score),
    Cell_Expression_Score = as.numeric(Cell_Expression_Score_Final),
    Aging_Score = as.numeric(Aging_Score)
  )

score_columns <- c(
  "Orthology_Score",
  "Structure_Score",
  "Tissue_Expression_Score",
  "Cell_Expression_Score",
  "Aging_Score"
)

if (any(is.na(score_table[score_columns]))) {
  warning(
    "One or more category scores are NA after numeric conversion. ",
    "Sensitivity scores will ignore missing weighted components; all-missing rows become NA."
  )
}

# ------------------------------------------------------------
# Define sensitivity models
# ------------------------------------------------------------

weight_models <- data.frame(
  Model = c(
    "Model 1: Primary model",
    "Model 2: Expression-heavy",
    "Model 3: Structure-heavy",
    "Model 4: Equal weights",
    "Model 5: Orthology-heavy",
    "Model 6: Biology-driven"
  ),
  Orthology_Weight = c(0.15, 0.10, 0.15, 0.20, 0.35, 0.10),
  Structure_Weight = c(0.25, 0.20, 0.35, 0.20, 0.20, 0.30),
  Tissue_Weight = c(0.25, 0.35, 0.20, 0.20, 0.20, 0.30),
  Cell_Weight = c(0.20, 0.20, 0.15, 0.20, 0.15, 0.20),
  Aging_Weight = c(0.15, 0.15, 0.15, 0.20, 0.10, 0.10),
  stringsAsFactors = FALSE
)

weight_models <- weight_models %>%
  mutate(
    Weight_Total =
      Orthology_Weight +
      Structure_Weight +
      Tissue_Weight +
      Cell_Weight +
      Aging_Weight
  )

invalid_weight_models <- weight_models %>%
  filter(abs(Weight_Total - 1) > 1e-8)

if (nrow(invalid_weight_models) > 0) {
  stop(
    "All sensitivity model weights must sum to 1. Invalid model(s): ",
    paste(invalid_weight_models$Model, collapse = ", ")
  )
}

print(weight_models)

# ------------------------------------------------------------
# Calculate rankings under each model
# ------------------------------------------------------------

calculate_model_score <- function(data, weights_row) {
  data %>%
    rowwise() %>%
    mutate(
      Model = weights_row$Model,
      Sensitivity_Score = sum_or_na(
        c(
          Orthology_Score * weights_row$Orthology_Weight,
          Structure_Score * weights_row$Structure_Weight,
          Tissue_Expression_Score * weights_row$Tissue_Weight,
          Cell_Expression_Score * weights_row$Cell_Weight,
          Aging_Score * weights_row$Aging_Weight
        )
      )
    ) %>%
    ungroup() %>%
    mutate(
      Sensitivity_Rank = rank(
        -Sensitivity_Score,
        ties.method = "min",
        na.last = "keep"
      )
    ) %>%
    arrange(Sensitivity_Rank, desc(Sensitivity_Score), Gene)
}

sensitivity_results <- bind_rows(
  lapply(
    seq_len(nrow(weight_models)),
    function(i) calculate_model_score(score_table, weight_models[i, ])
  )
)

if (any(is.na(sensitivity_results$Sensitivity_Score))) {
  warning("One or more sensitivity scores are missing.")
}

print(sensitivity_results)

write.csv(
  sensitivity_results,
  file.path(integrated_processed_dir, "Sensitivity_Analysis_All_Weighting_Models.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Top 3 table
# ------------------------------------------------------------

sensitivity_top3 <- sensitivity_results %>%
  filter(!is.na(Sensitivity_Rank), Sensitivity_Rank <= 3) %>%
  select(
    Model,
    Sensitivity_Rank,
    Gene,
    Sensitivity_Score
  ) %>%
  arrange(Model, Sensitivity_Rank, desc(Sensitivity_Score), Gene)

write.csv(
  sensitivity_top3,
  file.path(tables_dir, "Sensitivity_Analysis_Top3_By_Model.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Wide rank table
# ------------------------------------------------------------

sensitivity_rank_table <- sensitivity_results %>%
  select(Model, Gene, Sensitivity_Rank, Sensitivity_Score) %>%
  mutate(
    Score_Rank = if_else(
      is.na(Sensitivity_Score) | is.na(Sensitivity_Rank),
      NA_character_,
      paste0(
        round(Sensitivity_Score, 3),
        " (rank ",
        Sensitivity_Rank,
        ")"
      )
    )
  ) %>%
  select(Model, Gene, Score_Rank) %>%
  pivot_wider(
    names_from = Model,
    values_from = Score_Rank
  )

write.csv(
  sensitivity_rank_table,
  file.path(tables_dir, "Sensitivity_Analysis_Rank_Table.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot scores by model
# ------------------------------------------------------------

plot_data <- sensitivity_results %>%
  filter(!is.na(Sensitivity_Score)) %>%
  group_by(Model) %>%
  arrange(Sensitivity_Score, .by_group = TRUE) %>%
  mutate(
    Gene_Model = factor(
      paste(Gene, Model, sep = "___"),
      levels = paste(Gene, Model, sep = "___")
    )
  ) %>%
  ungroup()

if (nrow(plot_data) == 0) {
  warning("No non-missing sensitivity scores are available for plotting score panels.")
} else {
  x_min <- min(0, min(plot_data$Sensitivity_Score, na.rm = TRUE))
  x_max <- max(plot_data$Sensitivity_Score, na.rm = TRUE)

  if (x_max == x_min) {
    x_limits <- c(x_min - 0.05, x_max + 0.05)
  } else {
    x_padding <- (x_max - x_min) * 0.15
    x_limits <- c(x_min, x_max + x_padding)
  }

  p_sensitivity_scores <- ggplot(
    plot_data,
    aes(
      x = Sensitivity_Score,
      y = Gene_Model
    )
  ) +
    geom_col(width = 0.7, fill = "#2166AC") +
    facet_wrap(~ Model, ncol = 2, scales = "free_y") +
    geom_text(
      aes(label = round(Sensitivity_Score, 3)),
      hjust = -0.15,
      size = 3
    ) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    coord_cartesian(
      xlim = x_limits
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Sensitivity analysis of integrated candidate ranking",
      subtitle = "Candidate rankings under alternative evidence-weighting schemes",
      x = "Weighted integrated evidence score",
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  print(p_sensitivity_scores)

  ggsave(
    file.path(fig_dir, "Sensitivity_Analysis_All_Weighting_Models.png"),
    p_sensitivity_scores,
    width = 10,
    height = 10,
    dpi = 600,
    bg = "white"
  )
}

# ------------------------------------------------------------
# Rank stability plot
# ------------------------------------------------------------

rank_plot_data <- sensitivity_results %>%
  filter(!is.na(Sensitivity_Rank))

if (nrow(rank_plot_data) == 0) {
  warning("No non-missing sensitivity ranks are available for plotting rank stability.")
} else {
  p_rank_stability <- ggplot(
    rank_plot_data,
    aes(
      x = Model,
      y = Sensitivity_Rank,
      group = Gene,
      color = Gene
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_y_reverse(
      breaks = seq_len(max(rank_plot_data$Sensitivity_Rank, na.rm = TRUE))
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Rank stability across weighting models",
      x = NULL,
      y = "Rank"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )

  print(p_rank_stability)

  ggsave(
    file.path(fig_dir, "Sensitivity_Analysis_Rank_Stability.png"),
    p_rank_stability,
    width = 9,
    height = 5,
    dpi = 600,
    bg = "white"
  )
}
