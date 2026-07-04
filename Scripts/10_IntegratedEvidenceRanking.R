# ============================================================
# 10_IntegratedEvidenceRanking.R
# Final weighted integrated evidence ranking
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Prefer an environment variable for reproducibility across machines.
# If it is not set, fall back to the original local project path.
project_dir <- Sys.getenv(
  "xxxxxx_xxxxxxxx_PROJECT_DIR",
  unset = "D:/xxxxxx-xxxxxxxx lab/xxxxxx-xxxxxxxx lab"
)

chapter3_dir <- file.path(project_dir, "Bioinformatics chapter three")
expression_dir <- file.path(
  project_dir,
  "Bioinformatics_tissue expression collection",
  "Bioinformatics_tissue_expression"
)

chapter_processed_dir <- file.path(chapter3_dir, "Data", "Processed")
expression_processed_dir <- file.path(expression_dir, "Data", "Processed")

integrated_dir <- file.path(chapter3_dir, "Integrated_Evidence")
integrated_processed_dir <- file.path(integrated_dir, "Data", "Processed")
fig_dir <- file.path(integrated_dir, "Figures")
tables_dir <- file.path(integrated_dir, "Tables")

dir.create(integrated_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

selected_serpins <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)

if (length(selected_serpins) == 0) {
  stop("selected_serpins is empty. Add at least one candidate gene before ranking integrated evidence.")
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

mean_or_na <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

sum_or_na <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  sum(x, na.rm = TRUE)
}

check_required_columns <- function(data, required_columns, label) {
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

check_duplicate_genes <- function(data, label) {
  duplicated_genes <- data$Gene[duplicated(data$Gene)]

  if (length(duplicated_genes) > 0) {
    stop(
      label,
      " contains duplicate Gene values: ",
      paste(unique(duplicated_genes), collapse = ", ")
    )
  }
}

# ------------------------------------------------------------
# Load component files
# ------------------------------------------------------------

component_files <- c(
  DIOPT = file.path(chapter_processed_dir, "Table_DIOPT_candidate_serpin_orthology.csv"),
  Full_Length = file.path(chapter_processed_dir, "Full_Length_Alignment_Stats_vs_Human_PAI1.csv"),
  RCL = file.path(chapter_processed_dir, "RCL_Alignment_Stats_vs_Human_PAI1.csv"),
  Phylogeny = file.path(chapter_processed_dir, "Phylogenetic_Distances_to_Human_PAI1.csv"),
  Tissue = file.path(expression_processed_dir, "TissueExpression_Similarity_To_SERPINE1.csv"),
  Cell = file.path(expression_processed_dir, "FlyCellAtlas_Cell_Specificity_Score.csv"),
  Aging = file.path(integrated_processed_dir, "Matched_Tissue_CrossSpecies_Aging_Similarity_to_SERPINE1.csv")
)

missing_files <- component_files[!file.exists(component_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required component file(s): ",
    paste(names(missing_files), basename(missing_files), sep = " = ", collapse = "; ")
  )
}

diopt <- read.csv(component_files[["DIOPT"]], stringsAsFactors = FALSE)
full_length <- read.csv(component_files[["Full_Length"]], stringsAsFactors = FALSE)
rcl <- read.csv(component_files[["RCL"]], stringsAsFactors = FALSE)
phylo <- read.csv(component_files[["Phylogeny"]], stringsAsFactors = FALSE)
tissue <- read.csv(component_files[["Tissue"]], stringsAsFactors = FALSE)
cell <- read.csv(component_files[["Cell"]], stringsAsFactors = FALSE)
aging <- read.csv(component_files[["Aging"]], stringsAsFactors = FALSE)

required_diopt_columns <- c(
  "Gene",
  "Primary_Human_Ortholog",
  "Primary_DIOPT_Score",
  "Primary_Weighted_Score",
  "Supporting_DBs",
  "Reverse_Best",
  "SERPINE1_DIOPT_Score",
  "SERPINE1_Weighted_Score"
)

required_full_columns <- c(
  "Gene",
  "Full_Length_Identity",
  "Full_Length_Similarity"
)

required_rcl_columns <- c(
  "Gene",
  "RCL_Identity",
  "RCL_Similarity"
)

required_phylo_columns <- c(
  "Gene",
  "Full_Length_ML_Distance_to_PAI1"
)

required_tissue_columns <- c(
  "Gene",
  "Pearson",
  "Spearman",
  "Cosine",
  "Euclidean",
  "Integrated_Tissue_Expression_Similarity"
)

required_cell_columns <- c(
  "Gene",
  "Cell_Tau_Specificity",
  "Cell_Max_Class_Expression",
  "Cell_Strong_Class_Count",
  "Cell_Detected_Class_Count",
  "Cell_Expression_Score"
)

required_aging_columns <- c(
  "Gene",
  "Aging_Tissue_Pearson",
  "Aging_Tissue_Spearman",
  "Aging_Tissue_Cosine",
  "Aging_Tissue_Euclidean",
  "Matched_Tissue_Aging_Similarity"
)

check_required_columns(diopt, required_diopt_columns, "DIOPT component table")
check_required_columns(full_length, required_full_columns, "Full-length alignment component table")
check_required_columns(rcl, required_rcl_columns, "RCL alignment component table")
check_required_columns(phylo, required_phylo_columns, "Phylogeny component table")
check_required_columns(tissue, required_tissue_columns, "Tissue-expression component table")
check_required_columns(cell, required_cell_columns, "Cell-expression component table")
check_required_columns(aging, required_aging_columns, "Aging component table")

check_duplicate_genes(diopt, "DIOPT component table")
check_duplicate_genes(full_length, "Full-length alignment component table")
check_duplicate_genes(rcl, "RCL alignment component table")
check_duplicate_genes(phylo, "Phylogeny component table")
check_duplicate_genes(tissue, "Tissue-expression component table")
check_duplicate_genes(cell, "Cell-expression component table")
check_duplicate_genes(aging, "Aging component table")

# ------------------------------------------------------------
# Build master evidence table
# ------------------------------------------------------------

master <- data.frame(Gene = selected_serpins, stringsAsFactors = FALSE) %>%
  left_join(diopt, by = "Gene") %>%
  left_join(full_length, by = "Gene") %>%
  left_join(rcl, by = "Gene") %>%
  left_join(phylo, by = "Gene") %>%
  left_join(
    tissue %>%
      select(
        Gene,
        Tissue_Pearson = Pearson,
        Tissue_Spearman = Spearman,
        Tissue_Cosine = Cosine,
        Tissue_Euclidean = Euclidean,
        Tissue_Expression_Similarity =
          Integrated_Tissue_Expression_Similarity
      ),
    by = "Gene"
  ) %>%
  left_join(
    cell %>%
      select(
        Gene,
        Cell_Tau_Specificity,
        Cell_Max_Class_Expression,
        Cell_Strong_Class_Count,
        Cell_Detected_Class_Count,
        Cell_Expression_Score
      ),
    by = "Gene"
  ) %>%
  left_join(
    aging %>%
      select(
        Gene,
        Aging_Tissue_Pearson,
        Aging_Tissue_Spearman,
        Aging_Tissue_Cosine,
        Aging_Tissue_Euclidean,
        Matched_Tissue_Aging_Similarity
      ),
    by = "Gene"
  )

component_presence_columns <- c(
  Primary_DIOPT_Score = "DIOPT",
  Full_Length_Similarity = "Full-length alignment",
  RCL_Similarity = "RCL alignment",
  Full_Length_ML_Distance_to_PAI1 = "Phylogeny",
  Tissue_Expression_Similarity = "Tissue expression",
  Cell_Expression_Score = "Cell expression",
  Matched_Tissue_Aging_Similarity = "Aging"
)

for (column_name in names(component_presence_columns)) {
  missing_genes <- master %>%
    filter(is.na(.data[[column_name]])) %>%
    pull(Gene)

  if (length(missing_genes) > 0) {
    warning(
      component_presence_columns[[column_name]],
      " evidence is missing for gene(s): ",
      paste(missing_genes, collapse = ", ")
    )
  }
}

write.csv(
  master,
  file.path(integrated_processed_dir, "Master_Evidence_Table_All_Components.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Category scores
# ------------------------------------------------------------

score_numeric_columns <- c(
  "Primary_DIOPT_Score",
  "Primary_Weighted_Score",
  "Supporting_DBs",
  "Reverse_Best",
  "SERPINE1_DIOPT_Score",
  "SERPINE1_Weighted_Score",
  "Full_Length_Identity",
  "Full_Length_Similarity",
  "RCL_Identity",
  "RCL_Similarity",
  "Full_Length_ML_Distance_to_PAI1",
  "Tissue_Pearson",
  "Tissue_Spearman",
  "Tissue_Cosine",
  "Tissue_Euclidean",
  "Tissue_Expression_Similarity",
  "Cell_Expression_Score",
  "Matched_Tissue_Aging_Similarity"
)

master <- master %>%
  mutate(
    across(all_of(score_numeric_columns), as.numeric)
  )

if (any(is.na(master[score_numeric_columns]))) {
  warning(
    "One or more numeric evidence values are NA after conversion. ",
    "Category scores will ignore missing components; all-missing categories become NA."
  )
}

category_weights <- c(
  Orthology_Score = 0.15,
  Structure_Score = 0.25,
  Tissue_Expression_Score = 0.25,
  Cell_Expression_Score_Final = 0.20,
  Aging_Score = 0.15
)

if (abs(sum(category_weights) - 1) > 1e-8) {
  stop("Integrated evidence weights must sum to 1.")
}

scored <- master %>%
  mutate(
    # Orthology
    DIOPT_scaled = scale_01(Primary_DIOPT_Score),
    Weighted_DIOPT_scaled = scale_01(Primary_Weighted_Score),
    Supporting_DB_scaled = scale_01(Supporting_DBs),
    Reverse_Best_scaled = Reverse_Best,
    SERPINE1_DIOPT_scaled = scale_01(SERPINE1_DIOPT_Score),
    SERPINE1_Weighted_scaled = scale_01(SERPINE1_Weighted_Score),

    # Structure
    Full_Length_Identity_scaled = scale_01(Full_Length_Identity),
    Full_Length_Similarity_scaled = scale_01(Full_Length_Similarity),
    RCL_Identity_scaled = scale_01(RCL_Identity),
    RCL_Similarity_scaled = scale_01(RCL_Similarity),
    Phylo_Distance_scaled =
      1 - scale_01(Full_Length_ML_Distance_to_PAI1),

    # Tissue expression
    Tissue_Pearson_scaled = scale_01(Tissue_Pearson),
    Tissue_Spearman_scaled = scale_01(Tissue_Spearman),
    Tissue_Cosine_scaled = scale_01(Tissue_Cosine),
    Tissue_Euclidean_scaled = 1 - scale_01(Tissue_Euclidean),
    Tissue_Expression_Similarity_scaled =
      scale_01(Tissue_Expression_Similarity),

    # Cell
    Cell_Expression_scaled = scale_01(Cell_Expression_Score),

    # Aging
    Aging_scaled = scale_01(Matched_Tissue_Aging_Similarity)
  ) %>%
  rowwise() %>%
  mutate(
    Orthology_Score = mean_or_na(
      c(
        DIOPT_scaled,
        Weighted_DIOPT_scaled,
        Supporting_DB_scaled,
        Reverse_Best_scaled,
        SERPINE1_DIOPT_scaled,
        SERPINE1_Weighted_scaled
      )
    ),

    Structure_Score = mean_or_na(
      c(
        Full_Length_Identity_scaled,
        Full_Length_Similarity_scaled,
        RCL_Identity_scaled,
        RCL_Similarity_scaled,
        Phylo_Distance_scaled
      )
    ),

    Tissue_Expression_Score = mean_or_na(
      c(
        Tissue_Pearson_scaled,
        Tissue_Spearman_scaled,
        Tissue_Cosine_scaled,
        Tissue_Euclidean_scaled,
        Tissue_Expression_Similarity_scaled
      )
    ),

    Cell_Expression_Score_Final = Cell_Expression_scaled,

    Aging_Score = Aging_scaled,

    Weighted_Integrated_Evidence_Score = sum_or_na(
      c(
        Orthology_Score * category_weights[["Orthology_Score"]],
        Structure_Score * category_weights[["Structure_Score"]],
        Tissue_Expression_Score * category_weights[["Tissue_Expression_Score"]],
        Cell_Expression_Score_Final * category_weights[["Cell_Expression_Score_Final"]],
        Aging_Score * category_weights[["Aging_Score"]]
      )
    )
  ) %>%
  ungroup() %>%
  mutate(
    Weighted_Rank = rank(
      -Weighted_Integrated_Evidence_Score,
      ties.method = "min",
      na.last = "keep"
    )
  ) %>%
  arrange(Weighted_Rank, desc(Weighted_Integrated_Evidence_Score), Gene)

if (any(is.na(scored$Weighted_Integrated_Evidence_Score))) {
  warning("One or more genes have missing weighted integrated evidence scores.")
}

write.csv(
  scored,
  file.path(integrated_processed_dir, "Final_Weighted_Integrated_Evidence_Table.csv"),
  row.names = FALSE
)

final_matrix <- scored %>%
  transmute(
    Rank = Weighted_Rank,
    Gene,
    Primary_Human_Ortholog,
    Orthology_Score = round(Orthology_Score, 3),
    Structure_Score = round(Structure_Score, 3),
    Tissue_Expression_Score = round(Tissue_Expression_Score, 3),
    Cell_Expression_Score = round(Cell_Expression_Score_Final, 3),
    Aging_Score = round(Aging_Score, 3),
    Weighted_Integrated_Evidence_Score =
      round(Weighted_Integrated_Evidence_Score, 3)
  )

if (nrow(final_matrix) == 0) {
  warning("The final integrated evidence matrix has zero rows.")
}

print(final_matrix)

write.csv(
  final_matrix,
  file.path(tables_dir, "Table_3_2_Weighted_Integrated_Evidence_Matrix_Final_Ranking.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot final ranking
# ------------------------------------------------------------

plot_data <- final_matrix %>%
  filter(!is.na(Weighted_Integrated_Evidence_Score))

if (nrow(plot_data) == 0) {
  warning("No non-missing weighted integrated evidence scores are available for plotting.")
} else {
  x_min <- min(0, min(plot_data$Weighted_Integrated_Evidence_Score, na.rm = TRUE))
  x_max <- max(plot_data$Weighted_Integrated_Evidence_Score, na.rm = TRUE)

  if (x_max == x_min) {
    x_limits <- c(x_min - 0.05, x_max + 0.05)
  } else {
    x_padding <- (x_max - x_min) * 0.15
    x_limits <- c(x_min, x_max + x_padding)
  }

  p_rank <- ggplot(
    plot_data,
    aes(
      x = Weighted_Integrated_Evidence_Score,
      y = reorder(Gene, Weighted_Integrated_Evidence_Score)
    )
  ) +
    geom_col(width = 0.7, fill = "#2166AC") +
    geom_text(
      aes(label = Weighted_Integrated_Evidence_Score),
      hjust = -0.15,
      size = 3.5
    ) +
    coord_cartesian(
      xlim = x_limits
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Weighted final integrated candidate ranking",
      x = "Weighted integrated evidence score",
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  print(p_rank)

  ggsave(
    file.path(fig_dir, "Weighted_final_integrated_candidate_ranking.png"),
    p_rank,
    width = 7,
    height = 4.5,
    dpi = 600,
    bg = "white"
  )
}
