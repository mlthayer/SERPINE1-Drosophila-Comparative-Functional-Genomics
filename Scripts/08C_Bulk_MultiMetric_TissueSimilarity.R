# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Run from the repository root after step 08B.
# ==============================================================================
# MULTI-METRIC CROSS-SPECIES BULK TISSUE-EXPRESSION SIMILARITY
#
# Drosophila:
#   FlyAtlas2 adult bulk tissue expression
#
# Human:
#   GTEx v11 bulk tissue expression
#   Human Protein Atlas bulk tissue expression
#
# Goal:
#   Compare tissue-expression DISTRIBUTION patterns for:
#
#       6 Drosophila serpins
#                 vs
#       6 human SERPIN genes
#
# separately by:
#   - Female / Male
#   - GTEx / HPA
#
# Similarity metrics:
#   1. Spearman rho
#   2. Kendall tau
#   3. Cosine similarity of within-pair relative expression profiles
#   4. Top-3 functional-system Jaccard overlap
#
# Descriptive multi-metric score:
#
#   Spearman and Kendall are transformed from [-1,1] to [0,1].
#   Cosine and Jaccard are already [0,1].
#
#   Composite score =
#       mean of available similarity components
#
# IMPORTANT:
#   The composite score is DESCRIPTIVE.
#   It is not a statistical test and the component metrics are not
#   independent forms of evidence.
#
# Profiles with:
#   - fewer than 3 matched systems
#   - no variation in either gene
#   - zero total expression
#
# are flagged as not estimable.
# ==============================================================================


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

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

human_root <- "Expression/Human_SERPIN_reanalysis"


bulk_root <- file.path(
  human_root,
  "CrossSpecies_Bulk_BySex"
)


input_table_dir <- file.path(
  bulk_root,
  "Tables"
)


output_root <- file.path(
  human_root,
  "CrossSpecies_Bulk_MultiMetric"
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


# ==============================================================================
# 3. INPUT FILES
# ==============================================================================

fly_file <- file.path(
  input_table_dir,
  "FlyAtlas2_system_expression_by_sex.csv"
)


gtex_file <- file.path(
  input_table_dir,
  "GTEx_system_expression_by_sex.csv"
)


hpa_file <- file.path(
  input_table_dir,
  "HPA_system_expression_by_sex.csv"
)


stopifnot(
  file.exists(fly_file),
  file.exists(gtex_file),
  file.exists(hpa_file)
)


# ==============================================================================
# 4. READ SYSTEM-LEVEL TABLES
# ==============================================================================

fly_system <- read_csv(
  fly_file,
  show_col_types = FALSE
)


gtex_system <- read_csv(
  gtex_file,
  show_col_types = FALSE
)


hpa_system <- read_csv(
  hpa_file,
  show_col_types = FALSE
)


cat("\n============================================\n")
cat("INPUT TABLES LOADED\n")
cat("============================================\n")


cat(
  "\nFlyAtlas2 rows:",
  nrow(fly_system),
  "\n"
)


cat(
  "GTEx rows:",
  nrow(gtex_system),
  "\n"
)


cat(
  "HPA rows:",
  nrow(hpa_system),
  "\n"
)


# ==============================================================================
# 5. TARGET GENES
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
# 6. FUNCTIONAL SYSTEM ORDER
# ==============================================================================

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
# 7. HELPER FUNCTIONS
# ==============================================================================


# ------------------------------------------------------------------------------
# 7A. Convert an expression vector to relative abundance
#
# Example:
#   20, 30, 50
#
# becomes:
#   0.20, 0.30, 0.50
#
# This removes absolute FPKM/TPM scale differences.
# ------------------------------------------------------------------------------

relative_profile <- function(x) {

  if (
    all(is.na(x))
  ) {

    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }


  total <- sum(
    x,
    na.rm = TRUE
  )


  if (
    !is.finite(total) ||
      total <= 0
  ) {

    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }


  x / total
}


# ------------------------------------------------------------------------------
# 7B. Cosine similarity
#
# Range for these non-negative relative profiles:
#   approximately 0 to 1
#
# 1 = same expression-distribution shape
# ------------------------------------------------------------------------------

cosine_similarity <- function(x, y) {

  good <- complete.cases(
    x,
    y
  )


  x <- x[good]
  y <- y[good]


  if (
    length(x) < 2
  ) {

    return(
      NA_real_
    )
  }


  denominator <-
    sqrt(
      sum(
        x^2
      )
    ) *
    sqrt(
      sum(
        y^2
      )
    )


  if (
    !is.finite(denominator) ||
      denominator == 0
  ) {

    return(
      NA_real_
    )
  }


  sum(
    x * y
  ) / denominator
}


# ------------------------------------------------------------------------------
# 7C. Top-system set
#
# Uses min_rank() so tied systems at the top are retained rather than
# arbitrarily discarded.
# ------------------------------------------------------------------------------

get_top_systems <- function(
    systems,
    expression,
    k = 3
) {

  dat <- tibble(
    System = systems,
    Expression = expression
  ) %>%

    filter(
      !is.na(Expression)
    )


  if (
    nrow(dat) == 0
  ) {

    return(
      character(0)
    )
  }


  # Constant profiles have no meaningful "top tissues"
  if (
    n_distinct(
      dat$Expression
    ) < 2
  ) {

    return(
      character(0)
    )
  }


  dat %>%

    mutate(
      Rank =
        min_rank(
          desc(
            Expression
          )
        )
    ) %>%

    filter(
      Rank <= k
    ) %>%

    pull(
      System
    )
}


# ------------------------------------------------------------------------------
# 7D. Jaccard similarity
# ------------------------------------------------------------------------------

jaccard_similarity <- function(
    set_a,
    set_b
) {

  if (
    length(set_a) == 0 ||
      length(set_b) == 0
  ) {

    return(
      NA_real_
    )
  }


  union_set <- union(
    set_a,
    set_b
  )


  intersection_set <- intersect(
    set_a,
    set_b
  )


  if (
    length(union_set) == 0
  ) {

    return(
      NA_real_
    )
  }


  length(
    intersection_set
  ) /
    length(
      union_set
    )
}


# ==============================================================================
# 8. FUNCTION TO BUILD ONE MATCHED FLY-HUMAN PROFILE
# ==============================================================================

build_matched_profile <- function(
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


  # Fly ------------------------------------------------------------------------

  fly_one <- fly_system %>%

    filter(
      Sex == .env$sex_value,
      Gene == .env$fly_gene,
      System %in% .env$allowed_systems
    ) %>%

    transmute(
      System,
      Fly_Expression =
        Median_FPKM
    )


  # Human ----------------------------------------------------------------------

  human_source <- if (
    dataset_value == "GTEx"
  ) {

    gtex_system

  } else {

    hpa_system
  }


  human_one <- human_source %>%

    filter(
      Sex == .env$sex_value,
      Gene == .env$human_gene,
      System %in% .env$allowed_systems
    ) %>%

    transmute(
      System,
      Human_Expression =
        Median_TPM
    )


  # Match ----------------------------------------------------------------------

  inner_join(
    fly_one,
    human_one,
    by = "System",
    relationship = "one-to-one"
  ) %>%

    filter(
      !is.na(Fly_Expression),
      !is.na(Human_Expression)
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


# ==============================================================================
# 9. BUILD ALL 144 MATCHED COMPARISONS
# ==============================================================================

comparison_grid <- crossing(
  dataset_value = c(
    "GTEx",
    "HPA"
  ),
  sex_value = c(
    "Female",
    "Male"
  ),
  fly_gene =
    fly_genes,
  human_gene =
    human_genes
)


matched_profiles <- comparison_grid %>%

  purrr::pmap_dfr(
    build_matched_profile
  )


write_csv(
  matched_profiles,
  file.path(
    table_dir,
    "All_matched_bulk_expression_profiles.csv"
  )
)


# ==============================================================================
# 10. DUPLICATE CHECK
# ==============================================================================

duplicate_check <- matched_profiles %>%

  count(
    Dataset,
    Sex,
    Fly_Gene,
    Human_Gene,
    System
  ) %>%

  filter(
    n > 1
  )


cat("\n============================================\n")
cat("DUPLICATE CHECK\n")
cat("============================================\n")


print(
  duplicate_check,
  n = Inf
)


if (
  nrow(duplicate_check) > 0
) {

  stop(
    "Duplicate matched system rows detected."
  )
}


# ==============================================================================
# 11. FUNCTION TO CALCULATE MULTIPLE SIMILARITY METRICS
# ==============================================================================

calculate_similarity_metrics <- function(dat) {

  dat <- dat %>%

    filter(
      complete.cases(
        Fly_Expression,
        Human_Expression
      )
    )


  n_systems <- nrow(
    dat
  )


  n_unique_fly <-
    n_distinct(
      dat$Fly_Expression
    )


  n_unique_human <-
    n_distinct(
      dat$Human_Expression
    )


  fly_total <- sum(
    dat$Fly_Expression,
    na.rm = TRUE
  )


  human_total <- sum(
    dat$Human_Expression,
    na.rm = TRUE
  )


  # --------------------------------------------------------------------------
  # Determine whether comparison is estimable
  # --------------------------------------------------------------------------

  estimable <-
    n_systems >= 3 &&
    n_unique_fly >= 2 &&
    n_unique_human >= 2 &&
    fly_total > 0 &&
    human_total > 0


  if (
    !estimable
  ) {

    reason <- case_when(

      n_systems < 3 ~
        "Fewer than 3 matched systems",

      n_unique_fly < 2 ~
        "Fly profile has insufficient variation",

      n_unique_human < 2 ~
        "Human profile has insufficient variation",

      fly_total <= 0 ~
        "Fly profile has zero total expression",

      human_total <= 0 ~
        "Human profile has zero total expression",

      TRUE ~
        "Not estimable"
    )


    return(
      tibble(
        N_Systems =
          n_systems,

        N_Unique_Fly =
          n_unique_fly,

        N_Unique_Human =
          n_unique_human,

        Estimable =
          FALSE,

        Exclusion_Reason =
          reason,

        Spearman_Rho =
          NA_real_,

        Spearman_P =
          NA_real_,

        Kendall_Tau =
          NA_real_,

        Kendall_P =
          NA_real_,

        Cosine_Similarity =
          NA_real_,

        Top3_Jaccard =
          NA_real_,

        N_Top_Fly =
          NA_integer_,

        N_Top_Human =
          NA_integer_,

        Shared_Top_Systems =
          NA_character_,

        Composite_Similarity =
          NA_real_
      )
    )
  }


  # --------------------------------------------------------------------------
  # Spearman
  # --------------------------------------------------------------------------

  spearman_test <- suppressWarnings(
    cor.test(
      dat$Fly_Expression,
      dat$Human_Expression,
      method = "spearman",
      exact = FALSE
    )
  )


  spearman_rho <-
    unname(
      spearman_test$estimate
    )


  # --------------------------------------------------------------------------
  # Kendall
  # --------------------------------------------------------------------------

  kendall_test <- suppressWarnings(
    cor.test(
      dat$Fly_Expression,
      dat$Human_Expression,
      method = "kendall",
      exact = FALSE
    )
  )


  kendall_tau <-
    unname(
      kendall_test$estimate
    )


  # --------------------------------------------------------------------------
  # Relative abundance profiles
  # --------------------------------------------------------------------------

  fly_relative <-
    relative_profile(
      dat$Fly_Expression
    )


  human_relative <-
    relative_profile(
      dat$Human_Expression
    )


  cosine <-
    cosine_similarity(
      fly_relative,
      human_relative
    )


  # --------------------------------------------------------------------------
  # Top-3 system overlap
  #
  # If fewer than 3 systems exist, all systems would be used, but estimable
  # comparisons already require at least 3 systems.
  # --------------------------------------------------------------------------

  top_k <- min(
    3,
    n_systems
  )


  fly_top <- get_top_systems(
    systems =
      dat$System,
    expression =
      dat$Fly_Expression,
    k =
      top_k
  )


  human_top <- get_top_systems(
    systems =
      dat$System,
    expression =
      dat$Human_Expression,
    k =
      top_k
  )


  jaccard <-
    jaccard_similarity(
      fly_top,
      human_top
    )


  shared_top <- intersect(
    fly_top,
    human_top
  )


  # --------------------------------------------------------------------------
  # Transform rank metrics to common [0,1] descriptive scale
  # --------------------------------------------------------------------------

  spearman_scaled <-
    (
      spearman_rho + 1
    ) / 2


  kendall_scaled <-
    (
      kendall_tau + 1
    ) / 2


  components <- c(
    spearman_scaled,
    kendall_scaled,
    cosine,
    jaccard
  )


  composite <- if (
    all(
      is.na(
        components
      )
    )
  ) {

    NA_real_

  } else {

    mean(
      components,
      na.rm = TRUE
    )
  }


  tibble(

    N_Systems =
      n_systems,

    N_Unique_Fly =
      n_unique_fly,

    N_Unique_Human =
      n_unique_human,

    Estimable =
      TRUE,

    Exclusion_Reason =
      NA_character_,

    Spearman_Rho =
      spearman_rho,

    Spearman_P =
      spearman_test$p.value,

    Kendall_Tau =
      kendall_tau,

    Kendall_P =
      kendall_test$p.value,

    Cosine_Similarity =
      cosine,

    Top3_Jaccard =
      jaccard,

    N_Top_Fly =
      length(
        fly_top
      ),

    N_Top_Human =
      length(
        human_top
      ),

    Shared_Top_Systems =
      if (
        length(
          shared_top
        ) == 0
      ) {
        "None"
      } else {
        paste(
          shared_top,
          collapse = "; "
        )
      },

    Composite_Similarity =
      composite
  )
}


# ==============================================================================
# 12. CALCULATE METRICS FOR EVERY GENE PAIR
# ==============================================================================

similarity_results <- matched_profiles %>%

  group_by(
    Dataset,
    Sex,
    Fly_Gene,
    Human_Gene
  ) %>%

  group_modify(
    ~ calculate_similarity_metrics(.x)
  ) %>%

  ungroup()


# ==============================================================================
# 13. MULTIPLE-TEST CORRECTION FOR INFERENTIAL RANK METRICS
#
# Spearman and Kendall are corrected separately within Dataset × Sex.
#
# Composite similarity itself does NOT receive a p-value.
# ==============================================================================

similarity_results <- similarity_results %>%

  group_by(
    Dataset,
    Sex
  ) %>%

  mutate(

    Spearman_FDR =
      p.adjust(
        Spearman_P,
        method = "BH"
      ),

    Kendall_FDR =
      p.adjust(
        Kendall_P,
        method = "BH"
      )
  ) %>%

  ungroup()


# ==============================================================================
# 14. RANK FLY GENES FOR EACH HUMAN SERPIN
#
# Rank 1 = strongest multi-metric tissue-profile similarity.
# ==============================================================================

similarity_results <- similarity_results %>%

  group_by(
    Dataset,
    Sex,
    Human_Gene
  ) %>%

  mutate(

    Composite_Rank =
      if_else(
        Estimable,
        min_rank(
          desc(
            Composite_Similarity
          )
        ),
        NA_integer_
      ),

    Spearman_Rank =
      if_else(
        Estimable,
        min_rank(
          desc(
            Spearman_Rho
          )
        ),
        NA_integer_
      ),

    Kendall_Rank =
      if_else(
        Estimable,
        min_rank(
          desc(
            Kendall_Tau
          )
        ),
        NA_integer_
      ),

    Cosine_Rank =
      if_else(
        Estimable,
        min_rank(
          desc(
            Cosine_Similarity
          )
        ),
        NA_integer_
      ),

    Jaccard_Rank =
      if_else(
        Estimable,
        min_rank(
          desc(
            Top3_Jaccard
          )
        ),
        NA_integer_
      )
  ) %>%

  ungroup()


write_csv(
  similarity_results,
  file.path(
    table_dir,
    "Bulk_multimetric_similarity_all_pairs.csv"
  )
)


# ==============================================================================
# 15. ESTIMABILITY AUDIT
# ==============================================================================

estimability_audit <- similarity_results %>%

  count(
    Dataset,
    Sex,
    Human_Gene,
    Estimable,
    Exclusion_Reason
  ) %>%

  arrange(
    Dataset,
    Sex,
    Human_Gene,
    desc(
      Estimable
    )
  )


write_csv(
  estimability_audit,
  file.path(
    table_dir,
    "Bulk_multimetric_estimability_audit.csv"
  )
)


cat("\n============================================\n")
cat("ESTIMABILITY AUDIT\n")
cat("============================================\n")


print(
  estimability_audit,
  n = Inf
)


# ==============================================================================
# 16. CROSS-DATASET CONSENSUS WITHIN EACH SEX
#
# Rather than averaging raw metric values from GTEx and HPA, average the
# within-analysis fly-gene ranks.
#
# This prevents one dataset from dominating merely because its numeric
# similarity distribution differs from the other dataset.
# ==============================================================================

sex_consensus <- similarity_results %>%

  filter(
    Estimable,
    !is.na(
      Composite_Rank
    )
  ) %>%

  group_by(
    Sex,
    Human_Gene,
    Fly_Gene
  ) %>%

  summarise(

    N_Human_Datasets =
      n_distinct(
        Dataset
      ),

    Mean_Composite_Rank =
      mean(
        Composite_Rank,
        na.rm = TRUE
      ),

    Median_Composite_Rank =
      median(
        Composite_Rank,
        na.rm = TRUE
      ),

    Mean_Composite_Similarity =
      mean(
        Composite_Similarity,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%

  group_by(
    Sex,
    Human_Gene
  ) %>%

  mutate(

    Sex_Consensus_Rank =
      min_rank(
        Mean_Composite_Rank
      )
  ) %>%

  ungroup()


write_csv(
  sex_consensus,
  file.path(
    table_dir,
    "Bulk_cross_dataset_consensus_by_sex.csv"
  )
)


# ==============================================================================
# 17. OVERALL BULK CONSENSUS
#
# Integrates:
#   Female GTEx
#   Male GTEx
#   Female HPA
#   Male HPA
#
# Mean rank is used only across analyses in which that pair was estimable.
# Coverage is retained explicitly.
# ==============================================================================

overall_consensus <- similarity_results %>%

  filter(
    Estimable,
    !is.na(
      Composite_Rank
    )
  ) %>%

  group_by(
    Human_Gene,
    Fly_Gene
  ) %>%

  summarise(

    N_Analyses_Available =
      n(),

    N_Datasets =
      n_distinct(
        Dataset
      ),

    N_Sexes =
      n_distinct(
        Sex
      ),

    Mean_Composite_Rank =
      mean(
        Composite_Rank,
        na.rm = TRUE
      ),

    Median_Composite_Rank =
      median(
        Composite_Rank,
        na.rm = TRUE
      ),

    Mean_Composite_Similarity =
      mean(
        Composite_Similarity,
        na.rm = TRUE
      ),

    Mean_Spearman =
      mean(
        Spearman_Rho,
        na.rm = TRUE
      ),

    Mean_Kendall =
      mean(
        Kendall_Tau,
        na.rm = TRUE
      ),

    Mean_Cosine =
      mean(
        Cosine_Similarity,
        na.rm = TRUE
      ),

    Mean_Jaccard =
      mean(
        Top3_Jaccard,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%

  group_by(
    Human_Gene
  ) %>%

  mutate(

    Overall_Bulk_Rank =
      min_rank(
        Mean_Composite_Rank
      )
  ) %>%

  ungroup() %>%

  arrange(
    Human_Gene,
    Overall_Bulk_Rank,
    desc(
      N_Analyses_Available
    )
  )


write_csv(
  overall_consensus,
  file.path(
    table_dir,
    "Bulk_overall_consensus_all_SERPINs.csv"
  )
)


# ==============================================================================
# 18. TOP MATCH FOR EACH HUMAN SERPIN
# ==============================================================================

top_matches <- overall_consensus %>%

  group_by(
    Human_Gene
  ) %>%

  filter(
    Overall_Bulk_Rank ==
      min(
        Overall_Bulk_Rank,
        na.rm = TRUE
      )
  ) %>%

  ungroup()


write_csv(
  top_matches,
  file.path(
    table_dir,
    "Bulk_top_fly_match_for_each_human_SERPIN.csv"
  )
)


cat("\n============================================\n")
cat("TOP BULK-EXPRESSION MATCH FOR EACH HUMAN SERPIN\n")
cat("============================================\n")


print(
  top_matches,
  n = Inf
)


# ==============================================================================
# 19. SERPINE1-SPECIFIC SUMMARY
# ==============================================================================

serpine1_summary <- overall_consensus %>%

  filter(
    Human_Gene ==
      "SERPINE1"
  ) %>%

  arrange(
    Overall_Bulk_Rank
  )


write_csv(
  serpine1_summary,
  file.path(
    table_dir,
    "SERPINE1_bulk_multimetric_consensus.csv"
  )
)


cat("\n============================================\n")
cat("SERPINE1 MULTI-METRIC CONSENSUS\n")
cat("============================================\n")


print(
  serpine1_summary,
  n = Inf
)


# ==============================================================================
# 20. Spn42Dd ACROSS ALL HUMAN SERPINS
# ==============================================================================

spn42dd_summary <- overall_consensus %>%

  filter(
    Fly_Gene ==
      "Spn42Dd"
  ) %>%

  arrange(
    Overall_Bulk_Rank
  )


write_csv(
  spn42dd_summary,
  file.path(
    table_dir,
    "Spn42Dd_similarity_to_all_human_SERPINs.csv"
  )
)


cat("\n============================================\n")
cat("Spn42Dd ACROSS HUMAN SERPINs\n")
cat("============================================\n")


print(
  spn42dd_summary,
  n = Inf
)


# ==============================================================================
# 21. FIGURE 1
# DATASET- AND SEX-SPECIFIC COMPOSITE SIMILARITY
#
# One facet per human SERPIN.
# Each point represents one fly serpin.
# GTEx and HPA are shown separately.
#
# Generates separate Female and Male figures.
# ==============================================================================

make_dataset_similarity_plot <- function(
    dat,
    sex_value
) {

  plot_data <- dat %>%

    filter(
      Sex == .env$sex_value,
      Estimable
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
        )
    )


  ggplot(
    plot_data,
    aes(
      x =
        Composite_Similarity,
      y =
        Fly_Gene,
      shape =
        Dataset
    )
  ) +

    geom_point(
      size = 3,
      position =
        position_dodge(
          width = 0.45
        )
    ) +

    facet_wrap(
      ~ Human_Gene,
      ncol = 2
    ) +

    scale_x_continuous(
      limits = c(
        0,
        1
      ),
      breaks = seq(
        0,
        1,
        by = 0.2
      )
    ) +

    labs(

      title =
        paste0(
          sex_value,
          " cross-species bulk tissue-expression similarity"
        ),

      subtitle =
        "FlyAtlas2 compared with GTEx and HPA",

      x =
        "Multi-metric tissue-profile similarity",

      y =
        "Drosophila serpin",

      shape =
        "Human dataset",

      caption =
        paste0(
          "Similarity integrates Spearman rank correlation, Kendall rank ",
          "concordance, cosine similarity of relative system-expression ",
          "profiles, and top-3 functional-system overlap. ",
          "Higher values indicate more similar tissue-expression patterns. ",
          "The score is descriptive and is not a statistical significance test."
        )
    ) +

    theme_classic(
      base_size = 10
    ) +

    theme(

      axis.text.x =
        element_text(
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
        ),

      legend.title =
        element_text(
          size = 9
        )
    )
}


# Female -----------------------------------------------------------------------

female_dataset_plot <- make_dataset_similarity_plot(
  similarity_results,
  "Female"
)


female_file <- file.path(
  figure_dir,
  "Female_bulk_multimetric_similarity_all_SERPINs.png"
)


grDevices::png(
  filename =
    female_file,
  width =
    10,
  height =
    8.5,
  units =
    "in",
  res =
    600
)


print(
  female_dataset_plot
)


grDevices::dev.off()


# Male -------------------------------------------------------------------------

male_dataset_plot <- make_dataset_similarity_plot(
  similarity_results,
  "Male"
)


male_file <- file.path(
  figure_dir,
  "Male_bulk_multimetric_similarity_all_SERPINs.png"
)


grDevices::png(
  filename =
    male_file,
  width =
    10,
  height =
    8.5,
  units =
    "in",
  res =
    600
)


print(
  male_dataset_plot
)


grDevices::dev.off()


# ==============================================================================
# 22. FIGURE 2
# CONSENSUS RANK ACROSS GTEx AND HPA, SEPARATE BY SEX
#
# Rank 1 = best fly match for that human SERPIN.
# ==============================================================================

make_sex_consensus_plot <- function(
    dat,
    sex_value
) {

  plot_data <- dat %>%

    filter(
      Sex ==
        .env$sex_value
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
        )
    )


  ggplot(
    plot_data,
    aes(
      x =
        Mean_Composite_Rank,
      y =
        Fly_Gene
    )
  ) +

    geom_point(
      size = 3
    ) +

    facet_wrap(
      ~ Human_Gene,
      ncol = 2
    ) +

    scale_x_continuous(
      limits = c(
        1,
        6
      ),
      breaks = 1:6
    ) +

    labs(

      title =
        paste0(
          sex_value,
          " consensus tissue-expression ranking"
        ),

      subtitle =
        "Consensus across GTEx and HPA bulk tissue analyses",

      x =
        "Mean fly-gene rank (1 = strongest match)",

      y =
        "Drosophila serpin",

      caption =
        paste0(
          "Ranks are based on the descriptive multi-metric similarity score ",
          "within each human dataset and are averaged across GTEx and HPA. ",
          "Lower mean rank indicates a more consistently similar tissue profile."
        )
    ) +

    theme_classic(
      base_size = 10
    ) +

    theme(

      axis.text =
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
        )
    )
}


female_consensus_plot <- make_sex_consensus_plot(
  sex_consensus,
  "Female"
)


female_consensus_file <- file.path(
  figure_dir,
  "Female_GTEx_HPA_consensus_rank_all_SERPINs.png"
)


grDevices::png(
  filename =
    female_consensus_file,
  width =
    10,
  height =
    8.5,
  units =
    "in",
  res =
    600
)


print(
  female_consensus_plot
)


grDevices::dev.off()


male_consensus_plot <- make_sex_consensus_plot(
  sex_consensus,
  "Male"
)


male_consensus_file <- file.path(
  figure_dir,
  "Male_GTEx_HPA_consensus_rank_all_SERPINs.png"
)


grDevices::png(
  filename =
    male_consensus_file,
  width =
    10,
  height =
    8.5,
  units =
    "in",
  res =
    600
)


print(
  male_consensus_plot
)


grDevices::dev.off()


# ==============================================================================
# 23. FIGURE 3
# OVERALL BULK CONSENSUS
#
# Integrates:
#   GTEx female
#   GTEx male
#   HPA female
#   HPA male
#
# This is the broad bulk tissue-expression ranking.
# ==============================================================================

overall_plot_data <- overall_consensus %>%

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
      )
  )


overall_consensus_plot <- ggplot(
  overall_plot_data,
  aes(
    x =
      Mean_Composite_Rank,
    y =
      Fly_Gene
  )
) +

  geom_point(
    aes(
      size =
        N_Analyses_Available
    )
  ) +

  facet_wrap(
    ~ Human_Gene,
    ncol = 2
  ) +

  scale_x_continuous(
    limits = c(
      1,
      6
    ),
    breaks = 1:6
  ) +

  scale_size_continuous(
    range = c(
      2.5,
      4.5
    ),
    breaks = 1:4
  ) +

  labs(

    title =
      "Overall cross-species bulk tissue-expression consensus",

    subtitle =
      "FlyAtlas2 similarity to six human SERPIN genes across GTEx, HPA, and sex",

    x =
      "Mean fly-gene rank (1 = strongest match)",

    y =
      "Drosophila serpin",

    size =
      "Analyses\navailable",

    caption =
      paste0(
        "Mean rank integrates the available female and male GTEx and HPA ",
        "bulk tissue-profile comparisons. Point size indicates how many of ",
        "the four possible analyses were estimable. Lower mean rank indicates ",
        "a more consistently similar cross-species tissue-expression pattern."
      )
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(

    axis.text =
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
        size = 8
      ),

    legend.title =
      element_text(
        size = 9
      )
  )


overall_file <- file.path(
  figure_dir,
  "Overall_bulk_consensus_rank_all_SERPINs.png"
)


grDevices::png(
  filename =
    overall_file,
  width =
    10,
  height =
    8.5,
  units =
    "in",
  res =
    600
)


print(
  overall_consensus_plot
)


grDevices::dev.off()


# ==============================================================================
# 24. FIGURE 4
# SERPINE1-FOCUSED OVERALL RANKING
# ==============================================================================

serpine1_plot_data <- overall_consensus %>%

  filter(
    Human_Gene ==
      "SERPINE1"
  ) %>%

  mutate(

    Fly_Gene =
      reorder(
        Fly_Gene,
        Mean_Composite_Rank
      )
  )


serpine1_rank_plot <- ggplot(
  serpine1_plot_data,
  aes(
    x =
      Mean_Composite_Rank,
    y =
      Fly_Gene
  )
) +

  geom_point(
    aes(
      size =
        N_Analyses_Available
    ),
    shape = 16
  ) +

  scale_x_continuous(
    limits = c(
      1,
      6
    ),
    breaks = 1:6
  ) +

  scale_size_continuous(
    range = c(
      3,
      5
    ),
    breaks = 1:4
  ) +

  labs(

    title =
      "Bulk tissue-expression similarity to human SERPINE1",

    subtitle =
      "Consensus across GTEx, HPA, female, and male comparisons",

    x =
      "Mean fly-gene rank (1 = strongest match)",

    y =
      "Drosophila serpin",

    size =
      "Analyses\navailable",

    caption =
      paste0(
        "Ranking integrates Spearman, Kendall, relative-profile cosine ",
        "similarity, and top-system overlap within each bulk expression ",
        "analysis before combining ranks across datasets and sex."
      )
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(

    axis.text =
      element_text(
        size = 9
      ),

    axis.title =
      element_text(
        size = 10
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
      )
  )


serpine1_file <- file.path(
  figure_dir,
  "SERPINE1_overall_bulk_consensus_rank.png"
)


grDevices::png(
  filename =
    serpine1_file,
  width =
    7.5,
  height =
    5.5,
  units =
    "in",
  res =
    600
)


print(
  serpine1_rank_plot
)


grDevices::dev.off()


# ==============================================================================
# 25. VERIFY OUTPUT FILES
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
cat("FIGURES CREATED\n")
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
# 26. PRINT ALL OVERALL RANKINGS
# ==============================================================================

cat("\n============================================\n")
cat("OVERALL BULK CONSENSUS RANKINGS\n")
cat("============================================\n")


overall_consensus %>%

  select(
    Human_Gene,
    Fly_Gene,
    N_Analyses_Available,
    Mean_Composite_Rank,
    Mean_Composite_Similarity,
    Mean_Spearman,
    Mean_Kendall,
    Mean_Cosine,
    Mean_Jaccard,
    Overall_Bulk_Rank
  ) %>%

  arrange(
    Human_Gene,
    Overall_Bulk_Rank
  ) %>%

  print(
    n = Inf
  )


# ==============================================================================
# 27. FINAL MESSAGE
# ==============================================================================

cat("\n============================================\n")
cat("MULTI-METRIC BULK TISSUE ANALYSIS COMPLETE\n")
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