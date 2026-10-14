# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 10_IntegratedEvidenceRanking.R
# Run from the repository root.
# ============================================================
# 10_IntegratedEvidenceRanking_REVISED.R
# Comparative functional genomics final integrated evidence analysis
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

analysis_root <- "."

legacy_project_root <- "."
legacy_analysis_dir <- legacy_project_root
legacy_analysis_processed <- file.path(legacy_analysis_dir, "Data", "Processed")

default_output_dir <- file.path(
  analysis_root,
  "Integrated_Evidence",
  "Revised_Final_Integration"
)

# An optional output override also works in noninteractive sessions.
output_dir <- Sys.getenv("SERPINE1_INTEGRATED_OUTPUT", unset = default_output_dir)

tables_dir <- file.path(output_dir, "Tables")
figures_dir <- file.path(output_dir, "Figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Candidate genes
# ------------------------------------------------------------

selected_serpins <- c(
  "Spn42Da", "Spn42Db", "Spn42Dc",
  "Spn42Dd", "Spn28B", "Spn28F"
)

# Isoforms analyzed separately in the sequence-alignment analysis.
# Spn42Da-B is NOT combined with Spn42Da.
sequence_only_isoforms <- c("Spn42Da-B")

alignment_candidates <- c(
  selected_serpins,
  sequence_only_isoforms
)

target_human_gene <- "SERPINE1"

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

mean_or_na <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

clamp_01 <- function(x) pmin(1, pmax(0, as.numeric(x)))

proportion_01 <- function(x) {
  x <- as.numeric(x)
  out <- ifelse(is.na(x), NA_real_, ifelse(abs(x) > 1, x / 100, x))
  clamp_01(out)
}

correlation_to_01 <- function(x) clamp_01((as.numeric(x) + 1) / 2)

distance_to_similarity <- function(x) {
  x <- as.numeric(x)
  ifelse(is.na(x), NA_real_, exp(-pmax(x, 0)))
}

relative_scale_01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  if (isTRUE(all.equal(x_min, x_max))) return(ifelse(is.na(x), NA_real_, 0.5))
  (x - x_min) / (x_max - x_min)
}

check_required_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      "\nAvailable columns:\n", paste(names(data), collapse = ", ")
    )
  }
}

check_unique_gene <- function(data, gene_col, label) {
  dup <- data[[gene_col]][duplicated(data[[gene_col]])]
  if (length(dup) > 0) {
    stop(label, " contains duplicate candidate genes: ", paste(unique(dup), collapse = ", "))
  }
}

find_exact_file <- function(search_root, filename) {
  hits <- list.files(
    search_root,
    pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0) stop("Could not find file: ", filename, "\nSearch root: ", search_root)
  if (length(hits) > 1) {
    stop(
      "More than one copy of ", filename, " was found.\nMatches:\n",
      paste(hits, collapse = "\n")
    )
  }
  hits[[1]]
}

find_column <- function(data, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    if (required) {
      stop(
        "Could not identify ", label,
        ".\nTried: ", paste(candidates, collapse = ", "),
        "\nAvailable columns:\n", paste(names(data), collapse = ", ")
      )
    }
    return(NA_character_)
  }
  hit[[1]]
}

# Read an already-aligned FASTA file into a named character vector.
read_aligned_fasta <- function(path) {
  x <- readLines(path, warn = FALSE)
  header_idx <- grep("^>", x)

  if (length(header_idx) == 0) {
    stop("No FASTA headers were detected in: ", path)
  }

  end_idx <- c(header_idx[-1] - 1, length(x))
  seqs <- map2_chr(
    header_idx,
    end_idx,
    ~paste0(x[(.x + 1):.y], collapse = "")
  )

  names(seqs) <- sub("^>", "", x[header_idx])
  seqs <- toupper(gsub("\\s+", "", seqs))

  if (length(unique(nchar(seqs))) != 1) {
    stop("Sequences do not all have the same aligned length: ", path)
  }

  seqs
}

find_fasta_sequence <- function(seqs, patterns, label) {
  nm <- names(seqs)
  first_token <- sub("\\s.*$", "", nm)

  # 1. Prefer an exact full-header match.
  for (p in patterns) {
    exact_full <- nm[tolower(nm) == tolower(p)]

    if (length(exact_full) == 1) {
      return(exact_full[[1]])
    }

    if (length(exact_full) > 1) {
      stop(
        "More than one exact FASTA header matched ", label, ":\n",
        paste(exact_full, collapse = "\n")
      )
    }
  }

  # 2. Then prefer an exact first-token match.
  # This allows headers such as "Spn42Da description..." while
  # preventing "Spn42Da" from also matching "Spn42Da-B".
  for (p in patterns) {
    exact_token <- nm[tolower(first_token) == tolower(p)]

    if (length(exact_token) == 1) {
      return(exact_token[[1]])
    }

    if (length(exact_token) > 1) {
      stop(
        "More than one FASTA header had the exact first token ",
        p, " for ", label, ":\n",
        paste(exact_token, collapse = "\n")
      )
    }
  }

  # 3. Only if no exact match exists, use a unique partial match.
  hit <- unique(unlist(lapply(patterns, function(p) {
    grep(p, nm, ignore.case = TRUE, value = TRUE)
  })))

  if (length(hit) == 0) {
    stop(
      "Could not identify ", label, " in FASTA headers.\n",
      "Tried patterns: ", paste(patterns, collapse = ", "),
      "\nAvailable headers:\n", paste(nm, collapse = "\n")
    )
  }

  if (length(hit) > 1) {
    stop(
      "More than one FASTA header matched ", label, ":\n",
      paste(hit, collapse = "\n"),
      "\nPlease specify the desired FASTA header explicitly."
    )
  }

  hit[[1]]
}

aligned_percent_identity <- function(seq1, seq2) {
  a <- strsplit(seq1, "", fixed = TRUE)[[1]]
  b <- strsplit(seq2, "", fixed = TRUE)[[1]]

  if (length(a) != length(b)) {
    stop("Aligned sequences have unequal lengths.")
  }

  informative <- !(a == "-" & b == "-")
  denom <- sum(informative)
  if (denom == 0) return(NA_real_)

  100 * sum(a[informative] == b[informative]) / denom
}

# ------------------------------------------------------------
# 4. Configure input files
# ------------------------------------------------------------

# Set these variables to the original input paths. Interactive selection
# remains available when no configured path is supplied.
choose_input_file <- function(label, variable) {
  chosen <- Sys.getenv(variable, unset = "")
  if (!nzchar(chosen)) {
    if (!interactive()) {
      stop("Missing input configuration: set ", variable,
           " to the path for ", label, ".", call. = FALSE)
    }
    message("Select file for: ", label)
    chosen <- file.choose()
  }
  if (!file.exists(chosen) || dir.exists(chosen)) {
    stop("Input file not found for ", variable, ": ", chosen, call. = FALSE)
  }
  chosen
}

diopt_file <- choose_input_file(
  "DIOPT pairwise predictions CSV", "SERPINE1_DIOPT_INPUT")
full_length_alignment_file <- choose_input_file(
  "Full-length human-fly aligned FASTA", "SERPINE1_FULL_ALIGNMENT_INPUT")
rcl_alignment_file <- choose_input_file(
  "RCL human-fly aligned FASTA", "SERPINE1_RCL_ALIGNMENT_INPUT")
bulk_file <- choose_input_file(
  "Bulk cross-dataset consensus by sex CSV", "SERPINE1_BULK_INPUT")
single_cell_file <- choose_input_file(
  "FCA-GTEx snRNA Spearman by sex CSV", "SERPINE1_CELL_INPUT")
aging_file <- choose_input_file(
  "Cross-species aging overall similarity v2 CSV", "SERPINE1_AGING_INPUT")

input_files <- c(
  DIOPT = diopt_file,
  Full_Length_Alignment = full_length_alignment_file,
  RCL_Alignment = rcl_alignment_file,
  Bulk = bulk_file,
  Single_Cell = single_cell_file,
  Aging = aging_file
)

message("\nInput files selected successfully:")
walk2(
  names(input_files),
  input_files,
  ~message("  ", .x, ": ", .y)
)

# ------------------------------------------------------------
# 5. Orthology domain
# ------------------------------------------------------------

diopt <- read.csv(diopt_file, stringsAsFactors = FALSE)

check_required_columns(
  diopt,
  c(
    "Drosophila_Gene",
    "Human_Gene",
    "DIOPT_Score",
    "Weighted_Score",
    "Best_Score_Binary",
    "Reverse_Best_Binary",
    "Supporting_DB_Count"
  ),
  "DIOPT pairwise table"
)

orthology <- diopt %>%
  filter(
    Drosophila_Gene %in% selected_serpins,
    Human_Gene == target_human_gene
  ) %>%
  transmute(
    Fly_Gene = Drosophila_Gene,
    Human_Gene,

    DIOPT_Score =
      as.numeric(DIOPT_Score),

    Weighted_Score =
      as.numeric(Weighted_Score),

    Best_Score_Binary =
      as.numeric(Best_Score_Binary),

    Reverse_Best_Binary =
      as.numeric(Reverse_Best_Binary),

    Supporting_DB_Count =
      as.numeric(Supporting_DB_Count)
  )

check_unique_gene(
  orthology,
  "Fly_Gene",
  "DIOPT SERPINE1 pairwise subset"
)

if (nrow(orthology) != length(selected_serpins)) {
  warning(
    "Expected ",
    length(selected_serpins),
    " SERPINE1 DIOPT rows, but found ",
    nrow(orthology),
    "."
  )
}

orthology <- orthology %>%
  mutate(
    # DIOPT and weighted scores do not share the same natural 0-1
    # scale as the expression-domain similarities. They are therefore
    # normalized only relative to the six candidate fly serpins.
    DIOPT_Relative =
      relative_scale_01(DIOPT_Score),

    Weighted_DIOPT_Relative =
      relative_scale_01(Weighted_Score),

    Supporting_DB_Relative =
      relative_scale_01(Supporting_DB_Count)
  ) %>%
  rowwise() %>%
  mutate(
    Orthology_Score = mean_or_na(
      c(
        DIOPT_Relative,
        Weighted_DIOPT_Relative,
        Supporting_DB_Relative,
        Best_Score_Binary,
        Reverse_Best_Binary
      )
    )
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 6. Sequence-alignment domain
# ------------------------------------------------------------

# The supplied files are aligned FASTA files rather than precomputed
# statistics tables. Full-length identity is calculated directly against
# human SERPINE1/PAI-1. RCL identity uses validated step 03 scores
# with a denominator of 25 SERPINE1 residue positions. A phylogenetic-distance term is NOT included here because
# no validated phylogenetic-distance table was supplied with these paths.

full_alignment <- read_aligned_fasta(full_length_alignment_file)
rcl_alignment <- read_aligned_fasta(rcl_alignment_file)

full_serpine1_header <- find_fasta_sequence(
  full_alignment,
  c("SERPINE1", "PAI-1", "PAI1"),
  "human SERPINE1/PAI-1 full-length sequence"
)

rcl_serpine1_header <- find_fasta_sequence(
  rcl_alignment,
  c("SERPINE1", "PAI-1", "PAI1"),
  "human SERPINE1/PAI-1 RCL sequence"
)

# Use the validated step 03 SERPINE1-referenced 25-position identities.
# Candidate gaps count as nonmatches; SERPINE1-gap columns are excluded.
rcl_scores_file <- Sys.getenv("SERPINE1_RCL_SCORES_INPUT", unset =
  "Sequence/RCL/Processed/03_SERPINE1_referenced_RCL_region_scores.csv")
if (!file.exists(rcl_scores_file)) stop("Run step 03 first or set SERPINE1_RCL_SCORES_INPUT.")
rcl_scores <- read.csv(rcl_scores_file, stringsAsFactors = FALSE)
check_required_columns(rcl_scores,
  c("Gene", "Reference_Residues", "RCL_Region_Identity_Percent"), "Validated RCL scores")
validated_rcl_identity <- function(gene) {
  key <- switch(gene, "Spn42Da" = "Spn42Da_Q8T0M5",
    "Spn42Da-B" = "Spn42Da_B_A1Z6R1", gene)
  row <- rcl_scores[rcl_scores$Gene == key, , drop = FALSE]
  if (nrow(row) != 1L || row$Reference_Residues != 25L ||
      !is.finite(row$RCL_Region_Identity_Percent) ||
      row$RCL_Region_Identity_Percent < 0 || row$RCL_Region_Identity_Percent > 100) {
    stop("Expected one validated 25-position RCL score for ", gene)
  }
  row$RCL_Region_Identity_Percent
}

sequence_alignment_all <- map_dfr(
  alignment_candidates,
  function(gene) {
    full_fly_header <- find_fasta_sequence(
      full_alignment,
      c(gene),
      paste0(gene, " full-length sequence")
    )

    rcl_fly_header <- find_fasta_sequence(
      rcl_alignment,
      c(gene),
      paste0(gene, " RCL sequence")
    )

    tibble(
      Fly_Gene = gene,
      Full_Length_Header = full_fly_header,
      RCL_Header = rcl_fly_header,
      Full_Length_Identity = aligned_percent_identity(
        full_alignment[[full_fly_header]],
        full_alignment[[full_serpine1_header]]
      ),
      RCL_Identity = validated_rcl_identity(gene)
    )
  }
) %>%
  mutate(
    Full_Length_Identity_01 = proportion_01(Full_Length_Identity),
    RCL_Identity_01 = proportion_01(RCL_Identity)
  ) %>%
  rowwise() %>%
  mutate(
    Sequence_Alignment_Score = mean_or_na(
      c(Full_Length_Identity_01, RCL_Identity_01)
    )
  ) %>%
  ungroup()

# Export the isoform-resolved alignment results.
# Spn42Da and Spn42Da-B remain separate rows.
write.csv(
  sequence_alignment_all,
  file.path(
    tables_dir,
    "SERPINE1_SequenceAlignment_IsoformResolved.csv"
  ),
  row.names = FALSE
)

# The integrated evidence ranking remains gene-level because the
# orthology, bulk-expression, single-cell, and aging datasets do not
# distinguish Spn42Da from the Spn42Da-B isoform. Therefore Spn42Da-B
# is retained in the sequence analysis but is not merged into Spn42Da
# and is not given a separate integrated score.
sequence_alignment <- sequence_alignment_all %>%
  filter(Fly_Gene %in% selected_serpins)

# ------------------------------------------------------------
# 7. Revised bulk-expression domain
# ------------------------------------------------------------

bulk_raw <- read.csv(bulk_file, stringsAsFactors = FALSE)

check_required_columns(
  bulk_raw,
  c(
    "Sex",
    "Human_Gene",
    "Fly_Gene",
    "N_Human_Datasets",
    "Mean_Composite_Rank",
    "Median_Composite_Rank",
    "Mean_Composite_Similarity",
    "Sex_Consensus_Rank"
  ),
  "Bulk cross-dataset consensus table"
)

bulk_long <- bulk_raw %>%
  filter(
    Fly_Gene %in% selected_serpins,
    Human_Gene == target_human_gene
  ) %>%
  transmute(
    Sex = case_when(
      str_to_lower(as.character(Sex)) %in% c("female", "f") ~ "Female",
      str_to_lower(as.character(Sex)) %in% c("male", "m") ~ "Male",
      TRUE ~ as.character(Sex)
    ),
    Human_Gene,
    Fly_Gene,
    N_Human_Datasets = as.numeric(N_Human_Datasets),
    Mean_Composite_Rank = as.numeric(Mean_Composite_Rank),
    Median_Composite_Rank = as.numeric(Median_Composite_Rank),
    Mean_Composite_Similarity_Raw =
      as.numeric(Mean_Composite_Similarity),
    Sex_Consensus_Rank = as.numeric(Sex_Consensus_Rank),
    Bulk_Expression_Score =
      proportion_01(Mean_Composite_Similarity)
  ) %>%
  filter(Sex %in% c("Female", "Male"))

if (nrow(bulk_long) == 0) {
  stop("No sex-specific SERPINE1 bulk-expression rows were found.")
}

bulk_by_gene_sex <- bulk_long %>%
  group_by(Fly_Gene, Sex) %>%
  summarise(
    Bulk_Expression_Score =
      mean_or_na(Bulk_Expression_Score),
    Bulk_N_Human_Datasets =
      mean_or_na(N_Human_Datasets),
    Bulk_Mean_Composite_Rank =
      mean_or_na(Mean_Composite_Rank),
    Bulk_Median_Composite_Rank =
      mean_or_na(Median_Composite_Rank),
    Bulk_Sex_Consensus_Rank =
      mean_or_na(Sex_Consensus_Rank),
    .groups = "drop"
  )

if (any(duplicated(bulk_by_gene_sex[c("Fly_Gene", "Sex")]))) {
  stop("Bulk evidence contains duplicate Fly_Gene x Sex rows.")
}

# ------------------------------------------------------------
# 8. Revised FCA <-> GTEx snRNA cellular-expression domain
# ------------------------------------------------------------

single_cell_raw <- read.csv(
  single_cell_file,
  stringsAsFactors = FALSE
)

sc_fly_col <- find_column(
  single_cell_raw,
  c("Fly_Gene", "FlyGene"),
  "single-cell fly-gene column"
)

sc_human_col <- find_column(
  single_cell_raw,
  c("Human_Gene", "HumanGene"),
  "single-cell human-gene column"
)

sc_sex_col <- find_column(
  single_cell_raw,
  c("Sex", "sex"),
  "single-cell sex column"
)

sc_rho_col <- find_column(
  single_cell_raw,
  c("Spearman_Rho", "Spearman", "rho", "Rho"),
  "single-cell Spearman column"
)

sc_n_col <- find_column(
  single_cell_raw,
  c("N_Systems", "N", "n_systems"),
  "single-cell matched-system count",
  required = FALSE
)

single_cell_by_gene_sex <- single_cell_raw %>%
  transmute(
    Sex = case_when(
      str_to_lower(as.character(.data[[sc_sex_col]])) %in%
        c("female", "f") ~ "Female",
      str_to_lower(as.character(.data[[sc_sex_col]])) %in%
        c("male", "m") ~ "Male",
      TRUE ~ as.character(.data[[sc_sex_col]])
    ),
    Fly_Gene = .data[[sc_fly_col]],
    Human_Gene = .data[[sc_human_col]],
    Cellular_Spearman_Rho =
      as.numeric(.data[[sc_rho_col]]),
    Cellular_N_Systems =
      if (!is.na(sc_n_col)) {
        as.numeric(.data[[sc_n_col]])
      } else {
        NA_real_
      },
    Cellular_Expression_Score =
      correlation_to_01(.data[[sc_rho_col]])
  ) %>%
  filter(
    Fly_Gene %in% selected_serpins,
    Human_Gene == target_human_gene,
    Sex %in% c("Female", "Male")
  )

if (any(duplicated(single_cell_by_gene_sex[c("Fly_Gene", "Sex")]))) {
  stop("Single-cell evidence contains duplicate Fly_Gene x Sex rows.")
}

# ------------------------------------------------------------
# 9. Revised aging-v2 domain
# ------------------------------------------------------------

aging_raw <- read.csv(aging_file, stringsAsFactors = FALSE)

check_required_columns(
  aging_raw,
  c(
    "Fly_Gene", "Human_Gene", "Spearman_Rho", "Pearson_R_Overall",
    "Cosine_Overall", "Z_RMSE_Overall", "Z_MAE_Overall",
    "Direction_Concordance", "Integrated_Aging_Similarity"
  ),
  "Aging-v2 overall similarity table"
)

aging_by_gene <- aging_raw %>%
  filter(Fly_Gene %in% selected_serpins, Human_Gene == target_human_gene) %>%
  transmute(
    Fly_Gene,
    Aging_Spearman_Rho = as.numeric(Spearman_Rho),
    Aging_Pearson_R = as.numeric(Pearson_R_Overall),
    Aging_Cosine = as.numeric(Cosine_Overall),
    Aging_Z_RMSE = as.numeric(Z_RMSE_Overall),
    Aging_Z_MAE = as.numeric(Z_MAE_Overall),
    Aging_Direction_Concordance = as.numeric(Direction_Concordance),
    Aging_Integrated_Raw = as.numeric(Integrated_Aging_Similarity),
    Aging_Score = proportion_01(Integrated_Aging_Similarity)
  )

check_unique_gene(aging_by_gene, "Fly_Gene", "Aging-v2 SERPINE1 subset")

# ------------------------------------------------------------
# 10. Build sex-specific master evidence table
# ------------------------------------------------------------

sex_backbone <- tidyr::expand_grid(
  Fly_Gene = selected_serpins,
  Sex = c("Female", "Male")
) %>%
  mutate(Human_Gene = target_human_gene)

master <- sex_backbone %>%
  left_join(orthology, by = "Fly_Gene") %>%
  left_join(sequence_alignment, by = "Fly_Gene") %>%
  left_join(
    bulk_by_gene_sex,
    by = c("Fly_Gene", "Sex")
  ) %>%
  left_join(
    single_cell_by_gene_sex,
    by = c("Fly_Gene", "Sex")
  ) %>%
  left_join(
    aging_by_gene,
    by = "Fly_Gene"
  )

if (nrow(master) != length(selected_serpins) * 2) {
  stop(
    "Expected 12 rows (6 fly genes x 2 sexes) in the master table, but found ",
    nrow(master),
    "."
  )
}

# ------------------------------------------------------------
# 11. Equal-domain integration, separately by sex
# ------------------------------------------------------------

domain_columns <- c(
  "Orthology_Score",
  "Sequence_Alignment_Score",
  "Bulk_Expression_Score",
  "Cellular_Expression_Score",
  "Aging_Score"
)

final <- master %>%
  rowwise() %>%
  mutate(
    N_Available_Domains =
      sum(!is.na(c_across(all_of(domain_columns)))),
    Integrated_Evidence_Score =
      mean_or_na(c_across(all_of(domain_columns)))
  ) %>%
  ungroup() %>%
  group_by(Sex) %>%
  mutate(
    Integrated_Rank = rank(
      -Integrated_Evidence_Score,
      ties.method = "min",
      na.last = "keep"
    )
  ) %>%
  ungroup() %>%
  arrange(Sex, Integrated_Rank, Fly_Gene)

# ------------------------------------------------------------
# 12. Leave-one-domain-out sensitivity analysis by sex
# ------------------------------------------------------------

loo_results <- map_dfr(
  domain_columns,
  function(omit_domain) {

    keep_domains <- setdiff(domain_columns, omit_domain)

    tmp <- final %>%
      select(Fly_Gene, Sex, all_of(keep_domains))

    tmp$Leave_One_Domain_Out_Score <- apply(
      tmp[, keep_domains, drop = FALSE],
      1,
      mean_or_na
    )

    tmp %>%
      select(
        Fly_Gene,
        Sex,
        Leave_One_Domain_Out_Score
      ) %>%
      mutate(
        Omitted_Domain = omit_domain
      ) %>%
      group_by(Sex) %>%
      mutate(
        Leave_One_Domain_Out_Rank = rank(
          -Leave_One_Domain_Out_Score,
          ties.method = "min",
          na.last = "keep"
        )
      ) %>%
      ungroup()
  }
)

rank_stability <- loo_results %>%
  group_by(Fly_Gene, Sex) %>%
  summarise(
    Best_LOO_Rank =
      min(Leave_One_Domain_Out_Rank, na.rm = TRUE),
    Worst_LOO_Rank =
      max(Leave_One_Domain_Out_Rank, na.rm = TRUE),
    Mean_LOO_Rank =
      mean(Leave_One_Domain_Out_Rank, na.rm = TRUE),
    SD_LOO_Rank =
      sd(Leave_One_Domain_Out_Rank, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 13. Concise study evidence matrix
# ------------------------------------------------------------

final_matrix <- final %>%
  transmute(
    Sex,
    Rank = Integrated_Rank,
    Fly_Gene,
    Human_Gene,
    Orthology =
      round(Orthology_Score, 3),
    Sequence_Alignment =
      round(Sequence_Alignment_Score, 3),
    Bulk_Expression =
      round(Bulk_Expression_Score, 3),
    Cellular_Expression =
      round(Cellular_Expression_Score, 3),
    Aging =
      round(Aging_Score, 3),
    Available_Domains =
      N_Available_Domains,
    Integrated_Evidence_Score =
      round(Integrated_Evidence_Score, 3)
  )

print(final_matrix)

# ------------------------------------------------------------
# 14. Export tables
# ------------------------------------------------------------

write.csv(
  master,
  file.path(
    tables_dir,
    "SERPINE1_Master_Evidence_BySex_All_Metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  final,
  file.path(
    tables_dir,
    "SERPINE1_Integrated_Evidence_BySex_Full.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_matrix,
  file.path(
    tables_dir,
    "SERPINE1_Integrated_Evidence_BySex_Final_Ranking.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_matrix %>% filter(Sex == "Female"),
  file.path(
    tables_dir,
    "SERPINE1_Female_Integrated_Ranking.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_matrix %>% filter(Sex == "Male"),
  file.path(
    tables_dir,
    "SERPINE1_Male_Integrated_Ranking.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_results,
  file.path(
    tables_dir,
    "SERPINE1_LeaveOneDomainOut_BySex.csv"
  ),
  row.names = FALSE
)

write.csv(
  rank_stability,
  file.path(
    tables_dir,
    "SERPINE1_Rank_Stability_BySex.csv"
  ),
  row.names = FALSE
)

write.csv(
  bulk_long,
  file.path(
    tables_dir,
    "SERPINE1_Bulk_Evidence_Used_BySex.csv"
  ),
  row.names = FALSE
)

write.csv(
  single_cell_by_gene_sex,
  file.path(
    tables_dir,
    "SERPINE1_Cellular_Evidence_Used_BySex.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Figure A: integrated rankings separated by sex
# ------------------------------------------------------------

p_rank <- ggplot(
  final_matrix,
  aes(
    x = Integrated_Evidence_Score,
    y = reorder(Fly_Gene, Integrated_Evidence_Score)
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%.3f", Integrated_Evidence_Score)),
    hjust = -0.15,
    size = 3
  ) +
  facet_wrap(~Sex, ncol = 2, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw(base_size = 10) +
  labs(
    title = "Integrated cross-species evidence relative to human SERPINE1",
    subtitle = "Female and male evidence integrated separately",
    x = "Integrated evidence score",
    y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_rank)

ggsave(
  file.path(
    figures_dir,
    "SERPINE1_Integrated_Evidence_Ranking_BySex.png"
  ),
  p_rank,
  width = 9,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 16. Figure B: evidence-domain profiles separated by sex
# ------------------------------------------------------------

domain_plot_data <- final %>%
  select(Fly_Gene, Sex, all_of(domain_columns)) %>%
  pivot_longer(
    cols = all_of(domain_columns),
    names_to = "Evidence_Domain",
    values_to = "Score"
  ) %>%
  mutate(
    Evidence_Domain = recode(
      Evidence_Domain,
      Orthology_Score = "Orthology",
      Sequence_Alignment_Score = "Sequence alignment",
      Bulk_Expression_Score = "Bulk expression",
      Cellular_Expression_Score = "Cellular expression",
      Aging_Score = "Aging"
    ),
    Evidence_Domain = factor(
      Evidence_Domain,
      levels = c(
        "Orthology",
        "Sequence alignment",
        "Bulk expression",
        "Cellular expression",
        "Aging"
      )
    )
  )

p_domains <- ggplot(
  domain_plot_data,
  aes(
    x = Evidence_Domain,
    y = Score,
    group = Fly_Gene
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.9) +
  facet_grid(Sex ~ Fly_Gene) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  theme_bw(base_size = 8) +
  labs(
    title = "Evidence-domain profiles for Drosophila serpin candidates",
    subtitle = "Bulk and cellular expression remain sex-specific",
    x = NULL,
    y = "Domain score"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_domains)

ggsave(
  file.path(
    figures_dir,
    "SERPINE1_Evidence_Domain_Profiles_BySex.png"
  ),
  p_domains,
  width = 12,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 17. Console summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL SERPINE1 INTEGRATED EVIDENCE RANKING BY SEX\n")
cat("============================================================\n")

cat("\nFEMALE:\n")
print(
  final_matrix %>%
    filter(Sex == "Female") %>%
    arrange(Rank)
)

cat("\nMALE:\n")
print(
  final_matrix %>%
    filter(Sex == "Male") %>%
    arrange(Rank)
)

cat("\nInterpretation note:\n")
cat(
  "Female and male integrated rankings are calculated separately. ",
  "Bulk-expression and FCA-GTEx snRNA evidence are sex-specific. ",
  "Orthology and sequence-alignment evidence are not sex-specific ",
  "and therefore have the same values in the Female and Male rankings. ",
  "The current aging reanalysis is also not sex-stratified, so the same ",
  "aging score is used in both rankings. No female/male averaging is ",
  "performed for the integrated analysis.\n",
  sep = ""
)

cat("\nOutput directory:\n", output_dir, "\n")
