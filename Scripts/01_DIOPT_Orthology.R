# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 01_DIOPT_Orthology.R
# Run from the repository root.
# ============================================================
# 01_DIOPT_Orthology.R
# DIOPT orthology summary for candidate Drosophila serpins
# ============================================================

library(dplyr)
library(readr)
library(stringr)

# Run from the repository root; analysis paths are relative to it.
# No machine-specific project-directory configuration is required.
project_dir <- "."

analysis_dir <- project_dir

data_dir <- file.path(analysis_dir, "Data")
processed_dir <- file.path(data_dir, "Processed")
tables_dir <- file.path(analysis_dir, "Tables", "Orthology")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

diopt_file <- file.path(data_dir, "orthology_mapping.csv")

# Prepared score-tier summaries retain ties and use a different schema from
# the earlier primary-ortholog mapping. Preserve this supplied representation.
prepared_summary_file <- Sys.getenv(
  "SERPINE1_DIOPT_SUMMARY_INPUT",
  unset = file.path("Orthology", "Inputs", "Table_DIOPT_candidate_serpin_orthology.csv"))
use_prepared_summary <- nzchar(Sys.getenv("SERPINE1_DIOPT_SUMMARY_INPUT")) ||
  (!file.exists(diopt_file) && file.exists(prepared_summary_file))

if (use_prepared_summary) {
  diopt_summary <- read.csv(prepared_summary_file, stringsAsFactors = FALSE)
  selected_serpins <- c("Spn42Da", "Spn42Db", "Spn42Dc", "Spn42Dd", "Spn28B", "Spn28F")
  required_prepared_columns <- c(
    "Gene", "FlyBase_ID", "Tier1_DIOPT_Score", "Tier2_DIOPT_Score", "Tier3_DIOPT_Score",
    "Tier1_Human_Orthologs", "Tier2_Human_Orthologs", "Tier3_Human_Orthologs",
    "SERPINE1_Predicted", "SERPINE1_DIOPT_Score", "SERPINE1_Weighted_Score",
    "SERPINE1_Best_Binary", "SERPINE1_Reverse_Best", "SERPINE1_Supporting_DBs")
  missing_prepared_columns <- setdiff(required_prepared_columns, names(diopt_summary))
  if (length(missing_prepared_columns)) {
    stop("Prepared DIOPT summary is missing columns: ",
         paste(missing_prepared_columns, collapse = ", "), call. = FALSE)
  }
  if (nrow(diopt_summary) != 6L || anyDuplicated(diopt_summary$Gene) ||
      !setequal(diopt_summary$Gene, selected_serpins)) {
    stop("Prepared DIOPT summary must represent all six candidates exactly once.", call. = FALSE)
  }
  print(diopt_summary)
  write.csv(diopt_summary, file.path(processed_dir, "Table_DIOPT_candidate_serpin_orthology.csv"), row.names = FALSE)
  write.csv(diopt_summary, file.path(tables_dir, "Table_DIOPT_candidate_serpin_orthology.csv"), row.names = FALSE)
  message("DIOPT QC PASSED: supplied score-tier summary retained for all six candidates.")
} else {
if (!file.exists(diopt_file)) {
  stop("DIOPT input file not found: ", diopt_file)
}

diopt <- read.csv(diopt_file, stringsAsFactors = FALSE)

if (nrow(diopt) == 0) {
  warning("DIOPT input table has zero rows. Output tables will be empty.")
}

if (ncol(diopt) == 0) {
  stop("DIOPT input table has zero columns.")
}

selected_serpins <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)

if (length(selected_serpins) == 0) {
  stop("selected_serpins is empty. Add at least one candidate gene before summarizing DIOPT orthology.")
}

required_columns <- c(
  "Drosophila_Gene",
  "FlyBase_ID",
  "Primary_Human_Ortholog",
  "Primary_Human_Ensembl_ID",
  "Primary_DIOPT_Score",
  "Primary_Weighted_Score",
  "Primary_Confidence_Rank",
  "Primary_Best_Score",
  "Primary_Reverse_Best_Score",
  "Primary_Supporting_Database_Count",
  "SERPINE1_Predicted",
  "SERPINE1_DIOPT_Score",
  "SERPINE1_Weighted_Score",
  "SERPINE1_Confidence_Rank",
  "SERPINE1_Reverse_Best_Score"
)

missing_columns <- setdiff(required_columns, names(diopt))

if (length(missing_columns) > 0) {
  stop(
    "DIOPT input table is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

numeric_columns <- c(
  "Primary_DIOPT_Score",
  "Primary_Weighted_Score",
  "Primary_Confidence_Rank",
  "Primary_Supporting_Database_Count",
  "SERPINE1_DIOPT_Score",
  "SERPINE1_Weighted_Score",
  "SERPINE1_Confidence_Rank"
)

diopt_for_summary <- diopt %>%
  mutate(
    across(all_of(numeric_columns), as.numeric)
  )

if (any(is.na(diopt_for_summary[numeric_columns]))) {
  warning(
    "One or more DIOPT score/rank/count values are NA after numeric conversion. ",
    "Missing values will be retained in the output tables."
  )
}

diopt_candidates <- diopt_for_summary %>%
  filter(Drosophila_Gene %in% selected_serpins)

missing_serpins <- setdiff(selected_serpins, unique(diopt_candidates$Drosophila_Gene))

if (length(missing_serpins) > 0) {
  warning(
    "The following selected serpin(s) were not found in the DIOPT input table: ",
    paste(missing_serpins, collapse = ", ")
  )
}

duplicate_candidates <- diopt_candidates %>%
  count(Drosophila_Gene) %>%
  filter(n > 1)

if (nrow(duplicate_candidates) > 0) {
  warning(
    "The DIOPT input table contains multiple rows for selected candidate(s): ",
    paste(duplicate_candidates$Drosophila_Gene, collapse = ", "),
    ". All rows will be retained."
  )
}

yes_no_to_binary <- function(x) {
  case_when(
    str_to_lower(as.character(x)) %in% c("yes", "y", "true", "1") ~ 1,
    str_to_lower(as.character(x)) %in% c("no", "n", "false", "0") ~ 0,
    TRUE ~ NA_real_
  )
}

diopt_summary <- diopt_candidates %>%
  transmute(
    Gene = Drosophila_Gene,
    FlyBase_ID,
    Primary_Human_Ortholog,
    Primary_Human_Ensembl_ID,
    Primary_DIOPT_Score,
    Primary_Weighted_Score,
    Primary_Confidence_Rank,
    Primary_Best_Score,
    Primary_Reverse_Best_Score,
    Supporting_DBs = Primary_Supporting_Database_Count,
    SERPINE1_Predicted,
    SERPINE1_DIOPT_Score,
    SERPINE1_Weighted_Score,
    SERPINE1_Confidence_Rank,
    SERPINE1_Reverse_Best_Score
  ) %>%
  mutate(
    Reverse_Best = yes_no_to_binary(Primary_Reverse_Best_Score),
    SERPINE1_Predicted_Binary = yes_no_to_binary(SERPINE1_Predicted),
    SERPINE1_Reverse_Best = yes_no_to_binary(SERPINE1_Reverse_Best_Score)
  ) %>%
  arrange(
    desc(Primary_DIOPT_Score),
    desc(Primary_Weighted_Score),
    Gene
  )

# ============================================================
# QC: Require exactly one row per candidate
# ============================================================

expected_n <- length(selected_serpins)

if (nrow(diopt_summary) != expected_n) {
  stop(
    "Expected exactly ", expected_n,
    " candidate rows, but found ", nrow(diopt_summary),
    ". Resolve missing or duplicate candidates before continuing."
  )
}

if (anyDuplicated(diopt_summary$Gene)) {
  stop(
    "Duplicate candidate genes remain in the DIOPT summary. ",
    "Resolve duplicates before downstream analysis."
  )
}

if (!setequal(diopt_summary$Gene, selected_serpins)) {
  stop(
    "Candidate genes in the DIOPT summary do not exactly match ",
    "the predefined candidate set."
  )
}

message("DIOPT QC PASSED: all six candidate serpins are represented exactly once."
)

print(diopt_summary)

write.csv(
  diopt_summary,
  file.path(processed_dir, "Table_DIOPT_candidate_serpin_orthology.csv"),
  row.names = FALSE
)

write.csv(
  diopt_summary,
  file.path(tables_dir, "Table_DIOPT_candidate_serpin_orthology.csv"),
  row.names = FALSE
)

}
