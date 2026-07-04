# ============================================================
# 01_DIOPT_Orthology.R
# DIOPT orthology summary for candidate Drosophila serpins
# ============================================================

library(dplyr)
library(readr)
library(stringr)

# Prefer an environment variable for reproducibility across machines.
# If it is not set, fall back to the original local project path.
project_dir <- Sys.getenv(
  "filename",
  unset = "filepath"
)

chapter3_dir <- file.path(project_dir, "Bioinformatics chapter three")

data_dir <- file.path(chapter3_dir, "Data")
processed_dir <- file.path(data_dir, "Processed")
tables_dir <- file.path(chapter3_dir, "Tables", "Orthology")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

diopt_file <- file.path(data_dir, "orthology_mapping.csv")

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

if (nrow(diopt_summary) == 0) {
  warning("No DIOPT summary rows were generated for the selected serpins.")
}

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
