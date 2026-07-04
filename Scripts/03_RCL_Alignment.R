# ============================================================
# 03_RCL_Alignment.R
# Reactive center loop alignment identity/similarity
# ============================================================

library(dplyr)
library(stringr)
library(readr)

# Prefer an environment variable for reproducibility across machines.
# If it is not set, fall back to the original local project path.
project_dir <- Sys.getenv(
  "xxxxx_xxxxxxxx_PROJECT_DIR",
  unset = "D:/xxxxx_xxxxxxxx xxx/xxxxx_xxxxxxxx xxx"
)

chapter3_dir <- file.path(project_dir, "Bioinformatics chapter three")

data_dir <- file.path(chapter3_dir, "Data")
processed_dir <- file.path(data_dir, "Processed")
tables_dir <- file.path(chapter3_dir, "Tables", "Alignment")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

rcl_file <- file.path(data_dir, "RCL_sequences.csv")

if (!file.exists(rcl_file)) {
  stop("RCL input file not found: ", rcl_file)
}

rcl <- read.csv(rcl_file, stringsAsFactors = FALSE)

if (nrow(rcl) == 0) {
  warning("RCL input table has zero rows. Output tables will be empty.")
}

if (ncol(rcl) == 0) {
  stop("RCL input table has zero columns.")
}

required_columns <- c("Gene", "RCL")
missing_columns <- setdiff(required_columns, names(rcl))

if (length(missing_columns) > 0) {
  stop(
    "RCL input table is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

selected_serpins <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)

reference_gene <- "PAI1"

if (length(selected_serpins) == 0) {
  stop("selected_serpins is empty. Add at least one candidate gene before calculating RCL statistics.")
}

aa_similarity_groups <- list(
  c("A", "V", "L", "I", "M"),
  c("F", "Y", "W"),
  c("S", "T"),
  c("K", "R", "H"),
  c("D", "E"),
  c("N", "Q"),
  c("C"),
  c("G"),
  c("P")
)

is_similar_aa <- function(a, b) {
  if (is.na(a) || is.na(b)) {
    return(FALSE)
  }

  a <- str_to_upper(a)
  b <- str_to_upper(b)

  if (a == "-" || b == "-") {
    return(FALSE)
  }

  if (a == b) {
    return(TRUE)
  }

  any(vapply(
    aa_similarity_groups,
    function(group) a %in% group && b %in% group,
    logical(1)
  ))
}

empty_rcl_stats <- function(rcl_length = NA_real_, compared = NA_real_) {
  data.frame(
    RCL_Length = rcl_length,
    RCL_Compared_Positions = compared,
    RCL_Identity_Count = NA_real_,
    RCL_Similarity_Count = NA_real_,
    RCL_Identity = NA_real_,
    RCL_Similarity = NA_real_
  )
}

pairwise_rcl_stats <- function(query_seq, ref_seq) {
  if (length(query_seq) != 1 || length(ref_seq) != 1 || is.na(query_seq) || is.na(ref_seq)) {
    return(empty_rcl_stats())
  }

  q <- strsplit(str_to_upper(query_seq), "")[[1]]
  r <- strsplit(str_to_upper(ref_seq), "")[[1]]

  n <- min(length(q), length(r))

  if (n == 0) {
    return(empty_rcl_stats(rcl_length = 0, compared = 0))
  }

  q <- q[1:n]
  r <- r[1:n]

  aligned_positions <- q != "-" & r != "-"
  compared <- sum(aligned_positions)

  if (compared == 0) {
    return(empty_rcl_stats(rcl_length = n, compared = compared))
  }

  q_compared <- q[aligned_positions]
  r_compared <- r[aligned_positions]

  identity_count <- sum(q_compared == r_compared)
  similarity_count <- sum(vapply(
    seq_along(q_compared),
    function(i) is_similar_aa(q_compared[i], r_compared[i]),
    logical(1)
  ))

  data.frame(
    RCL_Length = n,
    RCL_Compared_Positions = compared,
    RCL_Identity_Count = identity_count,
    RCL_Similarity_Count = similarity_count,
    RCL_Identity = 100 * identity_count / compared,
    RCL_Similarity = 100 * similarity_count / compared
  )
}

reference_rows <- rcl %>%
  filter(Gene == reference_gene)

if (nrow(reference_rows) == 0) {
  stop("Reference gene not found in RCL input table: ", reference_gene)
}

if (nrow(reference_rows) > 1) {
  warning(
    "Multiple rows found for reference gene ",
    reference_gene,
    ". The first non-missing RCL sequence will be used."
  )
}

ref_rcl <- reference_rows %>%
  filter(!is.na(RCL), RCL != "") %>%
  pull(RCL)

if (length(ref_rcl) == 0) {
  stop("Reference gene has no usable RCL sequence: ", reference_gene)
}

ref_rcl <- ref_rcl[1]

missing_serpins <- setdiff(selected_serpins, unique(rcl$Gene))

if (length(missing_serpins) > 0) {
  warning(
    "The following selected serpin(s) were not found in the RCL input table: ",
    paste(missing_serpins, collapse = ", ")
  )
}

candidate_rcls <- rcl %>%
  filter(Gene %in% selected_serpins)

duplicate_candidates <- candidate_rcls %>%
  count(Gene) %>%
  filter(n > 1)

if (nrow(duplicate_candidates) > 0) {
  warning(
    "Multiple RCL rows were found for selected candidate(s): ",
    paste(duplicate_candidates$Gene, collapse = ", "),
    ". All rows will be retained."
  )
}

rcl_stats <- candidate_rcls %>%
  rowwise() %>%
  mutate(
    stats = list(pairwise_rcl_stats(RCL, ref_rcl))
  ) %>%
  ungroup() %>%
  select(Gene, RCL, stats) %>%
  tidyr::unnest(stats) %>%
  mutate(
    across(
      c(
        RCL_Length,
        RCL_Compared_Positions,
        RCL_Identity_Count,
        RCL_Similarity_Count,
        RCL_Identity,
        RCL_Similarity
      ),
      as.numeric
    )
  ) %>%
  arrange(desc(RCL_Similarity), desc(RCL_Identity), Gene)

if (any(is.na(rcl_stats$RCL_Similarity))) {
  warning("One or more RCL similarity values are NA, likely due to missing sequences or zero compared positions.")
}

if (nrow(rcl_stats) == 0) {
  warning("No RCL alignment statistics were generated for the selected serpins.")
}

print(rcl_stats)

write.csv(
  rcl_stats,
  file.path(processed_dir, "RCL_Alignment_Stats_vs_Human_PAI1.csv"),
  row.names = FALSE
)

write.csv(
  rcl_stats,
  file.path(tables_dir, "RCL_Alignment_Stats_vs_Human_PAI1.csv"),
  row.names = FALSE
)
