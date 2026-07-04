# ============================================================
# 02_Protein_Alignment_FullLength.R
# Full-length protein alignment identity/similarity
# ============================================================

library(dplyr)
library(stringr)
library(readr)

# Prefer an environment variable for reproducibility across machines.
# If it is not set, fall back to the original local project path.
project_dir <- Sys.getenv(
  "xxxxxx_xxxxxxxx_PROJECT_DIR",
  unset = "x:/xxxxxx-xxxxxxxx xxx/xxxxxx-xxxxxxxx-xxx"
)

chapter3_dir <- file.path(project_dir, "Bioinformatics chapter three")

data_dir <- file.path(chapter3_dir, "Data")
processed_dir <- file.path(data_dir, "Processed")
tables_dir <- file.path(chapter3_dir, "Tables", "Alignment")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

full_align_file <- file.path(data_dir, "full_alignment_text_output.csv")

if (!file.exists(full_align_file)) {
  stop("Full-length alignment input file not found: ", full_align_file)
}

full_align <- read.csv(full_align_file, stringsAsFactors = FALSE)

if (nrow(full_align) == 0) {
  warning("Full-length alignment input table has zero rows. Output tables will be empty.")
}

if (ncol(full_align) == 0) {
  stop("Full-length alignment input table has zero columns.")
}

required_columns <- c("Gene", "Alignment")
missing_columns <- setdiff(required_columns, names(full_align))

if (length(missing_columns) > 0) {
  stop(
    "Full-length alignment table is missing required column(s): ",
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
  stop("selected_serpins is empty. Add at least one candidate gene before calculating alignment statistics.")
}

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

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

pairwise_alignment_stats <- function(query_seq, ref_seq) {
  if (length(query_seq) != 1 || length(ref_seq) != 1 || is.na(query_seq) || is.na(ref_seq)) {
    return(data.frame(
      Alignment_Length = NA_real_,
      Compared_Positions = NA_real_,
      Identity_Count = NA_real_,
      Similarity_Count = NA_real_,
      Full_Length_Identity = NA_real_,
      Full_Length_Similarity = NA_real_,
      Query_Gaps = NA_real_,
      Reference_Gaps = NA_real_
    ))
  }

  q <- strsplit(str_to_upper(query_seq), "")[[1]]
  r <- strsplit(str_to_upper(ref_seq), "")[[1]]

  n <- min(length(q), length(r))

  if (n == 0) {
    return(data.frame(
      Alignment_Length = 0,
      Compared_Positions = 0,
      Identity_Count = 0,
      Similarity_Count = 0,
      Full_Length_Identity = NA_real_,
      Full_Length_Similarity = NA_real_,
      Query_Gaps = 0,
      Reference_Gaps = 0
    ))
  }

  q <- q[1:n]
  r <- r[1:n]

  aligned_positions <- q != "-" & r != "-"
  compared <- sum(aligned_positions)

  identity_count <- sum(q[aligned_positions] == r[aligned_positions])
  similarity_count <- sum(vapply(
    seq_along(q[aligned_positions]),
    function(i) is_similar_aa(q[aligned_positions][i], r[aligned_positions][i]),
    logical(1)
  ))

  gaps_query <- sum(q == "-")
  gaps_reference <- sum(r == "-")

  data.frame(
    Alignment_Length = n,
    Compared_Positions = compared,
    Identity_Count = identity_count,
    Similarity_Count = similarity_count,
    Full_Length_Identity = ifelse(compared > 0, 100 * identity_count / compared, NA_real_),
    Full_Length_Similarity = ifelse(compared > 0, 100 * similarity_count / compared, NA_real_),
    Query_Gaps = gaps_query,
    Reference_Gaps = gaps_reference
  )
}

reference_rows <- full_align %>%
  filter(Gene == reference_gene)

if (nrow(reference_rows) == 0) {
  stop("Reference gene not found in full-length alignment table: ", reference_gene)
}

if (nrow(reference_rows) > 1) {
  warning(
    "Multiple rows found for reference gene ",
    reference_gene,
    ". The first non-missing alignment will be used."
  )
}

ref_seq <- reference_rows %>%
  filter(!is.na(Alignment), Alignment != "") %>%
  pull(Alignment)

if (length(ref_seq) == 0) {
  stop("Reference gene has no usable alignment sequence: ", reference_gene)
}

ref_seq <- ref_seq[1]

missing_serpins <- setdiff(selected_serpins, unique(full_align$Gene))

if (length(missing_serpins) > 0) {
  warning(
    "The following selected serpin(s) were not found in the full-length alignment table: ",
    paste(missing_serpins, collapse = ", ")
  )
}

candidate_alignments <- full_align %>%
  filter(Gene %in% selected_serpins)

duplicate_candidates <- candidate_alignments %>%
  count(Gene) %>%
  filter(n > 1)

if (nrow(duplicate_candidates) > 0) {
  warning(
    "Multiple alignment rows were found for selected candidate(s): ",
    paste(duplicate_candidates$Gene, collapse = ", "),
    ". All rows will be retained."
  )
}

full_length_stats <- candidate_alignments %>%
  rowwise() %>%
  mutate(
    stats = list(pairwise_alignment_stats(Alignment, ref_seq))
  ) %>%
  ungroup() %>%
  select(Gene, stats) %>%
  tidyr::unnest(stats) %>%
  mutate(
    across(
      -Gene,
      as.numeric
    )
  ) %>%
  arrange(desc(Full_Length_Similarity), desc(Full_Length_Identity), Gene)

if (any(is.na(full_length_stats$Full_Length_Similarity))) {
  warning("One or more full-length similarity values are NA, likely due to missing sequences or zero compared positions.")
}

if (nrow(full_length_stats) == 0) {
  warning("No full-length alignment statistics were generated for the selected serpins.")
}

print(full_length_stats)

write.csv(
  full_length_stats,
  file.path(processed_dir, "Full_Length_Alignment_Stats_vs_Human_PAI1.csv"),
  row.names = FALSE
)

write.csv(
  full_length_stats,
  file.path(tables_dir, "Full_Length_Alignment_Stats_vs_Human_PAI1.csv"),
  row.names = FALSE
)
