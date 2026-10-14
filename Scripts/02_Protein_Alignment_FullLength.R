# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 02_Protein_Alignment_FullLength.R
# Run from the repository root.
# ==============================================================================
# 02_Protein_Alignment_FullLength.R
#
# Comparative functional genomics analysis
#
# PURPOSE
# -------
# Validate full-length amino-acid sequences, generate a new multiple-sequence
# alignment, and calculate pairwise full-length sequence identity and
# BLOSUM62-based amino-acid similarity between each Drosophila serpin
# candidate and human SERPINE1/PAI-1.
#
# INPUT
# -----
# Uniport FASTA FlySerpins_Human PAI-1.txt
#
# OUTPUTS
# -------
# QC/
#   02_Sequence_QC.csv
#   02_Alignment_QC.csv
#   02_sessionInfo.txt
#
# Processed/
#   02_Validated_FullLength_Sequences.fasta
#   02_FullLength_Multiple_Alignment.fasta
#   02_Full_Length_Alignment_Stats_vs_Human_PAI1.csv
#
# Tables/
#   02_Full_Length_Alignment_Stats_vs_Human_PAI1.csv
#
# IMPORTANT
# ---------
# This script does NOT use the previous full_alignment_text_output.csv.
# The alignment is regenerated directly from the source FASTA sequences.
#
# RCL analysis is intentionally NOT performed here.
# ==============================================================================


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_packages <- c(
  "Biostrings",
  "msa",
  "dplyr",
  "readr",
  "stringr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running Step 2."
  )
}

library(Biostrings)
library(msa)
library(dplyr)
library(readr)
library(stringr)
library(tibble)


# ==============================================================================
# 2. FILE PATHS
# ==============================================================================

step2_dir <- "Sequence/FullLength"

fasta_file <- file.path(
  step2_dir,
  "Uniport FASTA FlySerpins_Human PAI-1.txt"
)

if (!file.exists(fasta_file)) {
  stop(
    "FASTA input file not found:\n",
    fasta_file
  )
}

message("FASTA input file found:")
message(fasta_file)


# ==============================================================================
# 3. OUTPUT DIRECTORIES
# ==============================================================================

processed_dir <- file.path(step2_dir, "Processed")
tables_dir    <- file.path(step2_dir, "Tables")
qc_dir        <- file.path(step2_dir, "QC")

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 4. EXPECTED PROTEINS
# ==============================================================================

expected_genes <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F",
  "SERPINE1"
)

candidate_genes <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)

reference_gene <- "SERPINE1"


# ==============================================================================
# 5. READ RAW FASTA
# ==============================================================================

raw_sequences <- Biostrings::readAAStringSet(
  fasta_file,
  format = "fasta"
)

if (length(raw_sequences) == 0) {
  stop("No protein sequences were detected in the FASTA file.")
}

message(
  "Number of FASTA records detected: ",
  length(raw_sequences)
)


# ==============================================================================
# 6. PARSE FASTA HEADERS
# ==============================================================================

raw_headers <- names(raw_sequences)


parse_record_type <- function(header) {

  case_when(
    str_detect(header, "^sp\\|") ~ "Swiss-Prot",
    str_detect(header, "^tr\\|") ~ "TrEMBL",
    TRUE                         ~ "Other"
  )
}


parse_accession <- function(header) {

  pieces <- str_split(header, "\\|", simplify = TRUE)

  if (ncol(pieces) >= 2) {
    return(pieces[, 2])
  }

  return(NA_character_)
}


parse_gene <- function(header) {

  gene <- str_match(
    header,
    "GN=([^ ]+)"
  )[, 2]

  gene
}


sequence_metadata <- tibble(
  Original_Header = raw_headers,
  Record_Type = parse_record_type(raw_headers),
  Accession = parse_accession(raw_headers),
  Gene = parse_gene(raw_headers),
  Sequence = as.character(raw_sequences),
  Sequence_Length = nchar(as.character(raw_sequences))
)

# ==============================================================================
# 7. STANDARDIZE HUMAN REFERENCE NAME
# ==============================================================================

sequence_metadata <- sequence_metadata %>%
  mutate(
    Gene = case_when(
      Gene == "SERPINE1" ~ "SERPINE1",
      TRUE ~ Gene
    )
  )


# ==============================================================================
# 8. BASIC HEADER QC
# ==============================================================================

if (any(is.na(sequence_metadata$Gene))) {

  bad_headers <- sequence_metadata %>%
    filter(is.na(Gene)) %>%
    pull(Original_Header)

  stop(
    "Could not extract GN= gene symbol from the following FASTA header(s):\n",
    paste(bad_headers, collapse = "\n")
  )
}


if (any(is.na(sequence_metadata$Accession))) {

  stop(
    "One or more UniProt accession numbers could not be parsed from ",
    "the FASTA headers."
  )
}


# ==============================================================================
# 9. VERIFY EXPECTED GENES
# ==============================================================================

observed_genes <- unique(sequence_metadata$Gene)

missing_genes <- setdiff(
  expected_genes,
  observed_genes
)

unexpected_genes <- setdiff(
  observed_genes,
  expected_genes
)


if (length(missing_genes) > 0) {

  stop(
    "Missing expected protein(s): ",
    paste(missing_genes, collapse = ", ")
  )
}


if (length(unexpected_genes) > 0) {

  stop(
    "Unexpected protein(s) detected in FASTA: ",
    paste(unexpected_genes, collapse = ", ")
  )
}


# ==============================================================================
# 10. CHECK FOR DUPLICATE GENE RECORDS
# ==============================================================================

duplicate_genes <- sequence_metadata %>%
  count(Gene, name = "Number_of_Records") %>%
  filter(Number_of_Records != 1)


if (nrow(duplicate_genes) > 0) {

  print(duplicate_genes)

  stop(
    "Each expected gene must occur exactly once in the Step 2 FASTA file."
  )
}


# ==============================================================================
# 11. CHECK FOR DUPLICATE ACCESSIONS
# ==============================================================================

duplicate_accessions <- sequence_metadata %>%
  count(Accession, name = "Number_of_Records") %>%
  filter(Number_of_Records > 1)


if (nrow(duplicate_accessions) > 0) {

  print(duplicate_accessions)

  stop(
    "Duplicate UniProt accession(s) detected."
  )
}


# ==============================================================================
# 12. AMINO-ACID CHARACTER QC
# ==============================================================================

# Standard 20 amino acids.
standard_amino_acids <- strsplit(
  "ACDEFGHIKLMNPQRSTVWY",
  ""
)[[1]]


find_nonstandard_aa <- function(sequence) {

  aa <- strsplit(
    str_to_upper(sequence),
    ""
  )[[1]]

  bad <- setdiff(
    unique(aa),
    standard_amino_acids
  )

  if (length(bad) == 0) {
    return("")
  }

  paste(
    sort(bad),
    collapse = ";"
  )
}


sequence_metadata <- sequence_metadata %>%
  rowwise() %>%
  mutate(
    Nonstandard_AA = find_nonstandard_aa(Sequence)
  ) %>%
  ungroup()


if (any(sequence_metadata$Nonstandard_AA != "")) {

  bad_sequences <- sequence_metadata %>%
    filter(Nonstandard_AA != "") %>%
    select(
      Gene,
      Accession,
      Nonstandard_AA
    )

  print(bad_sequences)

  stop(
    "Non-standard amino-acid characters were detected. ",
    "Resolve these sequences before alignment."
  )
}


# ==============================================================================
# 13. CHECK FOR IDENTICAL PROTEIN SEQUENCES
# ==============================================================================

sequence_metadata <- sequence_metadata %>%
  mutate(
    Duplicate_Sequence = duplicated(Sequence) |
      duplicated(Sequence, fromLast = TRUE)
  )


if (any(sequence_metadata$Duplicate_Sequence)) {

  warning(
    "Two or more genes have identical full-length protein sequences. ",
    "Review the sequence QC table."
  )
}


# ==============================================================================
# 14. BASIC LENGTH QC
# ==============================================================================

if (any(sequence_metadata$Sequence_Length <= 0)) {

  stop(
    "One or more protein sequences have zero length."
  )
}


# ==============================================================================
# 15. FINAL SEQUENCE QC STATUS
# ==============================================================================

sequence_qc <- sequence_metadata %>%
  mutate(
    Sequence_QC = case_when(
      Sequence_Length <= 0 ~ "FAIL",
      Nonstandard_AA != "" ~ "FAIL",
      TRUE ~ "PASS"
    )
  ) %>%
  select(
    Gene,
    Accession,
    Record_Type,
    Sequence_Length,
    Nonstandard_AA,
    Duplicate_Sequence,
    Sequence_QC,
    Original_Header
  ) %>%
  arrange(
    factor(
      Gene,
      levels = expected_genes
    )
  )


print(sequence_qc)


write_csv(
  sequence_qc,
  file.path(
    qc_dir,
    "02_Sequence_QC.csv"
  )
)


if (any(sequence_qc$Sequence_QC != "PASS")) {

  stop(
    "Sequence QC failed. Alignment will not be performed."
  )
}


# ==============================================================================
# 16. CREATE VALIDATED AAStringSet
# ==============================================================================

validated_metadata <- sequence_metadata %>%
  arrange(
    factor(
      Gene,
      levels = expected_genes
    )
  )


validated_sequences <- AAStringSet(
  validated_metadata$Sequence
)

names(validated_sequences) <- validated_metadata$Gene


# ==============================================================================
# 17. WRITE VALIDATED FASTA
# ==============================================================================

validated_fasta_file <- file.path(
  processed_dir,
  "02_Validated_FullLength_Sequences.fasta"
)


Biostrings::writeXStringSet(
  validated_sequences,
  filepath = validated_fasta_file,
  format = "fasta"
)


message(
  "Validated FASTA written to:\n",
  validated_fasta_file
)


# ==============================================================================
# 18. GENERATE NEW MULTIPLE-SEQUENCE ALIGNMENT
# ==============================================================================

message("Generating full-length multiple-sequence alignment...")


full_alignment <- msa::msa(
  validated_sequences,
  method = "ClustalOmega",
  order = "input"
)


print(full_alignment)


# ==============================================================================
# 19. CONVERT ALIGNMENT TO ALIGNED SEQUENCES
# ==============================================================================

aligned_set <- Biostrings::AAStringSet(
  full_alignment
)

aligned_sequences <- as.character(
  aligned_set
)

aligned_names <- names(
  aligned_set
)


if (is.null(aligned_names)) {

  stop(
    "Aligned sequence names were lost during alignment."
  )
}


names(aligned_sequences) <- aligned_names


# ==============================================================================
# 20. ALIGNMENT LENGTH QC
# ==============================================================================

alignment_lengths <- nchar(
  aligned_sequences
)


alignment_qc <- tibble(
  Gene = names(aligned_sequences),
  Alignment_Length = alignment_lengths
)


if (length(unique(alignment_lengths)) != 1) {

  print(alignment_qc)

  stop(
    "Aligned sequences do not all have the same alignment length."
  )
}


alignment_length <- unique(
  alignment_lengths
)


alignment_qc <- alignment_qc %>%
  mutate(
    Alignment_Length_QC = "PASS"
  )


write_csv(
  alignment_qc,
  file.path(
    qc_dir,
    "02_Alignment_QC.csv"
  )
)


message(
  "Alignment length QC passed. Alignment contains ",
  alignment_length,
  " columns."
)


# ==============================================================================
# 21. WRITE NEW ALIGNMENT FASTA
# ==============================================================================

aligned_output <- AAStringSet(
  aligned_sequences
)

names(aligned_output) <- names(
  aligned_sequences
)


alignment_fasta_file <- file.path(
  processed_dir,
  "02_FullLength_Multiple_Alignment.fasta"
)


Biostrings::writeXStringSet(
  aligned_output,
  filepath = alignment_fasta_file,
  format = "fasta"
)


# ==============================================================================
# 22. LOAD BLOSUM62 SUBSTITUTION MATRIX
# ==============================================================================

data(
  "BLOSUM62",
  package = "Biostrings"
)


if (!exists("BLOSUM62")) {

  stop(
    "BLOSUM62 substitution matrix could not be loaded."
  )
}


# ==============================================================================
# 23. PAIRWISE ALIGNMENT STATISTICS
# ==============================================================================

pairwise_alignment_stats <- function(
    query_seq,
    reference_seq,
    query_gene,
    reference_gene
) {

  q <- strsplit(
    str_to_upper(query_seq),
    ""
  )[[1]]

  r <- strsplit(
    str_to_upper(reference_seq),
    ""
  )[[1]]


  # --------------------------------------------------------------------------
  # Strict alignment-length QC
  # --------------------------------------------------------------------------

  if (length(q) != length(r)) {

    stop(
      "Aligned sequences have different lengths for ",
      query_gene,
      " and ",
      reference_gene,
      ": ",
      length(q),
      " vs ",
      length(r)
    )
  }


  alignment_length <- length(q)


  # --------------------------------------------------------------------------
  # Identify alignment-column types
  # --------------------------------------------------------------------------

  both_residues <- q != "-" & r != "-"

  gap_columns <- q == "-" | r == "-"

  query_gap_columns <- q == "-"

  reference_gap_columns <- r == "-"


  compared_positions <- sum(
    both_residues
  )


  if (compared_positions == 0) {

    stop(
      "No mutually occupied alignment positions were found for ",
      query_gene,
      " versus ",
      reference_gene,
      "."
    )
  }


  # --------------------------------------------------------------------------
  # Identity
  # --------------------------------------------------------------------------

  q_compared <- q[both_residues]
  r_compared <- r[both_residues]


  identity_vector <- q_compared == r_compared

  identity_count <- sum(
    identity_vector
  )


  pairwise_nongap_identity <- 100 *
    identity_count /
    compared_positions


  # --------------------------------------------------------------------------
  # BLOSUM62 similarity
  #
  # A substitution is considered similar if its BLOSUM62 score is > 0.
  # Identical residues are therefore included as similar residues.
  # --------------------------------------------------------------------------

  blosum_scores <- vapply(
    seq_along(q_compared),
    function(i) {

      BLOSUM62[
        q_compared[i],
        r_compared[i]
      ]

    },
    numeric(1)
  )


  similarity_vector <- blosum_scores > 0

  similarity_count <- sum(
    similarity_vector
  )


  pairwise_nongap_similarity <- 100 *
    similarity_count /
    compared_positions


  # --------------------------------------------------------------------------
  # Alignment coverage
  #
  # Percentage of complete MSA columns in which both proteins contain
  # an amino-acid residue.
  # --------------------------------------------------------------------------

  compared_coverage <- 100 *
    compared_positions /
    alignment_length


  # --------------------------------------------------------------------------
  # Identity relative to complete alignment
  #
  # Retained as a QC/descriptive metric. Gap-containing columns contribute
  # to the denominator.
  # --------------------------------------------------------------------------

  identity_complete_alignment <- 100 *
    identity_count /
    alignment_length


  similarity_complete_alignment <- 100 *
    similarity_count /
    alignment_length


  tibble(
    Gene = query_gene,
    Reference = reference_gene,

    Alignment_Length = alignment_length,

    Compared_Positions = compared_positions,

    Gap_Columns = sum(gap_columns),

    Query_Gaps = sum(query_gap_columns),

    Reference_Gaps = sum(reference_gap_columns),

    Compared_Coverage_Percent = compared_coverage,

    Identity_Count = identity_count,

    Similarity_Count_BLOSUM62 = similarity_count,

    Full_Length_Identity = pairwise_nongap_identity,

    Full_Length_Similarity = pairwise_nongap_similarity,

    Identity_Complete_Alignment = identity_complete_alignment,

    Similarity_Complete_Alignment = similarity_complete_alignment
  )
}


# ==============================================================================
# 24. IDENTIFY HUMAN PAI-1 REFERENCE
# ==============================================================================

if (!reference_gene %in% names(aligned_sequences)) {

  stop(
    "Human reference ",
    reference_gene,
    " was not found in the final alignment."
  )
}


reference_alignment <- aligned_sequences[
  reference_gene
]


# ==============================================================================
# 25. CALCULATE STATISTICS FOR SIX FLY SERPINS
# ==============================================================================

full_length_stats <- lapply(
  candidate_genes,
  function(gene) {

    if (!gene %in% names(aligned_sequences)) {

      stop(
        "Candidate gene missing from alignment: ",
        gene
      )
    }


    pairwise_alignment_stats(
      query_seq = aligned_sequences[[gene]],
      reference_seq = reference_alignment[[1]],
      query_gene = gene,
      reference_gene = reference_gene
    )
  }
) %>%
  bind_rows()


# ==============================================================================
# 26. ADD ORIGINAL UNALIGNED SEQUENCE LENGTHS
# ==============================================================================

# Calculate unaligned protein lengths directly from the validated sequences.
sequence_lengths <- tibble(
  Gene = names(validated_sequences),
  Query_Sequence_Length = nchar(as.character(validated_sequences))
)

# Determine the length of the human SERPINE1 / PAI-1 reference sequence.
reference_length <- nchar(
  as.character(validated_sequences[[reference_gene]])
)

if (length(reference_length) != 1) {
  stop(
    "Could not uniquely determine human PAI-1 sequence length."
  )
}

if (is.na(reference_length) || reference_length <= 0) {
  stop(
    "Human PAI-1 reference sequence has an invalid sequence length."
  )
}

# Join sequence lengths to the alignment statistics.
full_length_stats <- full_length_stats %>%
  left_join(
    sequence_lengths,
    by = "Gene"
  ) %>%
  mutate(
    Query_Sequence_Length = as.numeric(Query_Sequence_Length),
    Reference_Sequence_Length = as.numeric(reference_length)
  ) %>%
  select(
    Gene,
    Reference,

    Query_Sequence_Length,
    Reference_Sequence_Length,

    Alignment_Length,
    Compared_Positions,
    Compared_Coverage_Percent,

    Gap_Columns,
    Query_Gaps,
    Reference_Gaps,

    Identity_Count,
    Similarity_Count_BLOSUM62,

    Full_Length_Identity,
    Full_Length_Similarity,

    Identity_Complete_Alignment,
    Similarity_Complete_Alignment
  ) %>%
  arrange(
    desc(Full_Length_Similarity),
    desc(Full_Length_Identity),
    Gene
  )

# ==============================================================================
# 27. NUMERIC QC
# ==============================================================================

numeric_columns <- c(
  "Compared_Coverage_Percent",
  "Full_Length_Identity",
  "Full_Length_Similarity",
  "Identity_Complete_Alignment",
  "Similarity_Complete_Alignment"
)


if (
  any(
    is.na(
      full_length_stats[
        numeric_columns
      ]
    )
  )
) {

  stop(
    "NA values were detected in required full-length alignment statistics."
  )
}


if (
  any(
    full_length_stats$Full_Length_Identity < 0 |
      full_length_stats$Full_Length_Identity > 100
  )
) {

  stop(
    "Full-length identity values outside 0-100% were detected."
  )
}


if (
  any(
    full_length_stats$Full_Length_Similarity < 0 |
      full_length_stats$Full_Length_Similarity > 100
  )
) {

  stop(
    "Full-length similarity values outside 0-100% were detected."
  )
}


if (
  any(
    full_length_stats$Full_Length_Similarity <
      full_length_stats$Full_Length_Identity
  )
) {

  stop(
    "BLOSUM62 similarity is lower than identity for at least one protein. ",
    "Review the similarity calculation."
  )
}


# ==============================================================================
# 28. VERIFY EXACTLY SIX CANDIDATES
# ==============================================================================

if (nrow(full_length_stats) != length(candidate_genes)) {

  stop(
    "Expected ",
    length(candidate_genes),
    " candidate comparisons but obtained ",
    nrow(full_length_stats),
    "."
  )
}


if (
  !setequal(
    full_length_stats$Gene,
    candidate_genes
  )
) {

  stop(
    "Final alignment statistics do not contain exactly the six ",
    "expected Drosophila candidate genes."
  )
}


# ==============================================================================
# 29. PRINT FINAL RESULTS
# ==============================================================================

message(
  "\n============================================================"
)

message(
  "FULL-LENGTH ALIGNMENT RESULTS"
)

message(
  "============================================================\n"
)


print(
  full_length_stats,
  n = Inf
)


# ==============================================================================
# 30. WRITE PROCESSED OUTPUT
# ==============================================================================

processed_output_file <- file.path(
  processed_dir,
  "02_Full_Length_Alignment_Stats_vs_Human_PAI1.csv"
)


write_csv(
  full_length_stats,
  processed_output_file
)


# ==============================================================================
# 31. WRITE study TABLE
# ==============================================================================

table_output_file <- file.path(
  tables_dir,
  "02_Full_Length_Alignment_Stats_vs_Human_PAI1.csv"
)


write_csv(
  full_length_stats,
  table_output_file
)


# ==============================================================================
# 32. SAVE SESSION INFORMATION
# ==============================================================================

session_file <- file.path(
  qc_dir,
  "02_sessionInfo.txt"
)


capture.output(
  sessionInfo(),
  file = session_file
)


# ==============================================================================
# 33. COMPLETION MESSAGE
# ==============================================================================

message(
  "\n============================================================"
)

message(
  "STEP 2 COMPLETED SUCCESSFULLY"
)

message(
  "============================================================"
)

message(
  "\nSequence QC:"
)

message(
  file.path(
    qc_dir,
    "02_Sequence_QC.csv"
  )
)

message(
  "\nValidated FASTA:"
)

message(
  validated_fasta_file
)

message(
  "\nNew multiple-sequence alignment:"
)

message(
  alignment_fasta_file
)

message(
  "\nFull-length alignment statistics:"
)

message(
  processed_output_file
)

message(
  "\nSession information:"
)

message(
  session_file
)

message(
  "\nAll required sequence and alignment QC checks passed."
)
