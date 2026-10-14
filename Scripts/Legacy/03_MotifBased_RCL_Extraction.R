# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 03_RCL_Alignment.R
# Run from the repository root.
# ==============================================================================
# 03_RCL_Domain.R
#
# Comparative functional genomics analysis
# REACTIVE CENTER LOOP (RCL) DOMAIN ANALYSIS
#
# PURPOSE
# -------
# This script:
#
# 1. Reads the selected full-length human and Drosophila serpin FASTA file.
#
# 2. Identifies proteins using their UniProt accession numbers.
#
# 3. Validates experimentally or curator-supported P1/P1' reactive sites
#    against the actual FASTA sequences.
#
# 4. Extracts standardized P17-P8' (25-aa) RCL windows for proteins for which
#    a defensible P1/P1' reactive bond has been established:
#
#       SERPINE1
#       SERPINB3
#       SERPINB4
#       SERPINB11
#       SERPINI1
#       Spn42Da isoform B (A1Z6R1)
#
# 5. Extracts broader candidate RCL-containing regions for:
#
#       Spn42Da Q8T0M5
#       Spn42Db
#       Spn42Dc
#       Spn42Dd
#       Spn28B
#       Spn28F
#
#    WITHOUT assigning predicted P1/P1' sites from Garrett et al. (2009).
#
# 6. Aligns the RCL/RCL-containing sequences independently using
#    Clustal Omega.
#
# 7. Locates the alignment columns corresponding specifically to the
#    experimentally established SERPINE1 P1/P1' reactive bond.
#
# 8. Produces QC tables, FASTA files, an alignment table, and an RCL
#    alignment figure.
#
#
# IMPORTANT INTERPRETIVE RULE
# ---------------------------
# The dashed boundaries in the final alignment figure identify the
# SERPINE1 P1/P1' alignment positions.
#
# They DO NOT indicate experimentally established P1/P1' sites for the
# other Drosophila serpins.
#
# Predicted Drosophila P1/P1' sites from Garrett et al. (2009) are NOT
# used to define the RCL boundaries in this analysis.
#
# ==============================================================================



# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(Biostrings)
library(msa)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(tibble)
library(ggplot2)
library(ggmsa)



# ==============================================================================
# 1. FILE PATHS
# ==============================================================================

step3_dir <- "Sequence/RCL"


# ------------------------------------------------------------------------------
# Input FASTA
# ------------------------------------------------------------------------------

fasta_file <- file.path(
  step3_dir,
  "selected_human_fly_Serpins_Fasta.txt"
)


# ------------------------------------------------------------------------------
# Output folders
# ------------------------------------------------------------------------------

processed_dir <- file.path(
  step3_dir,
  "Processed"
)

table_dir <- file.path(
  step3_dir,
  "Tables"
)

figure_dir <- file.path(
  step3_dir,
  "Figures"
)

qc_dir <- file.path(
  step3_dir,
  "QC"
)


dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
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

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# ==============================================================================
# 2. VERIFY INPUT FILE
# ==============================================================================

cat(
  "\n============================================================\n",
  "STEP 3: RCL DOMAIN ANALYSIS\n",
  "============================================================\n\n",
  sep = ""
)


cat(
  "Input FASTA:\n",
  fasta_file,
  "\n\n"
)


if (!file.exists(fasta_file)) {

  stop(
    "\nFASTA file was not found:\n",
    fasta_file,
    "\n\nCheck the file name and directory."
  )

}


cat(
  "FASTA file found successfully.\n\n"
)



# ==============================================================================
# 3. READ FULL-LENGTH PROTEIN FASTA
# ==============================================================================

seqs_raw <- readAAStringSet(
  fasta_file,
  format = "fasta"
)


cat(
  "Number of sequences read:",
  length(seqs_raw),
  "\n\n"
)


cat(
  "Original FASTA headers:\n\n"
)


print(
  names(seqs_raw)
)



# ==============================================================================
# 4. DEFINE ACCESSION KEY
# ==============================================================================

accession_key <- tribble(

  ~Accession, ~Gene, ~Display_Name, ~Species,

  "P05121",
  "SERPINE1",
  "SERPINE1",
  "Human",

  "P29508",
  "SERPINB3",
  "SERPINB3",
  "Human",

  "P48594",
  "SERPINB4",
  "SERPINB4",
  "Human",

  "Q96P15",
  "SERPINB11",
  "SERPINB11",
  "Human",

  "Q99574",
  "SERPINI1",
  "SERPINI1",
  "Human",

  "Q8T0M5",
  "Spn42Da_Q8T0M5",
  "Spn42Da (Q8T0M5)",
  "Drosophila",

  "A1Z6R1",
  "Spn42Da_B_A1Z6R1",
  "Spn42Da-B (A1Z6R1)",
  "Drosophila",

  "A1Z6R3",
  "Spn42Db",
  "Spn42Db",
  "Drosophila",

  "A1Z6R4",
  "Spn42Dc",
  "Spn42Dc",
  "Drosophila",

  "Q7YTY6",
  "Spn42Dd",
  "Spn42Dd",
  "Drosophila",

  "Q9VLZ8",
  "Spn28B",
  "Spn28B",
  "Drosophila",

  "Q9VLQ7",
  "Spn28F",
  "Spn28F",
  "Drosophila"

)



# ==============================================================================
# 5. EXTRACT UNIPROT ACCESSIONS FROM FASTA HEADERS
# ==============================================================================

extract_accession <- function(header) {

  parts <- str_split(
    header,
    "\\|"
  )[[1]]

  if (length(parts) >= 2) {

    return(parts[2])

  } else {

    return(NA_character_)

  }

}


fasta_accessions <- sapply(
  names(seqs_raw),
  extract_accession
)



# ==============================================================================
# 6. BUILD FASTA INVENTORY
# ==============================================================================

fasta_inventory <- tibble(

  Original_Header = names(seqs_raw),

  Accession = fasta_accessions,

  Protein_Length = width(seqs_raw)

) %>%

  left_join(
    accession_key,
    by = "Accession"
  )


cat(
  "\nFASTA inventory:\n\n"
)


print(
  fasta_inventory
)


write_csv(
  fasta_inventory,
  file.path(
    qc_dir,
    "03_FASTA_inventory.csv"
  )
)



# ==============================================================================
# 7. VERIFY ALL REQUIRED PROTEINS ARE PRESENT
# ==============================================================================

missing_accessions <- setdiff(
  accession_key$Accession,
  fasta_inventory$Accession
)


if (length(missing_accessions) > 0) {

  stop(
    "\nRequired UniProt accession(s) missing from FASTA:\n",
    paste(
      missing_accessions,
      collapse = "\n"
    )
  )

}


unexpected_accessions <- setdiff(
  fasta_inventory$Accession,
  accession_key$Accession
)


if (length(unexpected_accessions) > 0) {

  warning(
    "\nUnexpected sequence(s) found in FASTA:\n",
    paste(
      unexpected_accessions,
      collapse = "\n"
    )
  )

}


if (anyDuplicated(fasta_inventory$Accession) > 0) {

  stop(
    "\nDuplicate UniProt accession numbers were detected in the FASTA."
  )

}


cat(
  "\nAll 12 required UniProt accessions are present.\n"
)



# ==============================================================================
# 8. RENAME SEQUENCES
# ==============================================================================

new_names <- accession_key$Gene[

  match(
    fasta_accessions,
    accession_key$Accession
  )

]


names(seqs_raw) <- new_names



# ==============================================================================
# 9. SET ANALYSIS ORDER
# ==============================================================================

analysis_order <- c(

  "SERPINE1",

  "SERPINB3",

  "SERPINB4",

  "SERPINB11",

  "SERPINI1",

  "Spn42Da_Q8T0M5",

  "Spn42Da_B_A1Z6R1",

  "Spn42Db",

  "Spn42Dc",

  "Spn42Dd",

  "Spn28B",

  "Spn28F"

)


seqs <- seqs_raw[
  analysis_order
]



# ==============================================================================
# 10. FULL-LENGTH SEQUENCE INVENTORY
# ==============================================================================

sequence_qc <- tibble(

  Gene = names(seqs),

  Protein_Length = width(seqs)

) %>%

  left_join(
    accession_key %>%
      select(
        Gene,
        Accession,
        Display_Name,
        Species
      ),
    by = "Gene"
  ) %>%

  select(
    Gene,
    Display_Name,
    Accession,
    Species,
    Protein_Length
  )


cat(
  "\nFull-length sequence inventory:\n\n"
)


print(
  sequence_qc
)


write_csv(
  sequence_qc,
  file.path(
    qc_dir,
    "03_full_length_sequence_inventory.csv"
  )
)



# ==============================================================================
# 11. DEFINE VALIDATED / SUPPORTED REACTIVE SITES
# ==============================================================================
#
# These sites are NOT generated computationally by this script.
#
# Human serpins:
#
# SERPINE1
#   P1  = R369
#   P1' = M370
#
# SERPINB3
#   P1  = S354
#   P1' = S355
#
# SERPINB4
#   P1  = L354
#   P1' = S355
#
# SERPINB11
#   P1  = K357
#   P1' = S358
#
# SERPINI1
#   P1  = R362
#   P1' = M363
#
# Spn42Da isoform B:
#
#   P1  = R342
#   P1' = A343
#
# The Spn42Da-B assignment is supported by the experimental/structural
# characterization of the A/B RCL.
#
# Predicted sites for the other Drosophila proteins are deliberately omitted.
#
# ==============================================================================

reactive_sites <- tribble(

  ~Gene,
  ~P1_Position,
  ~P1,
  ~P1prime_Position,
  ~P1prime,
  ~Evidence_Class,

  "SERPINE1",
  369,
  "R",
  370,
  "M",
  "validated",

  "SERPINB3",
  354,
  "S",
  355,
  "S",
  "validated",

  "SERPINB4",
  354,
  "L",
  355,
  "S",
  "validated",

  "SERPINB11",
  357,
  "K",
  358,
  "S",
  "validated",

  "SERPINI1",
  362,
  "R",
  363,
  "M",
  "validated",

 "Spn42Da_B_A1Z6R1",
374,
"R",
375,
"A",
"experimentally_characterized_literature_R342_A343_full_FASTA_R374_A375"

)



# ==============================================================================
# 12. FUNCTION TO EXTRACT A SINGLE RESIDUE
# ==============================================================================

get_residue <- function(
  sequence_set,
  gene,
  position
) {

  as.character(

    subseq(
      sequence_set[[gene]],
      start = position,
      end = position
    )

  )

}



# ==============================================================================
# 13. VALIDATE P1/P1' RESIDUES AGAINST THE ACTUAL FASTA
# ==============================================================================

reactive_site_validation <- reactive_sites %>%

  rowwise() %>%

  mutate(

    Observed_P1 = get_residue(
      seqs,
      Gene,
      P1_Position
    ),

    Observed_P1prime = get_residue(
      seqs,
      Gene,
      P1prime_Position
    ),

    P1_OK =
      Observed_P1 == P1,

    P1prime_OK =
      Observed_P1prime == P1prime,

    Reactive_Site_OK =
      P1_OK & P1prime_OK

  ) %>%

  ungroup()


cat(
  "\nReactive-site validation:\n\n"
)


print(
  reactive_site_validation
)


write_csv(
  reactive_site_validation,
  file.path(
    qc_dir,
    "03_reactive_site_validation.csv"
  )
)



# ==============================================================================
# 14. HARD STOP IF REACTIVE-SITE VALIDATION FAILS
# ==============================================================================

if (
  any(
    !reactive_site_validation$Reactive_Site_OK
  )
) {

  failed_sites <- reactive_site_validation %>%
    filter(
      !Reactive_Site_OK
    )


  cat(
    "\nFAILED REACTIVE-SITE CHECKS:\n\n"
  )


  print(
    failed_sites
  )


  stop(
    "\nReactive-site validation FAILED.\n",
    "At least one stated P1/P1' residue does not match the supplied FASTA.\n",
    "Resolve the sequence/numbering discrepancy before continuing."
  )

}


cat(
  "\nAll supported reactive-site coordinates match the supplied FASTA.\n"
)



# ==============================================================================
# 15. EXTRACT STANDARDIZED P17-P8' RCL WINDOWS
# ==============================================================================
#
# Standardized window:
#
# P17 through P8'
#
# This gives:
#
#   16 residues before P1
#   P1
#   P1'
#   7 residues after P1'
#
# Total = 25 amino acids
#
# Therefore:
#
#   Start = P1 - 16
#   End   = P1' + 7
#
# Within the unaligned 25-aa sequence:
#
#   residue 17 = P1
#   residue 18 = P1'
#
# ==============================================================================

extract_standardized_rcl <- function(
  sequence_set,
  gene,
  p1_position,
  p1prime_position
) {

  start_position <- p1_position - 16

  end_position <- p1prime_position + 7


  protein_length <- length(
    sequence_set[[gene]]
  )


  if (
    start_position < 1 |
    end_position > protein_length
  ) {

    stop(
      "RCL extraction coordinates are outside the protein sequence for ",
      gene
    )

  }


  rcl <- subseq(
    sequence_set[[gene]],
    start = start_position,
    end = end_position
  )


  tibble(

    Gene = gene,

    RCL_Start = start_position,

    RCL_End = end_position,

    RCL_Length = length(rcl),

    RCL_Sequence = as.character(rcl)

  )

}



validated_rcls <- bind_rows(

  lapply(

    seq_len(
      nrow(reactive_sites)
    ),

    function(i) {

      extract_standardized_rcl(

        sequence_set = seqs,

        gene =
          reactive_sites$Gene[i],

        p1_position =
          reactive_sites$P1_Position[i],

        p1prime_position =
          reactive_sites$P1prime_Position[i]

      )

    }

  )

)



validated_rcls <- validated_rcls %>%

  left_join(
    reactive_sites,
    by = "Gene"
  ) %>%

  mutate(

    RCL_Definition =
      "P17-P8prime",

    Site_Status =
      Evidence_Class

  )


cat(
  "\nValidated standardized RCL windows:\n\n"
)


print(
  validated_rcls
)



# ==============================================================================
# 16. VERIFY ALL STANDARDIZED WINDOWS ARE 25 AA
# ==============================================================================

if (
  any(
    validated_rcls$RCL_Length != 25
  )
) {

  stop(
    "\nStandardized RCL extraction failed.\n",
    "Not all P17-P8' windows are exactly 25 amino acids."
  )

}



# ==============================================================================
# 17. SAVE VALIDATED RCL TABLE
# ==============================================================================

write_csv(
  validated_rcls,
  file.path(
    table_dir,
    "03_validated_standardized_RCL_sequences.csv"
  )
)



# ==============================================================================
# 18. HUMAN RCL SEQUENCE QC
# ==============================================================================
#
# These sequences were reconstructed from the supplied full-length FASTA
# using the stated reactive-site coordinates.
#
# This provides a second check against accidental coordinate/version changes.
#
# ==============================================================================

expected_human_rcls <- tribble(

  ~Gene,
  ~Expected_RCL,

  "SERPINE1",
  "ESGTVASSSTAVIVSARMAPEEIIM",

  "SERPINB3",
  "EEGAEAAAATAVVGFGSSPTSTNEE",

  "SERPINB4",
  "EEGVEAAAATAVVVVELSSPSTNEE",

  "SERPINB11",
  "EEGTEAAAATGDSIAVKSLPMRAQF",

  "SERPINI1",
  "EEGSEAAAVSGMIAISRMAVLYPQV"

)



human_rcl_check <- validated_rcls %>%

  inner_join(
    expected_human_rcls,
    by = "Gene"
  ) %>%

  mutate(

    Sequence_OK =
      RCL_Sequence == Expected_RCL

  )


cat(
  "\nHuman RCL sequence validation:\n\n"
)


print(
  human_rcl_check
)


write_csv(
  human_rcl_check,
  file.path(
    qc_dir,
    "03_human_RCL_sequence_validation.csv"
  )
)



if (
  any(
    !human_rcl_check$Sequence_OK
  )
) {

  stop(
    "\nHuman RCL sequence validation FAILED.\n",
    "At least one extracted RCL differs from the expected FASTA-derived sequence."
  )

}



# ==============================================================================
# 19. DISPLAY SPN42DA-B STANDARDIZED RCL
# ==============================================================================

spn42da_b_rcl <- validated_rcls %>%

  filter(
    Gene == "Spn42Da_B_A1Z6R1"
  )


cat(
  "\nSpn42Da-B experimentally supported standardized RCL:\n\n"
)


print(
  spn42da_b_rcl
)



# ==============================================================================
# 20. DEFINE FLY PROTEINS WITHOUT ASSIGNED P1/P1'
# ==============================================================================
#
# These proteins are intentionally NOT assigned predicted reactive sites:
#
#   Spn42Da Q8T0M5
#   Spn42Db
#   Spn42Dc
#   Spn42Dd
#   Spn28B
#   Spn28F
#
# Instead, a broader candidate RCL-containing C-terminal region is extracted.
#
# ==============================================================================

candidate_fly_genes <- c(

  "Spn42Da_Q8T0M5",

  "Spn42Db",

  "Spn42Dc",

  "Spn42Dd",

  "Spn28B",

  "Spn28F"

)



# ==============================================================================
# 21. EXTRACT CANDIDATE RCL-CONTAINING REGIONS
# ==============================================================================
#
# BIOLOGICAL RATIONALE
# --------------------
#
# The fly serpin sequences contain a conserved region preceding the highly
# variable C-terminal RCL-containing segment, followed downstream by a
# conserved DHPF-like serpin motif.
#
# Because the exact reactive bond is not independently established for these
# proteins in this analysis, we extract a BROADER region rather than forcing
# a P1/P1'-centered window.
#
# This extraction is used to visualize and inspect the RCL-containing region.
#
# It is NOT interpreted as an experimentally validated P17-P8' RCL.
#
# ==============================================================================

extract_candidate_region <- function(
  sequence_set,
  gene
) {

  x <- as.character(
    sequence_set[[gene]]
  )


  # ---------------------------------------------------------------------------
  # Locate upstream sequence anchor
  #
  # Several related variants occur among the selected fly serpins.
  # ---------------------------------------------------------------------------

  upstream_patterns <- c(
    "EVNEE",
    "DVNEA",
    "EINEV"
  )


  upstream_hits <- lapply(

    upstream_patterns,

    function(pattern) {

      str_locate_all(
        x,
        fixed(pattern)
      )[[1]]

    }

  )


  upstream_hits <- bind_rows(

    lapply(

      seq_along(upstream_hits),

      function(i) {

        hit <- upstream_hits[[i]]

        if (
          nrow(hit) == 0
        ) {

          return(NULL)

        }


        tibble(

          pattern =
            upstream_patterns[i],

          start =
            hit[, "start"],

          end =
            hit[, "end"]

        )

      }

    )

  )


  if (
    nrow(upstream_hits) == 0
  ) {

    stop(
      "\nCould not locate an upstream RCL-region anchor in ",
      gene,
      ".\nInspect the sequence manually."
    )

  }


  # Use the most C-terminal occurrence
  upstream_hit <- upstream_hits %>%

    arrange(
      desc(start)
    ) %>%

    slice(1)



  # ---------------------------------------------------------------------------
  # Locate downstream DHPF-like motif
  # ---------------------------------------------------------------------------

  downstream_match <- str_locate_all(
    x,
    fixed("DHPF")
  )[[1]]


  if (
    nrow(downstream_match) == 0
  ) {

    stop(
      "\nCould not locate downstream DHPF motif in ",
      gene,
      ".\nInspect the sequence manually."
    )

  }


  # Keep only downstream DHPF motifs occurring after upstream anchor
  downstream_candidates <- tibble(

    start =
      downstream_match[, "start"],

    end =
      downstream_match[, "end"]

  ) %>%

    filter(
      start > upstream_hit$end
    ) %>%

    arrange(
      start
    )


  if (
    nrow(downstream_candidates) == 0
  ) {

    stop(
      "\nNo DHPF motif occurs downstream of the selected upstream anchor in ",
      gene,
      "."
    )

  }


  # Use first downstream DHPF
  downstream_hit <-
    downstream_candidates %>%
    slice(1)



  # ---------------------------------------------------------------------------
  # Extract region
  # ---------------------------------------------------------------------------

  start_position <-
    upstream_hit$start


  end_position <-
    downstream_hit$end


  region <- substr(
    x,
    start_position,
    end_position
  )


  tibble(

    Gene = gene,

    Upstream_Anchor =
      upstream_hit$pattern,

    RCL_Start =
      start_position,

    RCL_End =
      end_position,

    RCL_Length =
      nchar(region),

    RCL_Sequence =
      region,

    RCL_Definition =
      "candidate_RCL_containing_region",

    Site_Status =
      "P1-P1prime_not_assigned"

  )

}



# ==============================================================================
# 22. APPLY CANDIDATE REGION EXTRACTION
# ==============================================================================

candidate_fly_rcls <- bind_rows(

  lapply(

    candidate_fly_genes,

    function(g) {

      extract_candidate_region(
        seqs,
        g
      )

    }

  )

)


cat(
  "\nCandidate Drosophila RCL-containing regions:\n\n"
)


print(
  candidate_fly_rcls
)



# ==============================================================================
# 23. SAVE CANDIDATE FLY REGIONS
# ==============================================================================

write_csv(
  candidate_fly_rcls,
  file.path(
    table_dir,
    "03_candidate_fly_RCL_regions.csv"
  )
)



# ==============================================================================
# 24. COMBINE VALIDATED AND CANDIDATE RCL REGIONS
# ==============================================================================

validated_for_alignment <- validated_rcls %>%

  select(
    Gene,
    RCL_Start,
    RCL_End,
    RCL_Length,
    RCL_Sequence,
    RCL_Definition,
    Site_Status
  )


candidate_for_alignment <- candidate_fly_rcls %>%

  select(
    Gene,
    RCL_Start,
    RCL_End,
    RCL_Length,
    RCL_Sequence,
    RCL_Definition,
    Site_Status
  )


rcl_table <- bind_rows(

  validated_for_alignment,

  candidate_for_alignment

)



# ==============================================================================
# 25. SET FIGURE / TABLE ORDER
# ==============================================================================

figure_order <- c(

  "SERPINE1",

  "SERPINB3",

  "SERPINB4",

  "SERPINB11",

  "SERPINI1",

  "Spn42Da_Q8T0M5",

  "Spn42Da_B_A1Z6R1",

  "Spn42Db",

  "Spn42Dc",

  "Spn42Dd",

  "Spn28B",

  "Spn28F"

)


rcl_table <- rcl_table %>%

  mutate(

    Gene = factor(
      Gene,
      levels = figure_order
    )

  ) %>%

  arrange(
    Gene
  ) %>%

  mutate(
    Gene = as.character(Gene)
  )



# ==============================================================================
# 26. ADD ACCESSION / DISPLAY INFORMATION
# ==============================================================================

rcl_table <- rcl_table %>%

  left_join(

    accession_key %>%
      select(
        Gene,
        Accession,
        Display_Name,
        Species
      ),

    by = "Gene"

  ) %>%

  select(
    Gene,
    Display_Name,
    Accession,
    Species,
    everything()
  )



# ==============================================================================
# 27. SAVE MASTER RCL TABLE
# ==============================================================================

write_csv(
  rcl_table,
  file.path(
    processed_dir,
    "03_RCL_sequences_for_alignment.csv"
  )
)


cat(
  "\nMaster RCL sequence table:\n\n"
)


print(
  rcl_table
)



# ==============================================================================
# 28. BUILD RCL AAStringSet
# ==============================================================================

rcl_set <- AAStringSet(

  setNames(

    rcl_table$RCL_Sequence,

    rcl_table$Gene

  )

)



# ==============================================================================
# 29. WRITE UNALIGNED RCL FASTA
# ==============================================================================

unaligned_fasta_file <- file.path(
  processed_dir,
  "03_RCL_sequences_unaligned.fasta"
)


writeXStringSet(
  rcl_set,
  filepath = unaligned_fasta_file,
  format = "fasta"
)


cat(
  "\nUnaligned RCL FASTA written to:\n",
  unaligned_fasta_file,
  "\n"
)



# ==============================================================================
# 30. ALIGN RCL REGIONS USING CLUSTAL OMEGA
# ==============================================================================

cat(
  "\nRunning Clustal Omega RCL alignment...\n"
)


rcl_alignment <- msa(

  rcl_set,

  method = "ClustalOmega",

  order = "input"

)


cat(
  "\nRCL alignment completed.\n\n"
)


print(
  rcl_alignment
)



# ==============================================================================
# 31. CONVERT ALIGNMENT TO AAStringSet
# ==============================================================================

aligned_rcl_set <- as(
  rcl_alignment,
  "AAStringSet"
)



# ==============================================================================
# 32. SAVE ALIGNED FASTA
# ==============================================================================

aligned_fasta_file <- file.path(
  processed_dir,
  "03_RCL_alignment_ClustalOmega.fasta"
)


writeXStringSet(
  aligned_rcl_set,
  filepath = aligned_fasta_file,
  format = "fasta"
)


cat(
  "\nAligned RCL FASTA written to:\n",
  aligned_fasta_file,
  "\n"
)



# ==============================================================================
# 33. CREATE HUMAN-READABLE ALIGNMENT TABLE
# ==============================================================================

aligned_strings <- as.character(
  aligned_rcl_set
)


alignment_table <- tibble(

  Gene = names(
    aligned_strings
  ),

  Aligned_RCL = unname(
    aligned_strings
  ),

  Alignment_Length = nchar(
    unname(aligned_strings)
  )

)


write_csv(
  alignment_table,
  file.path(
    processed_dir,
    "03_RCL_alignment_table.csv"
  )
)


cat(
  "\nAligned sequences:\n\n"
)


print(
  alignment_table
)



# ==============================================================================
# 34. LOCATE SERPINE1 P1/P1' IN THE ALIGNMENT
# ==============================================================================
#
# SERPINE1 standardized window:
#
#   P17-P8'
#
# Within the UNALIGNED 25-aa SERPINE1 RCL:
#
#   residue 17 = R369 = P1
#   residue 18 = M370 = P1'
#
# Clustal Omega may insert gaps.
#
# Therefore, we map residues 17 and 18 to their actual alignment columns.
#
# ==============================================================================

aligned_pai1 <- aligned_strings[
  "SERPINE1"
]


aligned_pai1_chars <- strsplit(
  aligned_pai1,
  ""
)[[1]]


pai1_nongap_columns <- which(
  aligned_pai1_chars != "-"
)


if (
  length(
    pai1_nongap_columns
  ) < 18
) {

  stop(
    "\nSERPINE1 alignment unexpectedly contains fewer than 18 non-gap residues."
  )

}


pai1_p1_alignment_column <-
  pai1_nongap_columns[17]


pai1_p1prime_alignment_column <-
  pai1_nongap_columns[18]



# ==============================================================================
# 35. VALIDATE SERPINE1 ALIGNMENT MAPPING
# ==============================================================================

observed_aligned_p1 <-
  aligned_pai1_chars[
    pai1_p1_alignment_column
  ]


observed_aligned_p1prime <-
  aligned_pai1_chars[
    pai1_p1prime_alignment_column
  ]


if (
  observed_aligned_p1 != "R" |
  observed_aligned_p1prime != "M"
) {

  stop(
    "\nSERPINE1 P1/P1' alignment mapping failed.\n",
    "Expected R/M at the mapped alignment columns."
  )

}



cat(
  "\nSERPINE1 alignment reference:\n",
  "P1 R369 alignment column = ",
  pai1_p1_alignment_column,
  "\n",
  "P1' M370 alignment column = ",
  pai1_p1prime_alignment_column,
  "\n\n",
  sep = ""
)



# ==============================================================================
# 36. SAVE SERPINE1 ALIGNMENT REFERENCE
# ==============================================================================

alignment_reference <- tibble(

  Reference_Gene =
    "SERPINE1",

  P1_Position_Full_Length =
    369,

  P1_Residue =
    "R",

  P1_Alignment_Column =
    pai1_p1_alignment_column,

  P1prime_Position_Full_Length =
    370,

  P1prime_Residue =
    "M",

  P1prime_Alignment_Column =
    pai1_p1prime_alignment_column

)


write_csv(
  alignment_reference,
  file.path(
    qc_dir,
    "03_SERPINE1_alignment_reference.csv"
  )
)



# ==============================================================================
# 37. OPTIONAL: LOCATE SPN42DA-B R342/A343 IN ALIGNMENT
# ==============================================================================
#
# # Spn42Da-B has an experimentally supported R-A reactive bond.
# Literature numbering: R342-A343.
# Full-length A1Z6R1 FASTA numbering used here: R374-A375.
#
# Its standardized P17-P8' window also places:
#
#   P1  at residue 17
#   P1' at residue 18
#
# This is recorded for QC but is NOT used to define the dashed boundaries
# in the primary figure.
#
# ==============================================================================

aligned_spn42da_b <- aligned_strings[
  "Spn42Da_B_A1Z6R1"
]


aligned_spn42da_b_chars <- strsplit(
  aligned_spn42da_b,
  ""
)[[1]]


spn42da_b_nongap_columns <- which(
  aligned_spn42da_b_chars != "-"
)


spn42da_b_p1_alignment_column <-
  spn42da_b_nongap_columns[17]


spn42da_b_p1prime_alignment_column <-
  spn42da_b_nongap_columns[18]


spn42da_b_reference <- tibble(

  Gene =
    "Spn42Da_B_A1Z6R1",

  P1_Position_Full_Length =
    374,

  P1_Residue =
    "R",

  P1_Alignment_Column =
    spn42da_b_p1_alignment_column,

  P1prime_Position_Full_Length =
    375,

  P1prime_Residue =
    "A",

  P1prime_Alignment_Column =
    spn42da_b_p1prime_alignment_column

)

write_csv(
  spn42da_b_reference,
  file.path(
    qc_dir,
    "03_Spn42DaB_alignment_reference.csv"
  )
)



# ==============================================================================
# 38. CREATE RCL ALIGNMENT FIGURE
# ==============================================================================
#
# Dashed boundaries surround the alignment columns corresponding to the
# SERPINE1 R369-M370 P1/P1' reactive bond.
#
# The lines are intentionally described as SERPINE1 alignment references.
#
# ==============================================================================

max_alignment_length <- max(
  width(
    aligned_rcl_set
  )
)


p_rcl <- ggmsa(

  aligned_rcl_set,

  start = 1,

  end = max_alignment_length,

  seq_name = TRUE,

  char_width = 0.5

) +

  geom_vline(

    xintercept =
      pai1_p1_alignment_column - 0.5,

    linetype = "dashed",

    linewidth = 0.6

  ) +

  geom_vline(

    xintercept =
      pai1_p1prime_alignment_column + 0.5,

    linetype = "dashed",

    linewidth = 0.6

  ) +

  labs(

    title =
      paste0(
        "Reactive center loop alignment of human ",
        "and Drosophila serpins"
      ),

    subtitle =
      paste0(
        "Dashed boundaries identify the SERPINE1 ",
        "P1-P1' alignment positions"
      ),

    x =
      "Multiple-sequence alignment position",

    y =
      NULL

  ) +

  theme(

    plot.title = element_text(
      face = "bold",
      size = 14
    ),

    plot.subtitle = element_text(
      size = 10
    )

  )



# ==============================================================================
# 39. DISPLAY FIGURE
# ==============================================================================

print(
  p_rcl
)



# ==============================================================================
# 40. SAVE FIGURE AS PNG
# ==============================================================================

ggsave(

  filename = file.path(
    figure_dir,
    "03_RCL_alignment_human_Drosophila.png"
  ),

  plot = p_rcl,

  width = 12,

  height = 8,

  dpi = 600

)



# ==============================================================================
# 41. SAVE FIGURE AS PDF
# ==============================================================================

ggsave(

  filename = file.path(
    figure_dir,
    "03_RCL_alignment_human_Drosophila.pdf"
  ),

  plot = p_rcl,

  width = 12,

  height = 8

)



# ==============================================================================
# 42. CREATE ALIGNMENT QC TABLE
# ==============================================================================

alignment_qc <- tibble(

  Gene = names(
    aligned_rcl_set
  ),

  Original_RCL_Length =
    width(
      rcl_set[
        names(aligned_rcl_set)
      ]
    ),

  Aligned_Length =
    width(
      aligned_rcl_set
    ),

  Number_of_Gaps =
    str_count(
      as.character(
        aligned_rcl_set
      ),
      fixed("-")
    )

)


write_csv(
  alignment_qc,
  file.path(
    qc_dir,
    "03_RCL_alignment_QC.csv"
  )
)


cat(
  "\nAlignment QC:\n\n"
)


print(
  alignment_qc
)



# ==============================================================================
# 43. DO NOT CALCULATE FINAL RCL SIMILARITY YET
# ==============================================================================
#
# IMPORTANT:
#
# The validated human sequences and Spn42Da-B are standardized P17-P8'
# windows.
#
# The other fly sequences currently represent broader candidate
# RCL-containing regions because their exact P1/P1' sites were intentionally
# not assigned.
#
# Therefore, direct percentage identity/similarity values between these
# unequal biological definitions should NOT yet be incorporated into the
# study's final quantitative structure/sequence score.
#
# The current alignment should first be inspected to determine whether
# equivalent structural RCL boundaries can be established for the remaining
# fly proteins without relying on unsupported predicted cleavage sites.
#
# ==============================================================================



# ==============================================================================
# 44. SAVE EXTRACTION STATUS TABLE
# ==============================================================================

extraction_status <- rcl_table %>%

  transmute(

    Gene,

    Display_Name,

    Accession,

    Species,

    RCL_Start,

    RCL_End,

    RCL_Length,

    RCL_Definition,

    Site_Status

  )


write_csv(
  extraction_status,
  file.path(
    table_dir,
    "03_RCL_extraction_status.csv"
  )
)



# ==============================================================================
# 45. SAVE SESSION INFORMATION
# ==============================================================================

capture.output(

  sessionInfo(),

  file = file.path(
    qc_dir,
    "03_sessionInfo.txt"
  )

)



# ==============================================================================
# 46. COMPLETION SUMMARY
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "03_RCL_Domain.R COMPLETED SUCCESSFULLY\n",
  "============================================================\n\n",
  sep = ""
)


cat(
  "Input:\n",
  "  ",
  fasta_file,
  "\n\n",
  sep = ""
)


cat(
  "Proteins analyzed: ",
  length(seqs),
  "\n\n",
  sep = ""
)


cat(
  "Validated P1/P1'-centered RCLs:\n",
  "  SERPINE1\n",
  "  SERPINB3\n",
  "  SERPINB4\n",
  "  SERPINB11\n",
  "  SERPINI1\n",
  "  Spn42Da-B (A1Z6R1)\n\n",
  sep = ""
)


cat(
  "Fly proteins retained without predicted P1/P1' assignment:\n",
  "  Spn42Da (Q8T0M5)\n",
  "  Spn42Db\n",
  "  Spn42Dc\n",
  "  Spn42Dd\n",
  "  Spn28B\n",
  "  Spn28F\n\n",
  sep = ""
)


cat(
  "SERPINE1 P1/P1' alignment reference:\n",
  "  R369 = alignment column ",
  pai1_p1_alignment_column,
  "\n",
  "  M370 = alignment column ",
  pai1_p1prime_alignment_column,
  "\n\n",
  sep = ""
)


cat(
  "Main outputs:\n",
  "  Processed/03_RCL_sequences_for_alignment.csv\n",
  "  Processed/03_RCL_sequences_unaligned.fasta\n",
  "  Processed/03_RCL_alignment_ClustalOmega.fasta\n",
  "  Processed/03_RCL_alignment_table.csv\n",
  "  Tables/03_validated_standardized_RCL_sequences.csv\n",
  "  Tables/03_candidate_fly_RCL_regions.csv\n",
  "  Tables/03_RCL_extraction_status.csv\n",
  "  Figures/03_RCL_alignment_human_Drosophila.png\n",
  "  Figures/03_RCL_alignment_human_Drosophila.pdf\n",
  "  QC/03_reactive_site_validation.csv\n",
  "  QC/03_RCL_alignment_QC.csv\n",
  "  QC/03_SERPINE1_alignment_reference.csv\n",
  "  QC/03_Spn42DaB_alignment_reference.csv\n",
  "  QC/03_sessionInfo.txt\n\n",
  sep = ""
)


cat(
  "IMPORTANT:\n",
  "No final RCL identity/similarity score was calculated for the\n",
  "integrated analysis because the remaining Drosophila RCL boundaries\n",
  "have not yet been standardized independently of predicted P1/P1' sites.\n\n",
  sep = ""
)


cat(
  "============================================================\n"
)
