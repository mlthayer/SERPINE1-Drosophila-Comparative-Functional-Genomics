# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 03B_RCL_ExtractedSequence_Sensitivity.R
# Run from the repository root.
# ==============================================================================
# 03B_RCL_EXTRACTED_SEQUENCE SENSITIVITY ANALYSIS
#
# PURPOSE
# -------
# Reconstruct the PREVIOUS RCL analysis in which RCL sequences were first
# extracted independently from each serpin and the isolated RCL sequences were
# subsequently aligned using Clustal Omega.
#
# This analysis is retained as a sensitivity/alternative analysis and is
# compared with the revised primary analysis:
#
#   PRIMARY:
#   Full-length MSA -> map SERPINE1 P17-P8' -> extract homologous region
#
#   SENSITIVITY:
#   Previously extracted RCLs -> align isolated RCL sequences
#
# IMPORTANT
# ---------
# The historical RCL sequences are reproduced exactly as reported in the
# previous study analysis. This script DOES NOT reinterpret or redefine
# their biological boundaries.
#
# Predicted Drosophila P1/P1' sites are NOT used for scoring.
#
# Identity and biochemical similarity are recalculated using explicit,
# reproducible scoring rules:
#
#   Identity:
#       exact amino-acid match to SERPINE1
#
#   Similarity:
#       exact match OR BLOSUM62 substitution score > 0
#
# Candidate gaps opposite SERPINE1 residues are counted as nonmatches.
#
# ==============================================================================


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

required_packages <- c(
  "Biostrings",
  "msa",
  "dplyr",
  "tidyr",
  "stringr",
  "readr",
  "tibble",
  "ggplot2"
)


for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    stop(
      "Required package is not installed: ",
      pkg
    )

  }

}


library(Biostrings)
library(msa)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(tibble)
library(ggplot2)


# ==============================================================================
# 1. DIRECTORY STRUCTURE
# ==============================================================================

step3_dir <- "Sequence/RCL"


sensitivity_dir <- file.path(
  step3_dir,
  "Extracted_RCL_Sensitivity_Analysis"
)


processed_dir <- file.path(
  sensitivity_dir,
  "Processed"
)


table_dir <- file.path(
  sensitivity_dir,
  "Tables"
)


figure_dir <- file.path(
  sensitivity_dir,
  "Figures"
)


qc_dir <- file.path(
  sensitivity_dir,
  "QC"
)


dirs_to_create <- c(
  sensitivity_dir,
  processed_dir,
  table_dir,
  figure_dir,
  qc_dir
)


for (dir_i in dirs_to_create) {

  if (!dir.exists(dir_i)) {

    dir.create(
      dir_i,
      recursive = TRUE
    )

  }

}


# ==============================================================================
# 2. HISTORICAL EXTRACTED RCL SEQUENCES
# ==============================================================================
#
# These sequences reproduce the previously reported study RCL dataset.
#
# IMPORTANT:
# These are NOT being redefined as experimentally validated RCL boundaries.
# They are retained specifically to reproduce the historical analytical
# approach.
#
# ==============================================================================

historical_rcl <- tibble(

  Gene = c(
    "SERPINE1",
    "Spn42Da",
    "Spn42Db",
    "Spn42Dc",
    "Spn42Dd",
    "Spn28B",
    "Spn28F"
  ),

  Accession = c(
    "P05121",
    "Q8T0M5",
    "A1Z6R3",
    "A1Z6R4",
    "Q7YTY6",
    "Q9VLZ8",
    "Q9VLQ7"
  ),

  Historical_RCL = c(
    "GTVASSSTAVIVSARMAPEEIIMDRP",
    "SPEPLKVSAIIHKAFIEVNEE",
    "PKVFHANRPFFYAIKDNTHGL",
    "IVPMMLNMNKKLFKADHPFVFY",
    "STFLMADHPFAFVIRDANTIYFQGR",
    "GGEASAATGVLTRRKK",
    "AAAATALLFVRYSMPMPSSQMVFNA"
  )

) %>%

  mutate(

    Historical_RCL_Length =
      nchar(Historical_RCL)

  )


print(
  historical_rcl,
  n = Inf
)


# ==============================================================================
# 3. EXPECTED LENGTH QC
# ==============================================================================

expected_lengths <- c(
  SERPINE1 = 26,
  Spn42Da = 21,
  Spn42Db = 21,
  Spn42Dc = 22,
  Spn42Dd = 25,
  Spn28B = 16,
  Spn28F = 25
)


historical_rcl <- historical_rcl %>%

  mutate(

    Expected_Length =
      expected_lengths[Gene],

    Length_QC =
      Historical_RCL_Length == Expected_Length

  )


if (!all(historical_rcl$Length_QC)) {

  stop(
    "Historical RCL sequence length QC failed."
  )

}


write_csv(
  historical_rcl,
  file.path(
    qc_dir,
    "03B_Historical_RCL_Input_QC.csv"
  )
)


# ==============================================================================
# 4. CREATE AAStringSet
# ==============================================================================

historical_rcl_set <- AAStringSet(
  historical_rcl$Historical_RCL
)


names(historical_rcl_set) <-
  historical_rcl$Gene


# ==============================================================================
# 5. SAVE INPUT FASTA
# ==============================================================================

writeXStringSet(

  historical_rcl_set,

  filepath = file.path(
    processed_dir,
    "03B_Historical_Extracted_RCL_Input.fasta"
  ),

  format = "fasta"

)


# ==============================================================================
# 6. ALIGN THE EXTRACTED RCL SEQUENCES
# ==============================================================================
#
# This reproduces the defining feature of the previous approach:
#
#       Extract RCL first
#               |
#               v
#       Align isolated RCLs
#
# rather than:
#
#       Full-length alignment
#               |
#               v
#       Map SERPINE1 RCL
#
# ==============================================================================

historical_alignment <- msa(

  historical_rcl_set,

  method = "ClustalOmega",

  order = "input"

)


print(
  historical_alignment
)


# ==============================================================================
# 7. CONVERT ALIGNMENT TO AAStringSet
# ==============================================================================

historical_alignment_set <-
  as(
    historical_alignment,
    "AAStringSet"
  )


# Restore expected order explicitly.

historical_alignment_set <-
  historical_alignment_set[
    historical_rcl$Gene
  ]


# ==============================================================================
# 8. SAVE ALIGNED FASTA
# ==============================================================================

writeXStringSet(

  historical_alignment_set,

  filepath = file.path(
    processed_dir,
    "03B_Historical_Extracted_RCL_ClustalOmega_Alignment.fasta"
  ),

  format = "fasta"

)


# ==============================================================================
# 9. CONVERT ALIGNMENT TO TABLE
# ==============================================================================

alignment_strings <-
  as.character(
    historical_alignment_set
  )


historical_alignment_table <- tibble(

  Gene =
    names(
      historical_alignment_set
    ),

  Aligned_Historical_RCL =
    alignment_strings,

  Alignment_Length =
    nchar(
      alignment_strings
    ),

  Ungapped_Length =
    nchar(
      gsub(
        "-",
        "",
        alignment_strings
      )
    )

)


print(
  historical_alignment_table,
  n = Inf,
  width = Inf
)


write_csv(

  historical_alignment_table,

  file.path(
    processed_dir,
    "03B_Historical_Extracted_RCL_Alignment.csv"
  )

)


# ==============================================================================
# 10. ALIGNMENT LENGTH QC
# ==============================================================================

alignment_lengths <-
  unique(
    historical_alignment_table$Alignment_Length
  )


if (length(alignment_lengths) != 1) {

  stop(
    "Aligned historical RCL sequences do not have equal alignment lengths."
  )

}


historical_alignment_length <-
  alignment_lengths[[1]]


message(
  "Historical extracted-RCL alignment length = ",
  historical_alignment_length
)


# ==============================================================================
# 11. BLOSUM62 MATRIX
# ==============================================================================

data(
  "BLOSUM62",
  package = "Biostrings"
)


blosum62_matrix <-
  BLOSUM62


# ==============================================================================
# 12. IDENTIFY SERPINE1 REFERENCE COLUMNS
# ==============================================================================
#
# Unlike the revised primary analysis, the historical sequence contains
# 26 SERPINE1 residues:
#
# GTVASSSTAVIVSARMAPEEIIMDRP
#
# Identity and similarity are therefore evaluated across the alignment
# positions occupied by those 26 SERPINE1 residues.
#
# ==============================================================================

serpine1_aligned <-
  alignment_strings[
    "SERPINE1"
  ]


serpine1_chars <-
  strsplit(
    serpine1_aligned,
    ""
  )[[1]]


serpine1_reference_columns <-
  which(
    serpine1_chars != "-"
  )


historical_reference_n <-
  length(
    serpine1_reference_columns
  )


message(
  "Number of SERPINE1 residues used as historical reference denominator = ",
  historical_reference_n
)


if (historical_reference_n != 26) {

  stop(
    paste0(
      "Expected 26 SERPINE1 residues in historical RCL, but found ",
      historical_reference_n,
      "."
    )
  )

}


# ==============================================================================
# 13. SCORING FUNCTION
# ==============================================================================

score_historical_rcl <- function(
  candidate_gene,
  aligned_sequences,
  reference_gene = "SERPINE1"
) {


  reference_string <-
    aligned_sequences[
      reference_gene
    ]


  candidate_string <-
    aligned_sequences[
      candidate_gene
    ]


  reference_chars <-
    strsplit(
      reference_string,
      ""
    )[[1]]


  candidate_chars <-
    strsplit(
      candidate_string,
      ""
    )[[1]]


  if (
    length(reference_chars) !=
    length(candidate_chars)
  ) {

    stop(
      "Alignment length mismatch for ",
      candidate_gene
    )

  }


  # ---------------------------------------------------------------------------
  # Only columns containing an actual SERPINE1 residue contribute to the
  # identity/similarity denominator.
  # ---------------------------------------------------------------------------

  reference_positions <-
    which(
      reference_chars != "-"
    )


  ref <-
    reference_chars[
      reference_positions
    ]


  cand <-
    candidate_chars[
      reference_positions
    ]


  denominator <-
    length(ref)


  # ---------------------------------------------------------------------------
  # Identity
  # ---------------------------------------------------------------------------

  identical_positions <-
    cand != "-" &
    cand == ref


  identical_n <-
    sum(
      identical_positions
    )


  identity_percent <-
    100 *
    identical_n /
    denominator


  # ---------------------------------------------------------------------------
  # BLOSUM62 similarity
  # ---------------------------------------------------------------------------

  similar_positions <-
    logical(
      denominator
    )


  blosum_scores <-
    rep(
      NA_real_,
      denominator
    )


  for (i in seq_len(denominator)) {


    ref_aa <-
      ref[i]


    cand_aa <-
      cand[i]


    if (cand_aa == "-") {

      similar_positions[i] <-
        FALSE

      blosum_scores[i] <-
        NA_real_

      next

    }


    if (
      !(ref_aa %in% rownames(blosum62_matrix)) ||
      !(cand_aa %in% colnames(blosum62_matrix))
    ) {

      stop(
        "Unknown amino-acid character while scoring ",
        candidate_gene,
        " at reference position ",
        i,
        "."
      )

    }


    score_i <-
      blosum62_matrix[
        ref_aa,
        cand_aa
      ]


    blosum_scores[i] <-
      score_i


    similar_positions[i] <-
      (
        ref_aa == cand_aa
      ) ||
      (
        score_i > 0
      )


  }


  similar_n <-
    sum(
      similar_positions
    )


  similarity_percent <-
    100 *
    similar_n /
    denominator


  # ---------------------------------------------------------------------------
  # Candidate gaps opposite SERPINE1 residues
  # ---------------------------------------------------------------------------

  candidate_gaps_at_reference <-
    sum(
      cand == "-"
    )


  # ---------------------------------------------------------------------------
  # Candidate residues opposite SERPINE1 gaps
  # ---------------------------------------------------------------------------

  reference_gap_columns <-
    which(
      reference_chars == "-"
    )


  candidate_insertions <-
    if (
      length(reference_gap_columns) == 0
    ) {

      0L

    } else {

      sum(
        candidate_chars[
          reference_gap_columns
        ] != "-"
      )

    }


  # ---------------------------------------------------------------------------
  # Mean BLOSUM62 score
  #
  # Descriptive QC only.
  # Not used as the similarity percentage.
  # ---------------------------------------------------------------------------

  mean_blosum <-
    if (
      all(
        is.na(
          blosum_scores
        )
      )
    ) {

      NA_real_

    } else {

      mean(
        blosum_scores,
        na.rm = TRUE
      )

    }


  tibble(

    Gene =
      candidate_gene,

    SERPINE1_Reference_Residues =
      denominator,

    Identical_Residues =
      identical_n,

    Historical_RCL_Identity_Percent =
      identity_percent,

    Similar_Residues_BLOSUM62_Positive =
      similar_n,

    Historical_RCL_Similarity_Percent =
      similarity_percent,

    Candidate_Gaps_at_SERPINE1_Residues =
      candidate_gaps_at_reference,

    Candidate_Residues_at_SERPINE1_Gap_Columns =
      candidate_insertions,

    Mean_BLOSUM62_NonGap =
      mean_blosum

  )

}


# ==============================================================================
# 14. SCORE ALL SEQUENCES
# ==============================================================================

historical_scores <-
  bind_rows(

    lapply(

      names(
        historical_alignment_set
      ),

      score_historical_rcl,

      aligned_sequences =
        alignment_strings

    )

  )


# ==============================================================================
# 15. ADD ALIGNMENT INFORMATION
# ==============================================================================

historical_scores <-
  historical_scores %>%

  left_join(

    historical_alignment_table %>%
      select(
        Gene,
        Aligned_Historical_RCL,
        Alignment_Length,
        Ungapped_Length
      ),

    by =
      "Gene"

  ) %>%

  left_join(

    historical_rcl %>%
      select(
        Gene,
        Accession,
        Historical_RCL,
        Historical_RCL_Length
      ),

    by =
      "Gene"

  ) %>%

  relocate(
    Gene,
    Accession,
    Historical_RCL,
    Aligned_Historical_RCL,
    Historical_RCL_Length,
    Alignment_Length,
    Ungapped_Length
  )


# ==============================================================================
# 16. PRINT COMPLETE RESULTS
# ==============================================================================

historical_scores %>%

  select(
    Gene,
    Accession,
    Historical_RCL,
    Aligned_Historical_RCL,
    SERPINE1_Reference_Residues,
    Identical_Residues,
    Historical_RCL_Identity_Percent,
    Similar_Residues_BLOSUM62_Positive,
    Historical_RCL_Similarity_Percent,
    Candidate_Gaps_at_SERPINE1_Residues,
    Candidate_Residues_at_SERPINE1_Gap_Columns
  ) %>%

  print(
    n = Inf,
    width = Inf
  )


# ==============================================================================
# 17. SAVE COMPLETE RESULTS
# ==============================================================================

write_csv(

  historical_scores,

  file.path(
    processed_dir,
    "03B_Historical_Extracted_RCL_Scores_vs_SERPINE1.csv"
  )

)


# ==============================================================================
# 18. PRIMARY SIX DROSOPHILA CANDIDATES
# ==============================================================================

primary_candidates <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)


historical_primary_scores <-
  historical_scores %>%

  filter(
    Gene %in%
      primary_candidates
  ) %>%

  mutate(

    Gene =
      factor(
        Gene,
        levels =
          primary_candidates
      )

  ) %>%

  arrange(
    Gene
  )


print(
  historical_primary_scores,
  n = Inf,
  width = Inf
)


write_csv(

  historical_primary_scores,

  file.path(
    table_dir,
    "03B_Historical_RCL_Scores_Primary_Candidates.csv"
  )

)


# ==============================================================================
# 19. RANK CANDIDATES
# ==============================================================================

historical_ranking <-
  historical_primary_scores %>%

  mutate(

    Historical_Identity_Rank =
      rank(
        -Historical_RCL_Identity_Percent,
        ties.method = "min"
      ),

    Historical_Similarity_Rank =
      rank(
        -Historical_RCL_Similarity_Percent,
        ties.method = "min"
      )

  ) %>%

  select(
    Gene,
    Historical_RCL_Identity_Percent,
    Historical_Identity_Rank,
    Historical_RCL_Similarity_Percent,
    Historical_Similarity_Rank
  )


print(
  historical_ranking,
  n = Inf
)


write_csv(

  historical_ranking,

  file.path(
    table_dir,
    "03B_Historical_RCL_Candidate_Ranking.csv"
  )

)


# ==============================================================================
# 20. CREATE LONG-FORM ALIGNMENT DATA FOR FIGURE
# ==============================================================================

historical_plot_df <-
  historical_alignment_table %>%

  select(
    Gene,
    Aligned_Historical_RCL
  ) %>%

  rowwise() %>%

  mutate(

    Residue_List =
      list(
        strsplit(
          Aligned_Historical_RCL,
          ""
        )[[1]]
      )

  ) %>%

  ungroup() %>%

  select(
    Gene,
    Residue_List
  ) %>%

  unnest_longer(
    Residue_List,
    values_to = "Residue"
  ) %>%

  group_by(
    Gene
  ) %>%

  mutate(
    Alignment_Position =
      row_number()
  ) %>%

  ungroup()


# ==============================================================================
# 21. AMINO-ACID CLASSES
# ==============================================================================
#
# These categories can be changed to exactly match the final primary RCL figure.
#
# ==============================================================================

historical_plot_df <-
  historical_plot_df %>%

  mutate(

    AA_Class =
      case_when(

        Residue %in% c(
          "D",
          "E"
        ) ~
          "Acidic",

        Residue %in% c(
          "R",
          "K",
          "H"
        ) ~
          "Basic",

        Residue %in% c(
          "A",
          "V",
          "L",
          "I",
          "M",
          "F",
          "W",
          "Y"
        ) ~
          "Hydrophobic",

        Residue %in% c(
          "S",
          "T",
          "N",
          "Q"
        ) ~
          "Polar",

        Residue ==
          "G" ~
          "Glycine",

        Residue ==
          "P" ~
          "Proline",

        Residue ==
          "C" ~
          "Cysteine",

        Residue ==
          "-" ~
          "Gap",

        TRUE ~
          "Other"

      )

  )


# ==============================================================================
# 22. ROW ORDER
# ==============================================================================

plot_order <-
  rev(
    c(
      "SERPINE1",
      "Spn42Da",
      "Spn42Db",
      "Spn42Dc",
      "Spn42Dd",
      "Spn28B",
      "Spn28F"
    )
  )


historical_plot_df <-
  historical_plot_df %>%

  mutate(

    Gene =
      factor(
        Gene,
        levels =
          plot_order
      )

  )


# ==============================================================================
# 23. COLOR PALETTE
# ==============================================================================

aa_colors <- c(

  "Acidic" =
    "#F8766D",

  "Basic" =
    "#F5A142",

  "Hydrophobic" =
    "#FFD84D",

  "Polar" =
    "#74C98B",

  "Glycine" =
    "#7FA2D6",

  "Proline" =
    "#B77CCB",

  "Cysteine" =
    "#9E9E9E",

  "Gap" =
    "white",

  "Other" =
    "#D9D9D9"

)


# ==============================================================================
# 24. HISTORICAL ALIGNMENT FIGURE
# ==============================================================================

historical_alignment_plot <-
  ggplot(

    historical_plot_df,

    aes(
      x = Alignment_Position,
      y = Gene
    )

  ) +

  geom_tile(

    aes(
      fill = AA_Class
    ),

    color =
      "white",

    linewidth =
      0.35

  ) +

  geom_text(

    aes(
      label = Residue
    ),

    size =
      4

  ) +

  scale_fill_manual(

    values =
      aa_colors,

    name =
      "Amino acid class"

  ) +

  scale_x_continuous(

    breaks =
      seq(
        1,
        historical_alignment_length,
        by = 5
      ),

    expand =
      c(
        0,
        0
      )

  ) +

  labs(

    title =
      "Alignment of previously extracted reactive center loop sequences",

    subtitle =
      paste0(
        "Sensitivity analysis reproducing the previous workflow: ",
        "RCL extraction followed by Clustal Omega alignment"
      ),

    x =
      "Extracted-RCL alignment position",

    y =
      NULL

  ) +

  theme_minimal(
    base_size = 12
  ) +

  theme(

    panel.grid =
      element_blank(),

    axis.text.y =
      element_text(
        color = "black",
        size = 11
      ),

    axis.text.x =
      element_text(
        color = "black"
      ),

    plot.title =
      element_text(
        face = "bold"
      ),

    legend.position =
      "bottom"

  )


print(
  historical_alignment_plot
)


ggsave(

  filename = file.path(
    figure_dir,
    "03B_Historical_Extracted_RCL_Alignment.png"
  ),

  plot =
    historical_alignment_plot,

  width =
    13,

  height =
    6.5,

  dpi =
    600,

  bg =
    "white"

)


# ==============================================================================
# 25. PREPARE CURRENT-APPROACH VALUES FOR COMPARISON
# ==============================================================================
#
# These are the audited values generated by the primary SERPINE1-referenced
# homologous-region analysis.
#
# They are included here so the historical approach can be compared directly
# with the current approach.
#
# ==============================================================================

current_scores <- tibble(

  Gene = c(
    "Spn42Da",
    "Spn42Db",
    "Spn42Dc",
    "Spn42Dd",
    "Spn28B",
    "Spn28F"
  ),

  Current_RCL_Identity_Percent = c(
    24,
    32,
    12,
    28,
    36,
    20
  ),

  Current_RCL_Similarity_Percent = c(
    56,
    56,
    44,
    52,
    56,
    40
  )

)


# ==============================================================================
# 26. BUILD METHOD-COMPARISON TABLE
# ==============================================================================

method_comparison <-
  historical_ranking %>%

  left_join(
    current_scores,
    by = "Gene"
  ) %>%

  mutate(

    Current_Identity_Rank =
      rank(
        -Current_RCL_Identity_Percent,
        ties.method = "min"
      ),

    Current_Similarity_Rank =
      rank(
        -Current_RCL_Similarity_Percent,
        ties.method = "min"
      ),


    Identity_Change_Percentage_Points =
      Current_RCL_Identity_Percent -
      Historical_RCL_Identity_Percent,


    Similarity_Change_Percentage_Points =
      Current_RCL_Similarity_Percent -
      Historical_RCL_Similarity_Percent,


    Identity_Rank_Change =
      Current_Identity_Rank -
      Historical_Identity_Rank,


    Similarity_Rank_Change =
      Current_Similarity_Rank -
      Historical_Similarity_Rank

  ) %>%

  select(
    Gene,

    Historical_RCL_Identity_Percent,
    Current_RCL_Identity_Percent,
    Identity_Change_Percentage_Points,

    Historical_Identity_Rank,
    Current_Identity_Rank,
    Identity_Rank_Change,

    Historical_RCL_Similarity_Percent,
    Current_RCL_Similarity_Percent,
    Similarity_Change_Percentage_Points,

    Historical_Similarity_Rank,
    Current_Similarity_Rank,
    Similarity_Rank_Change
  )


print(
  method_comparison,
  n = Inf,
  width = Inf
)


write_csv(

  method_comparison,

  file.path(
    table_dir,
    "03B_Historical_vs_Current_RCL_Method_Comparison.csv"
  )

)


# ==============================================================================
# 27. LONG-FORM METHOD COMPARISON
# ==============================================================================

comparison_long <-
  method_comparison %>%

  select(
    Gene,
    Historical_RCL_Identity_Percent,
    Current_RCL_Identity_Percent,
    Historical_RCL_Similarity_Percent,
    Current_RCL_Similarity_Percent
  ) %>%

  pivot_longer(

    cols =
      -Gene,

    names_to =
      "Metric_Method",

    values_to =
      "Percent"

  ) %>%

  mutate(

    Metric =
      if_else(
        str_detect(
          Metric_Method,
          "Identity"
        ),
        "Identity",
        "Similarity"
      ),

    Method =
      if_else(
        str_detect(
          Metric_Method,
          "Historical"
        ),
        "Extracted RCL first",
        "Full-length MSA first"
      )

  )


# ==============================================================================
# 28. COMPARISON FIGURE
# ==============================================================================

comparison_plot <-
  ggplot(

    comparison_long,

    aes(
      x = Gene,
      y = Percent,
      group = Method,
      shape = Method
    )

  ) +

  geom_point(
    size = 3
  ) +

  facet_wrap(

    ~Metric,

    ncol = 1

  ) +

  labs(

    title =
      "Effect of RCL alignment strategy on SERPINE1 sequence conservation",

    subtitle =
      paste0(
        "Comparison of independently extracted RCL alignment with ",
        "full-length-MSA-based homologous-region mapping"
      ),

    x =
      NULL,

    y =
      "SERPINE1-referenced conservation (%)",

    shape =
      "RCL analysis"

  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),

    plot.title =
      element_text(
        face = "bold"
      ),

    legend.position =
      "bottom"

  )


print(
  comparison_plot
)


ggsave(

  filename = file.path(
    figure_dir,
    "03B_Historical_vs_Current_RCL_Comparison.png"
  ),

  plot =
    comparison_plot,

  width =
    9,

  height =
    8,

  dpi =
    600,

  bg =
    "white"

)


# ==============================================================================
# 29. SESSION INFORMATION
# ==============================================================================

capture.output(

  sessionInfo(),

  file = file.path(
    qc_dir,
    "03B_sessionInfo.txt"
  )

)


# ==============================================================================
# 30. COMPLETION MESSAGE
# ==============================================================================

message(
  "\n============================================================"
)

message(
  "Historical extracted-RCL sensitivity analysis complete."
)

message(
  "Primary output:"
)

message(
  file.path(
    table_dir,
    "03B_Historical_vs_Current_RCL_Method_Comparison.csv"
  )
)

message(
  "============================================================\n"
)
