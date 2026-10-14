# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 04_Phylogenetic_Analysis.R
# Run from the repository root.
# ==============================================================================
# 04_PHYLOGENETIC_ANALYSIS.R
#
# REVISED PHYLOGENETIC ANALYSIS
#
# Analyses:
#
#   A. Full-length SERPINE1 + Drosophila
#      1. Neighbor Joining
#      2. Maximum Likelihood
#
#   B. SERPINE1-defined homologous RCL region + Drosophila
#      3. Neighbor Joining
#      4. Maximum Likelihood
#
#   C. Expanded human + Drosophila full-length comparison
#      5. Neighbor Joining
#      6. Maximum Likelihood
#
# Drosophila sequences:
#   Spn42Da          Q8T0M5
#   Spn42Da-B        A1Z6R1
#   Spn42Db          A1Z6R3
#   Spn42Dc          A1Z6R4
#   Spn42Dd          Q7YTY6
#   Spn28B           Q9VLZ8
#   Spn28F           Q9VLQ7
#
# Human sequences:
#   SERPINE1         P05121
#   SERPINB3         P29508
#   SERPINB4         P48594
#   SERPINB11        Q96P15
#   SERPINI1         Q99574
#
# IMPORTANT:
#   - Full-length proteins are aligned with Clustal Omega before tree inference.
#   - RCL phylogeny uses the REVISED SERPINE1-defined homologous regions.
#   - Historical independently extracted fly RCL regions are NOT used.
#   - JTT amino-acid substitution model is used for ML analysis.
#   - Bootstrap support is calculated for both NJ and ML.
#   - Trees are not interpreted as rooted evolutionary histories.
#   - Patristic distance to SERPINE1 is calculated separately.
#
# ==============================================================================


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

required_packages <- c(
  "Biostrings",
  "msa",
  "ape",
  "phangorn",
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
library(ape)
library(phangorn)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(tibble)
library(ggplot2)


set.seed(1)


# ==============================================================================
# 1. PROJECT PATHS
# ==============================================================================

step3_dir <- "Sequence/RCL"


# Master 12-protein FASTA used in the revised sequence/RCL analysis

master_fasta <- file.path(
  step3_dir,
  "selected_human_fly_Serpins_Fasta.txt"
)


# Revised SERPINE1-defined RCL-region table

rcl_file <- file.path(
  step3_dir,
  "Processed",
  "03_SERPINE1_defined_homologous_RCL_regions.csv"
)

# Reuse the supplied validated homologous-region table when generated input
# is absent. An explicit override takes precedence over both locations.
if (!file.exists(rcl_file)) {
  rcl_file <- file.path(step3_dir, "Reference_Uploads",
                        "03_SERPINE1_defined_homologous_RCL_regions.csv")
}
rcl_file <- Sys.getenv("SERPINE1_PHYLOGENY_RCL_INPUT", unset = rcl_file)



# Phylogeny output directory

phylo_dir <- file.path(
  dirname(step3_dir),
  "Phylogeny"
)


processed_dir <- file.path(
  phylo_dir,
  "Processed"
)

tree_dir <- file.path(
  phylo_dir,
  "Trees"
)

figure_dir <- file.path(
  phylo_dir,
  "Figures"
)

qc_dir <- file.path(
  phylo_dir,
  "QC"
)


for (dir_i in c(
  phylo_dir,
  processed_dir,
  tree_dir,
  figure_dir,
  qc_dir
)) {

  dir.create(
    dir_i,
    recursive = TRUE,
    showWarnings = FALSE
  )

}


# ==============================================================================
# 2. REQUIRED FILE QC
# ==============================================================================

required_files <- c(
  Master_FASTA = master_fasta,
  Revised_RCL = rcl_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {

  stop(
    paste0(
      "Missing required input file(s):\n",
      paste(
        names(missing_files),
        missing_files,
        sep = " = ",
        collapse = "\n"
      )
    )
  )

}


# ==============================================================================
# 3. ACCESSION / GENE KEY
# ==============================================================================

sequence_key <- tibble(

  Gene = c(
    "SERPINE1",
    "SERPINB3",
    "SERPINB4",
    "SERPINB11",
    "SERPINI1",
    "Spn42Da",
    "Spn42Da-B",
    "Spn42Db",
    "Spn42Dc",
    "Spn42Dd",
    "Spn28B",
    "Spn28F"
  ),

  Accession = c(
    "P05121",
    "P29508",
    "P48594",
    "Q96P15",
    "Q99574",
    "Q8T0M5",
    "A1Z6R1",
    "A1Z6R3",
    "A1Z6R4",
    "Q7YTY6",
    "Q9VLZ8",
    "Q9VLQ7"
  ),

  Species = c(
    rep("Human", 5),
    rep("Drosophila", 7)
  )

)


# ==============================================================================
# 4. ANALYSIS SETS
# ==============================================================================

fly_genes <- c(
  "Spn42Da",
  "Spn42Da-B",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)


serpine1_fly_set <- c(
  "SERPINE1",
  fly_genes
)


expanded_human_fly_set <- c(
  "SERPINE1",
  "SERPINB3",
  "SERPINB4",
  "SERPINB11",
  "SERPINI1",
  fly_genes
)


if (length(serpine1_fly_set) != 8) {
  stop("SERPINE1-fly analysis should contain exactly 8 proteins.")
}


if (length(expanded_human_fly_set) != 12) {
  stop("Expanded human-fly analysis should contain exactly 12 proteins.")
}


# ==============================================================================
# 5. READ MASTER FASTA
# ==============================================================================

master_sequences <- readAAStringSet(
  master_fasta,
  format = "fasta"
)


if (length(master_sequences) != 12) {

  warning(
    "Master FASTA contains ",
    length(master_sequences),
    " sequences; expected 12."
  )

}


# ==============================================================================
# 6. PARSE UNIPROT ACCESSIONS
# ==============================================================================

parse_accession <- function(header) {

  pieces <- strsplit(
    header,
    "\\|"
  )[[1]]

  if (length(pieces) >= 2) {

    pieces[2]

  } else {

    header

  }

}


fasta_accessions <- vapply(
  names(master_sequences),
  parse_accession,
  character(1)
)


# ==============================================================================
# 7. MAP FASTA TO GENE NAMES
# ==============================================================================

missing_accessions <- setdiff(
  sequence_key$Accession,
  fasta_accessions
)


if (length(missing_accessions) > 0) {

  stop(
    "Required accessions missing from master FASTA: ",
    paste(
      missing_accessions,
      collapse = ", "
    )
  )

}


if (anyDuplicated(fasta_accessions) > 0) {

  stop(
    "Duplicate UniProt accession(s) detected in master FASTA."
  )

}


ordered_sequences <- master_sequences[
  match(
    sequence_key$Accession,
    fasta_accessions
  )
]


names(ordered_sequences) <- sequence_key$Gene


# ==============================================================================
# 8. FULL-LENGTH SEQUENCE INVENTORY
# ==============================================================================

sequence_inventory <- sequence_key %>%

  mutate(
    Protein_Length =
      width(
        ordered_sequences[Gene]
      )
  )


print(
  sequence_inventory,
  n = Inf
)


write_csv(
  sequence_inventory,
  file.path(
    qc_dir,
    "04_sequence_inventory.csv"
  )
)


# ==============================================================================
# 9. HELPER: ALIGN FULL-LENGTH SEQUENCES
# ==============================================================================

align_proteins <- function(
  sequence_set,
  genes,
  analysis_name
) {

  missing_genes <- setdiff(
    genes,
    names(sequence_set)
  )


  if (length(missing_genes) > 0) {

    stop(
      analysis_name,
      " missing genes: ",
      paste(
        missing_genes,
        collapse = ", "
      )
    )

  }


  seq_subset <- sequence_set[
    genes
  ]


  cat(
    "\nRunning Clustal Omega for ",
    analysis_name,
    "...\n",
    sep = ""
  )


  alignment <- msa(
    seq_subset,
    method = "ClustalOmega",
    order = "input"
  )


  aligned_set <- as(
    alignment,
    "AAStringSet"
  )


  aligned_set <- aligned_set[
    genes
  ]


  alignment_lengths <- width(
    aligned_set
  )


  if (length(unique(alignment_lengths)) != 1) {

    stop(
      analysis_name,
      " alignment sequences have unequal lengths."
    )

  }


  cat(
    analysis_name,
    " alignment length = ",
    unique(alignment_lengths),
    "\n"
  )


  aligned_set

}


# ==============================================================================
# 10. FULL-LENGTH SERPINE1 + FLY ALIGNMENT
# ==============================================================================

full_serpine1_fly_alignment <- align_proteins(

  ordered_sequences,

  serpine1_fly_set,

  "Full-length SERPINE1 + Drosophila"

)


writeXStringSet(

  full_serpine1_fly_alignment,

  file.path(
    processed_dir,
    "04_full_length_SERPINE1_fly_alignment.fasta"
  ),

  format = "fasta"

)


# ==============================================================================
# 11. EXPANDED HUMAN + FLY ALIGNMENT
# ==============================================================================

expanded_alignment <- align_proteins(

  ordered_sequences,

  expanded_human_fly_set,

  "Expanded human + Drosophila full-length"

)


writeXStringSet(

  expanded_alignment,

  file.path(
    processed_dir,
    "04_full_length_expanded_human_fly_alignment.fasta"
  ),

  format = "fasta"

)


# ==============================================================================
# 12. READ REVISED RCL DATA
# ==============================================================================

rcl_data <- read_csv(
  rcl_file,
  show_col_types = FALSE
)


cat(
  "\nRCL input columns:\n"
)

print(
  names(rcl_data)
)


# ==============================================================================
# 13. IDENTIFY RCL GENE AND SEQUENCE COLUMNS
# ==============================================================================

gene_column <- "Gene"

rcl_sequence_column <- "Homologous_RCL_Alignment"


# Confirm required columns exist

required_rcl_columns <- c(
  gene_column,
  rcl_sequence_column
)

missing_rcl_columns <- setdiff(
  required_rcl_columns,
  names(rcl_data)
)

if (length(missing_rcl_columns) > 0) {

  stop(
    "Required RCL column(s) missing: ",
    paste(
      missing_rcl_columns,
      collapse = ", "
    )
  )
}


cat(
  "\nUsing RCL gene column: ",
  gene_column,
  "\n"
)

cat(
  "Using RCL sequence column: ",
  rcl_sequence_column,
  "\n"
)

# ==============================================================================
# 14. STANDARDIZE RCL GENE NAMES
# ==============================================================================

rcl_prepared <- rcl_data %>%

  transmute(

    Gene =
      as.character(
        .data[[gene_column]]
      ),

    Sequence =
      as.character(
        .data[[rcl_sequence_column]]
      )

  )


# Handle Step-3 internal names if necessary.

rcl_prepared <- rcl_prepared %>%

  mutate(

    Gene = case_when(

      Gene == "Spn42Da_Q8T0M5" ~
        "Spn42Da",

      Gene == "Spn42Da_B_A1Z6R1" ~
        "Spn42Da-B",

      TRUE ~
        Gene

    )

  )

# ==============================================================================
# 14A. RCL ALIGNMENT QC
# ==============================================================================

rcl_alignment_qc <- rcl_data %>%

  transmute(
    Gene = as.character(Gene),

    Alignment_Length =
      nchar(Homologous_RCL_Alignment),

    Calculated_Ungapped_Length =
      nchar(
        gsub(
          "-",
          "",
          Homologous_RCL_Alignment
        )
      ),

    Reported_Alignment_Length =
      Alignment_Interval_Length,

    Reported_Ungapped_Length =
      Ungapped_Length
  )


print(
  rcl_alignment_qc,
  n = Inf
)


# All sequences should occupy the same alignment interval

if (length(unique(rcl_alignment_qc$Alignment_Length)) != 1) {

  stop(
    "RCL QC FAILED: homologous RCL alignment strings have unequal lengths."
  )
}


# The revised SERPINE1-defined interval should contain 32 alignment columns

if (unique(rcl_alignment_qc$Alignment_Length) != 32) {

  stop(
    "RCL QC FAILED: expected a 32-column SERPINE1-defined alignment interval."
  )
}


# Calculated ungapped lengths should match the values saved by Step 3

if (
  !all(
    rcl_alignment_qc$Calculated_Ungapped_Length ==
      rcl_alignment_qc$Reported_Ungapped_Length
  )
) {

  stop(
    paste0(
      "RCL QC FAILED: calculated and reported ",
      "ungapped lengths do not agree."
    )
  )
}


message(
  "RCL alignment QC PASSED: all sequences span the same 32-column interval."
)

# ==============================================================================
# 15. RCL DATA QC
# ==============================================================================

missing_rcl_genes <- setdiff(
  serpine1_fly_set,
  rcl_prepared$Gene
)


if (length(missing_rcl_genes) > 0) {

  stop(
    "Revised RCL dataset missing required proteins: ",
    paste(
      missing_rcl_genes,
      collapse = ", "
    )
  )

}


if (
  anyDuplicated(
    rcl_prepared$Gene[
      rcl_prepared$Gene %in% serpine1_fly_set
    ]
  ) > 0
) {

  stop(
    "Duplicate genes detected in revised RCL dataset."
  )

}


# ==============================================================================
# 16. BUILD RCL AAStringSet
# ==============================================================================

rcl_subset <- rcl_prepared %>%

  filter(
    Gene %in% serpine1_fly_set
  ) %>%

  mutate(
    Gene = factor(
      Gene,
      levels = serpine1_fly_set
    )
  ) %>%

  arrange(
    Gene
  )


rcl_set <- AAStringSet(
  setNames(
    rcl_subset$Sequence,
    rcl_subset$Gene
  )
)

# ==============================================================================
# 16A. BUILD EXPANDED HUMAN + FLY RCL AAStringSet
# ==============================================================================
#
# Uses the same SERPINE1-defined homologous RCL alignment interval generated
# in Step 3, but retains all 12 human and Drosophila sequences.
#
# ==============================================================================

missing_expanded_rcl_genes <- setdiff(
  expanded_human_fly_set,
  rcl_prepared$Gene
)

if (length(missing_expanded_rcl_genes) > 0) {

  stop(
    "Expanded RCL dataset missing required proteins: ",
    paste(
      missing_expanded_rcl_genes,
      collapse = ", "
    )
  )
}


expanded_rcl_subset <- rcl_prepared %>%

  filter(
    Gene %in% expanded_human_fly_set
  ) %>%

  mutate(
    Gene = factor(
      Gene,
      levels = expanded_human_fly_set
    )
  ) %>%

  arrange(
    Gene
  )


expanded_rcl_set <- AAStringSet(
  setNames(
    expanded_rcl_subset$Sequence,
    expanded_rcl_subset$Gene
  )
)


expanded_rcl_widths <- width(
  expanded_rcl_set
)


if (length(unique(expanded_rcl_widths)) != 1) {

  stop(
    paste0(
      "Expanded RCL QC FAILED: sequences have unequal ",
      "alignment lengths."
    )
  )
}


if (unique(expanded_rcl_widths) != 32) {

  stop(
    paste0(
      "Expanded RCL QC FAILED: expected a 32-column ",
      "SERPINE1-defined homologous alignment interval."
    )
  )
}


message(
  "Expanded human-fly RCL QC PASSED: ",
  length(expanded_rcl_set),
  " proteins across ",
  unique(expanded_rcl_widths),
  " aligned columns."
)


writeXStringSet(
  expanded_rcl_set,

  file.path(
    processed_dir,
    "04_expanded_human_fly_RCL_alignment.fasta"
  ),

  format = "fasta"
)

# ==============================================================================
# 17. CHECK WHETHER RCL REGIONS ARE ALREADY ALIGNED
# ==============================================================================

rcl_widths <- width(
  rcl_set
)


if (length(unique(rcl_widths)) != 1) {

  stop(
    paste0(
      "Revised RCL sequences do not have equal aligned lengths. ",
      "The phylogenetic analysis requires the Step-3 ",
      "SERPINE1-defined homologous alignment regions."
    )
  )

}


cat(
  "\nRevised RCL alignment length = ",
  unique(rcl_widths),
  "\n"
)


writeXStringSet(

  rcl_set,

  file.path(
    processed_dir,
    "04_revised_RCL_SERPINE1_fly_alignment.fasta"
  ),

  format = "fasta"

)


# ==============================================================================
# 18. CONVERT ALIGNED AAStringSet TO PHYDAT
# ==============================================================================

alignment_to_phyDat <- function(aligned_set) {

  aligned_strings <- as.character(
    aligned_set
  )


  matrix_list <- strsplit(
    aligned_strings,
    ""
  )


  alignment_matrix <- do.call(
    rbind,
    matrix_list
  )


  rownames(alignment_matrix) <- names(
    aligned_set
  )


  phy_data <- phyDat(
    alignment_matrix,
    type = "AA"
  )


  list(
    phyDat = phy_data,
    matrix = alignment_matrix
  )

}


# ==============================================================================
# 19. BUILD NJ + ML TREES
# ==============================================================================

build_phylogenies <- function(
  aligned_set,
  analysis_name,
  bootstrap_replicates = 500
) {

  converted <- alignment_to_phyDat(
    aligned_set
  )


  phy_data <- converted$phyDat
  alignment_matrix <- converted$matrix


  # ---------------------------------------------------------------------------
  # JTT evolutionary distances
  # ---------------------------------------------------------------------------

  jtt_distance <- dist.ml(
    phy_data,
    model = "JTT"
  )


  # ---------------------------------------------------------------------------
  # Neighbor Joining
  # ---------------------------------------------------------------------------

  nj_tree <- NJ(
    jtt_distance
  )


  # ---------------------------------------------------------------------------
  # NJ bootstrap
  # ---------------------------------------------------------------------------

  nj_function <- function(x) {

    x_phy <- phyDat(
      x,
      type = "AA"
    )

    x_dist <- dist.ml(
      x_phy,
      model = "JTT"
    )

    NJ(
      x_dist
    )

  }


  cat(
    "\nRunning ",
    bootstrap_replicates,
    " NJ bootstrap replicates for ",
    analysis_name,
    "...\n",
    sep = ""
  )


  nj_bootstrap <- boot.phylo(

    phy = nj_tree,

    x = alignment_matrix,

    FUN = nj_function,

    B = bootstrap_replicates,

    quiet = TRUE

  )


  # ---------------------------------------------------------------------------
  # Maximum Likelihood
  # ---------------------------------------------------------------------------

  starting_fit <- pml(
  nj_tree,
  data = phy_data,
  k = 4
)


  cat(
    "\nOptimizing JTT maximum-likelihood tree for ",
    analysis_name,
    "...\n",
    sep = ""
  )


  ml_fit <- optim.pml(

    starting_fit,

    model = "JTT",

    optInv = TRUE,

    optGamma = TRUE,

    optNni = TRUE,

    rearrangement = "stochastic",

    control = pml.control(
      trace = 0
    )

  )


  ml_tree <- ml_fit$tree


  # ---------------------------------------------------------------------------
  # ML bootstrap
  # ---------------------------------------------------------------------------

  cat(
    "\nRunning ",
    bootstrap_replicates,
    " ML bootstrap replicates for ",
    analysis_name,
    "...\n",
    sep = ""
  )


  ml_bootstrap <- bootstrap.pml(

    ml_fit,

    bs = bootstrap_replicates,

    optNni = TRUE,

    multicore = FALSE,

    control = pml.control(
      trace = 0
    )

  )


  list(

    NJ_tree =
      nj_tree,

    NJ_bootstrap =
      nj_bootstrap,

    ML_tree =
      ml_tree,

    ML_fit =
      ml_fit,

    ML_bootstrap =
      ml_bootstrap,

    Alignment =
      aligned_set,

    PhyDat =
      phy_data,

    LogLik =
      ml_fit$logLik,

    Analysis =
      analysis_name

  )

}


# ==============================================================================
# 20. ANALYSIS A:
#     FULL-LENGTH SERPINE1 + DROSOPHILA
# ==============================================================================

trees_full_serpine1_fly <- build_phylogenies(

  full_serpine1_fly_alignment,

  analysis_name =
    "Full-length SERPINE1 + Drosophila",

  bootstrap_replicates =
    500

)


# ==============================================================================
# 21. ANALYSIS B:
#     REVISED RCL SERPINE1 + DROSOPHILA
# ==============================================================================

trees_rcl <- build_phylogenies(

  rcl_set,

  analysis_name =
    "SERPINE1-defined homologous RCL region",

  bootstrap_replicates =
    500

)


# ==============================================================================
# 22. ANALYSIS C:
#     EXPANDED HUMAN + DROSOPHILA
# ==============================================================================

trees_expanded <- build_phylogenies(

  expanded_alignment,

  analysis_name =
    "Expanded human + Drosophila full-length",

  bootstrap_replicates =
    500

)

# ==============================================================================
# 22A. ANALYSIS D:
#      EXPANDED HUMAN + DROSOPHILA RCL
# ==============================================================================

trees_expanded_rcl <- build_phylogenies(

  expanded_rcl_set,

  analysis_name =
    "Expanded human + Drosophila SERPINE1-defined homologous RCL region",

  bootstrap_replicates =
    500

)

# ==============================================================================
# 23. WRITE NEWICK TREES
# ==============================================================================

write.tree(

  trees_full_serpine1_fly$NJ_tree,

  file.path(
    tree_dir,
    "04A_full_length_SERPINE1_fly_NJ.nwk"
  )

)


write.tree(

  trees_full_serpine1_fly$ML_tree,

  file.path(
    tree_dir,
    "04A_full_length_SERPINE1_fly_ML.nwk"
  )

)


write.tree(

  trees_rcl$NJ_tree,

  file.path(
    tree_dir,
    "04B_revised_RCL_SERPINE1_fly_NJ.nwk"
  )

)


write.tree(

  trees_rcl$ML_tree,

  file.path(
    tree_dir,
    "04B_revised_RCL_SERPINE1_fly_ML.nwk"
  )

)


write.tree(

  trees_expanded$NJ_tree,

  file.path(
    tree_dir,
    "04C_expanded_human_fly_NJ.nwk"
  )

)


write.tree(

  trees_expanded$ML_tree,

  file.path(
    tree_dir,
    "04C_expanded_human_fly_ML.nwk"
  )

)

# ==============================================================================
# EXPANDED HUMAN + FLY RCL NEWICK TREES
# ==============================================================================

write.tree(

  trees_expanded_rcl$NJ_tree,

  file.path(
    tree_dir,
    "04D_expanded_human_fly_RCL_NJ.nwk"
  )

)


write.tree(

  trees_expanded_rcl$ML_tree,

  file.path(
    tree_dir,
    "04D_expanded_human_fly_RCL_ML.nwk"
  )

)

# ==============================================================================
# 24. NJ TREE PLOTTING FUNCTION
# ==============================================================================
#
# boot.phylo() returns the NUMBER of bootstrap replicates supporting each node.
# Because 500 replicates were performed, convert these counts to percentages
# before displaying them on the tree.
#
# ==============================================================================

plot_nj_tree <- function(
  tree,
  bootstrap_values,
  title,
  output_base,
  bootstrap_replicates = 500
) {

  # Convert bootstrap counts to percentages
  bootstrap_percent <- round(
    100 * bootstrap_values / bootstrap_replicates,
    1
  )


  # ---------------------------------------------------------------------------
  # PDF
  # ---------------------------------------------------------------------------

  pdf(
    paste0(
      output_base,
      ".pdf"
    ),
    width = 8,
    height = 7
  )


  plot(
    tree,
    main = title,
    cex = 0.9,
    no.margin = FALSE
  )


  nodelabels(
    text = bootstrap_percent,
    frame = "none",
    cex = 0.7,
    adj = c(
      1.1,
      -0.2
    )
  )


  dev.off()


  # ---------------------------------------------------------------------------
  # PNG
  # ---------------------------------------------------------------------------

  png(
    paste0(
      output_base,
      ".png"
    ),
    width = 2000,
    height = 1700,
    res = 250
  )


  plot(
    tree,
    main = title,
    cex = 0.9,
    no.margin = FALSE
  )


  nodelabels(
    text = bootstrap_percent,
    frame = "none",
    cex = 0.7,
    adj = c(
      1.1,
      -0.2
    )
  )


  dev.off()
}

# ==============================================================================
# 25. ML TREE PLOTTING FUNCTION
# ==============================================================================

plot_ml_tree <- function(
  tree,
  bootstrap_trees,
  title,
  output_base
) {

  # ---------------------------------------------------------------------------
  # Calculate bootstrap support for each internal node
  # ---------------------------------------------------------------------------

  bootstrap_support <- prop.clades(
    tree,
    bootstrap_trees
  )


  bootstrap_percent <- round(
    100 * bootstrap_support / length(bootstrap_trees),
    1
  )


  cat(
    "\nML bootstrap support for:",
    title,
    "\n"
  )

  print(
    bootstrap_percent
  )


  # ---------------------------------------------------------------------------
  # PDF
  # ---------------------------------------------------------------------------

  pdf(
    paste0(
      output_base,
      ".pdf"
    ),
    width = 8,
    height = 7
  )


  plot(
    tree,
    main = title,
    cex = 0.9,
    no.margin = FALSE
  )


  nodelabels(
    text = bootstrap_percent,
    frame = "none",
    cex = 0.7,
    adj = c(
      1.1,
      -0.2
    )
  )


  dev.off()


  # ---------------------------------------------------------------------------
  # PNG
  # ---------------------------------------------------------------------------

  png(
    paste0(
      output_base,
      ".png"
    ),
    width = 2000,
    height = 1700,
    res = 250
  )


  plot(
    tree,
    main = title,
    cex = 0.9,
    no.margin = FALSE
  )


  nodelabels(
    text = bootstrap_percent,
    frame = "none",
    cex = 0.7,
    adj = c(
      1.1,
      -0.2
    )
  )


  dev.off()


  # Return support values invisibly
  invisible(
    bootstrap_percent
  )
}

# ==============================================================================
# 26. SAVE FULL-LENGTH SERPINE1 + FLY TREES
# ==============================================================================

plot_nj_tree(

  trees_full_serpine1_fly$NJ_tree,

  trees_full_serpine1_fly$NJ_bootstrap,

  "Full-length serpin Neighbor-Joining tree",

  file.path(
    figure_dir,
    "04A_full_length_SERPINE1_fly_NJ"
  )

)


plot_ml_tree(

  trees_full_serpine1_fly$ML_tree,

  trees_full_serpine1_fly$ML_bootstrap,

  "Full-length serpin Maximum-Likelihood tree",

  file.path(
    figure_dir,
    "04A_full_length_SERPINE1_fly_ML"
  )

)


# ==============================================================================
# 27. SAVE RCL TREES
# ==============================================================================

plot_nj_tree(

  trees_rcl$NJ_tree,

  trees_rcl$NJ_bootstrap,

  "SERPINE1-defined RCL Neighbor-Joining tree",

  file.path(
    figure_dir,
    "04B_revised_RCL_SERPINE1_fly_NJ"
  )

)


plot_ml_tree(

  trees_rcl$ML_tree,

  trees_rcl$ML_bootstrap,

  "SERPINE1-defined RCL Maximum-Likelihood tree",

  file.path(
    figure_dir,
    "04B_revised_RCL_SERPINE1_fly_ML"
  )

)


# ==============================================================================
# 28. SAVE EXPANDED HUMAN + FLY TREES
# ==============================================================================

plot_nj_tree(

  trees_expanded$NJ_tree,

  trees_expanded$NJ_bootstrap,

  "Expanded human-Drosophila serpin Neighbor-Joining tree",

  file.path(
    figure_dir,
    "04C_expanded_human_fly_NJ"
  )

)


plot_ml_tree(

  trees_expanded$ML_tree,

  trees_expanded$ML_bootstrap,

  "Expanded human-Drosophila serpin Maximum-Likelihood tree",

  file.path(
    figure_dir,
    "04C_expanded_human_fly_ML"
  )

)

# ==============================================================================
# 28A. SAVE EXPANDED HUMAN + FLY RCL TREES
# ==============================================================================

plot_nj_tree(

  trees_expanded_rcl$NJ_tree,

  trees_expanded_rcl$NJ_bootstrap,

  "Expanded human-Drosophila homologous RCL Neighbor-Joining tree",

  file.path(
    figure_dir,
    "04D_expanded_human_fly_RCL_NJ"
  )

)


plot_ml_tree(

  trees_expanded_rcl$ML_tree,

  trees_expanded_rcl$ML_bootstrap,

  "Expanded human-Drosophila homologous RCL Maximum-Likelihood tree",

  file.path(
    figure_dir,
    "04D_expanded_human_fly_RCL_ML"
  )

)

# ==============================================================================
# 28B. FOUR-PANEL TREE DRAWING FUNCTIONS
# ==============================================================================


# ------------------------------------------------------------------------------
# Draw NJ tree inside an existing graphics device
# ------------------------------------------------------------------------------

draw_nj_panel <- function(
  tree,
  bootstrap_values,
  panel_label,
  panel_title,
  bootstrap_replicates = 500
) {

  bootstrap_percent <- round(
    100 * bootstrap_values / bootstrap_replicates,
    1
  )


  plot(
    tree,
    main = panel_title,
    cex = 0.75,
    no.margin = FALSE
  )


  nodelabels(
    text = bootstrap_percent,
    frame = "none",
    cex = 0.55,
    adj = c(
      1.1,
      -0.2
    )
  )


  mtext(
    panel_label,
    side = 3,
    adj = 0,
    line = 0.3,
    font = 2,
    cex = 1.1
  )
}


# ------------------------------------------------------------------------------
# Draw ML tree inside an existing graphics device
# ------------------------------------------------------------------------------

draw_ml_panel <- function(
  tree,
  bootstrap_trees,
  panel_label,
  panel_title
) {

  bootstrap_support <- prop.clades(
    tree,
    bootstrap_trees
  )


  bootstrap_percent <- round(
    100 * bootstrap_support /
      length(bootstrap_trees),
    1
  )


  plot(
    tree,
    main = panel_title,
    cex = 0.75,
    no.margin = FALSE
  )


  nodelabels(
    text = bootstrap_percent,
    frame = "none",
    cex = 0.55,
    adj = c(
      1.1,
      -0.2
    )
  )


  mtext(
    panel_label,
    side = 3,
    adj = 0,
    line = 0.3,
    font = 2,
    cex = 1.1
  )
}

# ==============================================================================
# 28C. FIGURE 1:
#      SERPINE1 + DROSOPHILA PHYLOGENETIC ANALYSIS
#
#      A = Full-length NJ
#      B = Full-length ML
#      C = Homologous RCL NJ
#      D = Homologous RCL ML
# ==============================================================================


figure1_pdf <- file.path(
  figure_dir,
  "04_Figure1_SERPINE1_Fly_FullLength_and_RCL_Phylogeny.pdf"
)


pdf(
  figure1_pdf,
  width = 14,
  height = 12
)


par(
  mfrow = c(2, 2),
  mar = c(2, 2, 4, 1),
  oma = c(1, 1, 3, 1)
)


# Panel A
draw_nj_panel(

  trees_full_serpine1_fly$NJ_tree,

  trees_full_serpine1_fly$NJ_bootstrap,

  panel_label = "A",

  panel_title =
    "Full-length Neighbor-Joining"

)


# Panel B
draw_ml_panel(

  trees_full_serpine1_fly$ML_tree,

  trees_full_serpine1_fly$ML_bootstrap,

  panel_label = "B",

  panel_title =
    "Full-length Maximum-Likelihood"

)


# Panel C
draw_nj_panel(

  trees_rcl$NJ_tree,

  trees_rcl$NJ_bootstrap,

  panel_label = "C",

  panel_title =
    "Homologous RCL-region Neighbor-Joining"

)


# Panel D
draw_ml_panel(

  trees_rcl$ML_tree,

  trees_rcl$ML_bootstrap,

  panel_label = "D",

  panel_title =
    "Homologous RCL-region Maximum-Likelihood"

)


mtext(
  "SERPINE1 and Drosophila serpin phylogenetic analysis",
  outer = TRUE,
  side = 3,
  line = 1,
  font = 2,
  cex = 1.4
)


dev.off()

figure1_png <- file.path(
  figure_dir,
  "04_Figure1_SERPINE1_Fly_FullLength_and_RCL_Phylogeny.png"
)


png(
  figure1_png,
  width = 3600,
  height = 3000,
  res = 250
)


par(
  mfrow = c(2, 2),
  mar = c(2, 2, 4, 1),
  oma = c(1, 1, 3, 1)
)


draw_nj_panel(
  trees_full_serpine1_fly$NJ_tree,
  trees_full_serpine1_fly$NJ_bootstrap,
  "A",
  "Full-length Neighbor-Joining"
)


draw_ml_panel(
  trees_full_serpine1_fly$ML_tree,
  trees_full_serpine1_fly$ML_bootstrap,
  "B",
  "Full-length Maximum-Likelihood"
)


draw_nj_panel(
  trees_rcl$NJ_tree,
  trees_rcl$NJ_bootstrap,
  "C",
  "Homologous RCL-region Neighbor-Joining"
)


draw_ml_panel(
  trees_rcl$ML_tree,
  trees_rcl$ML_bootstrap,
  "D",
  "Homologous RCL-region Maximum-Likelihood"
)


mtext(
  "SERPINE1 and Drosophila serpin phylogenetic analysis",
  outer = TRUE,
  side = 3,
  line = 1,
  font = 2,
  cex = 1.4
)


dev.off()

# ==============================================================================
# 28D. FIGURE 2:
#      EXPANDED HUMAN + DROSOPHILA PHYLOGENETIC ANALYSIS
#
#      A = Full-length NJ
#      B = Full-length ML
#      C = Homologous RCL NJ
#      D = Homologous RCL ML
# ==============================================================================


figure2_pdf <- file.path(
  figure_dir,
  "04_Figure2_Expanded_Human_Fly_FullLength_and_RCL_Phylogeny.pdf"
)


pdf(
  figure2_pdf,
  width = 15,
  height = 13
)


par(
  mfrow = c(2, 2),
  mar = c(2, 2, 4, 1),
  oma = c(1, 1, 3, 1)
)


# Panel A
draw_nj_panel(

  trees_expanded$NJ_tree,

  trees_expanded$NJ_bootstrap,

  panel_label = "A",

  panel_title =
    "Full-length Neighbor-Joining"

)


# Panel B
draw_ml_panel(

  trees_expanded$ML_tree,

  trees_expanded$ML_bootstrap,

  panel_label = "B",

  panel_title =
    "Full-length Maximum-Likelihood"

)


# Panel C
draw_nj_panel(

  trees_expanded_rcl$NJ_tree,

  trees_expanded_rcl$NJ_bootstrap,

  panel_label = "C",

  panel_title =
    "Homologous RCL-region Neighbor-Joining"

)


# Panel D
draw_ml_panel(

  trees_expanded_rcl$ML_tree,

  trees_expanded_rcl$ML_bootstrap,

  panel_label = "D",

  panel_title =
    "Homologous RCL-region Maximum-Likelihood"

)


mtext(
  "Expanded human and Drosophila serpin phylogenetic analysis",
  outer = TRUE,
  side = 3,
  line = 1,
  font = 2,
  cex = 1.4
)


dev.off()

figure2_png <- file.path(
  figure_dir,
  "04_Figure2_Expanded_Human_Fly_FullLength_and_RCL_Phylogeny.png"
)


png(
  figure2_png,
  width = 3800,
  height = 3200,
  res = 250
)


par(
  mfrow = c(2, 2),
  mar = c(2, 2, 4, 1),
  oma = c(1, 1, 3, 1)
)


draw_nj_panel(
  trees_expanded$NJ_tree,
  trees_expanded$NJ_bootstrap,
  "A",
  "Full-length Neighbor-Joining"
)


draw_ml_panel(
  trees_expanded$ML_tree,
  trees_expanded$ML_bootstrap,
  "B",
  "Full-length Maximum-Likelihood"
)


draw_nj_panel(
  trees_expanded_rcl$NJ_tree,
  trees_expanded_rcl$NJ_bootstrap,
  "C",
  "Homologous RCL-region Neighbor-Joining"
)


draw_ml_panel(
  trees_expanded_rcl$ML_tree,
  trees_expanded_rcl$ML_bootstrap,
  "D",
  "Homologous RCL-region Maximum-Likelihood"
)


mtext(
  "Expanded human and Drosophila serpin phylogenetic analysis",
  outer = TRUE,
  side = 3,
  line = 1,
  font = 2,
  cex = 1.4
)


dev.off()


# ==============================================================================
# 29. PATRISTIC DISTANCE HELPER
# ==============================================================================

get_patristic_distance <- function(
  tree,
  reference = "SERPINE1",
  distance_name
) {

  if (!reference %in% tree$tip.label) {

    stop(
      "Reference protein not present in tree: ",
      reference
    )

  }


  distance_matrix <- cophenetic.phylo(
    tree
  )


  tibble(

    Gene =
      rownames(
        distance_matrix
      ),

    Distance =
      as.numeric(
        distance_matrix[
          ,
          reference
        ]
      )

  ) %>%

    filter(
      Gene != reference
    ) %>%

    rename(
      !!distance_name :=
        Distance
    )

}


# ==============================================================================
# 30. FULL-LENGTH SERPINE1 + FLY DISTANCES
# ==============================================================================

full_nj_distances <- get_patristic_distance(

  trees_full_serpine1_fly$NJ_tree,

  distance_name =
    "Full_Length_NJ_Patristic_Distance_to_SERPINE1"

)


full_ml_distances <- get_patristic_distance(

  trees_full_serpine1_fly$ML_tree,

  distance_name =
    "Full_Length_ML_Patristic_Distance_to_SERPINE1"

)


# ==============================================================================
# 31. RCL DISTANCES
# ==============================================================================

rcl_nj_distances <- get_patristic_distance(

  trees_rcl$NJ_tree,

  distance_name =
    "RCL_NJ_Patristic_Distance_to_SERPINE1"

)


rcl_ml_distances <- get_patristic_distance(

  trees_rcl$ML_tree,

  distance_name =
    "RCL_ML_Patristic_Distance_to_SERPINE1"

)


# ==============================================================================
# 32. COMBINE SERPINE1-FOCUSED DISTANCES
# ==============================================================================

serpine1_phylogenetic_distances <- full_nj_distances %>%

  full_join(
    full_ml_distances,
    by = "Gene"
  ) %>%

  full_join(
    rcl_nj_distances,
    by = "Gene"
  ) %>%

  full_join(
    rcl_ml_distances,
    by = "Gene"
  ) %>%

  filter(
    Gene %in% fly_genes
  ) %>%

  arrange(
    Full_Length_ML_Patristic_Distance_to_SERPINE1
  )


print(
  serpine1_phylogenetic_distances,
  n = Inf,
  width = Inf
)


write_csv(

  serpine1_phylogenetic_distances,

  file.path(
    processed_dir,
    "04_SERPINE1_Fly_Phylogenetic_Distances.csv"
  )

)


# ==============================================================================
# 33. EXPANDED-TREE DISTANCES TO SERPINE1
# ==============================================================================

expanded_nj_distances <- get_patristic_distance(

  trees_expanded$NJ_tree,

  distance_name =
    "Expanded_NJ_Patristic_Distance_to_SERPINE1"

)


expanded_ml_distances <- get_patristic_distance(

  trees_expanded$ML_tree,

  distance_name =
    "Expanded_ML_Patristic_Distance_to_SERPINE1"

)


expanded_distances <- expanded_nj_distances %>%

  full_join(
    expanded_ml_distances,
    by = "Gene"
  ) %>%

  left_join(
    sequence_key,
    by = "Gene"
  ) %>%

  arrange(
    Expanded_ML_Patristic_Distance_to_SERPINE1
  )


print(
  expanded_distances,
  n = Inf,
  width = Inf
)


write_csv(

  expanded_distances,

  file.path(
    processed_dir,
    "04_Expanded_Human_Fly_Distances_to_SERPINE1.csv"
  )

)


# ==============================================================================
# 34. RANK SEVEN DROSOPHILA SEQUENCES
# ==============================================================================

fly_phylogenetic_ranking <- serpine1_phylogenetic_distances %>%

  mutate(

    Full_Length_NJ_Rank =
      rank(
        Full_Length_NJ_Patristic_Distance_to_SERPINE1,
        ties.method = "min"
      ),

    Full_Length_ML_Rank =
      rank(
        Full_Length_ML_Patristic_Distance_to_SERPINE1,
        ties.method = "min"
      ),

    RCL_NJ_Rank =
      rank(
        RCL_NJ_Patristic_Distance_to_SERPINE1,
        ties.method = "min"
      ),

    RCL_ML_Rank =
      rank(
        RCL_ML_Patristic_Distance_to_SERPINE1,
        ties.method = "min"
      )

  )


print(
  fly_phylogenetic_ranking,
  n = Inf,
  width = Inf
)


write_csv(

  fly_phylogenetic_ranking,

  file.path(
    processed_dir,
    "04_Fly_Phylogenetic_Ranking.csv"
  )

)


# ==============================================================================
# 35. ML MODEL SUMMARY
# ==============================================================================

ml_summary <- tibble(

  Analysis = c(
    "Full-length SERPINE1 + fly",
    "SERPINE1-defined homologous RCL + fly",
    "Expanded human + fly full-length",
    "Expanded human + fly homologous RCL"
  ),

  Number_of_Tips = c(
    length(
      trees_full_serpine1_fly$ML_tree$tip.label
    ),

    length(
      trees_rcl$ML_tree$tip.label
    ),

    length(
      trees_expanded$ML_tree$tip.label
    ),

    length(
      trees_expanded_rcl$ML_tree$tip.label
    )
  ),

  Log_Likelihood = c(
    trees_full_serpine1_fly$LogLik,
    trees_rcl$LogLik,
    trees_expanded$LogLik,
    trees_expanded_rcl$LogLik
  ),

  Model = "JTT",

  Bootstrap_Replicates = 500

)



print(
  ml_summary,
  n = Inf
)


write_csv(

  ml_summary,

  file.path(
    qc_dir,
    "04_ML_Model_Summary.csv"
  )

)


# ==============================================================================
# 36. TREE TIP QC
# ==============================================================================

tree_tip_qc <- bind_rows(

    tibble(
        Analysis = "Full_SERPINE1_Fly_NJ",
        Gene = trees_full_serpine1_fly$NJ_tree$tip.label
    ),

    tibble(
        Analysis = "Full_SERPINE1_Fly_ML",
        Gene = trees_full_serpine1_fly$ML_tree$tip.label
    ),

    tibble(
        Analysis = "RCL_SERPINE1_Fly_NJ",
        Gene = trees_rcl$NJ_tree$tip.label
    ),

    tibble(
        Analysis = "RCL_SERPINE1_Fly_ML",
        Gene = trees_rcl$ML_tree$tip.label
    ),

    tibble(
        Analysis = "Expanded_Human_Fly_NJ",
        Gene = trees_expanded$NJ_tree$tip.label
    ),

    tibble(
        Analysis = "Expanded_Human_Fly_ML",
        Gene = trees_expanded$ML_tree$tip.label
    ),

    tibble(
        Analysis = "Expanded_Human_Fly_RCL_NJ",
        Gene = trees_expanded_rcl$NJ_tree$tip.label
    ),

    tibble(
        Analysis = "Expanded_Human_Fly_RCL_ML",
        Gene = trees_expanded_rcl$ML_tree$tip.label
    )

)   # <-- closes bind_rows()


# Print QC table

print(
    tree_tip_qc,
    n = Inf
)


# Save QC table

write_csv(
    tree_tip_qc,
    file.path(
        qc_dir,
        "04_Tree_Tip_QC.csv"
    )
)

# ==============================================================================
# 37. SESSION INFORMATION
# ==============================================================================

capture.output(

  sessionInfo(),

  file = file.path(
    qc_dir,
    "04_sessionInfo.txt"
  )

)


# ==============================================================================
# 38. COMPLETION
# ==============================================================================

message(
    "\n============================================================"
)

message(
    "PHYLOGENETIC ANALYSIS COMPLETE"
)

message(
    ""
)

message(
    "Eight trees generated:"
)

message(
    "  1. Full-length SERPINE1 + fly NJ"
)

message(
    "  2. Full-length SERPINE1 + fly ML"
)

message(
    "  3. Homologous RCL SERPINE1 + fly NJ"
)

message(
    "  4. Homologous RCL SERPINE1 + fly ML"
)

message(
    "  5. Expanded human + fly full-length NJ"
)

message(
    "  6. Expanded human + fly full-length ML"
)

message(
    "  7. Expanded human + fly homologous RCL NJ"
)

message(
    "  8. Expanded human + fly homologous RCL ML"
)

message(
    "============================================================\n"
)
