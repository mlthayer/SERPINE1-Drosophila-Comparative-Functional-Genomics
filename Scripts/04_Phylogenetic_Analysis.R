# ============================================================
# 04_Phylogenetic_Analysis.R
# Full-length and RCL phylogenetic trees
# Neighbor joining and maximum likelihood
# ============================================================

library(dplyr)
library(ape)
library(phangorn)
library(ggplot2)

# Prefer an environment variable for reproducibility across machines.
# If it is not set, fall back to the original local project path.
project_dir <- Sys.getenv(
  "xxxxx_xxxxxxxx_PROJECT_DIR",
  unset = "D:/xxxxx-xxxxxxxx xxx/xxxxx-xxxxxxxx xxx"
)

chapter3_dir <- file.path(project_dir, "Bioinformatics chapter three")

data_dir <- file.path(chapter3_dir, "Data")
processed_dir <- file.path(data_dir, "Processed")
fig_dir <- file.path(chapter3_dir, "Figures", "Phylogeny")
tree_dir <- file.path(data_dir, "Trees")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)

full_align_file <- file.path(data_dir, "full_alignment_text_output.csv")
rcl_file <- file.path(data_dir, "RCL_sequences.csv")

required_files <- c(
  Full_Length = full_align_file,
  RCL = rcl_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required alignment file(s): ",
    paste(names(missing_files), basename(missing_files), sep = " = ", collapse = "; ")
  )
}

full_align <- read.csv(full_align_file, stringsAsFactors = FALSE)
rcl <- read.csv(rcl_file, stringsAsFactors = FALSE)

required_full_columns <- c("Gene", "Alignment")
required_rcl_columns <- c("Gene", "RCL")

missing_full_columns <- setdiff(required_full_columns, names(full_align))
missing_rcl_columns <- setdiff(required_rcl_columns, names(rcl))

if (length(missing_full_columns) > 0) {
  stop(
    "Full-length alignment table is missing required column(s): ",
    paste(missing_full_columns, collapse = ", ")
  )
}

if (length(missing_rcl_columns) > 0) {
  stop(
    "RCL table is missing required column(s): ",
    paste(missing_rcl_columns, collapse = ", ")
  )
}

selected_genes <- c(
  "PAI1",
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)

reference_gene <- "PAI1"

if (!reference_gene %in% selected_genes) {
  stop("reference_gene must be included in selected_genes.")
}

if (length(selected_genes) < 3) {
  stop("At least three selected genes are required to build phylogenetic trees.")
}

set.seed(1)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

prepare_alignment_table <- function(df, seq_col, selected_genes, dataset_label) {
  if (anyDuplicated(df$Gene) > 0) {
    duplicated_genes <- unique(df$Gene[duplicated(df$Gene)])
    warning(
      dataset_label,
      " table contains duplicate Gene values: ",
      paste(duplicated_genes, collapse = ", "),
      ". The first non-missing sequence per gene will be used."
    )
  }

  missing_genes <- setdiff(selected_genes, unique(df$Gene))

  if (length(missing_genes) > 0) {
    warning(
      dataset_label,
      " table is missing selected gene(s): ",
      paste(missing_genes, collapse = ", ")
    )
  }

  prepared <- df %>%
    filter(Gene %in% selected_genes) %>%
    mutate(
      Gene = as.character(Gene),
      Sequence = toupper(as.character(.data[[seq_col]]))
    ) %>%
    filter(!is.na(Sequence), Sequence != "") %>%
    group_by(Gene) %>%
    slice(1) %>%
    ungroup()

  if (!reference_gene %in% prepared$Gene) {
    stop(dataset_label, " alignment does not contain a usable sequence for reference gene: ", reference_gene)
  }

  if (nrow(prepared) < 3) {
    stop(dataset_label, " alignment has fewer than three usable selected sequences.")
  }

  invalid_sequences <- prepared %>%
    filter(grepl("[^A-Z\\-\\?]", Sequence))

  if (nrow(invalid_sequences) > 0) {
    warning(
      dataset_label,
      " alignment contains nonstandard sequence character(s) for gene(s): ",
      paste(invalid_sequences$Gene, collapse = ", "),
      ". Characters other than A-Z, '-', and '?' may not be interpreted as expected."
    )
  }

  prepared %>%
    select(Gene, Sequence)
}

alignment_to_phyDat <- function(df) {
  seqs <- df$Sequence
  names(seqs) <- df$Gene

  seqs_split <- strsplit(seqs, "")
  max_len <- max(lengths(seqs_split))

  if (!is.finite(max_len) || max_len == 0) {
    stop("Alignment contains no usable sequence characters.")
  }

  seqs_padded <- lapply(seqs_split, function(x) {
    length(x) <- max_len
    x[is.na(x)] <- "-"
    x
  })

  mat <- do.call(rbind, seqs_padded)
  rownames(mat) <- names(seqs)

  phyDat(mat, type = "AA")
}

build_phylogenies <- function(phy_data, dataset_label) {
  distance <- dist.ml(phy_data)
  nj_tree <- NJ(distance)

  fit <- pml(nj_tree, data = phy_data)
  fit_opt <- optim.pml(
    fit,
    model = "JTT",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    control = pml.control(trace = 0)
  )

  list(
    nj = nj_tree,
    ml = fit_opt$tree,
    logLik = fit_opt$logLik,
    dataset = dataset_label
  )
}

plot_tree_to_files <- function(tree, title, pdf_file, png_file) {
  pdf(pdf_file, width = 7, height = 6)
  plot(tree, main = title)
  dev.off()

  png(png_file, width = 1600, height = 1300, res = 200)
  plot(tree, main = title)
  dev.off()
}

get_distance_to_pai1 <- function(tree, label = "PAI1") {
  if (!label %in% tree$tip.label) {
    stop("Reference label not found in tree tips: ", label)
  }

  dist_mat <- cophenetic.phylo(tree)

  data.frame(
    Gene = rownames(dist_mat),
    Distance_to_PAI1 = as.numeric(dist_mat[, label]),
    stringsAsFactors = FALSE
  ) %>%
    filter(Gene != label)
}

# ------------------------------------------------------------
# Full-length tree
# ------------------------------------------------------------

full_prepared <- prepare_alignment_table(
  full_align,
  seq_col = "Alignment",
  selected_genes = selected_genes,
  dataset_label = "Full-length"
)

full_phy <- alignment_to_phyDat(full_prepared)
full_trees <- build_phylogenies(full_phy, "Full-length")

full_nj <- full_trees$nj
full_ml <- full_trees$ml

write.tree(
  full_nj,
  file.path(tree_dir, "full_length_serpin_NJ_tree.nwk")
)

write.tree(
  full_ml,
  file.path(tree_dir, "full_length_serpin_ML_tree.nwk")
)

plot_tree_to_files(
  full_nj,
  "Full-length serpin NJ tree",
  file.path(fig_dir, "full_length_NJ_tree.pdf"),
  file.path(fig_dir, "full_length_NJ_tree.png")
)

plot_tree_to_files(
  full_ml,
  "Full-length serpin ML tree",
  file.path(fig_dir, "full_length_ML_tree.pdf"),
  file.path(fig_dir, "full_length_ML_tree.png")
)

# ------------------------------------------------------------
# RCL tree
# ------------------------------------------------------------

rcl_prepared <- prepare_alignment_table(
  rcl,
  seq_col = "RCL",
  selected_genes = selected_genes,
  dataset_label = "RCL"
)

rcl_phy <- alignment_to_phyDat(rcl_prepared)
rcl_trees <- build_phylogenies(rcl_phy, "RCL")

rcl_nj <- rcl_trees$nj
rcl_ml <- rcl_trees$ml

write.tree(
  rcl_nj,
  file.path(tree_dir, "RCL_domain_serpin_NJ_tree.nwk")
)

write.tree(
  rcl_ml,
  file.path(tree_dir, "RCL_domain_serpin_ML_tree.nwk")
)

plot_tree_to_files(
  rcl_nj,
  "RCL domain NJ tree",
  file.path(fig_dir, "RCL_domain_NJ_tree.pdf"),
  file.path(fig_dir, "RCL_domain_NJ_tree.png")
)

plot_tree_to_files(
  rcl_ml,
  "RCL domain ML tree",
  file.path(fig_dir, "RCL_domain_ML_tree.pdf"),
  file.path(fig_dir, "RCL_domain_ML_tree.png")
)

# ------------------------------------------------------------
# Patristic distance to human PAI1
# ------------------------------------------------------------

full_ml_dist <- get_distance_to_pai1(full_ml, reference_gene) %>%
  rename(Full_Length_ML_Distance_to_PAI1 = Distance_to_PAI1)

full_nj_dist <- get_distance_to_pai1(full_nj, reference_gene) %>%
  rename(Full_Length_NJ_Distance_to_PAI1 = Distance_to_PAI1)

rcl_ml_dist <- get_distance_to_pai1(rcl_ml, reference_gene) %>%
  rename(RCL_ML_Distance_to_PAI1 = Distance_to_PAI1)

rcl_nj_dist <- get_distance_to_pai1(rcl_nj, reference_gene) %>%
  rename(RCL_NJ_Distance_to_PAI1 = Distance_to_PAI1)

phylo_distances <- full_ml_dist %>%
  full_join(full_nj_dist, by = "Gene") %>%
  full_join(rcl_ml_dist, by = "Gene") %>%
  full_join(rcl_nj_dist, by = "Gene") %>%
  arrange(Full_Length_ML_Distance_to_PAI1, Gene)

if (nrow(phylo_distances) == 0) {
  warning("No phylogenetic distance rows were generated.")
}

if (any(is.na(phylo_distances$Full_Length_ML_Distance_to_PAI1))) {
  warning("One or more genes are missing Full_Length_ML_Distance_to_PAI1 values.")
}

print(phylo_distances)

write.csv(
  phylo_distances,
  file.path(processed_dir, "Phylogenetic_Distances_to_Human_PAI1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Distance ranking plot
# ------------------------------------------------------------

plot_data <- phylo_distances %>%
  filter(!is.na(Full_Length_ML_Distance_to_PAI1))

if (nrow(plot_data) == 0) {
  warning("No non-missing full-length ML distances are available for plotting.")
} else {
  p_phylo <- ggplot(
    plot_data,
    aes(
      x = Full_Length_ML_Distance_to_PAI1,
      y = reorder(Gene, -Full_Length_ML_Distance_to_PAI1)
    )
  ) +
    geom_col(width = 0.7, fill = "#2166AC") +
    theme_bw(base_size = 12) +
    labs(
      title = "Phylogenetic distance to human PAI1",
      x = "ML patristic distance",
      y = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  print(p_phylo)

  ggsave(
    file.path(fig_dir, "Phylogenetic_Distance_to_Human_PAI1.png"),
    p_phylo,
    width = 7,
    height = 4.5,
    dpi = 600,
    bg = "white"
  )
}
