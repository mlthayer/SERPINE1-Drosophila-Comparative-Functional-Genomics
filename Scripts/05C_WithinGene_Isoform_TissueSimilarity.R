# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 05C_WithinGene_Isoform_TissueSimilarity.R
# Source after the adult isoform expression analysis in the same R session.
# ==============================================================================
# PART 3: WITHIN-GENE ISOFORM TISSUE-PROFILE SIMILARITY
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)

# Portable directories shared by the isoform analyses.
isoform_output_dir <- Sys.getenv(
  "SERPINE1_FLYATLAS_OUTPUT", unset = "Expression/FlyAtlas2_Analysis")
processed_dir <- file.path(isoform_output_dir, "Processed")
table_dir <- file.path(isoform_output_dir, "Tables")
figure_dir <- file.path(isoform_output_dir, "Figures")
qc_dir <- file.path(isoform_output_dir, "QC")
for (path in c(processed_dir, table_dir, figure_dir, qc_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# A standalone run uses the exact processed table written by step 05B.
if (!exists("adult_isoform_expression", inherits = FALSE)) {
  isoform_input <- file.path(processed_dir, "05_FlyAtlas2_adult_isoform_analysis.csv")
  if (!file.exists(isoform_input)) {
    stop("Run step 05B first; its processed adult isoform table is missing: ",
         isoform_input, call. = FALSE)
  }
  adult_isoform_expression <- readr::read_csv(isoform_input, show_col_types = FALSE)
}
required_isoform_columns <- c(
  "Gene", "Isoform_Label", "FlyBase_Transcript", "Sex", "Tissue", "log2_FPKM")
missing_isoform_columns <- setdiff(
  required_isoform_columns, names(adult_isoform_expression))
if (length(missing_isoform_columns)) {
  stop("Adult isoform table is missing columns: ",
       paste(missing_isoform_columns, collapse = ", "), call. = FALSE)
}

# ------------------------------------------------------------------------------
# Create matched adult tissue/sex identifier
# ------------------------------------------------------------------------------

isoform_similarity_data <- adult_isoform_expression %>%
    mutate(
        Sample = paste(
            as.character(Sex),
            as.character(Tissue),
            sep = " | "
        ),
        Isoform_Label = as.character(Isoform_Label)
    ) %>%
    select(
        Gene,
        Isoform_Label,
        FlyBase_Transcript,
        Sample,
        log2_FPKM
    )


# ==============================================================================
# FUNCTION: PAIRWISE CORRELATIONS WITHIN A GENE
# ==============================================================================

calculate_isoform_correlations <- function(data, gene_name) {

    gene_data <- data %>%
        filter(Gene == gene_name)

    wide <- gene_data %>%
        select(
            Sample,
            Isoform_Label,
            log2_FPKM
        ) %>%
        pivot_wider(
            names_from = Isoform_Label,
            values_from = log2_FPKM
        )

    expression_matrix <- wide %>%
        select(-Sample) %>%
        as.data.frame()

    # Pearson
    pearson <- cor(
        expression_matrix,
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    # Spearman
    spearman <- cor(
        expression_matrix,
        method = "spearman",
        use = "pairwise.complete.obs"
    )

    list(
        Pearson = pearson,
        Spearman = spearman
    )
}


# ==============================================================================
# MULTI-ISOFORM GENES
# ==============================================================================

multi_isoform_genes <- adult_isoform_expression %>%
    distinct(
        Gene,
        FlyBase_Transcript
    ) %>%
    count(Gene) %>%
    filter(n > 1) %>%
    pull(Gene)

print(multi_isoform_genes)


# ==============================================================================
# CALCULATE CORRELATIONS
# ==============================================================================

correlation_results <- setNames(
    lapply(
        multi_isoform_genes,
        function(g) {
            calculate_isoform_correlations(
                isoform_similarity_data,
                g
            )
        }
    ),
    multi_isoform_genes
)


# ==============================================================================
# CONVERT MATRICES TO LONG TABLE
# ==============================================================================

correlation_long <- map_dfr(
    names(correlation_results),
    function(g) {

        pearson <- correlation_results[[g]]$Pearson
        spearman <- correlation_results[[g]]$Spearman

        pearson_long <- as.data.frame(
            as.table(pearson)
        ) %>%
            rename(
                Isoform_1 = Var1,
                Isoform_2 = Var2,
                Pearson_r = Freq
            )

        spearman_long <- as.data.frame(
            as.table(spearman)
        ) %>%
            rename(
                Isoform_1 = Var1,
                Isoform_2 = Var2,
                Spearman_rho = Freq
            )

        left_join(
            pearson_long,
            spearman_long,
            by = c(
                "Isoform_1",
                "Isoform_2"
            )
        ) %>%
            mutate(
                Gene = g,
                .before = 1
            )
    }
)


# Remove self-comparisons and duplicate reciprocal comparisons
correlation_unique <- correlation_long %>%
    mutate(
        Isoform_1 = as.character(Isoform_1),
        Isoform_2 = as.character(Isoform_2)
    ) %>%
    filter(
        Isoform_1 != Isoform_2
    ) %>%
    rowwise() %>%
    mutate(
        Pair_ID = paste(
            sort(c(Isoform_1, Isoform_2)),
            collapse = "___"
        )
    ) %>%
    ungroup() %>%
    distinct(
        Gene,
        Pair_ID,
        .keep_all = TRUE
    ) %>%
    select(-Pair_ID)


# ==============================================================================
# SPN42DA-B COMPARISON WITH OTHER SPN42DA ISOFORMS
# ==============================================================================

spn42da_b_similarity <- correlation_unique %>%
    filter(
        Gene == "Spn42Da",
        Isoform_1 == "Spn42Da-B (RB)" |
        Isoform_2 == "Spn42Da-B (RB)"
    ) %>%
    mutate(
        Other_Isoform = if_else(
            Isoform_1 == "Spn42Da-B (RB)",
            Isoform_2,
            Isoform_1
        )
    ) %>%
    select(
        Other_Isoform,
        Pearson_r,
        Spearman_rho
    ) %>%
    arrange(
        desc(Pearson_r)
    )


cat("\n============================================================\n")
cat("SPN42DA-B ISOFORM PROFILE SIMILARITY\n")
cat("============================================================\n\n")

print(
    spn42da_b_similarity,
    n = Inf
)


# ==============================================================================
# SPN42DA PEARSON HEATMAP
# ==============================================================================

spn42da_pearson <- correlation_results[["Spn42Da"]]$Pearson

spn42da_pearson_long <- as.data.frame(
    as.table(spn42da_pearson)
) %>%
    rename(
        Isoform_1 = Var1,
        Isoform_2 = Var2,
        Pearson_r = Freq
    )

p_spn42da_correlation <- ggplot(
    spn42da_pearson_long,
    aes(
        x = Isoform_1,
        y = Isoform_2,
        fill = Pearson_r
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.4
    ) +
    geom_text(
        aes(
            label = sprintf("%.2f", Pearson_r)
        ),
        size = 3
    ) +
    scale_fill_gradient2(
        midpoint = 0,
        limits = c(-1, 1),
        name = "Pearson r"
    ) +
    labs(
        title = "Similarity among Spn42Da isoform tissue-expression profiles",
        subtitle = expression(
            "Pearson correlation of " *
            log[2] * "(FPKM + 1) across adult sex × tissue observations"
        ),
        x = NULL,
        y = NULL
    ) +
    coord_fixed() +
    theme_minimal(base_size = 10) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        plot.title = element_text(
            face = "bold"
        )
    )

ggsave(
    file.path(
        figure_dir,
        "04E_Spn42Da_Isoform_Correlation.png"
    ),
    p_spn42da_correlation,
    width = 9,
    height = 8,
    dpi = 600
)


# ==============================================================================
# SAVE RESULTS
# ==============================================================================

write_csv(
    correlation_unique,
    file.path(
        table_dir,
        "05_FlyAtlas2_within_gene_isoform_correlations.csv"
    )
)

write_csv(
    spn42da_b_similarity,
    file.path(
        table_dir,
        "05_Spn42Da-B_isoform_similarity.csv"
    )
)
