# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 05B_Adult_Isoform_TissueExpression.R
# Run from the repository root.
# ==============================================================================
# PART 2: ADULT ISOFORM-LEVEL TISSUE EXPRESSION ANALYSIS
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

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

# Reuse an upstream object when sourced, or load its original exported CSV.
if (!exists("flyatlas_adult_transcripts", inherits = FALSE)) {
  transcript_input <- Sys.getenv(
    "SERPINE1_ADULT_TRANSCRIPT_INPUT",
    unset = file.path(processed_dir, "05_FlyAtlas2_all_transcript_expression.csv"))
  if (!file.exists(transcript_input)) {
    stop("Adult transcript data are missing. Supply the original upstream CSV ",
         "using SERPINE1_ADULT_TRANSCRIPT_INPUT, or source its preparation first. ",
         "Expected columns: Gene, Transcript_Label, FlyBase_Transcript, Sex, Tissue, FPKM.",
         call. = FALSE)
  }
  flyatlas_adult_transcripts <- readr::read_csv(
    transcript_input, show_col_types = FALSE)
}
# The supplied export includes larval rows; the uploaded analysis targets adults.
if ("Developmental_Stage" %in% names(flyatlas_adult_transcripts)) {
  flyatlas_adult_transcripts <- flyatlas_adult_transcripts %>%
    filter(Developmental_Stage == "Adult")
}
required_transcript_columns <- c(
  "Gene", "Transcript_Label", "FlyBase_Transcript", "Sex", "Tissue", "FPKM")
missing_transcript_columns <- setdiff(
  required_transcript_columns, names(flyatlas_adult_transcripts))
if (length(missing_transcript_columns)) {
  stop("Adult transcript data are missing columns: ",
       paste(missing_transcript_columns, collapse = ", "), call. = FALSE)
}
if (!is.numeric(flyatlas_adult_transcripts$FPKM) ||
    any(!is.na(flyatlas_adult_transcripts$FPKM) &
        (!is.finite(flyatlas_adult_transcripts$FPKM) |
         flyatlas_adult_transcripts$FPKM < 0))) {
  stop("FPKM must contain nonnegative finite numbers or NA.", call. = FALSE)
}

# ==============================================================================
# 1. CREATE ANALYSIS LABELS
# ==============================================================================

adult_isoform_expression <- flyatlas_adult_transcripts %>%
    mutate(
        Isoform_Label = case_when(
            Gene == "Spn42Da" &
                FlyBase_Transcript == "FBtr0086196" ~
                "Spn42Da-B (RB)",

            TRUE ~ paste0(
                Gene,
                "-",
                Transcript_Label
            )
        ),

        log2_FPKM = log2(FPKM + 1)
    )


# ==============================================================================
# 2. QC: VERIFY EXPECTED ADULT STRUCTURE
# ==============================================================================

# Each transcript should have:
# 16 male tissues + 17 female tissues = 33 adult observations

adult_structure_qc <- adult_isoform_expression %>%
    group_by(
        Gene,
        Transcript_Label,
        FlyBase_Transcript,
        Isoform_Label
    ) %>%
    summarize(
        Male_Tissues = n_distinct(
            Tissue[Sex == "Male"]
        ),
        Female_Tissues = n_distinct(
            Tissue[Sex == "Female"]
        ),
        Total_Adult_Observations = n(),
        Missing_FPKM = sum(is.na(FPKM)),
        .groups = "drop"
    )

print(
    adult_structure_qc,
    n = Inf
)

if (
    any(adult_structure_qc$Male_Tissues != 16) ||
    any(adult_structure_qc$Female_Tissues != 17) ||
    any(adult_structure_qc$Total_Adult_Observations != 33)
) {
    stop(
        "Unexpected adult tissue structure detected. ",
        "Inspect adult_structure_qc before continuing."
    )
}


# ==============================================================================
# 3. ISOFORM EXPRESSION SUMMARY
# ==============================================================================

# FPKM > 1 is used here only as a descriptive expression threshold.
# It is NOT used to calculate the later cross-species similarity score.

isoform_summary <- adult_isoform_expression %>%
    group_by(
        Gene,
        Transcript_Label,
        FlyBase_Transcript,
        Isoform_Label
    ) %>%
    summarize(
        Mean_Adult_FPKM =
            mean(FPKM, na.rm = TRUE),

        Median_Adult_FPKM =
            median(FPKM, na.rm = TRUE),

        Max_Adult_FPKM =
            max(FPKM, na.rm = TRUE),

        Tissues_FPKM_gt_1 =
            sum(FPKM > 1, na.rm = TRUE),

        Percent_Adult_Observations_FPKM_gt_1 =
            100 * sum(FPKM > 1, na.rm = TRUE) /
            sum(!is.na(FPKM)),

        .groups = "drop"
    )


# ==============================================================================
# 4. IDENTIFY PEAK-EXPRESSION TISSUE
# ==============================================================================

peak_expression <- adult_isoform_expression %>%
    group_by(
        Gene,
        Transcript_Label,
        FlyBase_Transcript,
        Isoform_Label
    ) %>%
    slice_max(
        order_by = FPKM,
        n = 1,
        with_ties = TRUE
    ) %>%
    ungroup() %>%
    select(
        Gene,
        Transcript_Label,
        FlyBase_Transcript,
        Isoform_Label,
        Peak_Sex = Sex,
        Peak_Tissue = Tissue,
        Peak_FPKM = FPKM
    )


# Join peak information to summary
isoform_summary <- isoform_summary %>%
    left_join(
        peak_expression,
        by = c(
            "Gene",
            "Transcript_Label",
            "FlyBase_Transcript",
            "Isoform_Label"
        )
    )


# ==============================================================================
# 5. SEX-SPECIFIC SUMMARY
# ==============================================================================

isoform_sex_summary <- adult_isoform_expression %>%
    group_by(
        Gene,
        Transcript_Label,
        FlyBase_Transcript,
        Isoform_Label,
        Sex
    ) %>%
    summarize(
        Number_of_Tissues =
            n_distinct(Tissue),

        Mean_FPKM =
            mean(FPKM, na.rm = TRUE),

        Median_FPKM =
            median(FPKM, na.rm = TRUE),

        Max_FPKM =
            max(FPKM, na.rm = TRUE),

        Tissues_FPKM_gt_1 =
            sum(FPKM > 1, na.rm = TRUE),

        Percent_Tissues_FPKM_gt_1 =
            100 * sum(FPKM > 1, na.rm = TRUE) /
            sum(!is.na(FPKM)),

        .groups = "drop"
    )


# ==============================================================================
# 6. WITHIN-ISOFORM STANDARDIZATION
# ==============================================================================

# Standardize across all adult sex × tissue observations for each transcript.
#
# This answers:
# "Where is this transcript relatively enriched compared with its own
# adult expression profile?"
#
# This is for visualization only.

adult_isoform_expression <- adult_isoform_expression %>%
    group_by(
        Gene,
        FlyBase_Transcript
    ) %>%
    mutate(
        Within_Isoform_Z = case_when(
            sd(log2_FPKM, na.rm = TRUE) > 0 ~
                as.numeric(
                    (
                        log2_FPKM -
                        mean(log2_FPKM, na.rm = TRUE)
                    ) /
                    sd(log2_FPKM, na.rm = TRUE)
                ),
            TRUE ~ 0
        )
    ) %>%
    ungroup()


# ==============================================================================
# 7. ORDER TISSUES
# ==============================================================================

# Keep male and female reproductive tissues distinct.
# Common somatic tissues remain comparable between sexes.

tissue_order <- c(
    "Head",
    "Eye",
    "Brain",
    "Thoracicoabdominal ganglion",
    "Crop",
    "Midgut",
    "Hindgut",
    "Malpighian Tubules",
    "Rectal pad",
    "Salivary gland",
    "Fat body",
    "Heart",
    "Carcass",
    "Whole body",
    "Testis",
    "Accessory glands",
    "Ovary",
    "Virgin Spermatheca",
    "Mated Spermatheca"
)

adult_isoform_expression <- adult_isoform_expression %>%
    mutate(
        Tissue = factor(
            Tissue,
            levels = tissue_order
        ),

        Sex = factor(
            Sex,
            levels = c(
                "Male",
                "Female"
            )
        )
    )


# ==============================================================================
# 8. ORDER ISOFORMS BY GENE
# ==============================================================================

isoform_order <- c(
    "Spn42Da-RA",
    "Spn42Da-B (RB)",
    "Spn42Da-RD",
    "Spn42Da-RE",
    "Spn42Da-RF",
    "Spn42Da-RG",
    "Spn42Da-RH",
    "Spn42Da-RI",
    "Spn42Da-RJ",
    "Spn42Da-RK",
    "Spn42Da-RL",
    "Spn42Db-RC",
    "Spn42Db-RD",
    "Spn42Dc-RA",
    "Spn42Dc-RB",
    "Spn42Dd-RB",
    "Spn28B-RA",
    "Spn28B-RB",
    "Spn28F-RA"
)

adult_isoform_expression <- adult_isoform_expression %>%
    mutate(
        Isoform_Label = factor(
            Isoform_Label,
            levels = isoform_order
        )
    )


# ==============================================================================
# 9. FIGURE 1 — ALL ISOFORMS, ABSOLUTE EXPRESSION
# ==============================================================================

p_all_absolute <- ggplot(
    adult_isoform_expression,
    aes(
        x = Tissue,
        y = Isoform_Label,
        fill = log2_FPKM
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.25
    ) +
    facet_grid(
        . ~ Sex,
        scales = "free_x",
        space = "free_x"
    ) +
    scale_fill_viridis_c(
        option = "C",
        name = expression(
            log[2] * "(FPKM + 1)"
        )
    ) +
    labs(
        title = "Adult Drosophila serpin isoform expression",
        subtitle = "Absolute transcript-level expression from FlyAtlas2",
        x = NULL,
        y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
        ),
        strip.text = element_text(
            face = "bold"
        ),
        plot.title = element_text(
            face = "bold"
        )
    )

ggsave(
    file.path(
        figure_dir,
        "04A_All_Isoforms_Absolute_Expression.png"
    ),
    p_all_absolute,
    width = 15,
    height = 8,
    dpi = 600
)


# ==============================================================================
# 10. FIGURE 2 — ALL ISOFORMS, WITHIN-ISOFORM STANDARDIZED EXPRESSION
# ==============================================================================

p_all_standardized <- ggplot(
    adult_isoform_expression,
    aes(
        x = Tissue,
        y = Isoform_Label,
        fill = Within_Isoform_Z
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.25
    ) +
    facet_grid(
        . ~ Sex,
        scales = "free_x",
        space = "free_x"
    ) +
    scale_fill_gradient2(
        midpoint = 0,
        name = "Within-isoform\nz-score"
    ) +
    labs(
        title = "Relative adult tissue distribution of Drosophila serpin isoforms",
        subtitle = "Expression standardized within each transcript",
        x = NULL,
        y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
        ),
        strip.text = element_text(
            face = "bold"
        ),
        plot.title = element_text(
            face = "bold"
        )
    )

ggsave(
    file.path(
        figure_dir,
        "04B_All_Isoforms_Within_Isoform_Z.png"
    ),
    p_all_standardized,
    width = 15,
    height = 8,
    dpi = 600
)


# ==============================================================================
# 11. SPN42DA-ONLY DATASET
# ==============================================================================

spn42da_isoforms <- adult_isoform_expression %>%
    filter(
        Gene == "Spn42Da"
    )


# ==============================================================================
# 12. FIGURE 3 — SPN42DA ISOFORMS, ABSOLUTE EXPRESSION
# ==============================================================================

p_spn42da_absolute <- ggplot(
    spn42da_isoforms,
    aes(
        x = Tissue,
        y = Isoform_Label,
        fill = log2_FPKM
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.3
    ) +
    facet_grid(
        . ~ Sex,
        scales = "free_x",
        space = "free_x"
    ) +
    scale_fill_viridis_c(
        option = "C",
        name = expression(
            log[2] * "(FPKM + 1)"
        )
    ) +
    labs(
        title = "Spn42Da isoform-specific adult tissue expression",
        subtitle = "Spn42Da-B corresponds to RB / FBtr0086196",
        x = NULL,
        y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
        ),
        strip.text = element_text(
            face = "bold"
        ),
        plot.title = element_text(
            face = "bold"
        )
    )

ggsave(
    file.path(
        figure_dir,
        "04C_Spn42Da_Isoforms_Absolute_Expression.png"
    ),
    p_spn42da_absolute,
    width = 15,
    height = 6,
    dpi = 600
)


# ==============================================================================
# 13. FIGURE 4 — SPN42DA WITHIN-ISOFORM STANDARDIZATION
# ==============================================================================

p_spn42da_standardized <- ggplot(
    spn42da_isoforms,
    aes(
        x = Tissue,
        y = Isoform_Label,
        fill = Within_Isoform_Z
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.3
    ) +
    facet_grid(
        . ~ Sex,
        scales = "free_x",
        space = "free_x"
    ) +
    scale_fill_gradient2(
        midpoint = 0,
        name = "Within-isoform\nz-score"
    ) +
    labs(
        title = "Relative tissue distribution of Spn42Da isoforms",
        subtitle = "Expression standardized independently within each transcript",
        x = NULL,
        y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
        ),
        strip.text = element_text(
            face = "bold"
        ),
        plot.title = element_text(
            face = "bold"
        )
    )

ggsave(
    file.path(
        figure_dir,
        "04D_Spn42Da_Isoforms_Within_Isoform_Z.png"
    ),
    p_spn42da_standardized,
    width = 15,
    height = 6,
    dpi = 600
)


# ==============================================================================
# 14. SAVE ANALYSIS TABLES
# ==============================================================================

write_csv(
    adult_isoform_expression %>%
        mutate(
            Tissue = as.character(Tissue),
            Sex = as.character(Sex),
            Isoform_Label = as.character(Isoform_Label)
        ),
    file.path(
        processed_dir,
        "05_FlyAtlas2_adult_isoform_analysis.csv"
    )
)

write_csv(
    isoform_summary,
    file.path(
        table_dir,
        "05_FlyAtlas2_isoform_summary.csv"
    )
)

write_csv(
    isoform_sex_summary,
    file.path(
        table_dir,
        "05_FlyAtlas2_isoform_sex_summary.csv"
    )
)

write_csv(
    adult_structure_qc,
    file.path(
        qc_dir,
        "05_FlyAtlas2_adult_isoform_QC.csv"
    )
)


# ==============================================================================
# 15. DISPLAY KEY RESULTS
# ==============================================================================

cat("\n============================================================\n")
cat("ADULT ISOFORM EXPRESSION ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Number of isoforms analyzed: ",
    n_distinct(adult_isoform_expression$FlyBase_Transcript),
    "\n",
    sep = ""
)

cat(
    "Number of adult observations: ",
    nrow(adult_isoform_expression),
    "\n\n",
    sep = ""
)

cat("ISOFORM SUMMARY\n")
cat("---------------\n")

print(
    isoform_summary %>%
        arrange(
            desc(Max_Adult_FPKM)
        ),
    n = Inf
)


cat("\nSPN42DA-B PROFILE\n")
cat("-----------------\n")

print(
    adult_isoform_expression %>%
        filter(
            Gene == "Spn42Da",
            FlyBase_Transcript == "FBtr0086196"
        ) %>%
        arrange(
            desc(FPKM)
        ) %>%
        select(
            Gene,
            Isoform_Label,
            FlyBase_Transcript,
            Sex,
            Tissue,
            FPKM,
            any_of("SD"),
            log2_FPKM,
            Within_Isoform_Z
        ),
    n = Inf
)
