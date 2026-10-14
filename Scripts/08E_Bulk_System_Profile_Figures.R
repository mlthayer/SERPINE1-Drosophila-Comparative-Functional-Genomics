# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Run from the repository root.
# ==============================================================================
# COMPREHENSIVE CROSS-SPECIES BIOLOGICAL-SYSTEM PROFILE FIGURES
#
# PURPOSE:
#
# Create TWO comprehensive figures:
#
#   1. FEMALE:
#      FlyAtlas2 female vs GTEx female and HPA female
#
#   2. MALE:
#      FlyAtlas2 male vs GTEx male and HPA male
#
# FIGURE STRUCTURE:
#
#   Columns = six human SERPIN genes
#   Rows    = GTEx and HPA
#
# Within each panel:
#   - Human SERPIN = thick black solid line
#   - Spn42Da
#   - Spn42Db
#   - Spn42Dc
#   - Spn42Dd
#   - Spn28B
#   - Spn28F
#
# X-axis:
#   Matched biological systems
#
# Y-axis:
#   Relative expression across systems
#
# IMPORTANT:
#   FlyAtlas2 uses FPKM.
#   GTEx/HPA use TPM.
#
# Therefore absolute values are NOT directly compared.
#
# Each gene is normalized within itself across the biological systems
# available for that sex and dataset:
#
#     Relative expression =
#        expression in system / total expression across systems
#
# This allows comparison of tissue-distribution profile shape.
# ==============================================================================


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(scales)
library(tibble)


# ==============================================================================
# 2. PATHS
# ==============================================================================

human_root <- "Expression/Human_SERPIN_reanalysis"


bulk_root <- file.path(
  human_root,
  "CrossSpecies_Bulk_BySex"
)


input_table_dir <- file.path(
  bulk_root,
  "Tables"
)


fly_file <- file.path(
  input_table_dir,
  "FlyAtlas2_system_expression_by_sex.csv"
)


gtex_file <- file.path(
  input_table_dir,
  "GTEx_system_expression_by_sex.csv"
)


hpa_file <- file.path(
  input_table_dir,
  "HPA_system_expression_by_sex.csv"
)


output_root <- file.path(
  human_root,
  "CrossSpecies_Bulk_MultiMetric",
  "Comprehensive_System_Profile_Figures"
)


table_dir <- file.path(
  output_root,
  "Tables"
)


figure_dir <- file.path(
  output_root,
  "Figures"
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


stopifnot(
  file.exists(fly_file),
  file.exists(gtex_file),
  file.exists(hpa_file)
)


# ==============================================================================
# 3. READ SYSTEM-LEVEL DATA
# ==============================================================================

fly_system <- read_csv(
  fly_file,
  show_col_types = FALSE
)


gtex_system <- read_csv(
  gtex_file,
  show_col_types = FALSE
)


hpa_system <- read_csv(
  hpa_file,
  show_col_types = FALSE
)


# ==============================================================================
# 4. GENE ORDERS
# ==============================================================================

fly_genes <- c(
  "Spn42Da",
  "Spn42Db",
  "Spn42Dc",
  "Spn42Dd",
  "Spn28B",
  "Spn28F"
)


human_genes <- c(
  "SERPINE1",
  "SERPINI1",
  "SERPINB3",
  "SERPINB4",
  "SERPINB11",
  "SERPINB9"
)


# ==============================================================================
# 5. SYSTEM ORDERS
# ==============================================================================

female_systems <- c(
  "Heart",
  "Skeletal muscle",
  "Liver-metabolic",
  "Digestive",
  "Renal/excretory",
  "Respiratory",
  "Ovary",
  "Sensory/neural"
)


male_systems <- c(
  "Heart",
  "Skeletal muscle",
  "Liver-metabolic",
  "Digestive",
  "Renal/excretory",
  "Respiratory",
  "Male reproductive",
  "Sensory/neural"
)


# ==============================================================================
# 6. RELATIVE-EXPRESSION FUNCTION
# ==============================================================================

relative_expression <- function(x) {

  total <- sum(
    x,
    na.rm = TRUE
  )


  if (
    !is.finite(total) ||
    total <= 0
  ) {

    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }


  x / total
}


# ==============================================================================
# 7. FIND SYSTEMS THAT ACTUALLY EXIST IN EACH SEX × HUMAN DATASET
#
# We require that a system be represented in:
#
#   FlyAtlas2
#
# AND
#
#   the corresponding human dataset.
#
# This gives the biological systems that can actually be compared.
# ==============================================================================

get_common_systems <- function(
    sex_value,
    human_dataset
) {

  allowed_systems <- if (
    sex_value == "Female"
  ) {
    female_systems
  } else {
    male_systems
  }


  fly_available <- fly_system %>%

    filter(
      Sex == .env$sex_value,
      System %in% allowed_systems
    ) %>%

    distinct(
      System
    ) %>%

    pull(
      System
    )


  human_source <- if (
    human_dataset == "GTEx"
  ) {
    gtex_system
  } else {
    hpa_system
  }


  human_available <- human_source %>%

    filter(
      Sex == .env$sex_value,
      System %in% allowed_systems
    ) %>%

    distinct(
      System
    ) %>%

    pull(
      System
    )


  intersect(
    allowed_systems,
    intersect(
      fly_available,
      human_available
    )
  )
}


# ==============================================================================
# 8. PRINT COMMON SYSTEMS
# ==============================================================================

cat("\n============================================\n")
cat("COMMON BIOLOGICAL SYSTEMS\n")
cat("============================================\n")


for (
  sex_value in c(
    "Female",
    "Male"
  )
) {

  for (
    dataset_value in c(
      "GTEx",
      "HPA"
    )
  ) {

    systems_here <- get_common_systems(
      sex_value,
      dataset_value
    )


    cat(
      "\n",
      sex_value,
      " / ",
      dataset_value,
      ":\n",
      paste(
        systems_here,
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }
}


# ==============================================================================
# 9. PREPARE ONE SEX × DATASET
#
# Each fly gene is normalized once across the shared systems.
# Each human gene is normalized once across the same shared systems.
#
# This is different from pair-specific normalization and is important for this
# comprehensive figure because every fly gene must be directly comparable
# within the same panel.
# ==============================================================================

prepare_comprehensive_profiles <- function(
    sex_value,
    dataset_value
) {

  common_systems <- get_common_systems(
    sex_value,
    dataset_value
  )


  if (
    length(common_systems) < 2
  ) {

    stop(
      paste0(
        "Insufficient common systems for ",
        sex_value,
        " / ",
        dataset_value
      )
    )
  }


  # --------------------------------------------------------------------------
  # Fly profiles
  # --------------------------------------------------------------------------

  fly_profiles <- fly_system %>%

    filter(
      Sex == .env$sex_value,
      Gene %in% fly_genes,
      System %in% common_systems
    ) %>%

    select(
      Gene,
      System,
      Median_FPKM
    ) %>%

    group_by(
      Gene
    ) %>%

    mutate(

      Relative_Expression =
        relative_expression(
          Median_FPKM
        )
    ) %>%

    ungroup() %>%

    transmute(

      Dataset =
        dataset_value,

      Sex =
        sex_value,

      Gene,

      Gene_Type =
        "Drosophila",

      System,

      Relative_Expression
    )


  # --------------------------------------------------------------------------
  # Human profiles
  # --------------------------------------------------------------------------

  human_source <- if (
    dataset_value == "GTEx"
  ) {
    gtex_system
  } else {
    hpa_system
  }


  human_profiles <- human_source %>%

    filter(
      Sex == .env$sex_value,
      Gene %in% human_genes,
      System %in% common_systems
    ) %>%

    select(
      Gene,
      System,
      Median_TPM
    ) %>%

    group_by(
      Gene
    ) %>%

    mutate(

      Relative_Expression =
        relative_expression(
          Median_TPM
        )
    ) %>%

    ungroup() %>%

    transmute(

      Dataset =
        dataset_value,

      Sex =
        sex_value,

      Gene,

      Gene_Type =
        "Human",

      System,

      Relative_Expression
    )


  list(
    fly =
      fly_profiles,

    human =
      human_profiles,

    common_systems =
      common_systems
  )
}


# ==============================================================================
# 10. BUILD ALL PROFILE TABLES
# ==============================================================================

female_gtex <- prepare_comprehensive_profiles(
  sex_value = "Female",
  dataset_value = "GTEx"
)


female_hpa <- prepare_comprehensive_profiles(
  sex_value = "Female",
  dataset_value = "HPA"
)


male_gtex <- prepare_comprehensive_profiles(
  sex_value = "Male",
  dataset_value = "GTEx"
)


male_hpa <- prepare_comprehensive_profiles(
  sex_value = "Male",
  dataset_value = "HPA"
)


# ==============================================================================
# 11. BUILD DATA FOR COMPREHENSIVE FACET FIGURE
#
# Each human SERPIN gets its own panel.
#
# The same six fly profiles are repeated inside each human-gene panel so that
# the reader can directly compare all six candidates with that human SERPIN.
# ==============================================================================

build_figure_data <- function(
    prepared_data
) {

  fly_dat <- prepared_data$fly

  human_dat <- prepared_data$human


  figure_list <- lapply(
    human_genes,
    function(human_gene_value) {

      human_one <- human_dat %>%

        filter(
          Gene == .env$human_gene_value
        ) %>%

        mutate(

          Human_Gene =
            human_gene_value,

          Plot_Gene =
            "Human SERPIN"
        )


      fly_all <- fly_dat %>%

        mutate(

          Human_Gene =
            human_gene_value,

          Plot_Gene =
            Gene
        )


      bind_rows(
        human_one,
        fly_all
      )
    }
  )


  bind_rows(
    figure_list
  )
}


female_plot_data <- bind_rows(

  build_figure_data(
    female_gtex
  ),

  build_figure_data(
    female_hpa
  )
)


male_plot_data <- bind_rows(

  build_figure_data(
    male_gtex
  ),

  build_figure_data(
    male_hpa
  )
)


# ==============================================================================
# 12. FORMAT FACTORS
# ==============================================================================

format_plot_data <- function(
    dat,
    sex_value
) {

  system_levels <- if (
    sex_value == "Female"
  ) {
    female_systems
  } else {
    male_systems
  }


  dat %>%

    mutate(

      Dataset =
        factor(
          Dataset,
          levels = c(
            "GTEx",
            "HPA"
          )
        ),

      Human_Gene =
        factor(
          Human_Gene,
          levels = human_genes
        ),

      System =
        factor(
          System,
          levels = system_levels
        ),

      Plot_Gene =
        factor(
          Plot_Gene,
          levels = c(
            "Human SERPIN",
            fly_genes
          )
        )
    )
}


female_plot_data <- format_plot_data(
  female_plot_data,
  "Female"
)


male_plot_data <- format_plot_data(
  male_plot_data,
  "Male"
)


# ==============================================================================
# 13. SAVE EXACT FIGURE DATA
# ==============================================================================

write_csv(
  female_plot_data,
  file.path(
    table_dir,
    "Female_all_human_all_fly_SERPIN_system_profiles.csv"
  )
)


write_csv(
  male_plot_data,
  file.path(
    table_dir,
    "Male_all_human_all_fly_SERPIN_system_profiles.csv"
  )
)


# ==============================================================================
# 14. LINE/SHAPE DEFINITIONS
#
# Human SERPIN is visually emphasized.
# ==============================================================================

line_types <- c(

  "Human SERPIN" =
    "solid",

  "Spn42Da" =
    "dashed",

  "Spn42Db" =
    "dashed",

  "Spn42Dc" =
    "dashed",

  "Spn42Dd" =
    "dashed",

  "Spn28B" =
    "dashed",

  "Spn28F" =
    "dashed"
)


point_shapes <- c(

  "Human SERPIN" =
    16,

  "Spn42Da" =
    1,

  "Spn42Db" =
    2,

  "Spn42Dc" =
    0,

  "Spn42Dd" =
    5,

  "Spn28B" =
    6,

  "Spn28F" =
    4
)


# ==============================================================================
# 15. COLORS
#
# Human profile is black.
# Each fly serpin has its own consistent color.
# ==============================================================================

gene_colors <- c(

  "Human SERPIN" =
    "black",

  "Spn42Da" =
    "#0072B2",

  "Spn42Db" =
    "#D55E00",

  "Spn42Dc" =
    "#009E73",

  "Spn42Dd" =
    "#CC79A7",

  "Spn28B" =
    "#E69F00",

  "Spn28F" =
    "#56B4E9"
)


# ==============================================================================
# 16. COMPREHENSIVE FIGURE FUNCTION
# ==============================================================================

make_comprehensive_figure <- function(
    dat,
    sex_value
) {

  ggplot(
    dat,
    aes(
      x =
        System,
      y =
        Relative_Expression,
      group =
        Plot_Gene,
      color =
        Plot_Gene,
      linetype =
        Plot_Gene,
      shape =
        Plot_Gene
    )
  ) +

    # ------------------------------------------------------------------------
    # All profiles
    # ------------------------------------------------------------------------

    geom_line(
      aes(
        linewidth =
          Plot_Gene == "Human SERPIN"
      ),
      na.rm = TRUE
    ) +

    geom_point(
      aes(
        size =
          Plot_Gene == "Human SERPIN"
      ),
      na.rm = TRUE
    ) +

    # ------------------------------------------------------------------------
    # 2 datasets × 6 human SERPINs
    # ------------------------------------------------------------------------

    facet_grid(
      rows =
        vars(
          Dataset
        ),
      cols =
        vars(
          Human_Gene
        )
    ) +

    # ------------------------------------------------------------------------
    # Scales
    # ------------------------------------------------------------------------

    scale_color_manual(
      values =
        gene_colors,
      name =
        NULL
    ) +

    scale_linetype_manual(
      values =
        line_types,
      name =
        NULL
    ) +

    scale_shape_manual(
      values =
        point_shapes,
      name =
        NULL
    ) +

    scale_linewidth_manual(
      values = c(
        "TRUE" = 1.15,
        "FALSE" = 0.65
      ),
      guide = "none"
    ) +

    scale_size_manual(
      values = c(
        "TRUE" = 2.8,
        "FALSE" = 1.8
      ),
      guide = "none"
    ) +

    scale_y_continuous(
      labels =
        percent_format(
          accuracy = 1
        ),
      expand =
        expansion(
          mult = c(
            0.02,
            0.08
          )
        )
    ) +

    # ------------------------------------------------------------------------
    # Labels
    # ------------------------------------------------------------------------

    labs(

      title =
        paste0(
          sex_value,
          " cross-species SERPIN biological-system expression profiles"
        ),

      subtitle =
        paste0(
          "Six Drosophila serpins compared with six human SERPINs ",
          "across GTEx and HPA"
        ),

      x =
        "Biological system",

      y =
        "Relative expression across biological systems",

      caption =
        paste0(
          "Each gene was normalized independently across the biological systems ",
          "available for the corresponding sex and human dataset. Human SERPIN ",
          "profiles are shown as thick black lines; colored lines represent the ",
          "six Drosophila serpins. Relative expression describes tissue-profile ",
          "distribution and does not directly equate FlyAtlas2 FPKM with human TPM."
        )
    ) +

    # ------------------------------------------------------------------------
    # Theme
    # ------------------------------------------------------------------------

    theme_classic(
      base_size = 9
    ) +

    theme(

      axis.text.x =
        element_text(
          angle = 55,
          hjust = 1,
          size = 8
        ),

      axis.text.y =
        element_text(
          size = 8
        ),

      axis.title.x =
        element_text(
          size = 10
        ),

      axis.title.y =
        element_text(
          size = 10
        ),

      strip.text.x =
        element_text(
          size = 9,
          face = "bold"
        ),

      strip.text.y =
        element_text(
          size = 9,
          face = "bold"
        ),

      plot.title =
        element_text(
          size = 13,
          face = "bold"
        ),

      plot.subtitle =
        element_text(
          size = 10
        ),

      plot.caption =
        element_text(
          size = 8,
          hjust = 0,
          margin = margin(
            t = 8
          )
        ),

      legend.position =
        "bottom",

      legend.text =
        element_text(
          size = 9
        ),

      legend.key.width =
        unit(
          1.5,
          "cm"
        ),

      panel.spacing =
        unit(
          0.8,
          "lines"
        )
    ) +

    guides(

      color =
        guide_legend(
          nrow = 1,
          byrow = TRUE
        ),

      linetype =
        guide_legend(
          nrow = 1,
          byrow = TRUE
        ),

      shape =
        guide_legend(
          nrow = 1,
          byrow = TRUE
        )
    )
}


# ==============================================================================
# 17. FEMALE FIGURE
# ==============================================================================

female_plot <- make_comprehensive_figure(
  female_plot_data,
  "Female"
)


female_output <- file.path(
  figure_dir,
  "Female_all_human_all_fly_SERPIN_system_profiles.png"
)


grDevices::png(
  filename =
    female_output,
  width =
    18,
  height =
    8.5,
  units =
    "in",
  res =
    600
)


print(
  female_plot
)


grDevices::dev.off()


# ==============================================================================
# 18. MALE FIGURE
# ==============================================================================

male_plot <- make_comprehensive_figure(
  male_plot_data,
  "Male"
)


male_output <- file.path(
  figure_dir,
  "Male_all_human_all_fly_SERPIN_system_profiles.png"
)


grDevices::png(
  filename =
    male_output,
  width =
    18,
  height =
    8.5,
  units =
    "in",
  res =
    600
)


print(
  male_plot
)


grDevices::dev.off()


# ==============================================================================
# 19. OPTIONAL PDF VERSIONS
#
# Vector graphics retain text and line quality for publication.
# ==============================================================================

female_pdf <- file.path(
  figure_dir,
  "Female_all_human_all_fly_SERPIN_system_profiles.pdf"
)


pdf(
  female_pdf,
  width = 18,
  height = 8.5,
  useDingbats = FALSE
)


print(
  female_plot
)


dev.off()


male_pdf <- file.path(
  figure_dir,
  "Male_all_human_all_fly_SERPIN_system_profiles.pdf"
)


pdf(
  male_pdf,
  width = 18,
  height = 8.5,
  useDingbats = FALSE
)


print(
  male_plot
)


dev.off()


# ==============================================================================
# 20. VERIFY OUTPUT
# ==============================================================================

cat("\n============================================\n")
cat("COMPREHENSIVE FIGURES COMPLETE\n")
cat("============================================\n")


cat(
  "\nFemale PNG:\n",
  female_output,
  "\n"
)


cat(
  "\nMale PNG:\n",
  male_output,
  "\n"
)


cat(
  "\nFemale PDF:\n",
  female_pdf,
  "\n"
)


cat(
  "\nMale PDF:\n",
  male_pdf,
  "\n"
)


cat(
  "\nFemale PNG exists:",
  file.exists(
    female_output
  ),
  "\n"
)


cat(
  "Male PNG exists:",
  file.exists(
    male_output
  ),
  "\n"
)