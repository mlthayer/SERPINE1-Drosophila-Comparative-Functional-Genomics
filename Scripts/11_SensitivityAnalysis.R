# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# 11_SensitivityAnalysis.R
# Reproduce the per-sex leave-one-domain-out analysis from step 10.
#
# Orthology and sequence-alignment evidence are shared across sexes.
# Bulk-expression, cellular-expression, and aging evidence are sex-specific.
# Each sensitivity run removes one of the five evidence domains and
# recalculates the integrated score from the remaining available domains.
#
# The previous six weighting models are retained under Scripts/Legacy/.
library(dplyr)
library(purrr)
integration_dir <- Sys.getenv("SERPINE1_INTEGRATED_OUTPUT",
  unset = "Integrated_Evidence/Revised_Final_Integration")
sensitivity_input <- Sys.getenv("SERPINE1_SENSITIVITY_INPUT",
  unset = file.path(integration_dir, "Tables", "SERPINE1_Integrated_Evidence_BySex_Full.csv"))
if (!file.exists(sensitivity_input)) stop("Run step 10 first, or set SERPINE1_SENSITIVITY_INPUT: ", sensitivity_input, call. = FALSE)
final <- read.csv(sensitivity_input, stringsAsFactors = FALSE)
output_dir <- Sys.getenv("SERPINE1_SENSITIVITY_OUTPUT", unset = integration_dir)
dir.create(file.path(output_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
domain_columns <- c("Orthology_Score", "Sequence_Alignment_Score",
    "Bulk_Expression_Score", "Cellular_Expression_Score", "Aging_Score")
required_columns <- c("Fly_Gene", "Sex", domain_columns)
missing_columns <- setdiff(required_columns, names(final))
if (length(missing_columns)) stop("Sensitivity input is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
sex_counts <- final %>%
  count(Sex)

if (
  !setequal(final$Sex, c("Female", "Male")) ||
  any(sex_counts$n != 6) ||
  anyDuplicated(final[c("Fly_Gene", "Sex")])
) {
  stop(
    "Expected six unique fly candidates for each sex.",
    call. = FALSE
  )
}
if (!all(vapply(final[domain_columns], is.numeric, logical(1)))) stop("Evidence scores must be numeric.", call. = FALSE)
mean_or_na <- function(x) { if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE) }
loo_results <- map_dfr(domain_columns, function(omit_domain) {
    keep_domains <- setdiff(domain_columns, omit_domain)
    tmp <- final %>% select(Fly_Gene, Sex, all_of(keep_domains))
    tmp$Leave_One_Domain_Out_Score <- apply(tmp[, keep_domains,
        drop = FALSE], 1, mean_or_na)
    tmp %>% select(Fly_Gene, Sex, Leave_One_Domain_Out_Score) %>%
        mutate(Omitted_Domain = omit_domain) %>% group_by(Sex) %>%
        mutate(Leave_One_Domain_Out_Rank = rank(-Leave_One_Domain_Out_Score,
            ties.method = "min", na.last = "keep")) %>% ungroup()
})
rank_stability <- loo_results %>% group_by(Fly_Gene, Sex) %>%
    summarise(Best_LOO_Rank = min(Leave_One_Domain_Out_Rank,
        na.rm = TRUE), Worst_LOO_Rank = max(Leave_One_Domain_Out_Rank,
        na.rm = TRUE), Mean_LOO_Rank = mean(Leave_One_Domain_Out_Rank,
        na.rm = TRUE), SD_LOO_Rank = sd(Leave_One_Domain_Out_Rank,
        na.rm = TRUE), .groups = "drop")
write.csv(loo_results, file.path(output_dir, "Tables", "SERPINE1_LeaveOneDomainOut_BySex.csv"), row.names = FALSE)
write.csv(rank_stability, file.path(output_dir, "Tables", "SERPINE1_Rank_Stability_BySex.csv"), row.names = FALSE)
message("Completed existing per-sex sensitivity: ", nrow(loo_results), " domain-omission results.")
