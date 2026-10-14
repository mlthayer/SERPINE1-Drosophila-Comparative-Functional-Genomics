# Fly aging processed inputs

Eleven uploaded AFCA result tables are preserved in Reference_Uploads.
Processed/AFCA_serpin_aging_change_summary_table.csv is an unchanged copy
of the 78-row primary input for step 09.

Validated: log2(mean_expression+1), changes from the 5-day baseline,
correspondence of available change-summary values to the heatmap export,
and the exact 13-row Spn42Dd subset. The no-5d export has the same 228
nonbaseline observation keys and unmodified change values, but 16
change_clipped display values differ between heatmap versions. These are
retained and flagged; step 09 uses the unmodified change summary.

The uploaded single-nucleus long table exceeds the 32 MiB transfer limit.
Original cell-level aggregation and broad/detailed annotation summaries
have not been independently reproduced. No replacement cell data were invented.

Step 09 completed with these inputs and the validated adult GTEx change
summary. Generated 16 tables and figure outputs. All 30 overall gene-pair
integrated aging scores and six SERPINE1 rankings are finite. The existing
many-to-many join warning is expected when each fly trajectory is replicated
across the five human gene pairs; calculations and join behavior were preserved.
No inferential p-values are computed from interpolated trajectory points.
