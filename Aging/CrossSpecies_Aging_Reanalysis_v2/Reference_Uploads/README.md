# Step 09 reference-output comparison

All 16 uploaded original step 09 CSV outputs are preserved unchanged here.
Every table reproduces from the uploaded AFCA/GTEx inputs using the existing
step 09 script. Comparison covers column names, row order, missing values,
text fields, counts, numeric values and rankings; numeric tolerance 1e-10.
`comparison_report.csv` records each result. All 16 comparisons pass.

This includes mappings and coverage, observed and interpolated trajectories,
1,050 matched trajectory rows, 30 overall pair similarities/rankings,
210 system-specific similarities and the six SERPINE1 candidate rankings.
This verifies processed-input reproduction, not original cell-level extraction.
The signed-change methods, standardized age axes, interpolation parameters
and descriptive-only treatment of interpolated observations are unchanged.
