# GTEx single-nucleus processed references

Six author-supplied CSVs are preserved unchanged: 209,126 cell-metadata rows,
150 donor/tissue/gene summaries, 114 donor/system/gene summaries, 54 sex/system
summaries, nine donor-coverage rows and six gene checks.

Verified detection percentages and count bounds, donor-tissue to donor-system
aggregation (excluding unmapped systems as in step 07B), weighted means,
sex-specific donor statistics, coverage and agreement of every donor-system
nucleus count with cell metadata. Step 08D reran using the supplied FCA
sample-level table and GTEx donor-system table; all 72 correlation output rows
match the previously supplied cellular input including p-values/FDR and labels
(numeric tolerance 1e-10). No integrated evidence scores changed.

Original H5AD expression extraction is still not independently rerun. Cell
metadata alone cannot validate counts of expression-positive nuclei or raw
gene counts. Biological participant/sample identifiers remain for grouping;
no machine-specific paths are required.

To rerun step 08D, set SERPINE1_FCA_SUMMARY_INPUT to the preserved FCA
SampleLevel_MaleFemale_Only CSV and SERPINE1_GTEX_DONOR_INPUT to this folder's
GTEx_snRNA_SERPIN_donor_system_summary.csv, then run the script from the root.
