# Human sex-specific expression references

Preserved 23 author-supplied CSVs unchanged. SERPINE1 tissue summary statistics
reproduce from 8,054 GTEx and 186 HPA sample-level records, including counts,
means, medians, SD, IQR, range and log2(x+1) transforms.

Version differences requiring review:
- RAW combined system table: 228 rows; other combined system table: 208 rows.
  The 20 omitted rows are all HPA SERPINB4 (both sexes, ten systems each).
  This is a content difference, not a filename-only correction.
- The older `Human_SERPIN_system_expression_Zscore_by_sex.csv` does not reproduce
  within-dataset/sex/gene standardization (maximum absolute Z error 3.343521).
- `Human_SERPIN_system_expression_Zscore_BY_SEX_CORRECTED.csv` reproduces that
  standardization on its 208 rows; `Human_SERPIN_system_Zscores_by_sex.csv`
  reproduces it on 228 rows. Errors are below 7e-15. Neither is silently
  substituted for the other.

Step 08B accepts the supplied GTEx/HPA tissue tables via
SERPINE1_GTEX_BYSEX_INPUT and SERPINE1_HPA_BYSEX_INPUT, normalizing only the
median_TPM column label. Its existing mappings and equations are preserved.
The full run completes after correcting the named pmap argument adapter.
GTEx and FlyAtlas2 system summaries reproduce the earlier uploaded step 08B
results exactly. HPA matched values differ from those earlier results;
therefore the two sets must remain distinct until source/version differences
are resolved. In particular, six male HPA SERPINB4 correlations were NA in
the earlier results but are finite with the latest tissue tables. This is
not proof either version is the intended publication analysis.

Original source HPA/GTEx aggregation release provenance remains unverified.
These sex-specific bulk inputs do not supply GTEx single-nucleus inputs or
complete the separate four-dataset analysis. Study identifiers are retained
for joins and reproducibility.
