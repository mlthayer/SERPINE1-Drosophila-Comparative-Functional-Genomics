# Current human system-expression inputs

Author-selected current version: all 228 rows, including 20 HPA SERPINB4
system rows across both sexes. The raw and Z-score tables here are byte-identical
copies of the explicitly recorded reference sources in selection_manifest.json.
The Z-score table was verified within Dataset/Sex/Gene. The 208-row versions
remain historical references and are not the selected publication inputs.

Step 08B uses the latest supplied GTEx_SERPIN_tissue_expression_by_sex.csv
and HPA_SERPIN_tissue_expression_by_sex.csv from Reference_Uploads by default.
It maps and aggregates those tissue inputs with its existing definitions;
its bulk Spearman correlations use unstandardized expression, not these Z scores.
The broad-system 228-row tables and step 08B's comparative system tables
have different mapping definitions and must not be substituted for each other.
