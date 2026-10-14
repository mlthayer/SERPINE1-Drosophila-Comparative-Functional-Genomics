# Uploaded step 08B bulk expression results

Ten author-supplied CSVs are retained unchanged, including the original
abbreviated filenames `MALE_A~1.CSV` and `FEMALE~1.CSV`. Their columns identify
them as relative-expression plotting data, not original expression matrices.

Validated from the 828 matched system-profile rows:
- All 144 Spearman correlations and p-values reproduce using exact=FALSE.
- BH FDR values reproduce within each dataset/sex group.
- The 24 SERPINE1-only correlation rows match the full results.
- Matched fly values agree with the 84-row FlyAtlas2 system summary.
- Matched GTEx values agree with the 78-row GTEx system summary.

Comparisons use six matched systems for GTEx female, five for GTEx male,
and six for both HPA sexes. Thirty correlations are undefined and retained
as NA: 12 GTEx female, 12 GTEx male and six HPA male. Undefined correlations
are not biological zero correlation and must not be filled with zeros.

The mapping tables describe 18 human and 19 fly tissues. They do not by
themselves verify original file availability or upstream aggregation.
Full step 08B reproduction still requires the sex-specific human tissue
summary inputs named in the script and `HPA_system_expression_by_sex.csv`
for checking HPA system aggregation. The supplied matched HPA values permit
correlation validation but do not validate the source HPA aggregation.
These bulk results do not supply the separate step 08 four-dataset or
GTEx single-nucleus inputs. Upstream release provenance remains outstanding.
