# DIOPT source data and verification

All four CSV files are preserved byte-for-byte from the supplied uploads.

| File | Role | Rows |
| --- | --- | --- |
| `diopt_results_2026-06-13 174731.csv` | Original DIOPT export | 222 |
| `DIOPT_all_human_predictions_for_candidate_serpins.csv` | Complete pairwise predictions; integrated-ranking input | 222 |
| `DIOPT_top3_score_tiers_detailed.csv` | Detailed top-three score tiers, including ties | 67 |
| `Table_DIOPT_candidate_serpin_orthology.csv` | Candidate tier summary; prepared input for step 01 | 6 |

Verification matched raw and prepared prediction pairs, FlyBase and Ensembl IDs, DIOPT and weighted scores, best/reverse-best flags, and supporting-database lists. Supporting-database counts and binary flags agree with the source fields. All detailed top-tier memberships and candidate-summary tier scores, gene sets, and counts agree with the complete prediction table. Each candidate has one SERPINE1 prediction, and all summary SERPINE1 scores/counts/binary flags agree.

Prepared `Confidence_Rank` values are missing for all 222 predictions, while the original export retains categorical ranks in `Rank`. No categorical-to-numeric mapping was invented. This does not block current integrated orthology scoring, which does not use confidence rank.

Step 01 preserves the prepared six-row score-tier summary without treating it as the older primary-ortholog schema. The existing step-10 orthology equations executed successfully on all six SERPINE1 prediction rows. These checks establish the orthology component; other missing datasets still prevent full integrated ranking.

Additional author QC tables are preserved under `Reference_QC/`. All five matched the source data: candidate prediction counts (37 each), SERPINE1 row counts (one each), raw column names, top-three distinct score tiers, and selected-primary-ortholog scores/flags. The selected-primary summary describes the older representation; it does not replace the tie-preserving score-tier summary.
