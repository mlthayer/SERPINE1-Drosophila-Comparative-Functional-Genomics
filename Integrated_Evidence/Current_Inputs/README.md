# Current integration run inputs

The current bulk consensus was rebuilt with step 08C using the latest tissue
inputs processed by step 08B. The supplied original consensus remains under
../Reference_Uploads. Cellular and aging inputs remain unchanged. The aging table
matches the reproduced step 09 results. Bulk composite and cellular correlations now reproduce from the selected
processed tissue and donor/sample summaries. Original source expression
extraction remains separately unverified.
Bulk_Consensus_Update_Comparison.csv records all supplied/rebuilt differences.
All 12 SERPINE1 bulk composite scores change with the selected current inputs.

Step 10 now reads validated step 03 RCL identity scores with 25 SERPINE1
positions via SERPINE1_RCL_SCORES_INPUT (default step 03 processed score table).
It retains the full-length identity calculation, equal averaging of the two
sequence scores and equal averaging of the five evidence domains.
Cellular, orthology, aging and full-length values match historical integration;
RCL/sequence and bulk domains are updated. RCL_Update_Comparison.csv records
the earlier RCL-only update before the subsequent bulk update.
RCL_Update_Comparison.csv records old/new identities, scores and ranks.

Steps 10 and 11 completed. Rankings are recorded in Revised_Final_Integration/Tables after both updates.
Historical outputs remain under ../Reference_Uploads. Current and historical integration exports now have neutral filenames.
Original names and backups remain in private provenance outside the public tree.
This validates integration of supplied processed evidence, not all upstream
source datasets or original human tissue TPM/sex metadata joins.
