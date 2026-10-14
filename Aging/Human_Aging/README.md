# Human aging inputs for step 09

Original adult GTEx reference files are in Reference_Uploads. The required
GTEx_v11_SERPIN_change_from_20s.csv is preserved unchanged in Processed,
where step 09 reads it. Five aging genes are SERPINB1, SERPINB3, SERPINB9,
SERPINE1 and SERPINI1, matching the supplied aging script's gene set.

Verified 98,080 sample-level TPM values against the five-row GCT (19,616
samples). All 2,935 gene/tissue/age/sex summary rows reproduce: sample counts,
means, medians, SD, SE and mean log2(TPM+1). All 1,325 change rows retain
Mean_log2TPM minus the 20–29 baseline. Their age-level log2 means are the
unweighted mean of available sex-specific means, not sample-count-weighted
pooled means. This supplied method is preserved; Total_N sums sample counts.
Focused heatmap and supplementary tables are retained but not fully validated.

The original selected-gene GCT incorrectly declares 74,628 rows although it
contains five. Processed/GTEx_v11_Selected_SERPINs.gct corrects only that
row-count declaration; original sample headers and TPM values are preserved.

Developmental dGTEx metadata under Annotations are separate source materials,
not adult-aging input replacements. Its large TPM archive and the adult GTEx
SampleAttributesDS upload exceeded the 32 MiB transfer limit. The already
annotated adult sample table supports summary validation without those files,
but upstream metadata joins and developmental analysis were not verified.

Remaining step 09 input: the fly AFCA file
Aging/Fly_Aging_tissue_expression/Data/Processed/AFCA_serpin_aging_change_summary_table.csv,
containing Gene, cell_type_clean, 30d_vs_5d, 50d_vs_5d and 70d_vs_5d.
The cross-species run has not been performed without this file.
