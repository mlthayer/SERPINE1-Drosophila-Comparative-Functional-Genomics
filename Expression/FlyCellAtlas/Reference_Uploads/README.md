# Uploaded FlyCellAtlas result validation

Eight author-supplied CSVs are preserved byte-for-byte in this folder.
They provide processed evidence, not verification of the original H5AD extraction.

Verified:
- 28,629 sample-level rows; male/female subset exactly matches 26,733 rows.
- All reported detection percentages match counts; counts are within bounds.
- 2,553 pooled rows reproduce using the current step 07 pooling equations,
  including weighted means and expression-scale handling.
- All 83 overall file/gene detection counts and mean expressions reproduce
  from sample-level summaries.
- Numeric values in all three supplied top-five tables match their pooled rows.
  This checks row values, not every historical ranking/plot rule.
- Generated step 08 input has 1,317 rows using the unchanged filter:
  N_Cells >= 20 and N_Expressing > 0.

15 files have 90 expected gene/file combinations; 83 are present.
Missing combinations (do not convert to biological zeros):

- s_fca_biohub_haltere_10x.h5ad: Spn42Db
- s_fca_biohub_heart_10x.h5ad: Spn28F
- s_fca_biohub_leg_10x.h5ad: Spn28B
- s_fca_biohub_malpighian_tubule_10x.h5ad: Spn28F
- s_fca_biohub_ovary_10x.h5ad: Spn28F
- s_fca_biohub_proboscis_and_maxillary_palps_10x.h5ad: Spn28B
- s_fca_biohub_wing_10x.h5ad: Spn42Db

No local drive, OneDrive, Desktop, /Users or /home paths were found in these CSVs. Source filenames and biological sample IDs are retained for provenance and grouping. Sample labels include study/laboratory names; review their suitability for public release rather than deleting or rewriting them silently. Original H5AD files, exact release provenance and original matrix-scale verification remain outstanding.
