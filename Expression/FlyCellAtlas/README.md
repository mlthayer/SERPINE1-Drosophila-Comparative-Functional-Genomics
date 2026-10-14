# FlyCellAtlas data for step 07

Download the source datasets from the official FlyCellAtlas website:
**<https://flycellatlas.org/>**. Use its data/download links for the 10x H5AD
files. Source H5AD files are downloaded by users to their own computers;
they are not bundled with this repository and do not need to be uploaded here.

Place the files in `Expression/FlyCellAtlas/` and run from the repository root:

```sh
Rscript Scripts/07_FlyCellAtlas_CellExpression.R
```

For a different local folder, set `SERPINE1_FCA_INPUT` before running.
The script uses an explicit 15-file manifest:

- `s_fca_biohub_antenna_10x.h5ad`
- `s_fca_biohub_body_wall_10x.h5ad`
- `s_fca_biohub_fat_body_10x.h5ad`
- `s_fca_biohub_gut_10x.h5ad`
- `s_fca_biohub_haltere_10x.h5ad`
- `s_fca_biohub_heart_10x.h5ad`
- `s_fca_biohub_leg_10x.h5ad`
- `s_fca_biohub_male_reproductive_glands_10x.h5ad`
- `s_fca_biohub_malpighian_tubule_10x.h5ad`
- `s_fca_biohub_oenocyte_10x.h5ad`
- `s_fca_biohub_ovary_10x.h5ad`
- `s_fca_biohub_proboscis_and_maxillary_palps_10x.h5ad`
- `testis.h5ad`
- `s_fca_biohub_trachea_10x.h5ad`
- `s_fca_biohub_wing_10x.h5ad`

Match both the biological dataset and matrix variant, not only the filename.
The supplied analysis uses older flattened or nested `raw.X`/`raw/X` CSC
matrices for most datasets and a testis `X` CSR matrix. A similarly named
processed H5AD may have different gene coverage or expression values.
Do not rename a different matrix variant just to satisfy the manifest.
Exact historical release identifiers and checksums remain to be documented;
the website link alone does not prove release equivalence.

`Reference_Uploads/` preserves the author's eight processed result tables and
records successful consistency checks. These support downstream use without
uploading the large H5AD files. Reproducing source extraction requires users
to obtain the original matching inputs locally. Stored means are not compared
across incompatible expression scales; the primary metric remains percentage
of cells with expression greater than zero.

Large source files are download-only inputs. See [source data and verification scope](../../DATA_AVAILABILITY.md) for official portals and upstream steps not independently rerun here. Uploading large matrices to this chat is not required.
