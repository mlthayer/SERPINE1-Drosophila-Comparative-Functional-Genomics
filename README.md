# SERPINE1–Drosophila Comparative Functional Genomics Workflow

This workflow compares Drosophila serpins with human SERPINE1/PAI-1 using orthology, full-length and reactive-center-loop (RCL) sequences, phylogenetics, tissue and cell expression, aging profiles, and integrated evidence ranking. Expanded human SERPIN and isoform comparisons are retained as supplied.

## Requirements

Use R 4.5 or a compatible version. The scripts call these CRAN packages:

```r
install.packages(c("dplyr", "stringr", "readr", "tidyr", "ggplot2",
                   "tibble", "scales", "ape", "phangorn", "purrr", "hdf5r", "data.table"))
```

Sequence analyses also require Bioconductor packages:

```r
install.packages("BiocManager")
BiocManager::install(c("Biostrings", "msa", "ggmsa", "rhdf5"), ask = FALSE, update = FALSE)
```

Optional FCA inspection helpers additionally require `zellkonverter`,
`SingleCellExperiment`, and `SummarizedExperiment`.

Installation commands are guidance, not a claim that every dependency was installed or tested. `sessionInfo.txt` records a historical analysis environment, not a dependency lockfile. Record a new `sessionInfo()` with the final analysis.

## Portable data layout

Run from the repository root. Scripts use relative paths; no local project-directory environment variables are needed. The supplied `05_FlyAtlas2_all_transcript_expression.csv` is included beneath `Expression/FlyAtlas2_Analysis/Processed/`. Other original inputs are still needed. Retain original filenames and the column formats checked by each script.

| Directory | Purpose |
| --- | --- |
| `Orthology/Inputs/` | Supplied raw DIOPT export, complete pairwise predictions, detailed tiers, and candidate tier summary |
| `Data/` | Optional older `orthology_mapping.csv`; orthology outputs in `Processed/` |
| `Sequence/FullLength/` | `Uniport FASTA FlySerpins_Human PAI-1.txt`; full-length alignment outputs |
| `Sequence/RCL/` | `selected_human_fly_Serpins_Fasta.txt`; RCL and extracted-sequence sensitivity outputs |
| `Sequence/Phylogeny/` | Phylogenetic outputs; reads the master FASTA and processed homologous RCL table from `Sequence/RCL/` |
| `Expression/Data/Processed/` | Preprocessed `FlyAtlas2_serpin_tissue_summary.csv` for the supplied gene-level summary script |
| `Expression/FlyAtlas2_Analysis/` | Supplied transcript CSV, adult isoform outputs, and gene-level import outputs |
| `Expression/FlyAtlas2_TranscriptTables/` | Six supplied original transcript-level exports; verified against the all-transcript CSV |
| `Expression/FlyAtlas2_GeneTables/` | Supplied gene-level wide CSV for step 05A; original text exports can also be used |
| `Expression/GTEx/` | Original GTEx input data selected by the human-expression script |
| `Expression/HPA/` | Original HPA input data selected by the human-expression script |
| `Expression/Human_SERPIN_reanalysis/` | Human-expression outputs and both cross-species comparison output directories |
| `Expression/FlyCellAtlas/` | Expected H5AD inputs; single-cell outputs in `Reanalysis_Output/` |
| `Aging/` | Fly and human aging inputs in the subdirectories specified by the aging script; outputs in `CrossSpecies_Aging_Reanalysis_v2/` |
| `Integrated_Evidence/Revised_Final_Integration/` | Replacement integrated evidence outputs |

Analytical file contents are preserved. FlyAtlas2 filenames use step 05; integrated outputs and historical copies now use neutral publication-facing names. Current RCL identity follows the validated 25-position reference method; current bulk consensus uses the selected latest tissue inputs. Rerunning can overwrite generated outputs.

## Script inventory and execution

Run individual standalone steps from the repository root, for example:

```sh
Rscript Scripts/01_DIOPT_Orthology.R
```

Do not run all files blindly in a loop: the uploaded workflow has continuations, interactive steps, and unresolved prerequisites.

| Script | Analysis and execution notes |
| --- | --- |
| `01_DIOPT_Orthology.R` | DIOPT orthology summary; preserves the supplied tier summary, with the older mapping route retained |
| `02_Protein_Alignment_FullLength.R` | FASTA validation, multiple alignment, sequence identity and BLOSUM62 similarity |
| `03_RCL_Alignment.R` | Reconstructed SERPINE1 353–377 homologous-region extraction, supported-site QC, and identity/BLOSUM62 scoring |
| `03B_RCL_ExtractedSequence_Sensitivity.R` | Sensitivity analysis using the supplied extracted RCL sequences |
| `04_Phylogenetic_Analysis.R` | Full-length and RCL neighbor-joining and maximum-likelihood analyses |
| `05_FlyAtlas2_TissueExpression.R` | Supplied gene-level FlyAtlas2 summary; does not generate the transcript-level input needed below |
| `05A_FlyAtlas2_GeneTable_Import.R` | Imports and validates original gene-level text exports or the supplied wide CSV |
| `05B_Adult_Isoform_TissueExpression.R` | Loads the supplied transcript CSV and selects adult rows, or reuses an upstream adult object |
| `05C_WithinGene_Isoform_TissueSimilarity.R` | Loads the processed adult isoform table from 05B; also supports sourcing in the same session |
| `06_GTEx_HPA_SERPINE1_TissueExpression.R` | Expanded human SERPIN expression; `--gtex-only` supports partial bulk validation without HPA |
| `06B_GTEx_Bulk_BySex.R` | GTEx TPM/subject-sex join and summaries |
| `06C_HPA_Bulk_BySex.R` | Verified HPA sample-sex join and summaries |
| `07_FlyCellAtlas_CellExpression.R` | Multi-tissue H5AD single-cell reanalysis |
| `07B_GTEx_snRNA_SERPIN_Extraction.R` | Supplied human single-nucleus extraction; needs the original atlas H5AD |
| `08_CrossSpecies_TissueSimilarity.R` | Four-resource cross-species comparisons; input discovery requires unambiguous matching files |
| `08B_SexSpecific_Bulk_TissueSimilarity.R` | Sex-specific FlyAtlas2–GTEx/HPA comparisons |
| `08C_Bulk_MultiMetric_TissueSimilarity.R` | Descriptive bulk composite scores and cross-dataset consensus |
| `08D_FCA_GTEx_snRNA_BySex.R` | Sex-specific cellular comparisons from donor/sample summaries |
| `08E_Bulk_System_Profile_Figures.R` | Expanded bulk profile plotting workflow |
| `09_Aging_TissueExpression_Comparison.R` | Matched-system cross-species aging reanalysis |
| `10_IntegratedEvidenceRanking.R` | Replacement integrated evidence analysis; configure six input paths as described below, or select interactively |
| `11_SensitivityAnalysis.R` | Standalone rerun of the existing per-sex leave-one-domain-out method using the full table from step 10 |

Run the isoform steps as separate commands from the repository root:

```sh
Rscript Scripts/05B_Adult_Isoform_TissueExpression.R
Rscript Scripts/05C_WithinGene_Isoform_TissueSimilarity.R
```

Step 05B selects `Developmental_Stage == "Adult"` from the supplied all-stage table and retains the existing requirement of 16 male and 17 female observations per transcript. `SERPINE1_ADULT_TRANSCRIPT_INPUT` can select another original CSV. `SERPINE1_FLYATLAS_OUTPUT` overrides the shared output directory. Step 05A uses the supplied `Expression/FlyAtlas2_GeneTables/05_FlyAtlas2_gene_expression_wide.csv` when no raw text exports are present. `SERPINE1_FLYATLAS_GENE_WIDE_INPUT` selects a prepared CSV explicitly. It can also read exactly six original text exports from the same directory, overridable with `SERPINE1_FLYATLAS_GENE_INPUT`. Its parser and QC are preserved; no gene expression values are inferred from transcript expression.

## Noninteractive integration and sensitivity

Set these environment variables to the original input files, using paths relative to the repository root when possible:

| Variable | Input |
| --- | --- |
| `SERPINE1_DIOPT_INPUT` | `Orthology/Inputs/DIOPT_all_human_predictions_for_candidate_serpins.csv` |
| `SERPINE1_FULL_ALIGNMENT_INPUT` | Full-length aligned human–fly FASTA |
| `SERPINE1_RCL_ALIGNMENT_INPUT` | RCL aligned human–fly FASTA |
| `SERPINE1_RCL_SCORES_INPUT` | Validated step 03 scores; defaults to the processed score table |
| `SERPINE1_BULK_INPUT` | Bulk consensus by sex CSV |
| `SERPINE1_CELL_INPUT` | FCA–GTEx snRNA Spearman by sex CSV |
| `SERPINE1_AGING_INPUT` | Overall cross-species aging similarity CSV |

Then run `Rscript Scripts/10_IntegratedEvidenceRanking.R`. Missing configuration fails with the required variable name instead of opening a dialog in unattended runs. In interactive R, file selection remains available. `SERPINE1_INTEGRATED_OUTPUT` overrides the default `Integrated_Evidence/Revised_Final_Integration/`.

Run `Rscript Scripts/11_SensitivityAnalysis.R` after step 10. It reproduces the same per-sex leave-one-domain-out equations, ranks, and stability summaries from the exported full evidence table and writes the same sensitivity output basenames. Optional `SERPINE1_SENSITIVITY_INPUT` and `SERPINE1_SENSITIVITY_OUTPUT` select input and output paths. Recomputing from CSV can expose floating-point differences near numerical ties; step 10's in-memory results remain the primary sensitivity record.

The older six weighting models are retained unchanged in `Scripts/Legacy/11_AlternativeWeighting_Sensitivity.R` for use with the earlier weighted evidence table. They are not silently mapped to new evidence domains.

## Reproducibility status

The processed bulk, cellular, aging, integrated ranking and sensitivity routes
were rerun from an isolated source package. See
[reproduction instructions](Reproducibility/README.md) and
[data availability and verification limits](DATA_AVAILABILITY.md).
Human reactive sites are supplied sequence annotations without independent
primary-source verification. The historical expanded-RCL ML tree is selected
with the rerun difference disclosed in
[the selection note](Sequence/Phylogeny_Reference_Uploads/SELECTED_EXPANDED_RCL.md).
CITATION.cff retains attribution and the workflow title; release metadata will
be added when an actual release is created.

## Supplied DIOPT data

All four supplied DIOPT files are retained unchanged under `Orthology/Inputs/`. Step 01 uses the supplied score-tier summary when the older mapping input is absent; `SERPINE1_DIOPT_SUMMARY_INPUT` can select that summary explicitly. It preserves tier ties and does not infer a single primary ortholog. For integrated ranking, set `SERPINE1_DIOPT_INPUT=Orthology/Inputs/DIOPT_all_human_predictions_for_candidate_serpins.csv`; the six-row tier summary is not the pairwise input.

Raw and prepared predictions agree for all 222 pairs. Prepared `Confidence_Rank` values are missing; the categorical raw ranks remain in the original export. Integrated orthology scoring does not use confidence rank. See `Orthology/Inputs/README.md` for verification details.

## Full-length sequence reproduction

The original step-02 FASTA and seven accompanying alignment/statistics/QC files are retained under `Sequence/FullLength/`. `Rscript Scripts/02_Protein_Alignment_FullLength.R` completed and reproduced both uploaded FASTA outputs and all three uploaded CSV results. The older fly-only Clustal file includes an extra accession (`Q8MSJ3`), retained as a reference but excluded from the six-candidate analysis. See `Sequence/FullLength/README.md` for details. The expanded RCL/phylogeny input FASTA and homologous-region references have also been supplied; see `Sequence/RCL/README.md` for validation and extraction-method limitations.

## Expanded RCL references

All expanded sequence/RCL uploads are retained under `Sequence/RCL/Reference_Uploads/`. Step 03 reconstructs the validated SERPINE1-referenced homologous-region method and reproduces the supplied alignments, regions, scores, and inventories. The earlier motif-based script is retained under `Scripts/Legacy/`. Step 03B completed.

Step 04 uses the validated homologous-region table from `Reference_Uploads/` when generated input is absent. `SERPINE1_PHYLOGENY_RCL_INPUT` provides an explicit override. Model choices, bootstrap count, and random seed remain unchanged.

Step 04 completed all four comparisons with the original seed, JTT models, optimization settings, and 500 NJ/500 ML bootstrap replicates per comparison. Eight exported trees and four model summaries passed validation. The four uploaded phylogeny alignment references match the generated alignments; historical tree/support equivalence remains unverified without historical tree results.

## Reconstructed RCL extraction

Run `Rscript Scripts/03_RCL_Alignment.R` to extract regions from the validated alignment, or append `--realign` to regenerate the full-length ClustalOmega alignment from the original FASTA. Both modes reproduce the supplied results. Run `Rscript tests/test_rcl_reconstruction.R` for regression checks. See `METHODS.md` and `Sequence/RCL/README.md` for explicit boundaries, score definitions, supported-site QC, provenance, and path configuration. The original unavailable generator has been replaced by a clearly documented reconstruction; its calculations are validated against the supplied outputs.

## GTEx inputs received

Three of the 18 supplied bulk files transferred and passed validation. The remaining 15 exceed the upload transfer limit. Run step 06 with `--gtex-only` for the available tissues; HPA is not run in that mode. Compact six-gene GCT extracts and original metadata are retained for the three validated tissues. These extracts reproduce the same target-gene summaries. See `Expression/GTEx/README.md` for sample counts, source-data metadata flags, direct-download access, and single-nucleus prerequisites. Complete tissue coverage and the HPA archive are still needed for complete downstream comparisons.

HPA data are now downloaded and validated; see `Expression/HPA/README.md`. For a compact reproducible partial run, set `SERPINE1_GTEX_INPUT=Expression/GTEx/SERPIN_Subsets` and `SERPINE1_HPA_INPUT=Expression/HPA/SERPIN_Subsets` when running step 06. GTEx's remaining bulk tissues and snRNA atlas require outbound HTTPS access to the official download host `storage.googleapis.com`; exact URLs are recorded in the GTEx validation manifest.

### Download and select human expression inputs locally

Download GTEx v11 **gene-level TPM** tissue GCT files from
<https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression>.
Download HPA **transcript_rna_tissue.tsv.zip** from
<https://www.proteinatlas.org/about/download>; keep the ZIP compressed.
For the single-nucleus analysis, download
**GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad** from
<https://gtexportal.org/home/downloads/adult-gtex/single_cell>.

Run R/RStudio from the repository root:

```r
options(serpine1.choose_files = TRUE)
source("Scripts/06_GTEx_HPA_SERPINE1_TissueExpression.R")
# Run separately when the single-nucleus atlas is available:
source("Scripts/07B_GTEx_snRNA_SERPIN_Extraction.R")
options(serpine1.choose_files = FALSE)
```

Step 06 prompts for the tissue GCT files and then the HPA archive. Windows
supports selecting multiple tissue files together; other systems select one
file at a time and answer whether to add another. Step 07B prompts for one
atlas. Cancelling stops the run. Select every tissue needed for the study;
selection of fewer files produces partial coverage. Use the same data releases
as the analysis being reproduced. Downloads and selections stay on your computer.
No personal paths are written into source code.

For automated `Rscript` runs, place downloads in the documented relative input
folders or configure the existing input environment variables. Interactive
selection requires an R/RStudio session.

### FlyCellAtlas source downloads

Download step 07 H5AD inputs locally from **https://flycellatlas.org/**. They are not distributed with this repository and need not be uploaded to the chat. The [FlyCellAtlas input guide](Expression/FlyCellAtlas/README.md) lists all 15 required files, matrix variants, input-folder configuration and run instructions. Supplied processed tables remain available for downstream analysis.

### Current human expression version

The selected human system tables contain **228 rows**, including HPA SERPINB4. See [current input selection](Expression/Human_SERPIN_reanalysis/Current_Inputs/README.md). Step 08B defaults to the latest uploaded sex-specific tissue tables; earlier result tables remain references.

### Rebuild current bulk consensus and integration

Step `08B_SexSpecific_Bulk_TissueSimilarity.R` defaults to the selected uploaded
GTEx/HPA sex-specific tissue summaries. Then run
`08C_Bulk_MultiMetric_TissueSimilarity.R` to calculate the descriptive mean of
available Spearman, Kendall, cosine and top-three-system Jaccard components,
and the cross-dataset consensus. The component metrics are not independent
statistical evidence. Step 10 uses the rebuilt consensus in
`Integrated_Evidence/Current_Inputs/`; its supplied predecessor remains a
historical reference. `Bulk_Consensus_Update_Comparison.csv` records differences.

Additional source preparation scripts are `06B_GTEx_Bulk_BySex.R` and
`06C_HPA_Bulk_BySex.R`. Complete original GTEx tissue inputs are needed to verify all source joins.
The HPA sample-sex join has been executed and verified; its metadata are retained. `08D_FCA_GTEx_snRNA_BySex.R` provides the
cellular comparison workflow; it now reproduces the supplied correlations from the uploaded human
donor-system summary and FCA sample summaries. Original H5AD expression
extraction remains unverified. `08E_Bulk_System_Profile_Figures.R` contains the supplied
expanded plotting workflow and has passed syntax checks.

### Source-data availability and verification limits

Large expression matrices are **download-only inputs**; obtain them from the
linked official portals and run upstream extraction locally. No large-file
upload to this chat is required. See [data availability and verification
scope](DATA_AVAILABILITY.md) for source links, local inputs and the distinction
between reproduced processed analyses and source extraction or metadata joins
that were not independently rerun here.

## References

[Workflow-specific references](References/WORKFLOW_REFERENCES.md) cover the software, public datasets and sequence annotations used here. Installed software citation records are retained in References/software_citations.bib. General background and manuscript bibliography items are excluded.

## FlyAtlas 2 source access

FlyAtlas 2 source RNA-seq projects are **PRJEB22205** (main), **PRJEB11865**
(midgut), and **PRJEB48667** (garland cells). The [data availability guide](DATA_AVAILABILITY.md#flyatlas-2-data-availability)
provides the SQL/workbook snapshots; the [FlyAtlas 2 input guide](Expression/FlyAtlas2_GeneTables/README.md)
provides direct gene and transcript download links for all six candidate serpins.
