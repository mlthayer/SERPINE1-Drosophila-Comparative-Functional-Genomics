# SERPINE1-referenced homologous RCL workflow

Step 03 has been reconstructed from the validated alignment, region table, original FASTA, and author-provided methods. It uses the explicit SERPINE1 P17–P8′ interval, residues 353–377, containing R369–M370. In the current alignment this maps to columns 402–433: 25 reference residues across 32 columns.

From the repository root:

```sh
Rscript Scripts/03_RCL_Alignment.R
Rscript Scripts/03_RCL_Alignment.R --realign
Rscript tests/test_rcl_reconstruction.R
```

The default run uses the validated full-length alignment. `--realign` regenerates it from the original FASTA using the same ClustalOmega method and input order. Both modes reproduce all 12 regions and score rows and both supplied alignment FASTAs. All 12 supplied CSV references, including the additional reactive-site, sequence, interval, and gap-count QC tables, also match. Inventories match after ordering by accession. Regression checks verify repeatability, supported reactive-site alignment, and rejection of an inconsistent reference interval.

Inputs are retained unchanged under `Reference_Uploads/`, with the original FASTA at `selected_human_fly_Serpins_Fasta.txt`. Optional variables `SERPINE1_RCL_REFERENCE_DIR`, `SERPINE1_RCL_FASTA_INPUT`, and `SERPINE1_RCL_OUTPUT` configure paths. Generated outputs use `Processed/` and `QC/`; uploaded references are preserved separately. Output basenames match the supplied homologous-region results. Additional QC records document reactive-site checks, reconstruction provenance, and software versions.

Supported human reactive sites and Spn42Da-B are used only to check alignment accuracy; all six match within one column of SERPINE1's bond. Predicted fly sites and manual bold annotations are not used for extraction or scoring. See `METHODS.md` for the publication-facing method and citation/numbering items requiring author review.

The earlier motif-based workflow is retained under `Scripts/Legacy/03_MotifBased_RCL_Extraction.R`. Its DHPF requirement fails for Spn42Da and is not part of the primary reconstructed analysis. The marked older fly-only alignment is retained as annotation and includes Q8MSJ3, while lacking SERPINE1 and Spn42Da-B.

Step 03B's extracted-sequence sensitivity analysis completed, and its embedded current scores agree with the supplied homologous-region results. Step 04 previously completed all eight trees with the same validated region table and original seed, models, and 500-bootstrap settings. The reconstructed table matches that phylogeny input, so no phylogeny rebuild is required by this change.
