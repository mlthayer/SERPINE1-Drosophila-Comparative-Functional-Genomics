# FlyAtlas 2 gene and transcript downloads

Official interface: <https://motif.mvls.gla.ac.uk/FlyAtlas2/>.
Source accessions and SQL/workbook downloads are documented in
[DATA_AVAILABILITY.md](../../DATA_AVAILABILITY.md#flyatlas-2-data-availability).

## Study-specific direct links

| Gene | FlyBase ID | Gene FPKM export | Transcript FPKM export |
|---|---|---|---|
| Spn28F | FBgn0028987 | [gene](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0028987&tableOut=gene) | [transcripts](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0028987&tableOut=transcriptGene) |
| Spn42Dd | FBgn0028988 | [gene](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0028988&tableOut=gene) | [transcripts](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0028988&tableOut=transcriptGene) |
| Spn42Db | FBgn0033112 | [gene](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0033112&tableOut=gene) | [transcripts](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0033112&tableOut=transcriptGene) |
| Spn42Dc | FBgn0033113 | [gene](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0033113&tableOut=gene) | [transcripts](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0033113&tableOut=transcriptGene) |
| Spn28B | FBgn0083141 | [gene](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0083141&tableOut=gene) | [transcripts](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0083141&tableOut=transcriptGene) |
| Spn42Da | FBgn0265137 | [gene](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0265137&tableOut=gene) | [transcripts](https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn=FBgn0265137&tableOut=transcriptGene) |

The direct endpoint takes a FlyBase gene ID in `fbgn`. `tableOut=gene`
requests gene-level results; `tableOut=transcriptGene` requests transcript-level
results for protein-coding genes. The other supported types, `mir` and
`transcriptMir`, are not used for this six-serpin analysis.

For browser search, use:
`https://motif.mvls.gla.ac.uk/FlyAtlas2/index.html?search=gene&gene=GENE_ID&idtype=ID_TYPE`.
`ID_TYPE` can be `fbgn`, `cgnum`, `symbol`, or `name`.
For example, [Spn42Da by FlyBase ID](https://motif.mvls.gla.ac.uk/FlyAtlas2/index.html?search=gene&gene=FBgn0265137&idtype=fbgn).

## Local inputs and reproduction

Step 05A reads `05_FlyAtlas2_gene_expression_wide.csv` in this directory,
or exactly six original gene text exports. `SERPINE1_FLYATLAS_GENE_INPUT`
selects another raw-export directory; `SERPINE1_FLYATLAS_GENE_WIDE_INPUT`
selects another prepared CSV. Preserve the schema checked by the importer;
do not substitute transcript-level exports for gene totals.

Original transcript exports are retained in `../FlyAtlas2_TranscriptTables/`.
Step 05B reads the matching all-transcript CSV in
`../FlyAtlas2_Analysis/Processed/` and selects adult observations; step 05C
reads its processed adult isoform results. The supplied six transcript text
exports and combined table passed value-by-value verification.

Download links were supplied by the author; their current server responses
and release equivalence were not independently verified here. A new release
may change expression values, so preserve original inputs when reproducing
historical results. No analytical values or calculations changed in this update.
