# RCL conservation methods

RCL conservation was evaluated using the experimentally established human SERPINE1 RCL as a reference. The SERPINE1 P17–P8′ interval (residues 353–377), containing the R369–M370 P1–P1′ reactive bond, was mapped to the full-length Clustal Omega alignment. Corresponding alignment columns were extracted for each human and *Drosophila* serpin to define SERPINE1-referenced homologous RCL regions. Insertions and deletions expanded the 25-residue SERPINE1 interval to 32 alignment columns in the validated dataset.

Alignment consistency was checked using supplied reactive-site annotations for human serpins and Spn42Da-B, which mapped within one alignment position of the SERPINE1 reactive bond. Sequence-coordinate agreement does not independently establish experimental support for those annotations; the five human annotations were not independently checked against primary literature, and no such verification is claimed. RCL conservation was quantified relative to the 25 SERPINE1 residue positions. Sequence identity represented exact amino-acid matches, while biochemical similarity included identical residues and substitutions with positive BLOSUM62 scores. Candidate gaps opposite SERPINE1 residues were treated as nonmatches, whereas positions containing SERPINE1 gaps were excluded from the denominator and evaluated separately as insertion/deletion differences. Predicted *Drosophila* P1–P1′ sites were not used to define RCL regions or calculate similarity. The analysis represents SERPINE1-referenced homologous RCL conservation and does not assign reactive sites to the other fly candidates.

## Reconstruction and provenance

The original script generating the supplied homologous-region outputs was unavailable. The executable workflow was reconstructed from the validated full-length alignment, region table, original FASTA, and author-provided methods definition. Both extraction from the validated alignment and regeneration using `msa::msa(..., method = "ClustalOmega", order = "input")` reproduce the supplied regions and scores. The validated region table supplies record identities, accessions, display labels, species, and existing reactive-site annotations; computed sequences and scores are generated anew. Its interval is audited against the explicit SERPINE1 353–377 definition.

Supported-site QC uses the full-length FASTA coordinates retained from the supplied script: SERPINE1 R369–M370, SERPINB3 S354–S355, SERPINB4 L354–S355, SERPINB11 K357–S358, SERPINI1 R362–M363, and Spn42Da-B R374–A375. Ellisdon et al. (2014; doi:10.1186/1472-6807-14-14, pp. 2–4 and 6–7) report a 392-residue recombinant protein bearing the shared A/B RCL, identify P1 R342 and P1′ Ala, and describe structural evidence consistent with R–A cleavage in crystal form 1. Their paper/construct numbering maps to R374–A375 in the supplied 424-residue full-length B sequence A1Z6R1 with a +32 offset. The paper describes an alternative N-terminal signal-peptide exon in isoform B; this coordinate conversion is not a claim that signal-peptide cleavage at residue 32 was experimentally measured. Crystal form 2 is described as cleaved after P2 K341, so the two forms are not both evidence of the identical bond. See [the primary-source audit](Sequence/RCL/Citation_Audit/Spn42DaB_PRIMARY_SOURCE.md). The selected SERPINE1 P17–P8′ endpoints are derived from the annotated P1/P1′ coordinates (369−16=353; 370+7=377); they are the analysis window definition. See [the coordinate and citation audit](Sequence/RCL/Citation_Audit/README.md).

Metadata annotations are retained from the supplied inputs. Reconstruction and score agreement are not new experimental evidence for those annotations. Software versions and the mapped interval are written to the generated QC records.

## Expression-data verification scope

Selected processed expression summaries and reference outputs were supplied by
the author. Bulk consensus, aging comparison, integrated evidence and per-sex
sensitivity were reproduced from those processed inputs. The HPA sample-sex join was independently rerun and the selected tissue
summaries reproduced. Complete original GTEx bulk TPM-to-sex joins and
original single-cell/single-nucleus extractions were not independently rerun
in this environment. Source-data
portals, local reproduction requirements and precise validation limits are
recorded in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

## RCL primary-literature sources supplied by the author

The supplied bibliography identifies Ellisdon et al. (2014),
[doi:10.1186/1472-6807-14-14](https://doi.org/10.1186/1472-6807-14-14),
as the relevant primary Spn42Da structural study. The uploaded paper was directly inspected and its construct/shared-RCL
numbering mapped to A1Z6R1 as described above. Only the inspected Spn42Da source is retained here as a sequence-annotation reference. The additional human primary-source audit was waived, as documented below.

## Human reactive-site annotation scope

At the author's request, no additional primary-literature verification was
performed for the five human reactive-site annotations. Their residue identities
and alignment coordinates were checked against the supplied sequences. They
are described as supplied annotations, not newly or independently established
experimental cleavage sites. This does not change the SERPINE1-referenced
25-position identity/similarity calculations or the Spn42Da literature audit.
