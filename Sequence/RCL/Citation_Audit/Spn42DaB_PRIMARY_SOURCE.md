# Spn42Da A/B reactive-site primary-source verification

Ellisdon AM, Zhang Q, Henstridge MA, Johnson TK, Warr CG, Law RHP, Whisstock JC.
High resolution structure of cleaved Serpin 42 Da from Drosophila melanogaster.
BMC Structural Biology. 2014;14:14.
https://doi.org/10.1186/1472-6807-14-14

The author's uploaded eight-page article was inspected directly, including
rendered pages to verify text recovered by OCR. The PDF's printed page numbers
are used below; the plain PDF text layer was not readable.

| Evidence | Article location | Interpretation |
|---|---|---|
| A and B share the RCL; B has an alternative N-terminal secretion-signal exon | p. 2, Figure 1 and text | Shared RCL does not imply identical full-length numbering |
| P4–P1 is R-R-K-R and P1 is R342 | p. 3, cleaved-structure discussion | Primary paper explicitly identifies the RCL P1 reference |
| Crystal form 1 density ends at P1 Arg; P1prime is Ala | pp. 3–4 | Structural evidence consistent with cleavage at the R–A site |
| Crystal form 2 ends at P2 Lys, described as cleavage after K341 | pp. 4 and 6 | Not both crystal forms demonstrate the same cleavage position |
| Shared A/B RCL motif is shown | p. 6, Figure 5A | Sequence-context correspondence with the supplied B RCL |
| Recombinant untagged protein corresponds to residues 1–392 | p. 7, Methods | Paper's shorter construct coordinate context |

## Numbering used in this workflow

The supplied A1Z6R1 B sequence is 424 residues. Its R374–A375 and local
R-R-K-R-A motif correspond to the paper's R342–A343 in the 392-residue
construct/A–B shared-RCL coordinate context: **424−392=32**, **374−32=342**,
**375−32=343**. Full-length A1Z6R1 positions 342–343 are L–K, so literature
numbers must not be applied directly to that full sequence.

Use the phrase **paper construct/A–B shared-RCL numbering versus full-length
A1Z6R1 numbering**. The paper's Figure 1 supports alternative N-terminal
architecture, but this analysis does not establish that a signal peptide
is experimentally cleaved at exactly residue 32. Nor does the paper text supply
a complete accession-tagged recombinant sequence for a full-sequence identity
comparison. Do not describe the +32 conversion as a measured cleavage event.

The R–A site is a literature-supported reference annotation for alignment QC.
The source also cites earlier functional work (its references 10 and 12) for
cleavage after P1; these cited papers have not been independently inspected here.
Alignment-column agreement remains sequence QC, not new biochemical evidence.
No predicted sites for other fly candidates or numerical RCL scores changed.

The uploaded PDF contains local creator/author metadata; the original is kept
in session attachments, not copied into the public repository. Its SHA256 and
page-specific evidence are recorded in Spn42DaB_numbering_audit.json.
