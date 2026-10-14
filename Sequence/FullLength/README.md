# Full-length sequence inputs and reproduction

The eight files in `Reference_Uploads/` are preserved unchanged from the uploads. The original UniProt-format input is also placed at `Uniport FASTA FlySerpins_Human PAI-1.txt`, the path used by step 02. Generated outputs remain separate in `Processed/`, `Tables/`, and `QC/`.

Step 02 completed using the existing ClustalOmega alignment, identity, BLOSUM62 similarity, gap handling, and QC definitions. The seven validated sequences and seven aligned sequences match the uploaded FASTA records exactly, including headers and order. The six candidate statistics, seven sequence-QC records, and seven alignment-QC records match their uploaded CSV tables (numeric comparison tolerance 1e-12).

Both older Clustal text alignments retain their original sequence sets and alignments. All records shared with the input FASTA match after removing alignment gaps. The fly-only file has accession `Q8MSJ3` in addition to the six current candidates; that additional record was not introduced into the current analysis. The human–fly text alignment contains SERPINE1 and the six current fly candidates. Older text alignments were retained as references rather than substituted for the newly reproduced alignment.

This resolves the full-length analysis input and reproduction blocker. It does not provide the expanded `selected_human_fly_Serpins_Fasta.txt` required by the RCL and phylogenetic scripts, including their additional human SERPIN/isoform records.
