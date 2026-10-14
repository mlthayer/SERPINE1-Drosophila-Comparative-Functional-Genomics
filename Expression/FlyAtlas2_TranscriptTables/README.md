# FlyAtlas2 transcript source verification

The six original text exports in this directory are transcript-level tables, identified by the transcript label/accession headers and paired FPKM/SD columns. They are preserved byte-for-byte from the uploads. They are not the gene-level exports consumed by step 05A; no gene totals were calculated from them.

All 817 records match `Expression/FlyAtlas2_Analysis/Processed/05_FlyAtlas2_all_transcript_expression.csv` by gene, source filename, stage/sex, tissue, transcript label, and transcript accession. All FPKM and SD values match, including missing-value status and numerical zeros. There are 627 adult observations and 190 larval observations across 19 transcripts. Both adult isoform analyses were previously verified using this matching CSV.

| Source file | Gene | Transcripts | Observations |
| --- | --- | --- | --- |
| `FBgn0028987T.txt` | Spn28F | 1 | 43 |
| `FBgn0028988T.txt` | Spn42Dd | 1 | 43 |
| `FBgn0033112T.txt` | Spn42Db | 2 | 86 |
| `FBgn0033113T.txt` | Spn42Dc | 2 | 86 |
| `FBgn0083141T.txt` | Spn28B | 2 | 86 |
| `FBgn0265137T.txt` | Spn42Da | 11 | 473 |
