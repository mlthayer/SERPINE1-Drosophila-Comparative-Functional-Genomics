# Chapter 3 SERPINE1 Functional Orthology Pipeline

This repository contains the R scripts used for Chapter 3 analyses comparing selected Drosophila serpins with human SERPINE1/PAI-1.

## Project purpose

The workflow evaluates candidate Drosophila serpins using multiple evidence categories:

1. DIOPT orthology
2. Full-length protein alignment
3. Reactive center loop alignment
4. Phylogenetic analysis
5. FlyAtlas2 tissue expression
6. GTEx and HPA SERPINE1 tissue expression
7. Fly Cell Atlas cell expression
8. Cross-species tissue-expression similarity
9. Aging tissue-expression comparison
10. Integrated evidence ranking
11. Sensitivity analysis

## Script order

Run the scripts in numerical order:

```text
01_DIOPT_Orthology_revised.R
02_Protein_Alignment_FullLength_revised.R
03_RCL_Alignment_revised.R
04_Phylogenetic_Analysis_revised.R
05_FlyAtlas2_TissueExpression_revised.R
06_GTEx_HPA_SERPINE1_TissueExpression_revised.R
07_FlyCellAtlas_CellExpression_revised.R
08_CrossSpecies_TissueSimilarity_revised.R
09_Aging_TissueExpression_Comparison_revised.R
10_IntegratedEvidenceRanking_revised.R
11_SensitivityAnalysis_revised.R

cat > LICENSE <<'EOF'

MIT License



Copyright (c) 2026 Michelle Thayer 



Permission is hereby granted, free of charge, to any person obtaining a copy

of this software and associated documentation files, to deal in the Software

without restriction, including without limitation the rights to use, copy, modify,

merge, publish, distribute, sublicense, and/or sell copies of the Software,

subject to the following conditions:



The above copyright notice and this permission notice shall be included in all

copies or substantial portions of the Software.



THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR

IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,

FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE

AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER

LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,

OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE

SOFTWARE.

EOF~
cat > CITATION.cff <<'EOF'

cff-version: 1.2.0

message: "If you use this workflow, please cite it using the metadata below."

title: "Chapter 3 SERPINE1 Functional Orthology Pipeline"

authors:

  - family-names: Thayer 

    given-names: Michelle

year: 2026

version: "1.0.0"

license: MIT

repository-code: "https://github.com/YOUR-USERNAME/Chapter3_SERPINE1_Functional_Orthology"

keywords:

  - SERPINE1

  - PAI-1

  - Drosophila

  - serpin

  - functional orthology

  - FlyAtlas2

  - Fly Cell Atlas

  - GTEx

  - phylogenetics

  - aging

