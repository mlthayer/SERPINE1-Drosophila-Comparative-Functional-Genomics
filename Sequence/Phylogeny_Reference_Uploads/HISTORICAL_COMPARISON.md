# Historical phylogeny comparison

All eight original Newick files and the original `04_sessionInfo.txt` are preserved unchanged. All tip sets agree. Seven of eight topologies agree; only expanded human–fly homologous RCL ML differs. All four NJ trees match, including distances within floating-point precision. The expanded full-length ML tree also matches exactly; the other two matching ML topologies have small branch-distance differences.

| Tree | Tips | RF distance from original to current | Maximum patristic-distance difference |
| --- | --- | --- | --- |
| `04A_full_length_SERPINE1_fly_ML.nwk` | 8 | 0 | 7.339165e-05 |
| `04A_full_length_SERPINE1_fly_NJ.nwk` | 8 | 0 | 0 |
| `04B_revised_RCL_SERPINE1_fly_ML.nwk` | 8 | 0 | 0.000224556 |
| `04B_revised_RCL_SERPINE1_fly_NJ.nwk` | 8 | 0 | 0 |
| `04C_expanded_human_fly_ML.nwk` | 12 | 0 | 0 |
| `04C_expanded_human_fly_NJ.nwk` | 12 | 0 | 0 |
| `04D_expanded_human_fly_RCL_ML.nwk` | 12 | 8 | 4.474882 |
| `04D_expanded_human_fly_RCL_NJ.nwk` | 12 | 0 | 4.440892e-16 |

## Runtime and model diagnosis

The original run used R 4.5.1 on Windows, phangorn 2.12.1, ape 5.8-1, msa 1.40.0, and Biostrings 2.76.0. The current run uses R 4.5.3 on Linux, the same phangorn and ape versions, msa 1.44.1, and Biostrings 2.78.0. All four alignment sequences and gap positions match despite these runtime differences; runtime differences alone do not establish the cause of the different ML result.

The historical expanded RCL likelihood is −703.0974658653 versus the current −704.7610395782. A diagnostic fixed-topology optimization with identical JTT, four gamma categories, and optimized gamma/invariant-site parameters yields approximately −703.0981467390 for the historical topology and −704.7610399149 for the current topology. The historical topology therefore fits the same current data/model better. This supports a difference in the stochastic tree search outcome rather than an alignment inconsistency; it does not identify a specific software or platform cause.

The diagnostic fixes topology only for comparison and is not a replacement for the original stochastic-search pipeline. Primary model settings, tree-search settings, bootstrap count, original references, and generated results were not changed to force agreement. All four seven-candidate ranking columns agree. The expanded RCL disagreement remains material for the expanded-tree result and should be disclosed/resolved before selecting publication outputs.

Original ML Newick trees include internal labels, while original NJ trees do not. Node labels alone are not a complete saved 500-replicate bootstrap collection. For exact historical bootstrap-support reconciliation, provide the original saved bootstrap trees/model objects or node-support tables and their label definitions. No additional biological dataset is needed to execute step 04.
