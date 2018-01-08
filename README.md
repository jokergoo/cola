# cola
Unsupervised subgroup classification based on consensus clustering

## Features

1. Compare different methods of extracting top genes and different methods of clustering in parallel.
2. Hierarchical partitioning
3. Generate a detialed HTML report for the analysis

## Workflow:

```{r}
mat = adjust_matrix(mat)
res_list = run_all_consensus_partition_methods(mat, ...)
cola_report(res_list, output_dir = ...)

# or
res = hierarchical_partition(mat, ...)
cola_report(res, output_dir = ...)
```
