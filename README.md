# cola
A General Framework for Consensus and Hierarchical Partition

## Features

1. It modulizes the consensus clustering processes that variaty of methods can be easily integrated.
2. It provides rich visualizations for intepreting the results.
3. It allows running multiple methods at the same time and provides functionalities to compare results in a straightforward manner.
4. It provides a new method to extract features which are more efficient to separate subgroups.
5. It allows doing partition in a hierarchical way to detect subgroups with relatively smaller difference.
6. It generates detailed reports for the whole analysis.

## Consensus Partition

<img width="600" alt="consensus_partition_workflow" src="https://user-images.githubusercontent.com/449218/41279038-b1d73680-6e2b-11e8-94eb-d0d53eae5d60.png">


## Usage

```{r}
mat = adjust_matrix(mat)
rl = run_all_consensus_partition_methods(mat, 
    top_value_method = c("sd", "MAD", ...),
    partition_method = c("hclust", "kmeans", ...),
    mc.cores = ...)
cola_report(rl, output_dir = ...)
```

## Plots

Following plots compares consensus heatmap with k = 4 under all combination of methods.

<img src="https://user-images.githubusercontent.com/449218/41279146-f949b9f2-6e2b-11e8-902e-0715fe7b16fa.png">

## Hierarchical Partition

**cola** supports to look for subgroups in a hierarchical way, which is efficient to found secondary subgroups with small difference.

<img width="600", src="https://user-images.githubusercontent.com/449218/41279180-1060bb68-6e2c-11e8-8bb3-f6a197294f4d.png">

## Usage

```{r}
rh = hierarchical_partition(data, top_n = c(1000, 2000, 3000, 4000),
    top_value_method = "MAD", partition_method = "kmeans")
cola_report(rh, output_dir = ...)
```

## Plots

The hierarchy of the subgroups:

<img src="https://user-images.githubusercontent.com/449218/41279212-2362cb7a-6e2c-11e8-8afc-fea14013a95f.png">

## License

MIT

