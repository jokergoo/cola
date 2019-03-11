<img src="https://user-images.githubusercontent.com/449218/54158555-03e3af80-444b-11e9-9773-070823101263.png" width=250 align="right" style="border:4px solid black;" />


[![Build Status](https://travis-ci.org/jokergoo/cola.svg)](https://travis-ci.org/jokergoo/cola) 
[ ![bioc](http://bioconductor.org/shields/downloads/devel/cola.svg) ](http://bioconductor.org/packages/stats/bioc/cola)
[ ![bioc](http://bioconductor.org//shields/lastcommit/devel/bioc/cola.svg) ](http://bioconductor.org/checkResults/devel/bioc-LATEST/cola/)

# cola

A General Framework for Consensus and Hierarchical Partitioning

## Features

1. It modularizes the consensus clustering processes that various methods can
   be easily integrated in different steps of the analysis.
2. It provides rich visualizations for intepreting the results.
3. It allows running multiple methods at the same time and provides
   functionalities to compare results in a straightforward way.
4. It provides a new method to extract features which are more efficient to
   separate subgroups.
5. It allows doing partitioning in a hierarchical way to detect subgroups
   with relatively smaller difference.
6. It generates detailed HTML reports for the complete analysis.

## Install

`cola` is available on [Bioconductor](http://bioconductor.org/packages/devel/bioc/html/cola.html), you can install it by:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cola", version = "3.9")
```

If you want to try the latest version, install it directly from GitHub:

```r
library(devtools)
install_github("jokergoo/cola")
```

## Examples

Several examples for cola analysis can be found at https://jokergoo.github.io/cola_examples/.

## Vignettes

1. [A Quick Start of Using cola Package](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/a_quick_start.html)
2. [A Framework for Consensus and Hierarchical Partitioning](http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/cola.html)

## Consensus Partition

<img width="700" src="https://user-images.githubusercontent.com/449218/52628723-86af3400-2eb8-11e9-968d-b7f47a408818.png" />

The steps of consensus partition is:

1. Clean the input matrix. The processing are: adjusting outliers, imputing missing
   values and removing rows with very small variance. This step is optional.
2. Extract subset of rows with highest scores. Here "scores" are calculated by
   a certain method. For gene expression analysis or methylation data
   analysis, $n$ rows with highest variance are used in most cases, where
   the "method", or let's call it **"the top-value method"** is the variance (by
   `var()` or `sd()`). Note the choice of "the top-value method" can be
   general. It can be e.g. MAD (median absolute deviation) or any user-defined
   method.
3. Scale the rows in the sub-matrix (e.g. gene expression) or not (e.g. methylation data).
   This step is optional.
4. Randomly sample a subset of rows from the sub-matrix with probability $p$ and
   perform partition on the columns of the matrix by a certain partition
   method, with trying different numbers of subgroups.
5. Repeat step 4 several times and collect all the partitions.
6. Perform consensus partitioning analysis and determine the best number of
   subgroups which gives the most stable subgrouping.
7. Apply statistical tests to find rows that show significant difference
   between the predicted subgroups. E.g. to extract subgroup specific genes.
8. If rows in the matrix can be associated to genes, downstream analysis such
   as function enrichment analysis can be performed.

### Usage

Three lines of code to perfrom cola analysis:

```r
mat = adjust_matrix(mat)
rl = run_all_consensus_partition_methods(
    mat, 
    top_value_method = c("sd", "MAD", ...),
    partition_method = c("hclust", "kmeans", ...),
    mc.cores = ...)
cola_report(rl, output_dir = ...)
```

### Plots

Following plots compare consensus heatmaps with k = 4 under all combinations of methods.

<img src="https://user-images.githubusercontent.com/449218/52631118-3a66f280-2ebe-11e9-8dea-0172d9beab91.png" />


## Hierarchical Partition

**cola** supports to look for subgroups in a hierarchical way, which is efficient to found secondary subgroups with small difference.

<img width="500" src="https://user-images.githubusercontent.com/449218/52628768-9a5a9a80-2eb8-11e9-9de1-974afd10bfc2.png" />

### Usage

```r
rh = hierarchical_partition(
    mat,
    top_value_method = "MAD", 
    partition_method = "kmeans",
    mc.cores = ...
)
cola_report(rh, output_dir = ...)
```

### Plots

The hierarchy of the subgroups:

<img src="https://user-images.githubusercontent.com/449218/48982041-8be0bc00-f0dd-11e8-80ab-2732daedad94.png">

## License

MIT @ Zuguang Gu

