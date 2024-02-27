# cola: A General Framework for Consensus Partitioning <img src="https://user-images.githubusercontent.com/449218/54158555-03e3af80-444b-11e9-9773-070823101263.png" width=250 align="right" style="border:4px solid black;" />


[![R-CMD-check](https://github.com/jokergoo/cola/workflows/R-CMD-check/badge.svg)](https://github.com/jokergoo/cola/actions)
[ ![bioc](https://bioconductor.org/shields/downloads/devel/cola.svg) ](http://bioconductor.org/packages/stats/bioc/cola)
[ ![bioc](http://bioconductor.org//shields/lastcommit/devel/bioc/cola.svg) ](http://bioconductor.org/checkResults/devel/bioc-LATEST/cola/)



## Citation

Zuguang Gu, et al., cola: an R/Bioconductor package for consensus partitioning through a general framework, Nucleic Acids Research, 2021. https://doi.org/10.1093/nar/gkaa1146

Zuguang Gu, et al., Improve consensus partitioning via a hierarchical procedure. Briefings in bioinformatics 2022. https://doi.org/10.1093/bib/bbac048 



## Install

*cola* is available on [Bioconductor](http://bioconductor.org/packages/devel/bioc/html/cola.html), you can install it by:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cola")
```

The latest version can be installed directly from GitHub:

```r
library(devtools)
install_github("jokergoo/cola")
```

## Methods

The **cola** supports two types of consensus partitioning.

### Standard consensus partitioning

#### Features

1. It modularizes the consensus clustering processes that various methods can
   be easily integrated in different steps of the analysis.
2. It provides rich visualizations for intepreting the results.
3. It allows running multiple methods at the same time and provides
   functionalities to compare results in a straightforward way.
4. It provides a new method to extract features which are more efficient to
   separate subgroups.
5. It generates detailed HTML reports for the complete analysis.



#### Workflow

<img width="700" src="https://user-images.githubusercontent.com/449218/52628723-86af3400-2eb8-11e9-968d-b7f47a408818.png" />

The steps of consensus partitioning is:

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

#### Usage

Three lines of code to perfrom *cola* analysis:

```r
mat = adjust_matrix(mat) # optional
rl = run_all_consensus_partition_methods(
    mat, 
    top_value_method = c("SD", "MAD", ...),
    partition_method = c("hclust", "kmeans", ...),
    cores = ...)
cola_report(rl, output_dir = ...)
```

#### Plots

Following plots compare consensus heatmaps with k = 4 under all combinations of methods.

<img src="https://user-images.githubusercontent.com/449218/52631118-3a66f280-2ebe-11e9-8dea-0172d9beab91.png" />


### Hierarchical consensus partitioning


#### Features

1. It can detect subgroups which show major differences and also moderate differences.
2. It can detect subgroups with large sizes as well as with tiny sizes.
3. It generates detailed HTML reports for the complete analysis.


#### Hierarchical Consensus Partitioning


<img src="https://user-images.githubusercontent.com/449218/126491482-31a9496f-cc4d-4c4f-80b7-7b752d8d8d06.png" width="400" />



#### Usage

Three lines of code to perfrom hierarchical consensus partitioning analysis:

```r
mat = adjust_matrix(mat) # optional
rh = hierarchical_partition(mat, mc.cores = ...)
cola_report(rh, output_dir = ...)
```

#### Plots

Following figure shows the hierarchy of the subgroups.

<img src="https://user-images.githubusercontent.com/449218/100014572-d7b2c280-2dd6-11eb-9265-a84d324122f2.png" width="300" />

Following figure shows the signature genes.

<img src="https://user-images.githubusercontent.com/449218/100014657-f913ae80-2dd6-11eb-9bf7-53f733e9f8f0.png" width="600" />



## License

MIT @ Zuguang Gu

