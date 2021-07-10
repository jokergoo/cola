# cola: A General Framework for Consensus Partitioning <img src="https://user-images.githubusercontent.com/449218/54158555-03e3af80-444b-11e9-9773-070823101263.png" width=250 align="right" style="border:4px solid black;" />


[![R-CMD-check](https://github.com/jokergoo/cola/workflows/R-CMD-check/badge.svg)](https://github.com/jokergoo/cola/actions)
[ ![bioc](https://bioconductor.org/shields/downloads/devel/cola.svg) ](http://bioconductor.org/packages/stats/bioc/cola)
[ ![bioc](http://bioconductor.org//shields/lastcommit/devel/bioc/cola.svg) ](http://bioconductor.org/checkResults/devel/bioc-LATEST/cola/)


## Features

1. It modularizes the consensus clustering processes that various methods can
   be easily integrated in different steps of the analysis.
2. It provides rich visualizations for intepreting the results.
3. It allows running multiple methods at the same time and provides
   functionalities to compare results in a straightforward way.
4. It provides a new method to extract features which are more efficient to
   separate subgroups.
5. It generates detailed HTML reports for the complete analysis.

## Citation

Zuguang Gu, et al., cola: an R/Bioconductor package for consensus partitioning through a general framework, Nucleic Acids Research, 2020. https://doi.org/10.1093/nar/gkaa1146

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

## Links

Examples for *cola* analysis can be found at https://jokergoo.github.io/cola_examples/ and https://jokergoo.github.io/cola_collection/.

Online documentation is at https://jokergoo.github.io/cola.

Supplementary for the *cola* manuscript is at https://github.com/jokergoo/cola_supplementary and the scripts are at https://github.com/jokergoo/cola_manuscript.

## Vignettes

There are the following vignettes:

1. [A Quick Start of Using cola Package](https://jokergoo.github.io/cola_vignettes/cola_quick.html)
2. [A Framework for Consensus Partitioning](https://jokergoo.github.io/cola_vignettes/cola.html)
3. [Automatic Functional Enrichment on Signature genes](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html)
4. [Predict Classes for New Samples](https://jokergoo.github.io/cola_vignettes/predict.html)
5. [Work with Big Datasets](https://jokergoo.github.io/cola_vignettes/work_with_big_datasets.html)
6. [Hierarchical consensus partitioning](https://jokergoo.github.io/cola_vignettes/hierarchical.html)

## Consensus Partitioning

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

### Usage

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

### Plots

Following plots compare consensus heatmaps with k = 4 under all combinations of methods.

<img src="https://user-images.githubusercontent.com/449218/52631118-3a66f280-2ebe-11e9-8dea-0172d9beab91.png" />


## License

MIT @ Zuguang Gu

#### Acknowledgement

The cola icon in logo is made by <a href="https://www.flaticon.com/authors/photo3idea-studio" title="photo3idea_studio">photo3idea_studio</a> from <a href="https://www.flaticon.com/" title="Flaticon">www.flaticon.com</a>.

