<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{A Framework for Consensus and Hierarchical Partition}
-->

cola: A General Framework for Consensus and Hierarchical Partition
=============================================================

**Author**: Zuguang Gu ( z.gu@dkfz.de )

**Date**: 2018-06-11

**Package version**: 0.99.0

-------------------------------------------------------------



## Introduction

Subgroup classification is a basic task in high throughput data analysis,
especially for gene expression and methylation data analysis. Mostly
unsupervised clustering methods are applied and it predicts new subgroups or
tests the consistency with known clinical annotations.

To get a stable classification of subgroups, [consensus
clustering](https://en.wikipedia.org/wiki/Consensus_clustering) or consensus
partition is always performed. It repeatedly clusters samples with a randomly
sampled subset of data and checks the robustness of the clustering, finally
gives a consensus classification of all samples.

Here we present the **cola** package which provides a general framework for
consensus partition. It has following advantages:

1. It modulizes the consensus clustering processes that variaty of methods can
   be easily integrated.
2. It provides rich visualizations for intepreting the results.
3. It allows running multiple methods at the same time and provides
   functionalities to compare results in a straightforward manner.
4. It provides a new method to extract features which are more efficient to
   separate subgroups.
5. It allows doing partition in a hierarchical way to detect subgroups
   with relatively smaller difference.
6. It generates detailed reports for the whole analysis.

Following flowchart lists the general steps of consensus partition implemented by **cola**:

<img src="consensus_partition_workflow.png", width="600"/>

The steps are:

1. Clean the input matrix, such as adjusting outliers and imputing missing
   values.
2. Extract subset of rows with highest scores. Here "scores" are calculated by a
   certain method. For gene expression analysis or methylation data analysis,
   $n$ rows with highest variance are always used. Here the "method", or let's
   call it "the top-value method" is the variance. However the choice of "the
   top-value method" can be general. It can be e.g. MAD (median absolute
   deviation) or any user-defined method.
3. Rows in the matrix with top rows are scaled (e.g. gene expression) or not
   scaled (e.g. methylation data).
4. Random sample a subset of rows from the matrix with probability $p$ and perform partition on
   the columns of the matrix by a certain partition method, with trying
   different numbers of subgroups.
5. Repeat step 4 several times and collect all the partitions.
6. Perfrom consensus partition analysis and determine the best number of
   subgroups which gives the most stable subgrouping.
7. Apply statistical tests to find rows that show significant difference between
   the predicted subgroups. E.g. to extract subgroup specific genes.
8. If rows in the matrix can be associated to genes, downstream analysis such
   as function enrichment analysis can be performed.

All the steps will be discussed in detail in following sections.

## Clean the matrix

In following part of this vignette, we call the columns of the matrix as
"samples".

Before doing consensus partition, a simple but important step is to clean the
input matrix:


```r
data = adjust_matrix(data)
```

`adjust_matrix()` does following preprocessing:

1. Rows where more than 25% of the samples having `NA` values are removed;
2. Use `impute::impute.knn()` to impute missing data if there is any;
3. For each row in the matrix, it uses `adjust_outlier()` (also provided by
   **cola** package) to adjust outliers. Values larger than the 95^th
   percentle or less than the 5^th percentle are replaced by corresponding
   percentiles.
4. Rows with zero variance are removed.
5. Rows with variance less than 5^th percentile of all row variance (after
   removing rows with zero variance) are removed.

Some of the above steps are optional. For example, methylation matrix does not
need to be adjusted for outliers because all the methylation values are
already in a fixed data scale (0 ~ 1).

## Basic usage

`consensus_partition()` performs consensus partition for a single top-value
method and a single partition method. The major arguments for
`consensus_partition()` are:


```r
res = consensus_partition(data,
    top_value_method = "MAD",
    top_n = c(1000, 2000, 3000, 4000, 5000),
    partition_method = "kmeans",
    max_k = 6,
    p_sampling = 0.8,
    partition_repeat = 50,
    anno = NULL)
```

- `data`: a data matrix where subgroups are found by columns.
- `top_value_method`: name of the method to assign values to rows in the
  matrix. Later these values are used to order and extract rows with top
  values.
- `top_n`: number of rows with top values used for partition. Normally we try
  a list of different ones.
- `partition_method`: name of the method for partition.
- `max_k`: max number of subgroups to try. It will try from 2 to `max_k`.
- `p_sampling`: proportion of the `top_n` rows to sample. The sub-matrix with
  `p_sample * top_n` rows is used for partition.
- `partition_repeats`: times of random sampling and partitions to perform.
- `anno`: a vector or a data frame which contains known annotations of
  samples. If it is provided, it will drawn along side with the predicted
  subgroups in the plot generated by downstream functions and it can also be
  tested for the correlation to predicted subgroups.

Other arguments can be found in the on-line documentation of
`consensus_partition()`.

In most of the cases, we want to try different top-value methods and different
partition methods to see which combination of methods gives the best prediction.
The helper function `run_all_consensus_partition_methods()` is a convinient
way for doing this:


```r
rl = run_all_consensus_partition_methods(data, 
	top_value_method = c("sd", "MAD", ...),
	partition_method = c("hclust", "kmeans", ...),
	mc.cores = ...)
```

There are functions in **cola** package that can visualize and compare the
results for all combinations of methods simutanuously.

## Top-value methods

Top-value methods are used to assign scores to rows in the matrix, later the
scores are ordered and only the top $n$ rows with the highest scores are used for
consensus partition. The default top-value methods provided in the package
are:


```r
all_top_value_methods()
```

```
## [1] "sd"  "cv"  "MAD" "AAC"
```

These top methods are:

- `sd`: Standard deviation.
- `cv`: [Coefficient of
  variance](https://en.wikipedia.org/wiki/Coefficient_of_variation), defined
  as `sd/(mean + s0)` where `s0` is a penalty term which is the 10^th
  percentile of all row means to avoid small values dividing small values
  giving large values.
- `MAD`: [Median absolute
  deviation](https://en.wikipedia.org/wiki/Median_absolute_deviation).
- `AAC`: A new method proposed in **cola** package and it will be explained later
  in this section.

These methods can be used in consensus partition by providing the name.

You can register a new top-value method by `register_top_value_method()`. The
value should be functions. For each function, it should only have one
argument which is the matrix for analysis and it must return a vector with
scores for rows. In following example, the "max" method uses the row maximum
as the row score and we also add the "QCD" ([quartile coefficient of
dispersion](https://en.wikipedia.org/wiki/Quartile_coefficient_of_dispersion))
method as a second method here.


```r
register_top_value_method(
	max = function(mat) rowMaxs(mat),
	QCD = function(mat) {
		qa = matrixStats::rowQuantile(mat, probs = c(0.25, 0.75))
		(qa[, 2] - qa[, 1])/(qa[, 2] + qa[, 1])
	})
all_top_value_methods()
```

```
## [1] "sd"  "cv"  "MAD" "AAC" "max" "QCD"
```

By default the consensus partition functions use all registered top-value
methods, but still you can explicitly specify a subset of top-value methods.
To remove registered top-value methods, simply use `remove_top_value_method()`
by providing a vector of names.


```r
remove_top_value_method(c("max", "QCD"))
```

### The AAC method

Choosing the top rows in the matrix is important for the subgroup
classification. In most of the cases, we extract the most variable rows which
is defined by row variance. However, sometimes it won't give you meaningful
rows which are efficient for subgroup classification. When random noise in the
data increases, e.g. for single cell RNASeq data, the most variable genes are
too weak to detect any stable subgroups.

If we think reversely, assuming there exist stable subgroups in the data,
there must be groups of rows showing similar pattern to support the
subgrouping, in other words, rows in the same groups should have high
correlations to each other. Thus, if we can get rows that more correlate to
others, they are more strong to form a stable subgroup for the samples.
According to this thought, we designed the AAC method.

For row $i$ in a matrix, $X$ is a vector of the absolute correlation to all other rows,
the AAC (area above the correlation CDF curve) for row $i$ is defined as:

$$AAC_i = 1 - \int_0^1F(x)$$

where $F(x)$ is the empirical CDF (cumulative distribution function) of $X$.

In following plot, the line is the CDF curve. AAC is the area above the CDF
curve. It can be imagined that when row $i$ correlates with more other rows,
the CDF curve shifts more to the right, thus with higher AAC scores.

<img src="figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" />

Next we perform a simulation test to show the attributes of AAC. A matrix with
160 rows with random values are generated as follows:

1. 100 rows with mean of 0, covariance of 0 with 1 on the diagnal;
2. 10 rows with mean of 0, covariance of 0.8 with 1 on the diagnal;
3. 50 rows wiht mean of 0, covaraince of 0.5 with 1 on the diagnal.

In the left figure in following, they are ECDF curves of the correlation when
calculating AAC scores. The middle figure is the AAC for all 160 rows and the
right figure is the standard deviation for the 160 rows.

<img src="figure/unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" style="display: block; margin: auto;" />

All the 160 rows have theretical variance of 1 that they cannot be
distinguished very well by using variance. As a contrast, the rows with
covariance have higher AAC values, even higher when the number of correlated
rows increases (although the correlation value itself is small). This shows
AAC method can assign higher values for rows which have wide-range correlation
patterns.

## Partition methods

Partition methods are used to separate samples into $k$ subgroups where $k$ is
a known parameter for the partition. The default partition methods are:


```r
all_partition_methods()
```

```
## [1] "hclust"  "kmeans"  "skmeans" "pam"     "mclust"
```

These partition methods are:

- `hclust`: hierarchical clustering + cutree.
- `kmeans`: k-means clustering.
- `skmeans`: spherical k-means clustering, from **skmeans** package.
- `pam`: partitioning around medoids, from **cluster** package.
- `mclust`: model-based clustering, from **mclust** package. The clustering is based on
  the first three principle dimension from the original matrix.

Similarly, you can register a new partition method by
`register_partition_method()`. The value is should be functions with two
arguments which are the input matrix and number of partitions. There can be a
third argument for the function which is  `...` used for passing more
arguments from the main partition functions. The function should only return a
vector of subgroup/class labels. **Please note the partition is applied on
columns of the matrix and the number of unique levels of subgroup levels which
are predicted by the partition method are not necessarily to be the same as
$k$.**

Following example registers a partition method which randomly assign subgroup
labels to samples:


```r
register_partition_method(random = function(mat, k) {
	sample(letters[1:k], ncol(mat), replace = TRUE)
})
```

Here the subgroup labels can be in any types (numbers, characters). They only
need to be different for different classes. These labels will be re-encoded
with numeric indices internally.

To remove a partition method, use `remove_partition_method()`:


```r
remove_partition_method("random")
```

### The skmeans method

The skmeans method ([the spherical k-means
clustering](https://www.jstatsoft.org/article/view/v050i10)) is powerful to
detect subgroups where samples in a same subgroup show strong correlations.
skmeans clustering uses cosine similarity and projects data points onto a unit hyper-sphere.
As we have tested for many datasets, skmeans is very efficient to detect stable subgroups.

<img src="skmeans.png" width="400" />

## Consensus clustering

For a given number of top rows $n_i$, the corresponding matrix with top rows
denoted as $M_i$, a subset of rows with probability of $p$ are randomly
sampled from $M_i$ and a certain partition method is applied on it, generating
a partition $P_a$. In most of cases, we have no prior knowledge of which $n_i$
gives better results, thus, **cola** allows to try multiple $n_i$ and put
partitions from all $n_i$ together to find a consensus subgrouping. Let's
assume top rows are tried for $n_1$, $n_2$, ..., $n_m$ and the randomly
sampling is performed for $N_s$ times, then, for a given number of subgroups
for trying, the total number of partitions is $N_P = m*N_s$.

### Consensus matrix

The consensus matrix measures how consistently two samples are in a same
subgroup and it can be used to visualize or analysis the stability of the
subgrouping. The value $c_{ij}$ in the consensus matrix is the probability of
sample $i$ and sample $j$ in a same subgroup in all $N_P$ partitions. It is
calculated as:

$$c_{ij} = \sum_a^{N_p}I(s_{ia}, s_{ja})/N_P$$

where $s_{ia}$ is the subgroup label for sample $i$ in partition $a$ and $I()$
is the identity function there $I(x = y) = 1$ and $I(x \neq y) = 0$.

Following two heatmaps visualize two consensus matrices. The left one shows
less stability of subgrouping than the right one, while for the right one,
there are three very obvious blocks in the diagnal that in each block, the
corresponding samples are likely all in a same subgroup.

<img src="figure/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" style="display: block; margin: auto;" />

### Consensus subgroup labels

As long as we have a list of $N_P$ partitions for a given subgroup number $k$,
we need to find a consensus partition based on all $N_P$ partitions. Simply
speaking, the consensus subgroup label for sample $i$ should be the one which
have the maximal occurency among all partitions.

Internally, **cola** package uses the **clue** package to construct the
"partition ensemble" and predict the consensus subgroups. The "SE" method from
`clue::cl_consensus()` function (please check the online documentation of this
function) are used to calculate the consensus subgroup labels.

### Adjust subgroup labels

The subgroup labels are assigned with numeric indices (1, 2, ...), however, in
each partition, the assignment of the labels can be random, e.g. one same
subgroup can be assigned with 1 in one partition, while in the other
partition, it can be 2, but they are all identical for the sense of
subgrouping. E.g. following partitions are all identical:

```
1 1 1 1 1 1 1 2 2 2 2 2 2
2 2 2 2 2 2 2 1 1 1 1 1 1
a a a a a a a b b b b b b
```

The subgroups are indentical if switching the group label. This group label
adjustment is called the linear sum assignment problem, which is solved by the
`solve_LSAP()` function in **clue** package. The aim is to find a mapping `m()`
between two sets of labels to minimize $\sum I(s_{1i}, m(s_{2i}))$ where $s_1$ is
the first label set and $s_2$ is the second label set.

In following example, if the mapping is `1 -> 2, 2 -> 1`, the second partition
in following

```
1 1 1 1 1 1 1 2 2 2 2 2
2 2 2 2 2 1 1 1 1 1 1 1
```

is adjusted to

```
1 1 1 1 1 1 1 2 2 2 2 2
1 1 1 1 1 2 2 2 2 2 2 2   # switch 1 <-> 2
```

For the subgroup predicted by `clue::cl_consensus()`, the labels are
additionally adjusted by the average distance in each subgroup (calculated
from the original matrix), which means, the subgroup with label 1 always has
the smallest intra-distance.

This subgroup label adjustment is frequently used in **cola** to help the
visualization as well as downstream analysis.

### Membership matrix 

The $N_P$ partitions are stored as a membership matrix where rows are
partitions (grouped by the number of top rows) and all the subgroup labels are
adjusted according to the consensus subgroups. Following heatmap is a
visualization of all partitions and correspondance to the consensus partition.
The "p\*" annotation on top of the heatmap is the probability of being in
subgroup $i$ across all partitions.

<img src="figure/unnamed-chunk-14-1.png" title="plot of chunk unnamed-chunk-14" alt="plot of chunk unnamed-chunk-14" style="display: block; margin: auto;" />

## Determine the best number of subgroups

Consensus partition is applied with a specified number of subgroups (we termed
it as $k$). Normally, a list of $k$ are tried to find the best $k$. **cola**
provides following metrics to help to determine the best number of subgroups:

The `get_stat()` function returns statistics for all metrics mentioned below
and `select_partition_number()` plots the statistics with the number of
subgroups.

### Cophenetic correlation

It measures if hierarchical clustering is applied on the consensus matrix, how
good it correlates to the consensus matrix itself
(https://en.wikipedia.org/wiki/Cophenetic_correlation). With higher the value,
better the $k$.

### Silhouette score

[The silhouette
scores](https://en.wikipedia.org/wiki/Silhouette_%28clustering%29) measures
how close one sample is in its own subgroup compared to the closest
neighbouring subgroup. For sample $i$, the mean distance to all subgroups are
calculated, denoted as $d_1$, $d_2$, ..., $d_k$. The distance to the subgroup
where sample $i$ stays is denoted as $d_a$ and the sihouette score is defined
as:

$$silhouette_i = 1 - d_a/d_b$$

where $d_b$ is the minimal distance excluding $d_a$:

$$d_b = min_{j \neq a}^k d_j$$

Following plot illustrates how silhouette score is calculated.

<img src="figure/unnamed-chunk-15-1.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" style="display: block; margin: auto;" />

The mean silhouette score from all samples is used to choose the best $k$ where
higher the mean silhouette score, better the $k$.

### PAC score

[The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering) measures the proportion of
the ambiguous subgrouping. If the subgrouping is stable, in $N_P$ partitions,
sample $i$ and sample $j$, in most of the cases, are either always in a same
subgroup, or always in different subgroups, which results in that, in the
consensus matrix, the values are, in most of the cases, close to 1 or 0. Then
in the CDF of the consensus matrix, the curve will be very flatterned between
$x_1$ and $x_2$ where $x_1$ is very close to 0 and $x_2$ is very close to 1
because there are very few values between $x_1$ and $x_2$. And then $F(x_2) -
F(x_1)$ can be used to measure how much the ambiguous subgrouping is.

<img src="figure/unnamed-chunk-16-1.png" title="plot of chunk unnamed-chunk-16" alt="plot of chunk unnamed-chunk-16" style="display: block; margin: auto;" />

In the original implementation of PAC, $x_1$ and $x_2$ are fixed (e.g. 0.1 and
0.9). But in many cases, the choice of $x_1$ and $x_2$ is quite sensitive,
thus, in **cola** package, the PAC implementation is a slightly changed to its
orignal implementation. For following form of PAC:

$$PAC_a = F(x_2) - F(x_1)$$

$x_1$ takes a list of values from $[0.1, 0.3]$ and $x_2$ takes a list of values
in $[0.7, 0.9]$. The final PAC score is the mean of $PAC_a$ by removing the
top 10^th and bottom 10^th values. This method improves the robustness of
choosing $x_1$ and $x_2$.

Smaller the PAC score, better the $k$.

### Concordance

The concordance of partitions to the consensus partition is calculated as,
for each partition $a$, the probability that it fits the consensus partition:

$$c_{a} = \sum_i^{N_s}I(s_{ia} = sc_i)/N_s$$

where $N_s$ is the number of samples and $sc$ is the consensus subgroup label.

The final concordance score is the mean value of $c_a$. Higher the concordance score, better the $k$.

### Area increased 

It is the increased area under CDF compared to the previous $k$.

$$A_k = \int F_k(x) - \int F_{k-1}(x)$$

and when $k = 2$ or for the minimal $k$:

$$A_k = \int F_k(x)$$

In follow example, there are five consensus heatmaps corresponding to $k = 2..6$:

<img src="figure/unnamed-chunk-17-1.png" title="plot of chunk unnamed-chunk-17" alt="plot of chunk unnamed-chunk-17" style="display: block; margin: auto;" />

The corresponding CDF curves and the area increased are:

<img src="figure/unnamed-chunk-18-1.png" title="plot of chunk unnamed-chunk-18" alt="plot of chunk unnamed-chunk-18" style="display: block; margin: auto;" />

The $k$ before the elbow is taken as the best $k$ (in above example it is 3).
Basically when $k$ reaches a stable subgrouping, increasing $k$ won't change
the consens matrix too much, which results in less change of the difference of
area under the CDF curve.

### Rand index 

In some cases, when number of subgroups changes from $k-1$ to $k$, all the
statistics imply $k$ is a better choisce than $k-1$. However, when observing
the consensus heatmap, basically it is because a very small set of samples are
separated to form a new subgroup and in many cases, they are outlier samples.
In this case, it is better to still keep $k-1$ subgroups. Or in other words,
the subgrouping with $k$ is similar as $k-1$ and it is not worth to increase
$k$ from $k-1$. In **cola** package, there are two metrics Rand index and
Jaccard index to measure the similarity of two partitions for $k-1$ and
$k$. The two metrics are calculated by `clue::cl_agreement(..., method =
"Rand")` or `clue::cl_agreement(..., method = "Jaccard")`.

For all pairs of samples, denote following symbols
(https://en.wikipedia.org/wiki/Rand_index#Definition):

- $a$: the number of pairs of samples that are in the same subgroup in $k$ and
  in the same subgroup in $k-1$
- $b$: the number of pairs of samples that are in the different subgroup in
  $k$ and in the different subgroup in $k-1$
- $c$: the number of pairs of samples that are in the same subgroup in $k$ and
  in the different subgroup in $k-1$
- $d$: the number of pairs of samples that are in the different subgroup in
  $k$ and in the same subgroup in $k-1$

the Rand index which is the percent of pairs of samples that are both in a
same cluster or both are not in a same cluster in the partition of k and k-1.

$$Rand = \frac{a+b}{a+b+c+d}$$

If Rand index is too high, it means the subgrouping is very similar and not
sufficient to increase from $k-1$ to $k$.

### Jaccard index

The Jaccard index is the ratio of pairs of samples that are both in a same subgroup
in the partition of $k$ and $k-1$ and the pairs of samples are both in a same
subgroup in the partition of $k$ or $k-1$.

$$Jaccard = \frac{a}{a+c+d}$$

### Rule

**cola** provides a `guess_best_k()` function which determines the best $k$.
It is based on following rules:

- All $k$ with Rand index larger than 0.95 are removed because the increase of
  the partition number does not provides enough extra information.
- For $k$ with cophenetic correlation coefficient larger than 0.99 or PAC
  score less than 0.1 or concordance score larger than 0.95, the maximal $k$
  is taken as the “best k”.
- If it does not fit the second rule. The $k$ with maximal vote of highest
  cophenetic correlation coefficient, lowest PAC score, highest mean
  silhouette and highest concordance is taken as the “best k”.

## Find signatures

As long as there are stable subgroups, we can look for rows which show
distinct difference in one subgroup compared to others. They can be called
signature genes or signature probes if the corresponding dataset is gene
expression data or methylation data.

By default, samples with silhouette scores less than 0.5 are removed. **cola** provides
following methods:

- `ttest`: First it looks for the subgroup with highest mean value, compare to
  each of the other subgroups with t-test and take the maximum p-value. Second
  it looks for the subgroup with lowest mean value, compare to each of the
  other subgroups again with t-test and take the maximum p-values. Later for
  these two list of p-values take the minimal p-value as the final p-value.
- `samr` and `pamr`: use SAM/PAM method to find significantly different rows between
  subgroups.
- `Ftest` use F-test to find significantly different rows between subgroups.

Users can also provide their own method by providing a function with the matrix and subgroup labels
as input, 

## Compare multiple methods

`consensus_partition()` is the core function for doing consensus partition.
But it can only perform analysis with a single top-value method and a single
partition method. In most of the cases, we have no idea of which combination
of top-value method and partition method gives better results. Here
`run_all_consensus_partition_methods()` can perform anlaysis with multiple
methods:


```r
rl = run_all_consensus_partition_methods(...)
```

By default it runs analysis for all combinations of top-value methods in
`all_top_value_methods()` and partition methods in `all_partition_methods()`.

**cola** package provides functions to collect plots from all combinations of
methods for straightforward comparisons.


```r
collect_plots(rl, fun = consensus_heatmap, k = ...)
collect_plots(rl, fun = membership_heatmap, k = ...)
collect_plots(rl, fun = get_signatures, k = ...)
```

And `collect_classes()` compares consensus partition from all methos:


```r
collect_classes(rl, k = ...)
```

When runnig `run_all_consensus_partition_methods()`, there are consensus
partition for each combination of the methods. To compare partitions from
different method, there is a ultimate consensus partition which is calcualted from all
methods by weighting the mean silhouette scores. And for each combination of
methods, the subgroup labels are adjusted by the ultimate consensus subgroup.

## Visualizations

**cola** package provides rich visualizations for the results generated by a
single method or multiple methods.

### On the ConsensusPartition object

The object which is generated with a single top-value method and a single
partition method belongs to the class `ConsensusPartition`. There are several
visualization functions that can be applied to it.
 `select_partition_number()` makes several plots to show
different statistics along with different $k$, which helps to determine the
"best k".


```r
rl = readRDS("~/cola_test/TCGA_subgroup.rds")
res = rl["MAD:kmeans"] # the ConsensusPartition object
select_partition_number(res)
```

<img src="figure/unnamed-chunk-22-1.png" title="plot of chunk unnamed-chunk-22" alt="plot of chunk unnamed-chunk-22" style="display: block; margin: auto;" />

The heatmap for the consensus matrix with a certain $k$:


```r
consensus_heatmap(res, k = 4, show_row_names = FALSE)
```

<img src="figure/unnamed-chunk-23-1.png" title="plot of chunk unnamed-chunk-23" alt="plot of chunk unnamed-chunk-23" style="display: block; margin: auto;" />

The heatmap for the membership matrix with a certain $k$:


```r
membership_heatmap(res, k =4)
```

<img src="figure/unnamed-chunk-24-1.png" title="plot of chunk unnamed-chunk-24" alt="plot of chunk unnamed-chunk-24" style="display: block; margin: auto;" />

The PCA plot under the subgroups with a certain $k$:


```r
dimension_reduction(res, k = 4)
```

<img src="figure/unnamed-chunk-25-1.png" title="plot of chunk unnamed-chunk-25" alt="plot of chunk unnamed-chunk-25" style="display: block; margin: auto;" />

The heatmap for the signature rows with a certain $k$. The heatmap is split
into two parts by columns. The left heatmap contains samples with silhouette
scores larger than 0.5 and the right heatmap contains samples with silhouette
scores less than 0.5.


```r
get_signatures(res, k = 4, show_column_names = FALSE)
```

<img src="figure/unnamed-chunk-26-1.png" title="plot of chunk unnamed-chunk-26" alt="plot of chunk unnamed-chunk-26" style="display: block; margin: auto;" />

`collect_classes()` which is applied on the `ConsensusPartition` object
visualizes how subgroups are corresponded with increasing $k$:


```r
collect_classes(res, show_row_names = FALSE)
```

<img src="figure/unnamed-chunk-27-1.png" title="plot of chunk unnamed-chunk-27" alt="plot of chunk unnamed-chunk-27" style="display: block; margin: auto;" />

`collect_plots()` which is applied on the `ConsensusPartition` object puts
all the plots from all $k$ into one single page.


```r
collect_plots(res)
```

<img src="figure/unnamed-chunk-28-1.png" title="plot of chunk unnamed-chunk-28" alt="plot of chunk unnamed-chunk-28" style="display: block; margin: auto;" />

### On the ConsensusPartitionList object

`run_all_consensus_partition_methods()` returns a `ConsensusPartitionList` object.
There are two main functions which can visualize results from all combinations
of methods and compare directly.


```r
collect_plots(rl, fun = consensus_heatmap, k = 4)
```

<img src="figure/unnamed-chunk-29-1.png" title="plot of chunk unnamed-chunk-29" alt="plot of chunk unnamed-chunk-29" style="display: block; margin: auto;" />

`fun` can also be `membership_heatmap` or `get_signatures` that membership heatmap
and signature heatmap for each method will be plotted.

`collect_classes()` which is applied on the `ConsensusPartitionList` object plots 
the partition for each combination of methods and the lightness correspond to the 
silhouette scores for samples in each method. Rows are clustered by the dissimilarity
measurement from `clue::cl_dissimilarity(..., method = "comembership")`. On top the
consensus subgroup is inferred from all methods by taking the mean silhouette scores
as weight.


```r
collect_classes(rl, k = 4)
```

<img src="figure/unnamed-chunk-30-1.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" style="display: block; margin: auto;" />

## A Real example

In this example, we use [the TCGA GBM microarray data set](https://www.ncbi.nlm.nih.gov/pubmed/20129251)
where four subtypes is predicted.
The two files (`unifiedScaled.txt` and `TCGA_unified_CORE_ClaNC840.txt`) for use is from 
[here](https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/).

Following code is used to perform analysis of consensus partition.


```r
data = read.table("~/cola_test/unifiedScaled.txt", header = TRUE, 
    row.names = 1, check.names = FALSE)
data = as.matrix(data)

subtype = read.table("~/cola_test/TCGA_unified_CORE_ClaNC840.txt", 
    sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])

data = data[, names(subtype)]

data = adjust_matrix(data)
rl = run_all_consensus_partition_methods(data, top_n = c(1000, 2000, 3000, 4000), 
    max_k = 6, mc.cores = ncore, anno = data.frame(subtype = subtype), 
    anno_col = list(subtype = structure(seq_len(4), names = unique(subtype))))
```

Simply type `rl` gives a summary of the analysis:


```r
rl
```

```
## A 'ConsensusPartitionList' object with 20 methods.
##   On a matrix with 11268 rows and 173 columns.
##   Top rows are extracted by 'sd, cv, MAD, AAC' methods.
##   Subgroups are detected by 'hclust, kmeans, skmeans, pam, mclust' method.
##   Number of partitions are tried for k = 2, 3, 4, 5, 6.
##   Performed in total 4000 partitions.
## 
## Following methods can be applied to this 'ConsensusPartitionList' object:
##  [1] "cola_report"           "collect_classes"       "collect_plots"         "get_anno_col"         
##  [5] "get_anno"              "get_classes"           "get_matrix"            "get_membership"       
##  [9] "get_stat"              "guess_best_k"          "show"                  "test_to_known_factors"
## [13] "top_rows_heatmap"      "top_rows_overlap"     
## 
## You can get result for a single method by its node id, e.g. object["sd", "hclust"] or object["sd:hclust"]
## or a subset of methods by object[c("sd", "cv")], c("hclust", "kmeans")]
```

`cola_report()` function generates a detailed report (with all the tables and
plots) for this analysis which can be visited at https://jokergoo.github.io/cola_test/tcga_cola_rl_report/cola_report.html.


```r
cola_report(rl, output_dir = "~/cola_test/tcga_cola_rl_report")
```

## Hierarchical partition

Normal consensus partition methods aim to find $k$ subgroups at the same time. However, when
1. there are dominant subgroups, or 2. the number of potential subgroups are large, it is 
difficult to find secondary subgroups which show less difference with normal consensus partition
process. To solve this problem, the subgroups can be found in a hierarchical way, where dominant
subgroups are found first, later secondary subgroups are detected.

**cola** package implements hierarchical partition and the flowchart is as follows.

<img src="hierarchical_partition_workflow.png" width="600" />

Generally, at each recursive step, the consensus partition is performed for a
small set of $k$ (because the final larger number of subgroups will be
detedted hierarchically) and a subset of samples, If the PAC score of the
"best k" is less than 0.15, the samples are split into two subgroups where the
first subgroup contains samples with smallest intra-group mean distance and
the second subgroups are all other samples. For each of the two subgropus, if
the number of samples is less than 6, the hierarchial partition stops, or it
repeatedly perform the hierarchial partition on the subset of samples in the
corresponding subgroups.

Hierarchical partition is performed by `hierarchical_partition()` function.
Here you can only use one single top-value method and a single partition
methods.

In following example, we still use the TCGA GBM microarray datasets. The
consensus partition which is summarized from all methods are added as an
annotation to compare to the subgroups predicted by hierarchical partition.


```r
rh = hierarchical_partition(data, top_n = c(1000, 2000, 3000, 4000),
    top_value_method = "MAD", partition_method = "kmeans",
    anno = data.frame(subtype = subtype,
        consensus = get_classes(rl, k = 4)$class), 
    anno_col = list(subtype = structure(seq_len(4), names = unique(subtype)),
        consensus = structure(brewer.pal(4, "Set1"), names = 1:4)))
```



Simply typing `rh` gives summary of the analysis.


```r
rh
```

```
## A 'HierarchicalPartition' object with 'MAD:kmeans' method.
##   On a matrix with 11268 rows and 173 columns.
##   Performed in total 5400 partitions.
##   There are 5 groups.
## 
## Hierarchy of the partition:
##   0, 173 cols
##   |-- 02, 38 cols
##   `-- 00, 135 cols
##       |-- 001, 52 cols
##       |   |-- 0011, 37 cols
##       |   `-- 0010, 15 cols
##       `-- 000, 83 cols
##           |-- 0001, 46 cols
##           `-- 0000, 37 cols
## 
## Following methods can be applied to this 'HierarchicalPartition' object:
##  [1] "all_leaves"            "all_nodes"             "cola_report"           "collect_classes"      
##  [5] "dimension_reduction"   "get_anno_col"          "get_anno"              "get_classes"          
##  [9] "get_matrix"            "get_signatures"        "max_depth"             "show"                 
## [13] "test_to_known_factors"
## 
## You can get result for a single node by e.g. object["01"]
```

**cola** uses a special way to encode the node in the hierarchy. The length of the node name
is the depth of the node in the hierarchy and the substring excluding the last digit is the name
node of the parent node. E.g. for the node `0011`, the depth is 4 and the parent node is `001`.

`collect_classes()` plots the hierarchial of subgroups as well as the annotations set before.


```r
collect_classes(rh)
```

<img src="figure/unnamed-chunk-37-1.png" title="plot of chunk unnamed-chunk-37" alt="plot of chunk unnamed-chunk-37" style="display: block; margin: auto;" />

In above plot, generally hierarchical partition found similar subgroups as
consensus partition, but hierarchical partition additionally found two
subgroups for the Mesenchymal subtype samples, or subgroup 1 in consensus
partition (the `rl` object), which makes totally 5 subgroups. However, if we
directly check the 5 subgroups in `rl`, actually the subgrouping is not
stable. This means the difference between the two subgroups in Mesenchymal is
so small that it cannot be distinguished if we consider all samples in the
analysis.


```r
consensus_heatmap(rl["MAD:kmeans"], k = 5, show_row_names = FALSE)
```

<img src="figure/unnamed-chunk-38-1.png" title="plot of chunk unnamed-chunk-38" alt="plot of chunk unnamed-chunk-38" style="display: block; margin: auto;" />

In the hierarchical partition, the two subgroups of Mesenchymal subtype are
under the node `001`. If only applying consensus partition on node `001`,
actually there are two obvious subgroups and there are quite a lot signature
genes which are significantly different between the two subgroups. This proves
it is very meanningful that there are two secondary subgroups in Mesenchymal subtype.


```r
get_signatures(rh["001"], k = 2, show_column_names = FALSE)
```

<img src="figure/unnamed-chunk-39-1.png" title="plot of chunk unnamed-chunk-39" alt="plot of chunk unnamed-chunk-39" style="display: block; margin: auto;" />

Similar as the `ConsensusPartitionList` object, `cola_report()` function
can also be applied to this `HierarchicalPartition` object. The full report
for `rh` can be found at https://jokergoo.github.io/cola_test/tcga_cola_rh_report/cola_report.html.


```r
cola_report(rh, output_dir = "~/cola_test/tcga_cola_rh_report")
```

## Session info


```r
sessionInfo()
```

```
## R version 3.4.4 (2018-03-15)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.2
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] GetoptLong_0.1.6      mvtnorm_1.0-8         matrixStats_0.53.1    gridBase_0.4-7       
##  [5] cola_0.99.0           circlize_0.4.4        ComplexHeatmap_1.17.1 markdown_0.8         
##  [9] knitr_1.20            colorout_1.2-0       
## 
## loaded via a namespace (and not attached):
##   [1] bitops_1.0-6         bit64_0.9-7          RColorBrewer_1.1-2   httr_1.3.1          
##   [5] prabclus_2.2-6       data.tree_0.7.5      tools_3.4.4          R6_2.2.2            
##   [9] KernSmooth_2.23-15   DBI_1.0.0            BiocGenerics_0.24.0  lazyeval_0.2.1      
##  [13] colorspace_1.3-2     trimcluster_0.1-2    nnet_7.3-12          tidyselect_0.2.4    
##  [17] gridExtra_2.3        bit_1.1-14           compiler_3.4.4       Biobase_2.38.0      
##  [21] xml2_1.2.0           influenceR_0.1.0     microbenchmark_1.4-4 slam_0.1-43         
##  [25] samr_2.0             diptest_0.75-7       caTools_1.17.1       scales_0.5.0        
##  [29] DEoptimR_1.0-8       robustbase_0.93-0    genefilter_1.60.0    readr_1.1.1         
##  [33] stringr_1.3.1        digest_0.6.15        pkgconfig_2.0.1      htmltools_0.3.6     
##  [37] highr_0.7            htmlwidgets_1.2      rlang_0.2.1          GlobalOptions_0.1.0 
##  [41] RSQLite_2.1.1        rstudioapi_0.7       impute_1.52.0        shape_1.4.4         
##  [45] bindr_0.1.1          visNetwork_2.0.3     jsonlite_1.5         mclust_5.4          
##  [49] gtools_3.5.0         dendextend_1.8.0     dplyr_0.7.5          rgexf_0.15.3        
##  [53] RCurl_1.95-4.10      magrittr_1.5         modeltools_0.2-21    Matrix_1.2-14       
##  [57] S4Vectors_0.16.0     Rcpp_0.12.17         munsell_0.4.3        viridis_0.5.1       
##  [61] stringi_1.2.2        whisker_0.3-2        MASS_7.3-50          pamr_1.55           
##  [65] flexmix_2.3-14       gplots_3.0.1         Rtsne_0.13           plyr_1.8.4          
##  [69] blob_1.1.1           parallel_3.4.4       gdata_2.18.0         crayon_1.3.4        
##  [73] lattice_0.20-35      splines_3.4.4        annotate_1.56.2      hms_0.4.2           
##  [77] pillar_1.2.3         igraph_1.2.1         rjson_0.2.20         fpc_2.1-11          
##  [81] stats4_3.4.4         XML_3.98-1.11        glue_1.2.0           evaluate_0.10.1     
##  [85] downloader_0.4       png_0.1-7            gtable_0.2.0         purrr_0.2.5         
##  [89] tidyr_0.8.1          clue_0.3-55          kernlab_0.9-26       assertthat_0.2.0    
##  [93] ggplot2_2.2.1        xtable_1.8-2         skmeans_0.2-11       class_7.3-14        
##  [97] survival_2.42-3      viridisLite_0.3.0    tibble_1.4.2         IRanges_2.12.0      
## [101] memoise_1.1.0        AnnotationDbi_1.40.0 bindrcpp_0.2.2       cluster_2.0.7-1     
## [105] Rook_1.1-1           DiagrammeR_1.0.0     brew_1.0-6
```
