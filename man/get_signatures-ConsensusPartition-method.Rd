\name{get_signatures-ConsensusPartition-method}
\alias{get_signatures,ConsensusPartition-method}
\title{
Get signature rows
}
\description{
Get signature rows
}
\usage{
\S4method{get_signatures}{ConsensusPartition}(object, k,
    silhouette_cutoff = 0.5,
    fdr_cutoff = cola_opt$fdr_cutoff,
    top_signatures = NULL,
    group_diff = cola_opt$group_diff,
    scale_rows = object@scale_rows,
    row_km = NULL,
    diff_method = c("Ftest", "ttest", "samr", "pamr", "one_vs_others"),
    anno = get_anno(object),
    anno_col = get_anno_col(object),
    internal = FALSE,
    show_row_dend = FALSE,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    use_raster = TRUE,
    plot = TRUE, verbose = TRUE, seed = 888,
    left_annotation = NULL, right_annotation = NULL,
    col = if(scale_rows) c("green", "white", "red") else c("blue", "white", "red"),
    simplify = FALSE, prefix = "", enforce = FALSE, hash = NULL,
    ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{silhouette_cutoff}{Cutoff for silhouette scores. Samples with values  less than it are not used for finding signature rows. For selecting a  proper silhouette cutoff, please refer to \url{https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.}}
  \item{fdr_cutoff}{Cutoff for FDR of the difference test between subgroups.}
  \item{top_signatures}{Top signatures with most significant fdr. Note since fdr might be same for multiple rows, the final number of signatures might not be exactly the same as the one that has been set.}
  \item{group_diff}{Cutoff for the maximal difference between group means.}
  \item{scale_rows}{Whether apply row scaling when making the heatmap.}
  \item{row_km}{Number of groups for performing k-means clustering on rows. By default it is automatically selected.}
  \item{diff_method}{Methods to get rows which are significantly different between subgroups, see 'Details' section.}
  \item{anno}{A data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{internal}{Used internally.}
  \item{show_row_dend}{Whether show row dendrogram.}
  \item{show_column_names}{Whether show column names in the heatmap.}
  \item{column_names_gp}{Graphics parameters for column names.}
  \item{use_raster}{Internally used.}
  \item{plot}{Whether to make the plot.}
  \item{verbose}{Whether to print messages.}
  \item{seed}{Random seed.}
  \item{left_annotation}{Annotation put on the left of the heatmap. It should be a \code{\link[ComplexHeatmap]{HeatmapAnnotation-class}} object.  The number of items should be the same as the number of the original matrix rows. The subsetting to the significant  rows are automatically performed on the annotation object.}
  \item{right_annotation}{Annotation put on the right of the heatmap. Same format as \code{left_annotation}.}
  \item{col}{Colors.}
  \item{simplify}{Only used internally.}
  \item{prefix}{Only used internally.}
  \item{enforce}{The analysis is cached by default, so that the analysis with the same input will be automatically extracted without rerunning them. Set \code{enforce} to \code{TRUE} to enforce the funtion to re-perform the analysis.}
  \item{hash}{Userd internally.}
  \item{...}{Other arguments.}

}
\details{
Basically the function applies statistical test for the difference in subgroups for every
row. There are following methods which test significance of the difference:

\describe{
  \item{ttest}{First it looks for the subgroup with highest mean value, compare to each of the  other subgroups with t-test and take the maximum p-value. Second it looks for the subgroup with lowest mean value, compare to each of the other subgroups again with t-test and take the maximum p-values. Later for these two list of p-values take the minimal p-value as the final p-value. }
  \item{samr/pamr}{use SAM (from samr package)/PAM (from pamr package) method to find significantly different rows between subgroups.}
  \item{Ftest}{use F-test to find significantly different rows between subgroups.}
  \item{one_vs_others}{For each subgroup i in each row, it uses t-test to compare samples in current  subgroup to all other samples, denoted as p_i. The p-value for current row is selected as min(p_i).}
}

\code{diff_method} can also be a self-defined function. The function needs two arguments which are the matrix for the analysis
and the predicted classes. The function should returns a vector of FDR from the difference test.
}
\value{
A data frame with more than two columns:

\describe{
  \item{\code{which_row}:}{row index corresponding to the original matrix.}
  \item{\code{fdr}:}{the FDR.}
  \item{\code{km}:}{the k-means groups if \code{row_km} is set.}
  \item{other_columns:}{the mean value (depending rows are scaled or not) in each subgroup.}
}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
res = golub_cola["ATC", "skmeans"]
tb = get_signatures(res, k = 3)
head(tb)
get_signatures(res, k = 3, top_signatures = 100)
}
