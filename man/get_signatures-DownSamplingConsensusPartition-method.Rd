\name{get_signatures-DownSamplingConsensusPartition-method}
\alias{get_signatures,DownSamplingConsensusPartition-method}
\title{
Get signature rows
}
\description{
Get signature rows
}
\usage{
\S4method{get_signatures}{DownSamplingConsensusPartition}(object, k,
    p_cutoff = 0.05,
    fdr_cutoff = cola_opt$fdr_cutoff,
    group_diff = cola_opt$group_diff,
    scale_rows = object@scale_rows,
    row_km = NULL,
    diff_method = c("Ftest", "ttest", "samr", "pamr", "one_vs_others", "uniquely_high_in_one_group"),
    anno = get_anno(object),
    anno_col = get_anno_col(object),
    internal = FALSE,
    show_row_dend = FALSE,
    show_column_names = FALSE, use_raster = TRUE,
    plot = TRUE, verbose = TRUE, seed = 888,
    left_annotation = NULL, right_annotation = NULL,
    col = if(scale_rows) c("green", "white", "red") else c("blue", "white", "red"),
    simplify = FALSE, prefix = "",
    ...)
}
\arguments{

  \item{object}{A \code{\link{DownSamplingConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{p_cutoff}{Cutoff for p-values of class label prediction. Samples with values  higher than it are not used for finding signature rows.}
  \item{fdr_cutoff}{Cutoff for FDR of the difference test between subgroups.}
  \item{group_diff}{Cutoff for the maximal difference between group means.}
  \item{scale_rows}{Whether apply row scaling when making the heatmap.}
  \item{row_km}{Number of groups for performing k-means clustering on rows. By default it is automatically selected.}
  \item{diff_method}{Methods to get rows which are significantly different between subgroups, see 'Details' section.}
  \item{anno}{A data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{internal}{Used internally.}
  \item{show_row_dend}{Whether show row dendrogram.}
  \item{show_column_names}{Whether show column names in the heatmap.}
  \item{use_raster}{Internally used.}
  \item{plot}{Whether to make the plot.}
  \item{verbose}{Whether to print messages.}
  \item{seed}{Random seed.}
  \item{left_annotation}{Annotation put on the left of the heatmap. It should be a \code{\link[ComplexHeatmap]{HeatmapAnnotation-class}} object.  The number of items should be the same as the number of the original matrix rows. The subsetting to the significant  rows are automatically performed on the annotation object.}
  \item{right_annotation}{Annotation put on the right of the heatmap. Same format as \code{left_annotation}.}
  \item{col}{Colors.}
  \item{simplify}{Only use internally.}
  \item{prefix}{Only use internally.}
  \item{...}{Other arguments.}

}
\details{
This function is very similar as \code{\link{get_signatures,ConsensusPartition-method}}.
}
\examples{
\donttest{
data(golub_cola_ds)
get_signatures(golub_cola_ds, k = 2)
get_signatures(golub_cola_ds, k = 3)
}
}
