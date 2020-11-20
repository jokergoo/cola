\name{get_signatures-HierarchicalPartition-method}
\alias{get_signatures,HierarchicalPartition-method}
\title{
Get signatures rows
}
\description{
Get signatures rows
}
\usage{
\S4method{get_signatures}{HierarchicalPartition}(object, merge_node = merge_node_param(),
    group_diff = cola_opt$group_diff,
    row_km = NULL,
    scale_rows = object[1]@scale_rows,
    anno = get_anno(object),
    anno_col = get_anno_col(object),
    show_column_names = FALSE, column_names_gp = gpar(fontsize = 8),
    verbose = TRUE, plot = TRUE, seed = 888,
    ...)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{merge_node}{Parameters to merge sub-dendrograms, see \code{\link{merge_node_param}}.}
  \item{group_diff}{Cutoff for the maximal difference between group means.}
  \item{row_km}{Number of groups for performing k-means clustering on rows. By default it is automatically selected.}
  \item{scale_rows}{whether apply row scaling when making the heatmap.}
  \item{anno}{a data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{hierarchical_partition}}.}
  \item{anno_col}{a list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{show_column_names}{whether show column names in the heatmap.}
  \item{column_names_gp}{Graphic parameters for column names.}
  \item{verbose}{whether to print messages.}
  \item{plot}{whether to make the plot.}
  \item{seed}{Random seed.}
  \item{...}{other arguments pass to \code{\link{get_signatures,ConsensusPartition-method}}.}

}
\details{
The function calls \code{\link{get_signatures,ConsensusPartition-method}} to find signatures at
each node of the partition hierarchy.
}
\value{
A data frame with more than two columns:

\describe{
  \item{\code{which_row}:}{row index corresponding to the original matrix.}
  \item{\code{km}:}{the k-means groups if \code{row_km} is set.}
  \item{other_columns:}{the mean value (depending rows are scaled or not) in each subgroup.}
}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\donttest{
data(golub_cola_rh)
tb = get_signatures(golub_cola_rh)
head(tb)
}
}
