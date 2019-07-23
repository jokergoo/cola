\name{get_signatures-HierarchicalPartition-method}
\alias{get_signatures,HierarchicalPartition-method}
\title{
Get signatures rows
}
\description{
Get signatures rows
}
\usage{
\S4method{get_signatures}{HierarchicalPartition}(object, depth = max_depth(object),
    scale_rows = object[1]@scale_rows,
    anno = get_anno(object),
    anno_col = get_anno_col(object),
    show_column_names = FALSE,
    verbose = TRUE, plot = TRUE,
    silhouette_cutoff = 0.5,
    ...)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{Depth of the hierarchy.}
  \item{scale_rows}{Whether apply row scaling when making the heatmap.}
  \item{anno}{A data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{hierarchical_partition}}.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{show_column_names}{Whether show column names in the heatmap.}
  \item{verbose}{Whether to print messages.}
  \item{plot}{Whether to make the plot.}
  \item{silhouette_cutoff}{Cutoff for silhouette scores. Samples with values  less than it are not used for finding signature rows. For selecting a  proper silhouette cutoff, please refer to \url{https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.}}
  \item{...}{Other arguments}

}
\details{
The function calls \code{\link{get_signatures,ConsensusPartition-method}} to find signatures at
each node of the partition hierarchy. The final signatures are the union of all signatures
at all nodes.
}
\value{
A list of row indices where rows are significantly different between subgroups in at least one node.
Other columns in the returned data frames are whether the rows are significantly different in the node.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
data(cola_rh)
get_signatures(cola_rh)
}
}
