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
    ...)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{depth of the hierarchy.}
  \item{scale_rows}{whether apply row scaling when making the heatmap.}
  \item{anno}{a data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{hierarchical_partition}}.}
  \item{anno_col}{a list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{show_column_names}{whether show column names in the heatmap.}
  \item{verbose}{whether to print messages.}
  \item{plot}{whether to make the plot.}
  \item{...}{other arguments pass to \code{\link{get_signatures,ConsensusPartition-method}}.}

}
\details{
The function calls \code{\link{get_signatures,ConsensusPartition-method}} to find signatures at
each node of the partition hierarchy.
}
\value{
A list of row indices where rows are significantly different between subgroups in at least one node.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
