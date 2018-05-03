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
    anno = get_anno(object[1]),
    anno_col = get_anno_col(object[1]),
    show_column_names = TRUE,
    verbose = TRUE, plot = TRUE,
    silhouette_cutoff = -Inf,
    ...)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{minimal depth of the hierarchy}
  \item{scale_rows}{whether scale rows}
  \item{anno}{annotation}
  \item{anno_col}{annotation color}
  \item{show_column_names}{whether show column names}
  \item{verbose}{whether print messages}
  \item{plot}{whether make the heatmap}
  \item{silhouette_cutoff}{cutoff for silhouette scores.}
  \item{...}{other arguments}

}
\details{
The function called \code{\link{get_signatures,ConsensusPartition-method}} to find signatures at
each level of the partition hierarchy.
}
\value{
A list of signature names (row names of the original data matrix)
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
