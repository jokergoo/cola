\name{get_signatures-HierarchicalPartition-method}
\alias{get_signatures,HierarchicalPartition-method}
\title{
Get signatures rows
}
\description{
Get signatures rows
}
\usage{
\S4method{get_signatures}{HierarchicalPartition}(object, depth = max_depth(object), node = "0",
    scale_rows = object[1]@scale_rows,
    anno = get_anno(object),
    anno_col = get_anno_col(object),
    show_column_names = FALSE,
    verbose = TRUE, plot = TRUE,
    silhouette_cutoff = 0.5,
    ...)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{depth of the hierarchy.}
  \item{scale_rows}{whether apply row scaling when making the heatmap.}
  \item{anno}{a data frame with known annotation of samples.}
  \item{anno_col}{a list of colors for the annotations in \code{anno}.}
  \item{show_column_names}{whether show column names in the heatmap.}
  \item{verbose}{whether to print messages.}
  \item{plot}{whether to make the plot.}
  \item{silhouette_cutoff}{cutoff for silhouette scores. Samples with values  less than it are not used for finding signature rows. For selecting a  proper silhouette cutoff, please refer to \url{https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.}}
  \item{...}{other arguments}

}
\details{
The function called \code{\link{get_signatures,ConsensusPartition-method}} to find signatures at
each level of the partition hierarchy.
}
\value{
A list of signature indices (row indices of the original data matrix).
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
