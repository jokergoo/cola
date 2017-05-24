\name{dimension_reduction-ConsensusPartition-method}
\alias{dimension_reduction,ConsensusPartition-method}
\alias{dimension_reduction}
\title{
Visualize columns after dimension reduction
}
\description{
Visualize columns after dimension reduction
}
\usage{
\S4method{dimension_reduction}{ConsensusPartition}(object, k, top_n = NULL,
    method = c("mds", "pca", "tsne"),
    silhouette_cutoff = 0.5, remove = FALSE,
    tsne_param = list(), ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{k}{number of partitions}
  \item{top_n}{top n rows to use.}
  \item{method}{which method to reduce the dimension of the data. \code{mds} uses \code{\link[stats]{cmdscale}}, \code{pca} uses \code{stats::prcomp} and \code{tsne} uses \code{\link[Rtsne]{Rtsne}}.}
  \item{silhouette_cutoff}{cutoff of silhouette. Data points with values less than it will be mapped to small points.}
  \item{remove}{whether to remove columns which have less silhouette values than the cutoff.}
  \item{tsne_param}{parameters pass to \code{\link[Rtsne]{Rtsne}}}
  \item{...}{other arguments}

}
\examples{
# There is no example
NULL

}
