\name{dimension_reduction-ConsensusPartition-method}
\alias{dimension_reduction,ConsensusPartition-method}
\title{
Visualize samples after dimension reduction
}
\description{
Visualize samples after dimension reduction
}
\usage{
\S4method{dimension_reduction}{ConsensusPartition}(object, k, top_n = NULL,
    method = c("pca", "mds", "tsne"),
    silhouette_cutoff = 0.5, remove = FALSE,
    tsne_param = list(), ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions.}
  \item{top_n}{top n rows to use. By default it uses all rows in the original matrix.}
  \item{method}{which method to reduce the dimension of the data. \code{mds} uses \code{\link[stats]{cmdscale}}, \code{pca} uses \code{\link[stats]{prcomp}} and \code{tsne} uses \code{\link[Rtsne]{Rtsne}}.}
  \item{silhouette_cutoff}{cutoff of silhouette. Data points with values less than it will be mapped to small points.}
  \item{remove}{whether to remove columns which have less silhouette values than the cutoff.}
  \item{tsne_param}{parameters pass to \code{\link[Rtsne]{Rtsne}}}
  \item{...}{other arguments}

}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
dimension_reduction(cola_rl["sd", "kmeans"], k = 3)
}
