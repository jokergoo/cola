\name{dimension_reduction-ConsensusPartition-method}
\alias{dimension_reduction,ConsensusPartition-method}
\title{
Visualize column after dimension reduction
}
\description{
Visualize samples (the matrix columns) after dimension reduction
}
\usage{
\S4method{dimension_reduction}{ConsensusPartition}(object, k, top_n = NULL,
    method = c("PCA", "MDS"), internal = FALSE,
    silhouette_cutoff = 0.5, remove = FALSE,
    scale_rows = TRUE, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions.}
  \item{top_n}{top n rows to use. By default it uses all rows in the original matrix.}
  \item{method}{which method to reduce the dimension of the data. \code{MDS} uses \code{\link[stats]{cmdscale}}, \code{PCA} uses \code{\link[stats]{prcomp}}.}
  \item{internal}{internally used.}
  \item{silhouette_cutoff}{cutoff of silhouette score. Data points with values less than it will be mapped with cross symbols.}
  \item{remove}{whether to remove columns which have less silhouette scores than the cutoff.}
  \item{scale_rows}{whether perform scaling on matrix rows.}
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
