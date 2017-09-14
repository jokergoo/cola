\name{dimension_reduction-HierarchicalPartition-method}
\alias{dimension_reduction,HierarchicalPartition-method}
\title{
Visualize columns after dimension reduction
}
\description{
Visualize columns after dimension reduction
}
\usage{
\S4method{dimension_reduction}{HierarchicalPartition}(object, merge = FALSE, depth = NULL,
    top_n = NULL, method = c("mds", "pca", "tsne"),
    silhouette_cutoff = 0.5, tsne_param = list())
}
\arguments{

  \item{object}{a numeric matrix}
  \item{merge}{whether merge all samples into one single plot or make plots for every hierarchy}
  \item{depth}{depth of the hierarchy}
  \item{top_n}{top n genes to use}
  \item{method}{which method to reduce the dimension of the data. \code{mds} uses \code{\link[stats]{cmdscale}}, \code{pca} uses \code{\link[stats]{prcomp}} and \code{tsne} uses \code{\link[Rtsne]{Rtsne}}.}
  \item{silhouette_cutoff}{silhouette cutoff}
  \item{tsne_param}{parameters pass to \code{\link[Rtsne]{Rtsne}}}

}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
