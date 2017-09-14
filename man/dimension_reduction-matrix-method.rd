\name{dimension_reduction-matrix-method}
\alias{dimension_reduction,matrix-method}
\title{
Visualize columns after dimension reduction
}
\description{
Visualize columns after dimension reduction
}
\usage{
\S4method{dimension_reduction}{matrix}(object,
    pch = 16, col = "black", cex = 1, main = "",
    method = c("mds", "pca", "tsne"),
    tsne_param = list())
}
\arguments{

  \item{object}{a numeric matrix}
  \item{method}{which method to reduce the dimension of the data. \code{mds} uses \code{\link[stats]{cmdscale}}, \code{pca} uses \code{\link[stats]{prcomp}} and \code{tsne} uses \code{\link[Rtsne]{Rtsne}}.}
  \item{pch}{shape of points}
  \item{col}{color of points}
  \item{cex}{size of points}
  \item{main}{title of the plot}
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
