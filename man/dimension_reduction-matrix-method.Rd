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
    method = c("PCA", "MDS", "t-SNE", "UMAP"),
    pc = NULL, control = list(),
    scale_rows = TRUE, nr = 5000,
    internal = FALSE, verbose = TRUE)
}
\arguments{

  \item{object}{A numeric matrix.}
  \item{method}{Which method to reduce the dimension of the data. \code{MDS} uses \code{\link[stats]{cmdscale}}, \code{PCA} uses \code{\link[stats]{prcomp}}. \code{t-SNE} uses \code{\link[Rtsne]{Rtsne}}. \code{UMAP} uses \code{\link[umap]{umap}}.}
  \item{pc}{Which two principle components to visualize}
  \item{control}{A list of parameters for \code{\link[Rtsne]{Rtsne}} or \code{\link[umap]{umap}}.}
  \item{pch}{Ahape of points.}
  \item{col}{Color of points.}
  \item{cex}{Aize of points.}
  \item{main}{Title of the plot.}
  \item{scale_rows}{Whether perform scaling on matrix rows.}
  \item{nr}{If number of matrix rows is larger than this value, random \code{nr} rows are used.}
  \item{internal}{Internally used.}
  \item{verbose}{Whether print messages.}

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
