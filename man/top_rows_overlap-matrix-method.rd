\name{top_rows_overlap-matrix-method}
\alias{top_rows_overlap,matrix-method}
\title{
Overlap of top rows from different top methods
}
\description{
Overlap of top rows from different top methods
}
\usage{
\S4method{top_rows_overlap}{matrix}(object, top_method = ALL_TOP_VALUE_METHOD(), top_n = round(0.25*nrow(object)),
    type = c("venn", "correspondance"))
}
\arguments{

  \item{object}{a numeric matrix}
  \item{top_method}{methods defined in \code{\link{ALL_TOP_VALUE_METHOD}}.}
  \item{top_n}{number of top rows}
  \item{type}{\code{venn}: use venn euler plots; \code{correspondance}: use \code{\link{correspond_between_rankings}}.}

}
\examples{
mat = matrix(rnorm(1000), nrow = 100)
top_rows_overlap(mat, top_n = 25)
}
