\name{top_rows_overlap-matrix-method}
\alias{top_rows_overlap,matrix-method}
\title{
Overlap of top rows from different top methods
}
\description{
Overlap of top rows from different top methods
}
\usage{
\S4method{top_rows_overlap}{matrix}(object, top_method = all_top_value_methods(), top_n = round(0.25*nrow(object)),
    type = c("venn", "correspondance"))
}
\arguments{

  \item{object}{a numeric matrix}
  \item{top_method}{methods defined in \code{\link{all_top_value_methods}}.}
  \item{top_n}{number of top rows}
  \item{type}{\code{venn}: use Venn Euler diagram; \code{correspondance}: use \code{\link{correspond_between_rankings}}.}

}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
mat = matrix(rnorm(1000), nrow = 100)
top_rows_overlap(mat, top_n = 25)
}
