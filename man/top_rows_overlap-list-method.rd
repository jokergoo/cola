\name{top_rows_overlap-list-method}
\alias{top_rows_overlap,list-method}
\title{
Overlap of top rows from different top value methods
}
\description{
Overlap of top rows from different top value methods
}
\usage{
\S4method{top_rows_overlap}{list}(object, top_n = round(0.25*length(object[[1]])),
    method = c("venn", "venneuler", "correspondance"), ...)
}
\arguments{

  \item{object}{a list which contains rankings from different metrics.}
  \item{top_n}{number of top rows.}
  \item{method}{\code{venn}: use Venn diagram; \code{venneuler}: use Venn Euler diagram; \code{correspondance}: use \code{\link{correspond_between_rankings}}.}
  \item{...}{additional arguments passed to \code{\link{venn_euler}} or \code{\link{correspond_between_rankings}}.}

}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(matrixStats)
set.seed(123)
mat = matrix(rnorm(1000), nrow = 100)
lt = list(sd = rowSds(mat), mad = rowMads(mat))
top_rows_overlap(lt, top_n = 25)
}
