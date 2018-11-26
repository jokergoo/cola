\name{top_rows_overlap-matrix-method}
\alias{top_rows_overlap,matrix-method}
\title{
Overlap of top rows from different top-value methods
}
\description{
Overlap of top rows from different top-value methods
}
\usage{
\S4method{top_rows_overlap}{matrix}(object, top_value_method = all_top_value_methods(),
    top_n = round(0.25*nrow(object)), method = c("venn", "venneuler", "correspondance"), ...)
}
\arguments{

  \item{object}{a numeric matrix.}
  \item{top_value_method}{methods defined in \code{\link{all_top_value_methods}}.}
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
\seealso{
\code{\link{top_rows_overlap,list-method}}
}
\examples{
set.seed(123)
mat = matrix(rnorm(1000), nrow = 100)
top_rows_overlap(mat, top_n = 25)
}
