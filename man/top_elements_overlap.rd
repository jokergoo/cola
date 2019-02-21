\name{top_elements_overlap}
\alias{top_elements_overlap}
\title{
Overlap of top elements from different metrics
}
\description{
Overlap of top elements from different metrics
}
\usage{
top_elements_overlap(object, top_n = round(0.25*length(object[[1]])),
    method = c("venn", "venneuler", "correspondance"), ...)
}
\arguments{

  \item{object}{a list which contains values from different metrics.}
  \item{top_n}{number of top rows.}
  \item{method}{\code{venn}: use Venn diagram; \code{venneuler}: use Venn Euler diagram; \code{correspondance}: use \code{\link{correspond_between_rankings}}.}
  \item{...}{additional arguments passed to \code{\link{venn_euler}} or \code{\link{correspond_between_rankings}}.}

}
\details{
The i^th value in all vectors in the input should correspond to a same element from the original data.
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
top_element_overlap(lt, top_n = 25, method = "venn")
top_element_overlap(lt, top_n = 25, method = "correspondance")
}
