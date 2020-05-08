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
    method = c("euler", "venn", "correspondance"), fill = NULL, ...)
}
\arguments{

  \item{object}{A list which contains values from different metrics.}
  \item{top_n}{Number of top rows.}
  \item{method}{\code{euler}: plot Euler diagram by \code{\link[eulerr]{euler}}; \code{venn}: plot Venn diagram by \code{\link[gplots]{venn}};  \code{correspondance}: use \code{\link{correspond_between_rankings}}.}
  \item{fill}{Filled color for the Euler diagram. The value should be a color vector. Transparency of 0.5 are added internally.}
  \item{...}{Additional arguments passed to \code{\link[eulerr]{plot.euler}} or \code{\link{correspond_between_rankings}}.}

}
\details{
The i^th value in every vectors in \code{object} should correspond to the same element from the original data.
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
top_elements_overlap(lt, top_n = 25, method = "venn")
top_elements_overlap(lt, top_n = 25, method = "correspondance")
}
