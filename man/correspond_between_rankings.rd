\name{correspond_between_rankings}
\alias{correspond_between_rankings}
\title{
Correspond between a list of rankings
}
\description{
Correspond between a list of rankings
}
\usage{
correspond_between_rankings(lt, top_n = length(lt[[1]]),
    col = brewer_pal_set1_col[1:length(lt)], ...)
}
\arguments{

  \item{lt}{A list of scores under different metrics.}
  \item{top_n}{Top n elements to show correspondance.}
  \item{col}{A vector of colors for \code{lt}.}
  \item{...}{Pass to \code{\link{correspond_between_two_rankings}}.}

}
\details{
It makes plots for every pairwise comparisons in \code{lt}.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(matrixStats)
mat = matrix(runif(1000), ncol = 10)
x1 = rowSds(mat)
x2 = rowMads(mat)
x3 = rowSds(mat)/rowMeans(mat)
correspond_between_rankings(lt = list(sd = x1, mad = x2, vc = x3), 
    top_n = 20, col = c("red", "blue", "green"))
}
