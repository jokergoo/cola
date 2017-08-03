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

  \item{lt}{a list of scores under different metrics.}
  \item{top_n}{top n elements to show correspondance.}
  \item{col}{a vector of colors for \code{lt}.}
  \item{...}{pass to \code{\link{correspond_between_two_rankings}}.}

}
\details{
It makes plots for pairwise comparisons between different rankings.
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
    top_n = 20)
}
