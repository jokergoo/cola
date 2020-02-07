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
    col = cola_opt$color_set_1[1:length(lt)], ...)
}
\arguments{

  \item{lt}{A list of scores under different metrics.}
  \item{top_n}{Top n elements to show the correspondance.}
  \item{col}{A vector of colors for \code{lt}.}
  \item{...}{Pass to \code{\link{correspond_between_two_rankings}}.}

}
\details{
It makes plots for every pairwise comparison in \code{lt}.
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
correspond_between_rankings(lt = list(SD = x1, MAD = x2, CV = x3), 
    top_n = 20, col = c("red", "blue", "green"))
}
