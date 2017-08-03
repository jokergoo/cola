\name{correspond_between_two_rankings}
\alias{correspond_between_two_rankings}
\title{
Correspond two rankings
}
\description{
Correspond two rankings
}
\usage{
correspond_between_two_rankings(x1, x2, name1 = "", name2 = "",
    col1 = 1, col2 = 2, top_n = round(0.25*length(x1)), transparency = 0.9,
    pt_size = unit(1, "mm"), newpage = FALSE, ratio = c(1, 1, 1))
}
\arguments{

  \item{x1}{a vector of scores calculated by one metric.}
  \item{x2}{a vector of scores calculated by another metric.}
  \item{name1}{name of the first metric.}
  \item{name2}{name of the second metric.}
  \item{col1}{color of the lines for the first metric.}
  \item{col2}{color of the lines for the second metric.}
  \item{top_n}{top n elements to show correspondance.}
  \item{transparency}{transparency of the connection lines.}
  \item{pt_size}{size of the points, must be a \code{\link[grid]{unit}} object}
  \item{newpage}{whether to plot in a new graphic page.}
  \item{ratio}{ratio of width of the left plot, connection lines and right plot.}

}
\details{
In \code{x1} and \code{x2}, the i^{th} element is the same object but with different 
scores under different metrics.

\code{x1} and \code{x2} are sorted in the plot. Lines are connecting the same object
for those with top values. Shadows on the two sides are the sorted values for \code{x1}
and \code{x2}, and points corresponds to positions of the top objects but are randomly shifted
on x-directions.
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
correspond_between_two_rankings(x1, x2, name1 = "sd", name2 = "mad", top_n = 20)
}
