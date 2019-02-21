\name{correspond_between_two_rankings}
\alias{correspond_between_two_rankings}
\title{
Correspond two rankings
}
\description{
Correspond two rankings
}
\usage{
correspond_between_two_rankings(x1, x2, name1, name2,
    col1 = 2, col2 = 3, top_n = round(0.25*length(x1)), transparency = 0.9,
    pt_size = unit(1, "mm"), newpage = TRUE, ratio = c(1, 1, 1))
}
\arguments{

  \item{x1}{a vector of scores calculated by one metric.}
  \item{x2}{a vector of scores calculated by another metric.}
  \item{name1}{name of the first metric.}
  \item{name2}{name of the second metric.}
  \item{col1}{color for the first metric.}
  \item{col2}{color for the second metric.}
  \item{top_n}{top n elements to show correspondance.}
  \item{transparency}{transparency of the connection lines.}
  \item{pt_size}{size of the points, must be a \code{\link[grid]{unit}} object}
  \item{newpage}{whether to plot in a new graphic page.}
  \item{ratio}{ratio of width of the left barplot, connection lines and right barplot. The three values will be scaled to a sum of 1.}

}
\details{
In \code{x1} and \code{x2}, the i^{th} element is the same object (e.g. same row if they are calculated from a matrix) but with different 
scores under different metrics.

\code{x1} and \code{x2} are sorted in the left panel and right panel. The top n elements
under corresponding metric are highlighted by vertical color lines in both panels.
The left and right panels also show as barplots of the scores in the two metrics.
Between the left and right panels, there are lines connecting the same element (e.g. i^th element in \code{x1} and \code{x2})
in the two ordered vectors so that you can see how a same element has two different ranks in the two metrics.

Under the plot is a simple Venn diagram showing the overlaps of the top n elements 
by the two metrics.
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
