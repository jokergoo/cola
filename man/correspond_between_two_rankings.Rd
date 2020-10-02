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

  \item{x1}{A vector of scores calculated by one metric.}
  \item{x2}{A vector of scores calculated by another metric.}
  \item{name1}{Name of the first metric.}
  \item{name2}{Name of the second metric.}
  \item{col1}{Color for the first metric.}
  \item{col2}{Color for the second metric.}
  \item{top_n}{Top n elements to show the correspondance.}
  \item{transparency}{Transparency of the connecting lines.}
  \item{pt_size}{Size of the points, must be a \code{\link[grid]{unit}} object.}
  \item{newpage}{Whether to plot in a new graphic page.}
  \item{ratio}{Ratio of width of the left barplot, connection lines and right barplot. The three values will be scaled to a sum of 1.}

}
\details{
In \code{x1} and \code{x2}, the i^th element in both vectors corresponds to the same object (e.g. same row if they are calculated from a matrix) but with different 
scores under different metrics.

\code{x1} and \code{x2} are sorted in the left panel and right panel respectively. The top n elements
under corresponding metric are highlighted by vertical colored lines in both panels.
The left and right panels also shown as barplots of the scores in the two metrics.
Between the left and right panels, there are lines connecting the same element (e.g. i^th element in \code{x1} and \code{x2})
in the two ordered vectors so that you can see how a same element has two different ranks in the two metrics.

Under the plot is a simple Venn diagram showing the overlaps of the top n elements 
by the two metrics.
}
\value{
No value is returned.
}
\seealso{
\code{\link{correspond_between_rankings}} draws for more than 2 sets of rankings.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(matrixStats)
mat = matrix(runif(1000), ncol = 10)
x1 = rowSds(mat)
x2 = rowMads(mat)
correspond_between_two_rankings(x1, x2, name1 = "SD", name2 = "MAD", top_n = 20)
}
