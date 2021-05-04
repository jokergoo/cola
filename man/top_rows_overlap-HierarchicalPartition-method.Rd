\name{top_rows_overlap-HierarchicalPartition-method}
\alias{top_rows_overlap,HierarchicalPartition-method}
\title{
Overlap of top rows on different nodes
}
\description{
Overlap of top rows on different nodes
}
\usage{
\S4method{top_rows_overlap}{HierarchicalPartition}(object, method = c("euler", "upset", "venn"), fill = NULL, ...)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{method}{\code{euler}: plot Euler diagram by \code{\link[eulerr]{euler}};  \code{upset}: draw the Upset plot by \code{\link[ComplexHeatmap]{UpSet}}; \code{venn}: plot Venn diagram by \code{\link[gplots]{venn}};  \code{correspondance}: use \code{\link{correspond_between_rankings}}.}
  \item{fill}{Filled color for the Euler diagram. The value should be a color vector. Transparency of 0.5 are added internally.}
  \item{...}{Additional arguments passed to \code{\link[eulerr]{plot.euler}}, \code{\link[ComplexHeatmap]{UpSet}} or \code{\link{correspond_between_rankings}}.}

}
\value{
No value is returned.
}
\seealso{
\code{\link{top_elements_overlap}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_rh)
top_rows_overlap(golub_cola_rh, method = "euler")
top_rows_overlap(golub_cola_rh, method = "upset")
top_rows_overlap(golub_cola_rh, method = "venn")
}
