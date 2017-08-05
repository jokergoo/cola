\name{top_rows_heatmap-matrix-method}
\alias{top_rows_heatmap,matrix-method}
\title{
Heatmap of top rows from different top methods
}
\description{
Heatmap of top rows from different top methods
}
\usage{
\S4method{top_rows_heatmap}{matrix}(object, all_value_list = NULL, top_method = all_top_value_methods(),
    top_n = round(0.25*nrow(object)), layout_nr = 1)
}
\arguments{

  \item{object}{a numeric matrix}
  \item{all_value_list}{scores that have already been calculated from the matrix. If it is \code{NULL} the values are calculated by methods in \code{top_method}.}
  \item{top_method}{methods defined in \code{\link{all_top_value_methods}}.}
  \item{top_n}{number of top rows}
  \item{layout_nr}{number of rows in the layout}

}
\details{
The function makes heatmaps where the rows are scaled for the top k rows
from different top methods.

Top k rows are used to subgroup classification. The heatmaps show which top
method gives best candidate rows for the classification.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
mat = matrix(rnorm(1000), nrow = 100)
top_rows_heatmap(mat, top_n = 25)
}
