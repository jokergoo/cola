\name{top_rows_heatmap-matrix-method}
\alias{top_rows_heatmap,matrix-method}
\title{
Heatmap of top rows from different top methods
}
\description{
Heatmap of top rows from different top methods
}
\usage{
\S4method{top_rows_heatmap}{matrix}(object, all_value_list = NULL, top_method = ALL_TOP_VALUE_METHOD(),
    top_n = round(0.25*nrow(object)))
}
\arguments{

  \item{object}{a numeric matrix}
  \item{all_value_list}{scores that have already been calculated from the matrix}
  \item{top_method}{methods defined in \code{\link{ALL_TOP_VALUE_METHOD}}.}
  \item{top_n}{number of top rows}

}
\examples{
mat = matrix(rnorm(1000), nrow = 100)
top_rows_heatmap(mat, top_n = 25)
}
