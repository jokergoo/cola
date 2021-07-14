\name{top_rows_heatmap-matrix-method}
\alias{top_rows_heatmap,matrix-method}
\title{
Heatmap of top rows from different top-value methods
}
\description{
Heatmap of top rows from different top-value methods
}
\usage{
\S4method{top_rows_heatmap}{matrix}(object, all_top_value_list = NULL,
    top_value_method = all_top_value_methods(),
    bottom_annotation = NULL,
    top_n = round(0.25*nrow(object)), scale_rows = TRUE, ...)
}
\arguments{

  \item{object}{A numeric matrix.}
  \item{all_top_value_list}{Top-values that have already been calculated from the matrix. If it is \code{NULL} the values are calculated by methods in \code{top_value_method} argument.}
  \item{top_value_method}{Methods defined in \code{\link{all_top_value_methods}}.}
  \item{bottom_annotation}{A \code{\link[ComplexHeatmap]{HeatmapAnnotation-class}} object.}
  \item{top_n}{Number of top rows to show in the heatmap.}
  \item{scale_rows}{Whether to scale rows.}
  \item{...}{Pass to \code{\link[ComplexHeatmap]{Heatmap}}.}

}
\details{
The function makes heatmaps where the rows are scaled (or not scaled) for the top n rows
from different top-value methods.

The top n rows are used for subgroup classification in cola analysis, so the heatmaps show which 
top-value method gives better candidate rows for the classification.
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
