\name{dimension_reduction-HierarchicalPartition-method}
\alias{dimension_reduction,HierarchicalPartition-method}
\title{
Visualize columns after dimension reduction
}
\description{
Visualize columns after dimension reduction
}
\usage{
\S4method{dimension_reduction}{HierarchicalPartition}(object,
    depth = max_depth(object), parent_node,
    top_n = NULL, method = c("PCA", "MDS", "t-SNE", "UMAP"),
    silhouette_cutoff = 0.5, scale_rows = TRUE, verbose = TRUE)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{Depth of the hierarchy.}
  \item{top_n}{Top n rows to use. By default it uses all rows in the original matrix.}
  \item{parent_node}{Parent node. If it is set, the function call is identical to \code{dimension_reduction(object[parent_node])}}
  \item{method}{Which method to reduce the dimension of the data. \code{MDS} uses \code{\link[stats]{cmdscale}}, \code{PCA} uses \code{\link[stats]{prcomp}}. \code{t-SNE} uses \code{\link[Rtsne]{Rtsne}}. \code{UMAP} uses \code{\link[umap]{umap}}.}
  \item{silhouette_cutoff}{Cutoff of silhouette score. Data points with values less than it will be mapped to small points.}
  \item{scale_rows}{Whether perform scaling on matrix rows.}
  \item{verbose}{Whether print messages.}

}
\details{
The class IDs are extract at \code{depth}.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
dimension_reduction(cola_rh)
dimension_reduction(cola_rh, parent_node = "00")
}
