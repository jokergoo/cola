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
    top_n = NULL, method = c("PCA", "MDS"),
    silhouette_cutoff = 0.5, scale = TRUE)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{depth of the hierarchy.}
  \item{top_n}{top n rows to use. By default it uses all rows in the original matrix.}
  \item{parent_node}{parent node. If it is set, the function call is identical to \code{dimension_reduction(object[parent_node])}}
  \item{method}{which method to reduce the dimension of the data. \code{mds} uses \code{\link[stats]{cmdscale}}, \code{pca} uses \code{\link[stats]{prcomp}}.}
  \item{silhouette_cutoff}{cutoff of silhouette score. Data points with values less than it will be mapped to small points.}
  \item{scale}{whether perform scaling on matrix rows.}

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
