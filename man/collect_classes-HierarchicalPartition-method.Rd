\name{collect_classes-HierarchicalPartition-method}
\alias{collect_classes,HierarchicalPartition-method}
\title{
Collect classes from HierarchicalPartition object
}
\description{
Collect classes from HierarchicalPartition object
}
\usage{
\S4method{collect_classes}{HierarchicalPartition}(object, depth = max_depth(object),
    anno = get_anno(object[1]), anno_col = get_anno_col(object[1]))
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{Depth of the hierarchy.}
  \item{anno}{A data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{hierarchical_partition}}.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}

}
\details{
The function plots the hierarchy of the classes.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_rh)
collect_classes(golub_cola_rh)
collect_classes(golub_cola_rh, depth = 2)
}
