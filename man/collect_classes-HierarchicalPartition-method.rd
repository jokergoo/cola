\name{collect_classes-HierarchicalPartition-method}
\alias{collect_classes,HierarchicalPartition-method}
\title{
Collect classes from hierarchical_partition object
}
\description{
Collect classes from hierarchical_partition object
}
\usage{
\S4method{collect_classes}{HierarchicalPartition}(object, depth = max_depth(object),
    anno = get_anno(object[1]), anno_col = get_anno_col(object[1]), ...)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{minimal depth of the hierarchy}
  \item{anno}{a data frame with known annotation of samples.}
  \item{anno_col}{a list of colors for the annotations in \code{anno}.}
  \item{...}{other arguments.}

}
\details{
The function plots the hierarchy of the subgroups.
}
\value{
No value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
