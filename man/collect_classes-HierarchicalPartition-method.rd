\name{collect_classes-HierarchicalPartition-method}
\alias{collect_classes,HierarchicalPartition-method}
\title{
Collect classes from hierarchical_partition object
}
\description{
Collect classes from hierarchical_partition object
}
\usage{
\S4method{collect_classes}{HierarchicalPartition}(object, anno = object@list[[1]]@known_anno,
    anno_col = if(missing(anno)) object@list[[1]]@known_col else NULL,
    ...)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
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
