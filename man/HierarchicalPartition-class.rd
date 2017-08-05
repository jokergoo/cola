\name{HierarchicalPartition-class}
\docType{class}
\alias{HierarchicalPartition}
\alias{HierarchicalPartition-class}
\title{
The HierarchicalPartition class
}
\description{
The HierarchicalPartition class
}
\section{Methods}{
The \code{\link{HierarchicalPartition-class}} has following methods:

\describe{
  \item{\code{\link{hierarchical_partition}}:}{constructor method.}
  \item{\code{\link{collect_classes,HierarchicalPartition-method}}:}{plot the hierarchy of subgroups predicted.}
  \item{\code{\link{get_class,HierarchicalPartition-method}}:}{get the class IDs of subgroups.}
  \item{\code{\link{get_signatures,HierarchicalPartition-method}}:}{get the signatures for each subgroup.}
  \item{\code{\link{get_single_run,HierarchicalPartition-method}}:}{get a \code{\link{ConsensusPartition-class}} object at a specified hierarchy level.}
  \item{\code{\link{test_to_known_factors,HierarchicalPartition-method}}:}{test correlation between predicted subgrouping and known annotations, if available.}
}}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
