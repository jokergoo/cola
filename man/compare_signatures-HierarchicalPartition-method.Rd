\name{compare_signatures-HierarchicalPartition-method}
\alias{compare_signatures,HierarchicalPartition-method}
\title{
Compare Signatures from Different Nodes
}
\description{
Compare Signatures from Different Nodes
}
\usage{
\S4method{compare_signatures}{HierarchicalPartition}(object, depth = max_depth(object),
    method = c("euler", "upset"), verbose = interactive(), ...)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object. }
  \item{depth}{Depth of the hierarchy.}
  \item{method}{Method to visualize.}
  \item{verbose}{Whether to print message.}
  \item{...}{Other arguments passed to \code{\link{get_signatures,HierarchicalPartition-method}}.}

}
\details{
It plots an Euler diagram or a UpSet plot showing the overlap of signatures from different nodes.
On each node, the number of subgroups is inferred by \code{\link{suggest_best_k,ConsensusPartition-method}}.
}
\examples{
data(golub_cola_rh)
compare_signatures(golub_cola_rh)
}
