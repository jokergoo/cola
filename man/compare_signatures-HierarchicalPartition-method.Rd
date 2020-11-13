\name{compare_signatures-HierarchicalPartition-method}
\alias{compare_signatures,HierarchicalPartition-method}
\title{
Compare Signatures from Different Nodes
}
\description{
Compare Signatures from Different Nodes
}
\usage{
\S4method{compare_signatures}{HierarchicalPartition}(object, filter_node = filter_node_param(),
    method = c("euler", "upset"), upset_max_comb_sets = 20,
    verbose = interactive(), ...)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object. }
  \item{depth}{Depth of the hierarchy.}
  \item{min_n_signatures}{Minimal number of signatures on the node.}
  \item{method}{Method to visualize.}
  \item{upset_max_comb_sets}{Maximal number of combination sets to show.}
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
