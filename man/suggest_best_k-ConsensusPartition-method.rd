\name{suggest_best_k-ConsensusPartition-method}
\alias{suggest_best_k,ConsensusPartition-method}
\title{
Suggest the best number of partitions
}
\description{
Suggest the best number of partitions
}
\usage{
\S4method{suggest_best_k}{ConsensusPartition}(object, rand_index_cutoff = 0.95)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{rand_index_cutoff}{The cutoff for Rand index compared to previous k.}

}
\details{
The best k is selected according to following rules:

1. k with rand index larger than \code{rand_index_cutoff} are removed. If all k are removed, the best k is defined as \code{NA}.
2. If there are some k having \code{1-PAC} larger than 0.9, the largest k is selected as the best k.
3. If it does not fit rule 2, the k with highest vote of highest 1-PAC, mean_silhouette and concordance scores is
   selected as the best k.

\code{\link{suggest_best_k}} function only gives suggestion of selecting
a reasonable best k. Users still need to look at the plots (e.g. by \code{\link{select_partition_number}} or \code{\link{consensus_heatmap}} functions), or even
by checking whether the subgrouping gives a reasonable signatures by \code{\link{get_signatures}}, to pick a reasonable k that best explains their study.

The best k with 1-PAC larger than 0.9 is treated as a stable partition.
}
\value{
The best k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["sd", "kmeans"]
suggest_best_k(obj)
}
