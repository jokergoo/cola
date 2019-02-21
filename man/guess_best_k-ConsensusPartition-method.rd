\name{guess_best_k-ConsensusPartition-method}
\alias{guess_best_k,ConsensusPartition-method}
\title{
Guess the best number of partitions
}
\description{
Guess the best number of partitions
}
\usage{
\S4method{guess_best_k}{ConsensusPartition}(object, rand_index_cutoff = 0.95)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{rand_index_cutoff}{the Rand index compared to previous k is larger than this value, it is filtered out.}

}
\details{
The best k is voted from 1) the k with the maximal cophcor value, 2) the k with the minimal PAC value,
3) the k with the maximal mean silhouette value and 4) the k with the maximal concordance value.

There are scenarios that a better partition with k groups than k - 1 groups (e.g. for the sense of better sihouette score) 
is only because of one tiny group of samples are separated and it is better to still put them back to the original group
to improve the robustness of the subgrouping. For this, users can set the cutoff of Rand index by \code{rand_index_cutoff} to
get rid of or reduce the effect of such cirsumstances.

Honestly, it is hard or maybe impossible to say which k is the best one. \code{\link{guess_best_k}} function only gives suggestion of selecting
a reasonable k. Users still need to look at the plots (e.g. by \code{\link{select_partition_number}} or \code{\link{consensus_heatmap}} functions), or even
by checking whether the subgrouping gives a reasonable signatures by \code{\link{get_signatures}}, to pick a reasonable k that best explains thieir study.
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
guess_best_k(obj)
}
