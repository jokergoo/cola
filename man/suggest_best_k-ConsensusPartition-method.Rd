\name{suggest_best_k-ConsensusPartition-method}
\alias{suggest_best_k,ConsensusPartition-method}
\title{
Suggest the best number of subgroups
}
\description{
Suggest the best number of subgroups
}
\usage{
\S4method{suggest_best_k}{ConsensusPartition}(object, jaccard_index_cutoff = 0.95)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{jaccard_index_cutoff}{The cutoff for Jaccard index for comparing to previous k.}

}
\details{
The best k is selected according to following rules:

\itemize{
  \item All k with Jaccard index larger than 0.95 are removed because increasing k does not provide enough extra information. If all k are removed, it is marked as no subgroup is detected. 
  \item For all k with 1-PAC score larger than 0.9, the maximal k is taken as the best k, and other k are marked as optional k. 
  \item If it does not fit the second rule. The k with the maximal vote of the highest 1-PAC score, highest mean silhouette, and highest concordance is taken as the best k.
}

Additionally, if 1-PAC for the best k is larger than 0.9 (10% ambiguity for
the partition), cola marks it as a stable partition. It should be noted that
it is difficult to find the best k deterministically, we encourage users to
compare results for all k and determine a proper one which best explain
their studies.
}
\value{
The best k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["SD", "kmeans"]
suggest_best_k(obj)
}
