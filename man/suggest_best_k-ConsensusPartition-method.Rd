\name{suggest_best_k-ConsensusPartition-method}
\alias{suggest_best_k,ConsensusPartition-method}
\title{
Suggest the best number of subgroups
}
\description{
Suggest the best number of subgroups
}
\usage{
\S4method{suggest_best_k}{ConsensusPartition}(object,
    jaccard_index_cutoff = select_jaccard_cutoff(ncol(object)),
    mean_silhouette_cutoff = NULL,
    stable_PAC = 0.1, help = TRUE)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{jaccard_index_cutoff}{The cutoff for Jaccard index for comparing to previous k.}
  \item{mean_silhouette_cutoff}{Cutoff for mean silhourtte scores.}
  \item{stable_PAC}{Cutoff for stable PAC. This argument only take effect when \code{mean_silhouette_cutoff} is set to \code{NULL}.}
  \item{help}{Whether to print help message.}

}
\details{
The best k is selected according to following rules:

\itemize{
  \item All k with Jaccard index larger than \code{jaccard_index_cutoff} are removed because increasing k does not provide enough extra information. If all k are removed, it is marked as no subgroup is detected. 
  \item If all k with Jaccard index larger than 0.75, k with the highest mean silhourtte score is taken as the best k.
  \item For all k with mean silhouette score larger than \code{mean_silhouette_cutoff}, the maximal k is taken as the best k, and other k are marked as optional best k. 
  \item If argument \code{mean_silhouette_cutoff} is set to NULL, which means we do not filter by mean silhouette scores while by 1-PAC scores. Similarly, k with the highest 1-PAC is taken the best k and other k are marked as optional best k.
  \item If it does not fit the second rule. The k with the maximal vote of the highest 1-PAC score, highest mean silhouette, and highest concordance is taken as the best k.
}

It should be noted that
it is difficult to find the best k deterministically, we encourage users to
compare results for all k and determine a proper one which best explain
their studies.
}
\section{See}{
The selection of the best k can be visualized by \code{\link{select_partition_number}}.}
\value{
The best k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
obj = golub_cola["ATC", "skmeans"]
suggest_best_k(obj)
}
