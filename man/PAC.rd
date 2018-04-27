\name{PAC}
\alias{PAC}
\title{
PAC score
}
\description{
PAC score
}
\usage{
PAC(consensus_mat, x1 = seq(0.1, 0.3, by = 0.02),
    x2 = seq(0.7, 0.9, by = 0.02), trim = 0.2)
}
\arguments{

  \item{consensus_mat}{a consensus matrix.}
  \item{x1}{lower bound to define "ambiguous clustering". The value can be a vector.}
  \item{x2}{upper bound to define "ambihuous clustering". The value can be a vector.}
  \item{trim}{percent of extreme values to trim if combinations of \code{x1} and \code{x2} are more than 10.}

}
\details{
This a variant of the orignial PAC (proportion of ambiguous clustering) method.

For each \code{x_1i} in \code{x1} and \code{x_2j} in \code{x2}, \code{PAC_k = F(x_2j) - F(x_1i)}
where \code{F(x)} is the ecdf of the consensus matrix (the lower triangle matrix without diagnals). 
The final PAC is the mean of all \code{PAC_k} by removing top \code{trim/2} percent and bottom \code{trim/2} percent of all values.
}
\value{
A single numeric score.
}
\section{See}{
See \url{https://www.nature.com/articles/srep06207} for explanation of PAC score.}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
PAC(get_consensus(cola_rl[1, 1], k = 2))
PAC(get_consensus(cola_rl[1, 1], k = 3))
PAC(get_consensus(cola_rl[1, 1], k = 4))
PAC(get_consensus(cola_rl[1, 1], k = 5))
PAC(get_consensus(cola_rl[1, 1], k = 6))
}
