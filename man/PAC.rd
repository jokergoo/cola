\name{PAC}
\alias{PAC}
\title{
The proportion of ambiguous clustering (PAC score)
}
\description{
The proportion of ambiguous clustering (PAC score)
}
\usage{
PAC(consensus_mat, x1 = 0.1, x2 = 0.9)
}
\arguments{

  \item{consensus_mat}{A consensus matrix.}
  \item{x1}{Lower bound to define "ambiguous clustering".}
  \item{x2}{Upper bound to define "ambihuous clustering".}

}
\details{
The PAC score is defined as F(x2) - F(x1) where F(x) is the CDF of the consensus matrix.
}
\value{
A single numeric vaule.
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
