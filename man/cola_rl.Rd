\name{cola_rl}
\docType{data}
\alias{cola_rl}
\title{
Example ConsensusPartitionList object
}
\description{
Example ConsensusPartitionList object
}
\usage{
data(cola_rl)
}
\details{
Following code was used to generate \code{cola_rl}:

  \preformatted{
  set.seed(123)
  m = cbind(rbind(matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
                  matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
                  matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
            rbind(matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
                  matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
                  matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
            rbind(matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
                  matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
                  matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20))
           ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
  cola_rl = run_all_consensus_partition_methods(data = m, cores = 6)  }
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
cola_rl
}
