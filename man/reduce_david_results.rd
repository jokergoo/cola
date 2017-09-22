\name{reduce_david_results}
\alias{reduce_david_results}
\title{
Reduce DAVID results
}
\description{
Reduce DAVID results
}
\usage{
reduce_david_results(tb, fdr_cutoff = 0.05, hit_cutoff = 5)
}
\arguments{

  \item{tb}{object from \code{\link{submit_to_david}}}
  \item{fdr_cutoff}{cutoff for fdr}
  \item{hit_cutoff}{cutoff for number of genes in a term}

}
\details{
For each cluster and each functional category, the function picks the most
significant term.
}
\value{
A subset of rows in \code{tb}.
}
\examples{
# There is no example
NULL
}
