\name{gene_cooccurrence}
\alias{gene_cooccurrence}
\title{
Co-occurrence of genes in gene sets
}
\description{
Co-occurrence of genes in gene sets
}
\usage{
gene_cooccurrence(x, map = NULL, genesets, min_count = 50, max_count = 5000)
}
\arguments{

  \item{x}{the object returned from \code{\link{get_signatures}}}
  \item{map}{mapping between rows of \code{x$mat} and genes in \code{genesets}}
  \item{genesets}{a list}
  \item{min_count}{minimal number of genes in genesets}
  \item{max_count}{maximal number of genes in genesets}

}
\details{
For genes in each row group, the co-occurence of every gene pair to be in a same gene set
is calculated. The mean co-occurence of all genes is used as the final statistic which can
be understanded as the mean number of gene sets that two genes co-exist.
}
\examples{
# There is no example
NULL

}
