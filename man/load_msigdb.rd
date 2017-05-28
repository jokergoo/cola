\name{load_msigdb}
\alias{load_msigdb}
\title{
Load from MSigDB
}
\description{
Load from MSigDB
}
\usage{
load_msigdb(f)
}
\arguments{

  \item{f}{path of the xml file.}

}
\details{
The xml file can be downloaded from \url{http://software.broadinstitute.org/gsea/downloads.jsp} .
}
\value{
A \code{msigdb} class object with two elements:

\describe{
  \item{\code{meta}}{a data frame with geneset ids, organisms, description of the gene sets and categories.}
  \item{\code{list}}{a list of gene sets of gene symbols.}
}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
