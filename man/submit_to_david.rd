\name{submit_to_david}
\alias{submit_to_david}
\title{
Perform DAVID analysis
}
\description{
Perform DAVID analysis
}
\usage{
submit_to_david(genes, email,
    catalog = c("GOTERM_CC_FAT", "GOTERM_BP_FAT", "GOTERM_MF_FAT", "KEGG_PATHWAY"),
    idtype = "ENSEMBL_GENE_ID", species = "Homo sapiens")
}
\arguments{

  \item{genes}{a vector of gene identifiers.}
  \item{email}{the email that user registered on DAVID web service (\url{https://david.ncifcrf.gov/content.jsp?file=WS.html} ).}
  \item{catalog}{a vector of function catalogs. Valid values are in \code{cola:::DAVID_ALL_CATALOGS}.}
  \item{idtype}{id types for the input gene list. Valid values are in \code{cola:::DAVID_ALL_ID_TYPES}.}
  \item{species}{full species name if the id type is not uniquely mapped to one single species}

}
\details{
This function just performed a quick DAVID analysis. Users are encouraged to use other
gene set enrichment tools such as **clusterProfiler** (\url{http://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html} ),
}
\section{Fgsea}{
If you want to run this function multiple times, please set time intervals between runs.}
\value{
A data frame with functional enrichment results.
}
\seealso{
\url{https://david.ncifcrf.gov}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
