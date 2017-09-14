\name{submit_to_david}
\alias{submit_to_david}
\title{
Doing DAVID analysis
}
\description{
Doing DAVID analysis
}
\usage{
submit_to_david(genes, email, catalog = c("GOTERM_CC_FAT", "GOTERM_BP_FAT", "GOTERM_MF_FAT", "KEGG_PATHWAY"),
    idtype = "ENSEMBL_GENE_ID", species = "Homo sapiens")
}
\arguments{

  \item{genes}{a vector of genes}
  \item{email}{the email that user registered on DAVID web service (\url{https://david.ncifcrf.gov/content.jsp?file=WS.html} )}
  \item{catalog}{a vector of functional catalogs. Valid values are in \code{cola:::DAVID_ALL_ID_TYPES}.}
  \item{idtype}{id types for the input gene list. Valid  vaules are in \code{cola:::DAVID_ALL_CATALOGS}.}
  \item{species}{full species name if the id type is not unique to the species}

}
\details{
If you want to run this function multiple times, please set time intervals.
}
\value{
A data frame with functional enrichment results
}
\examples{
# There is no example
NULL

}
