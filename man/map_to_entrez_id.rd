\name{map_to_entrez_id}
\alias{map_to_entrez_id}
\title{
Map to Entrez IDs
}
\description{
Map to Entrez IDs
}
\usage{
map_to_entrez_id(from, org_db = "org.Hs.eg.db")
}
\arguments{

  \item{from}{The input gene ID type. Valid values should be in, e.g. \code{columns(org.Hs.eg.db)}.}
  \item{org_db}{The annotation database.}

}
\details{
If there are multiple mappings from the input ID type to an unique Entrez ID, just one is randomly picked.
}
\value{
A named vectors where names are IDs with input ID type and values are the Entrez IDs.

The returned object normally is used in \code{\link{GO_enrichment}}.
}
\examples{
# There is no example
NULL

}
