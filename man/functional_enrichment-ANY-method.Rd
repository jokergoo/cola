\name{functional_enrichment-ANY-method}
\alias{functional_enrichment,ANY-method}
\title{
Perform functional enrichment on signature genes
}
\description{
Perform functional enrichment on signature genes
}
\usage{
\S4method{functional_enrichment}{ANY}(object,
    id_mapping = guess_id_mapping(object, org_db, verbose),
    org_db = "org.Hs.eg.db", ontology = "BP",
    min_set_size = 10, max_set_size = 1000,
    verbose = TRUE, prefix = "", ...)
}
\arguments{

  \item{object}{A vector of gene IDs.}
  \item{id_mapping}{If the gene IDs which are row names of the original matrix are not Entrez IDs, a named vector should be provided where the names are the gene IDs in the matrix and values are correspoinding Entrez IDs. The value can also be a function that converts gene IDs.}
  \item{org_db}{Annotation database.}
  \item{ontology}{Following ontologies are allowed: \code{BP}, \code{CC}, \code{MF}, \code{KEGG}, \code{Reactome}. \code{MSigDb} with the gmt file set by \code{gmt_file} argument, or \code{gmt} for general gmt gene sets.}
  \item{min_set_size}{The minimal size of the gene sets.}
  \item{max_set_size}{The maximal size of the gene sets.}
  \item{verbose}{Whether to print messages.}
  \item{prefix}{Used internally.}
  \item{...}{Pass to \code{\link[clusterProfiler]{enrichGO}}, \code{\link[clusterProfiler]{enrichKEGG}}, \code{\link[clusterProfiler]{enricher}}, \code{\link[DOSE]{enrichDO}} or \code{\link[ReactomePA]{enrichPathway}}.}

}
\details{
The function enrichment is applied by clusterProfiler, DOSE or ReactomePA packages.
}
\value{
A data frame.
}
\examples{
# There is no example
NULL

}
