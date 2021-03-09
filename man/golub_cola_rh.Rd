\name{golub_cola_rh}
\docType{data}
\alias{golub_cola_rh}
\title{
Example HierarchicalPartition object from Golub dataset
}
\description{
Example HierarchicalPartition object from Golub dataset
}
\usage{
data(golub_cola_rh)
}
\details{
Following code was used to generate \code{golub_cola_rh}:

  \preformatted{
    library(cola)
    data(golub_cola)
    m = get_matrix(golub_cola)
    set.seed(123)
    golub_cola_rh = hierarchical_partition(
        m, cores = 6, 
        anno = get_anno(golub_cola), 
        anno_col = get_anno_col(golub_cola)
    )  }
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_rh)
golub_cola_rh
}
