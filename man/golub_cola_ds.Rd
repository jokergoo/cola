\name{golub_cola_ds}
\docType{data}
\alias{golub_cola_ds}
\title{
Example DownSamplingConsensusPartition object from Golub dataset
}
\description{
Example DownSamplingConsensusPartition object from Golub dataset
}
\usage{
data(golub_cola_ds)
}
\details{
Following code was used to generate \code{golub_cola_ds}:

  \preformatted{
    library(cola)
    data(golub_cola)
    m = get_matrix(golub_cola)
    set.seed(123)
    golub_cola_ds = consensus_partition_by_down_sampling(m, subset = 50,
        anno = get_anno(golub_cola), anno_col = get_anno_col(golub_cola),
        top_value_method = "SD", partition_method = "kmeans")  }
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_ds)
golub_cola_ds
}
