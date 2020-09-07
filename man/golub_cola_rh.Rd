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

    library(golubEsets)  # from bioc
    data(Golub_Merge)
    m = exprs(Golub_Merge)
    colnames(m) = paste0("sample_", colnames(m))
    anno = pData(Golub_Merge)

    m[m <= 1] = NA
    m = log10(m)

    m = adjust_matrix(m)

    library(preprocessCore)  # from bioc
    cn = colnames(m)
    rn = rownames(m)
    m = normalize.quantiles(m)
    colnames(m) = cn
    rownames(m) = rn

    register_NMF()

    set.seed(123)
    golub_cola_rh = hierarchical_partition(
        m,
        PAC_cutoff = 0.3,
        mc.cores = 4, 
        anno = anno[, c("ALL.AML"), drop = FALSE],
        anno_col = c("ALL" = "red", "AML" = "blue")
    )  }
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_rh)
golub_cola_rh
}
