\name{register_NMF}
\alias{register_NMF}
\title{
Register NMF partition method
}
\description{
Register NMF partition method
}
\usage{
register_NMF()
}
\details{
It actually runs following code:

  \preformatted{
    register_partition_methods(
        NMF = function(mat, k, ...) \{
            fit = NNLM::nnmf(A = mat, k = k, verbose = FALSE, ...)
            apply(fit$H, 2, which.max)
        \}, scale_method = "rescale"
    )  }
}
\examples{
# There is no example
NULL

}
