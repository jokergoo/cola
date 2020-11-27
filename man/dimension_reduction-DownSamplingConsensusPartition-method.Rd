\name{dimension_reduction-DownSamplingConsensusPartition-method}
\alias{dimension_reduction,DownSamplingConsensusPartition-method}
\title{
Visualize column after dimension reduction
}
\description{
Visualize samples (the matrix columns) after dimension reduction
}
\usage{
\S4method{dimension_reduction}{DownSamplingConsensusPartition}(object, k, top_n = NULL,
    method = c("PCA", "MDS", "t-SNE", "UMAP"),
    control = list(), color_by = NULL,
    internal = FALSE, nr = 5000,
    p_cutoff = 0.05, remove = FALSE,
    scale_rows = TRUE, verbose = TRUE, ...)
}
\arguments{

  \item{object}{A \code{\link{DownSamplingConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{top_n}{Top n rows to use. By default it uses all rows in the original matrix.}
  \item{method}{Which method to reduce the dimension of the data. \code{MDS} uses \code{\link[stats]{cmdscale}}, \code{PCA} uses \code{\link[stats]{prcomp}}. \code{t-SNE} uses \code{\link[Rtsne]{Rtsne}}. \code{UMAP} uses \code{\link[umap]{umap}}.}
  \item{color_by}{If annotation table is set, an annotation name can be set here.}
  \item{control}{A list of parameters for \code{\link[Rtsne]{Rtsne}} or \code{\link[umap]{umap}}.}
  \item{internal}{Internally used.}
  \item{nr}{If number of matrix rows is larger than this value, random \code{nr} rows are used.}
  \item{p_cutoff}{Cutoff of p-value of class label prediction. Data points with values higher than it will be mapped with cross symbols.}
  \item{remove}{Whether to remove columns which have high p-values than the cutoff.}
  \item{scale_rows}{Whether to perform scaling on matrix rows.}
  \item{verbose}{Whether print messages.}
  \item{...}{Other arguments.}

}
\details{
This function is basically very similar as \code{\link{dimension_reduction,ConsensusPartition-method}}.
}
\value{
No value is returned.
}
\examples{
data(golub_cola_ds)
dimension_reduction(golub_cola_ds, k = 2)
dimension_reduction(golub_cola_ds, k = 3)
}
