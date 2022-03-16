\name{predict_classes-ConsensusPartition-method}
\alias{predict_classes,ConsensusPartition-method}
\title{
Predict classes for new samples based on cola classification
}
\description{
Predict classes for new samples based on cola classification
}
\usage{
\S4method{predict_classes}{ConsensusPartition}(object, k, mat,
    silhouette_cutoff = 0.5,
    fdr_cutoff = cola_opt$fdr_cutoff,
    group_diff = cola_opt$group_diff,
    scale_rows = object@scale_rows,
    diff_method = "Ftest",
    method = "centroid",
    dist_method = c("euclidean", "correlation", "cosine"), nperm = 1000,
    p_cutoff = 0.05, plot = TRUE, col_fun = NULL,
    split_by_sigatures = FALSE, force = FALSE,
    verbose = TRUE, help = TRUE, prefix = "",
    mc.cores = 1, cores = mc.cores)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups to get the classifications.}
  \item{mat}{The new matrix where the sample classes are going to be predicted. The number of rows should be the same as the original matrix for cola analysis (also make sure the row orders are the same). Be careful that the scaling of \code{mat} should be the same as that applied in cola analysis.}
  \item{silhouette_cutoff}{Send to \code{\link{get_signatures,ConsensusPartition-method}} for determining signatures.}
  \item{fdr_cutoff}{Send to \code{\link{get_signatures,ConsensusPartition-method}} for determining signatures.}
  \item{group_diff}{Send to \code{\link{get_signatures,ConsensusPartition-method}} for determining signatures.}
  \item{scale_rows}{Send to \code{\link{get_signatures,ConsensusPartition-method}} for determining signatures.}
  \item{diff_method}{Send to \code{\link{get_signatures,ConsensusPartition-method}} for determining signatures.}
  \item{method}{Method for predicting class labels. Possible values are "centroid", "svm" and "randomForest".}
  \item{dist_method}{Distance method. Value should be "euclidean", "correlation" or "cosine". Send to \code{\link{predict_classes,matrix-method}}.}
  \item{nperm}{Number of permutatinos. It is used when \code{dist_method} is set to "euclidean" or "cosine". Send to \code{\link{predict_classes,matrix-method}}.}
  \item{p_cutoff}{Cutoff for the p-values for determining class assignment. Send to \code{\link{predict_classes,matrix-method}}.}
  \item{plot}{Whether to draw the plot that visualizes the process of prediction. Send to \code{\link{predict_classes,matrix-method}}.}
  \item{col_fun}{A color mapping function generated from \code{\link{colorRamp2}}. It is set to both heatmaps.}
  \item{split_by_sigatures}{Should the heatmaps be split based on k-means on the main heatmap, or on the patterns of the signature heatmap.}
  \item{force}{If the value is \code{TRUE} and when \code{\link{get_signatures,ConsensusPartition-method}} internally failed, top 1000 rows with the highest between-group mean difference are used for constructing the signature centroid matrix. It is basically used internally.}
  \item{verbose}{Whether to print messages. Send to \code{\link{predict_classes,matrix-method}}.}
  \item{help}{Whether to print help messages.}
  \item{prefix}{Used internally.}
  \item{mc.cores}{Number of cores. This argument will be removed in future versions.}
  \item{cores}{Number of cores, or a \code{cluster} object returned by \code{\link[parallel]{makeCluster}}.}

}
\details{
The prediction is based on the signature centroid matrix from cola classification. The processes are as follows:

1. For the provided \code{\link{ConsensusPartition-class}} object and a selected k, the signatures that discriminate classes
   are extracted by \code{\link{get_signatures,ConsensusPartition-method}}. If number of signatures is more than 2000, only 2000 signatures are randomly sampled.
2. The signature centroid matrix is a k-column matrix where each column is the centroid of samples in the corresponding
   class, i.e. the mean across samples. If rows were scaled in cola analysis, the signature centroid matrix is the mean of scaled
   values and vise versa. Please note the samples with silhouette score less than \code{silhouette_cutoff} are removed
   for calculating the centroids.
3. With the signature centroid matrix and the new matrix, it calls \code{\link{predict_classes,matrix-method}} to perform the prediction.
   Please see more details of the prediction on that help page. Please note, the scales of the new matrix should be the same as the matrix
   used for cola analysis.
}
\value{
A data frame with two columns: the class labels (in numeric) and the corresponding p-values.
}
\seealso{
\code{\link{predict_classes,matrix-method}} that predicts the classes for new samples.
}
\examples{
\donttest{
data(golub_cola)
res = golub_cola["ATC:skmeans"]
mat = get_matrix(res)
# note scaling should be applied here because the matrix was scaled in the cola analysis
mat2 = t(scale(t(mat)))
cl = predict_classes(res, k = 3, mat2)
# compare the real classification and the predicted classification
data.frame(cola_class = get_classes(res, k = 3)[, "class"],
           predicted = cl[, "class"])
# change to correlation method
cl = predict_classes(res, k = 3, mat2, dist_method = "correlation")
# compare the real classification and the predicted classification
data.frame(cola_class = get_classes(res, k = 3)[, "class"],
           predicted = cl[, "class"]) 
}
}
