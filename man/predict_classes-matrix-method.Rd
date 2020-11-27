\name{predict_classes-matrix-method}
\alias{predict_classes,matrix-method}
\title{
Predict classes for new samples based on signature centroid matrix
}
\description{
Predict classes for new samples based on signature centroid matrix
}
\usage{
\S4method{predict_classes}{matrix}(object, mat, dist_method = c("euclidean", "correlation", "cosine"),
    nperm = 1000, p_cutoff = 0.05, plot = TRUE, col_fun = NULL,
    verbose = TRUE, prefix = "", mc.cores = 1)
}
\arguments{

  \item{object}{The signature centroid matrix. See the Details section.}
  \item{mat}{The new matrix where the classes are going to be predicted. The number of rows should be the same as the signature centroid matrix (also make sure the row orders are the same). Be careful that \code{mat} should be in the same scale as the centroid matrix.}
  \item{dist_method}{Distance method. Value should be "euclidean", "correlation" or "cosine".}
  \item{nperm}{Number of permutatinos. It is used when \code{dist_method} is set to "euclidean" or "cosine".}
  \item{p_cutoff}{Cutoff for the p-values for determining class assignment.}
  \item{plot}{Whether to draw the plot that visualizes the process of prediction.}
  \item{col_fun}{A color mapping function generated from \code{\link{colorRamp2}}. It is set to both heatmaps.}
  \item{verbose}{Whether to print messages.}
  \item{prefix}{Used internally.}
  \item{mc.cores}{Number of cores.}

}
\details{
The signature centroid matrix is a k-column matrix where each column is the centroid of samples 
in the corresponding class (k-group classification).

For each sample in the new matrix, the task is basically to test which signature centroid the 
current sample is the closest to. There are two methods: the Euclidean distance and the 
correlation (Spearman) distance.

For the Euclidean/cosine distance method, for the vector denoted as x which corresponds to sample i 
in the new matrix, to test which class should be assigned to sample i, the distance between 
sample i and all k signature centroids are calculated and denoted as d_1, d_2, ..., d_k. The class with the smallest distance is assigned to sample i.
The distances for k centroids are sorted increasingly, and we design a statistic named "difference ratio", denoted as r
and calculated as: (|d_(1) - d_(2)|)/mean(d), which is the difference between the smallest distance and the second
smallest distance, normalized by the mean distance. 
To test the statistical significance of r, we randomly permute rows of the signature centroid matrix and calculate r_rand. 
The random permutation is performed \code{n_perm} times and the p-value is calculated as the proportion of r_rand being
larger than r.

For the correlation method, the distance is calculated as the Spearman correlation between sample i and signature
centroid k. The label for the class with the maximal correlation value is assigned to sample i. The 
p-value is simply calculated by \code{\link[stats]{cor.test}} between sample i and centroid k.

If a sample is tested with a p-value higher than \code{p_cutoff}, the corresponding class label is set to \code{NA}.
}
\value{
A data frame with two columns: the class labels (the column names of the signature centroid matrix are treated as class labels) and the corresponding p-values.
}
\examples{
\donttest{
data(golub_cola)
res = golub_cola["ATC:skmeans"]
mat = get_matrix(res)
# note scaling should be applied here because the matrix was scaled in the cola analysis
mat2 = t(scale(t(mat)))

tb = get_signatures(res, k = 3, plot = FALSE)
sig_mat = tb[, grepl("scaled_mean", colnames(tb))]
sig_mat = as.matrix(sig_mat)
colnames(sig_mat) = paste0("class", seq_len(ncol(sig_mat)))
# this is how the signature centroid matrix looks like:
head(sig_mat)

mat2 = mat2[tb$which_row, , drop = FALSE]

# now we predict the class for `mat2` based on `sig_mat`
predict_classes(sig_mat, mat2)
}
}
