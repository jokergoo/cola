\name{get_signatures-ConsensusPartition-method}
\alias{get_signatures,ConsensusPartition-method}
\title{
Get signature rows
}
\description{
Get signature rows
}
\usage{
\S4method{get_signatures}{ConsensusPartition}(object, k,
    silhouette_cutoff = 0.5,
    fdr_cutoff = ifelse(diff_method == "samr", 0.05, 0.1),
    scale_rows = object@scale_rows,
    diff_method = c("ttest", "Ftest", "samr", "pamr"),
    anno = get_anno(object),
    anno_col = get_anno_col(object),
    internal = FALSE,
    show_column_names = TRUE, use_raster = TRUE,
    plot = TRUE, verbose = TRUE,
    top_k_genes = 5000,
    ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object. The object can be returned from \code{\link{get_single_run}} or directly from \code{\link{consensus_partition}}.}
  \item{k}{number of partitions}
  \item{silhouette_cutoff}{cutoff for silhouette values. Columns with values  less than it are not used for finding signature rows. For selecting a  proper silhouette value, please refer to \url{https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.}}
  \item{fdr_cutoff}{cutoff for FDR of the difference between subgroups.}
  \item{scale_rows}{whether apply row scaling when making the heatmap.}
  \item{diff_method}{methods to get rows which are significantly different between subgroups, see 'Details' section.}
  \item{anno}{a data frame with known annotation of samples.}
  \item{anno_col}{a list of colors for the annotations in \code{anno}.}
  \item{show_legend}{whether draw the legends on the heatmap.}
  \item{show_column_names}{whether show column names on the heatmap.}
  \item{use_raster}{internally used}
  \item{mat_other}{other matrix you want to attach to the heatmap list. The matrix should have row names so that rows can be subsetted and matched to the main heatmap}
  \item{plot}{whether to make the plot}
  \item{verbose}{whether to print messages}
  \item{...}{other arguments}

}
\details{
Basically the function applies test for the difference of subgroups for every
row. There are three methods which test significance of the difference:

\describe{
  \item{ttest}{it first extracts the subgroup with higest value, then use t-test to test to  all the other subgroups. }
  \item{samr}{use SAM method to find significantly different rows between subgroups}
  \item{Ftest}{use F-test to find significantly different rows between subgroups}
}

Also, to call it a signature for a given subgroup, the values in the
corresponding subgroup should have the highest mean value compared to all
other subgroups. The minimal p-value compared to all other subgroups is taken
as the p-value of the row and used for FDR calculation.
}
\value{
A list of three elements:

\describe{
  \item{\code{mat}}{the matrix for the signatures}
  \item{\code{group}}{subgroups that the rows are significantly high}
}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
