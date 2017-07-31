\name{get_signatures-ConsensusPartition-method}
\alias{get_signatures,ConsensusPartition-method}
\alias{get_signatures}
\title{
Get signature rows
}
\description{
Get signature rows
}
\usage{
\S4method{get_signatures}{ConsensusPartition}(object, k,
    silhouette_cutoff = 0.5,
    fdr_cutoff = ifelse(row_diff_by == "samr", 0.1, 0.05),
    scale_rows = object@scale_rows,
    row_diff_by = c("compare_to_highest_subgroup", "Ftest", "samr"),
    anno = object@known_anno,
    anno_col = if(missing(anno)) object@known_col else NULL,
    show_legend = TRUE,
    show_column_names = TRUE, use_raster = TRUE,
    plot = TRUE, mat_other = NULL,
    ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object. The object can be returned from \code{\link{get_single_run}}.}
  \item{k}{number of partitions}
  \item{silhouette_cutoff}{cutoff for silhouette values. Columns with values  less than it are not used for finding signature rows. For selecting a  proper silhouette value, please refer to \url{https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.}}
  \item{fdr_cutoff}{cutoff for FDR of the difference between subgroups.}
  \item{scale_rows}{whether apply row scaling when making the heatmap.}
  \item{row_diff_by}{methods to get rows which are significantly different between subgroups.}
  \item{anno}{a data frame with known annotation of columns.}
  \item{anno_col}{a list of colors for the annotations in \code{anno}.}
  \item{show_legend}{whether draw the legends on the heatmap.}
  \item{show_column_names}{whether show column names on the heatmap.}
  \item{use_raster}{internally used}
  \item{mat_other}{other matrix you want to attach to the heatmap list. The matrix should have row names so that rows can be subsetted and matched to the main heatmap}
  \item{plot}{whether to make the plot}
  \item{...}{other arguments}

}
\details{
Basically the function apply test for the difference of subgroups in every
row. Also, to call it a signature for a given subgroup, the values in the
corresponding subgroup should have the highest mean value compared to all
other subgroups. The minimal p-value compared to all other subgroups is taken
as the p-value of the row and used for FDR calculation.
}
\value{
A list of three elements:

\describe{
  \item{\code{mat}}{the matrix for the signatures}
  \item{\code{fdr}}{FDR for rows}
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
