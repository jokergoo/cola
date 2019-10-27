\name{consensus_heatmap-ConsensusPartition-method}
\alias{consensus_heatmap,ConsensusPartition-method}
\alias{consensus_heatmap}
\title{
Heatmap for the consensus matrix
}
\description{
Heatmap for the consensus matrix
}
\usage{
\S4method{consensus_heatmap}{ConsensusPartition}(object, k, internal = FALSE,
    anno = get_anno(object), anno_col = get_anno_col(object),
    show_row_names = FALSE, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of partitions.}
  \item{internal}{Used internally.}
  \item{anno}{A data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{show_row_names}{Whether plot row names on the consensus heatmap (which are the column names in the original matrix)}
  \item{...}{other arguments}

}
\details{
For row i and column j in the consensus matrix, the value of corresponding x_ij
is the probability of sample i and sample j being in a same group from all partitions.

There are following heatmaps from left to right:

\itemize{
  \item probability of the sample to stay in the corresponding group
  \item silhouette scores which measure the distance of an item to the second closest subgroups.
  \item predicted classes.
  \item consensus matrix.
  \item more annotations if provided as \code{anno}
}

One thing that is very important to note is that since we already know the consensus classes from consensus
partition, in the heatmap, only rows or columns within the group is clustered.
}
\value{
No value is returned.
}
\seealso{
\code{\link{membership_heatmap,ConsensusPartition-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
consensus_heatmap(cola_rl["sd", "hclust"], k = 3)
}
