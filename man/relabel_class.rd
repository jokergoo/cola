\name{relabel_class}
\alias{relabel_class}
\title{
Relabel class IDs according to the reference ID
}
\description{
Relabel class IDs according to the reference ID
}
\usage{
relabel_class(class, ref, full_set = union(class, ref), return_map = TRUE)
}
\arguments{

  \item{class}{a vector of class IDs.}
  \item{ref}{a vector of reference IDs.}
  \item{full_set}{the full set of levels. }
  \item{return_map}{whether return the mapping or the adjusted labels.}

}
\details{
In partition, the exact value of the class ID is not of importance. E.g. for two partitions
\code{a, a, a, b, b, b, b} and \code{b, b, b, a, a, a, a}, they are the same partitions although the labels
of \code{a} and \code{b} are switched in the two partitions. Here \code{\link{relabel_class}} function switches the labels
in \code{class} vector accoring to the labels in \code{ref} vector to maximize \code{sum(class == ref)}.

Mathematically, this is called linear sum assignment problem and is solved by \code{\link[clue]{solve_LSAP}}.
}
\value{
A data frame with three columns:

\itemize{
  \item original IDs
  \item adjusted IDs
  \item reference IDs
}

The mapping between adjusted IDs and original IDs are stored as the \code{map} attribute of the data frame.
}
\examples{
class = c(rep("a", 10), rep("b", 3))
ref = c(rep("b", 4), rep("a", 9))
relabel_class(class, ref)
}
