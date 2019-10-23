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

  \item{class}{A vector of class IDs.}
  \item{ref}{A vector of reference IDs.}
  \item{full_set}{The full set of ID levels. }
  \item{return_map}{Whether return the mapping or the adjusted labels.}

}
\details{
In partition, the exact value of the class ID is not of importance. E.g. for two partitions
\code{a, a, a, b, b, b, b} and \code{b, b, b, a, a, a, a}, they are the same partitions although the labels
of \code{a} and \code{b} are switched in the two partitions. Here \code{\link{relabel_class}} function switches the labels
in \code{class} vector according to the labels in \code{ref} vector to maximize \code{sum(class == ref)}.

Mathematically, this is called linear sum assignment problem and it is solved by \code{\link[clue]{solve_LSAP}}.
}
\value{
A named vector where names correspond to the IDs in \code{class} and values correspond to \code{ref},
which means \code{map = relabel_class(class, ref); map[class]} returns the relabelled IDs.

The returned object attaches a data frame with three columns:

\itemize{
  \item original IDs in \code{class}
  \item adjusted IDs according to \code{ref}
  \item reference IDs in \code{ref}
}

If \code{return_map} in the \code{\link{relabel_class}} is set to \code{\link{FALSE}}, the function simply returns
a vector of adjusted class IDs.

If the function returns the mapping vector (when \code{return_map = TRUE}), the mapping variable
is always character, which means, if your \code{class} and \code{ref} are numeric, you need to convert
them back to numeric explicitely. If \code{return_map = FALSE}, the returned relabelled vector has
the same mode as \code{class}.
}
\examples{
class = c(rep("a", 10), rep("b", 3))
ref = c(rep("b", 4), rep("a", 9))
relabel_class(class, ref)
relabel_class(class, ref, return_map = FALSE)
}
