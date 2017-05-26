\name{test_between_factors}
\alias{test_between_factors}
\title{
Test whether known factors are correlated
}
\description{
Test whether known factors are correlated
}
\usage{
test_between_factors(x, y = NULL, all_factors = FALSE, verbose = TRUE)
}
\arguments{

  \item{x}{a data frame or a vector which contains discrete or continuous variables.}
  \item{y}{a data frame or a vector which contains discrete or continuous variables.}
  \item{all_factors}{are all columns in \code{x} and \code{y} are enforced to be factors.}
  \item{verbose}{whether print messages.}

}
\details{
Pairwise test is applied to every two columns in the data frame. Methods are:

\itemize{
  \item two numeric variables: correlation test by \code{\link[stats]{cor.test}} is applied;
  \item two character or factor variables: Chi-squared test by \code{\link[stats]{fisher.test}} is applied;
  \item one numeric variable and one character/factor variable: oneway ANOVA test by \code{\link[stats]{oneway.test}} is applied.
}

This function can be used to test the correlation between the predicted classes and other known factors.
}
\value{
A matrix of p-values.
}
\examples{
df = data.frame(v1 = rnorm(100), v2 = sample(letters[1:3], 100, replace = TRUE), 
    v3 = sample(LETTERS[5:6], 100, replace = TRUE))
test_between_factors(df)
x = runif(100)
test_between_factors(x, df)
}
