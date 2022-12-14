% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Match&Permutation_Functions.R
\name{permmatRank_joint}
\alias{permmatRank_joint}
\title{Permutation test to get joint components ranks}
\usage{
permmatRank_joint(matchedResults, nperms = 100)
}
\arguments{
\item{matchedResults}{results generated by angleMatchICA}

\item{nperms}{the number of permutation}
}
\value{
a list of matrixes
 ## pvalues: pvalues for the matched colunmns don't have correlation.
 ## corrperm: correlation value for original Mx with each random permutation of My.
 ## corrmatched: the correlation for each pair of matched columns.
}
\description{
Permutation test to get joint components ranks
}
