% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Match&Permutation_Functions.R
\name{greedymatch}
\alias{greedymatch}
\title{Greedy Match}
\usage{
greedymatch(Mx, My, Ux, Uy)
}
\arguments{
\item{Mx}{Subject Score for X with n x n.comp.X matrix}

\item{My}{Subject Score for Y with n x n.comp.Y matrix}

\item{Ux}{Matrix with n.comp x n, Mx = Lx^-1 \%*\% t Ux, Lx is the whitener matrix of dX.}

\item{Uy}{Matrix with n.comp x n, My = Ly^-1 \%*\% t Uy, Ly is the whitener matrix of dY.}
}
\value{
a list of matrices:
\describe{
       \item{\code{Mx}}{Columns of original Mx reordered from highest to lowest correlation with matched component in My}
       \item{\code{My}}{Columns of original My reordered from highest to lowest correlation with matched component in Mx}
       \item{\code{Ux}}{Permuted rows of original Ux corresponds to MapX}
       \item{\code{Uy}}{Permuted rows of original Uy corresponds to MapY}
       \item{\code{correlations}}{a vector of correlations for each pair of columns in permuted Mx and M}
       \item{\code{mapX}}{the sequence of the columns in original Mx.}
       \item{\code{mapY}}{the sequence of the columns in original MY.}
}
}
\description{
\code{Greedy Match} matches a column of Mx and My by minimizing chordal distance between vectors,
removes the matched columns and then finds the next pair.
This equivalent to maximizing absolute correlation for data in which each column has mean equal to zero.
Returns permuted columns of Mx and My. This function does not do any scaling or sign flipping.
For this matching to coincide with angle matching, the columns must have zero mean.
}
