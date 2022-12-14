% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.R
\name{create.graph.long}
\alias{create.graph.long}
\title{create graph dataset with netmat and mmp_order
a data.frame called with vectorization of reordered netmat by mmp_order.}
\usage{
create.graph.long(gmatrix, sort_indices = NULL)
}
\arguments{
\item{gmatrix}{netmat}

\item{sort_indices}{mmp_order}
}
\value{
a data.frame with vectors:
 ## X1: vector of numerics.
 ## X2: vector of numerics.
 ## value: vectorization of reordered netmat by mmp_order.
}
\description{
create graph dataset with netmat and mmp_order
a data.frame called with vectorization of reordered netmat by mmp_order.
}
