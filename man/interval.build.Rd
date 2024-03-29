% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.interval.R
\name{interval.build}
\alias{interval.build}
\title{JOIN MULTIPLE INTERVALS}
\usage{
interval.build(yyy, B, tr)
}
\arguments{
\item{yyy}{column vector of B lower bounds and B upper bounds}

\item{B}{number of replications}

\item{tr}{truncation threshold for the algorithm}
}
\value{
A single interval.
}
\description{
The function build the multisplit interval for a given x0) in each point of the evaluation
grid in t, by combining multiple splits according to the Multi Split algorithm.
}
