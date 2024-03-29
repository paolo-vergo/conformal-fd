% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{convert2data}
\alias{convert2data}
\title{CONVERSION OF FUNCTIONAL DATA OBJECTS TO GRID EVALUATIONS OF FUNCTIONS}
\usage{
convert2data(t, y, x, x0)
}
\arguments{
\item{t}{Either NULL or a list of vectors (grid points in the various p dimensions)}

\item{y}{An object of class "fData","mfData" or "fd".}

\item{x}{The input variable, a list of n elements. Each element is composed by a list}

\item{x0}{The input variable, a list of n0 elements. Each element is composed by a list}
}
\value{
t A list of vectors (grid points in the various p dimensions)
y A grid of functional evaluated points (a list of list of vectors)
}
\description{
CONVERSION OF FUNCTIONAL DATA OBJECTS TO GRID EVALUATIONS OF FUNCTIONS
}
