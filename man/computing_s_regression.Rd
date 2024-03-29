% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s_regression.R
\name{computing_s_regression}
\alias{computing_s_regression}
\title{COMPUTING THE MODULATION FUNCTION S}
\usage{
computing_s_regression(vec_residual, type, alpha, tau, grid_size)
}
\arguments{
\item{vec_residual}{A vector of the residuals obtained via functional modeling.}

\item{type}{A string indicating the type of modulation function chosen.
The alternatives are "identity","st-dev","alpha-max".}

\item{alpha}{The value of the confidence interval.}

\item{tau}{A number between 0 and 1 used for the randomized version of the algorithm.}

\item{grid_size}{A vector containing the number of grid points in each dimension.}
}
\value{
It returns a the values of a modulation function in each dimension of the response.
}
\description{
It computes modulation functions which allows local scaling of the prediction bands.
}
\details{
More details can be found in the help of conformal.fun.split function.
}
\references{
The function structure is taken from "Conformal Prediction Bands
for Multivariate Functional Data" by Diquigiovanni, Fontana, Vantini (2021) and,
also, from "The Importance of Being a Band: Finite-Sample Exact Distribution-Free
Prediction Sets for Functional Data" by Diquigiovanni, Fontana, Vantini (2021).
}
