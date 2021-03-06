% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fungible.R
\encoding{UTF-8}
\name{.fungible_extrema}
\alias{.fungible_extrema}
\title{Locate extrema of fungible OLS regression weights}
\usage{
.fungible_extrema(
  theta,
  Rxx,
  rxy,
  Nstarts = 1000,
  MaxMin = c("min", "max"),
  silent = FALSE
)
}
\arguments{
\item{theta}{The value of the R-squared decrement used to generate a family of fungible coefficients.}

\item{Rxx}{An intercorrelation matrix among the predictor variables}

\item{rxy}{A vector of predictor–criterion correlations}

\item{Nstarts}{The maximum number}

\item{MaxMin}{Should the cosine between the OLS and alternative weights be maximized ("max") to find the maximally similar coefficients or minimized ("min") to find the maximally dissimilar coefficients?}

\item{silent}{Should current optimization values be printed to the console (\code{FALSE}) or suppressed (\code{TRUE})?}
}
\value{
A list containing the alternative weights and other fungible weights estimation parameters
}
\description{
Locate extrema of fungible OLS regression weights
}
\author{
Adapted from \code{fungible::fungibleExtrema()} by Niels Waller and Jeff Jones
}
\keyword{internal}
