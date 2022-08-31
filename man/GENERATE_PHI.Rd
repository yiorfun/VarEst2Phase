% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GENERATE_PHI.R
\name{GENERATE_PHI}
\alias{GENERATE_PHI}
\title{Calculate the values of the polynomial spline at a given vector}
\usage{
GENERATE_PHI(x, knot, degree, boundary, coeff_spline)
}
\arguments{
\item{x}{The predictor variable (or vector). Missing values are allowed.}

\item{knot}{The internal breakpoints (i.e., the breakpoints without two ending points) that define the spline. The default is \code{NULL}, which results in a basis for the ordinary polynomial regression. Typical values are the mean or median for one knot, quantiles for more knots. See also \code{Boundary.knots}.}

\item{degree}{A non-negative integer degree of the piece-wise polynomial. The default value is 3 for cubic splines. Zero degree is allowed for this function, which is the only difference compared with \code{bs} function in package \code{splines}. Note that the order of a polynomial equals to the degree plus one.}

\item{boundary}{A 2-dimensional boundary points at which to anchor the B-spline basis. By default, they are the rangle of the non-\code{NA} data. If both \code{knots} and \code{Boundary.knots} are supplied, the basis parameters do not depend on \code{x}. Data can extend beyound \code{Boundary.knots}.}

\item{coeff_spline}{The coefficients variables that are evaluated corresponding to the B-spline basis matrix included an intercept. The length of \code{coeff_spline} equals \code{df} or \code{degre + length(knots) + 1} (due to the presence of an intercept).}
}
\value{
A matrix of dimension \code{length(x)} by \code{df = degree + length(knots)} (plus one if the intercept is included). Attributes that correspond to the arguments specificed are returned for useage of other functions in this package.
}
\description{
Calculate the values of the polynomial spline at a given vector
}
\examples{
library(VarEst2Phase)
GENERATE_PHI(x = c(1, 2, 3), knot = c(0.25, 0.75), degree = 2, boundary = c(0, 1), coeff_spline = c(0.1, 0.2, 0.3))
}