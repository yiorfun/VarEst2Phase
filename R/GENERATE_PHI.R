#' Calculate the values of the polynomial spline at a given vector
#'
#' @param x The predictor variable (or vector). Missing values are allowed.
#' @param knot The internal breakpoints (i.e., the breakpoints without two ending points) that define the spline. The default is `NULL`, which results in a basis for the ordinary polynomial regression. Typical values are the mean or median for one knot, quantiles for more knots. See also `Boundary.knots`.
#' @param degree A non-negative integer degree of the piece-wise polynomial. The default value is 3 for cubic splines. Zero degree is allowed for this function, which is the only difference compared with `bs` function in package `splines`. Note that the order of a polynomial equals to the degree plus one.
#' @param boundary A 2-dimensional boundary points at which to anchor the B-spline basis. By default, they are the rangle of the non-`NA` data. If both `knots` and `Boundary.knots` are supplied, the basis parameters do not depend on `x`. Data can extend beyound `Boundary.knots`.
#' @param coeff_spline The coefficients variables that are evaluated corresponding to the B-spline basis matrix included an intercept. The length of `coeff_spline` equals `df` or `degre + length(knots) + 1` (due to the presence of an intercept).
#'
#' @return A matrix of dimension `length(x)` by `df = degree + length(knots)` (plus one if the intercept is included). Attributes that correspond to the arguments specificed are returned for useage of other functions in this package.
#' @export
#'
#' @examples
#' library(VarEst2Phase)
#' GENERATE_PHI(x = c(1, 2, 3), knot = c(0.25, 0.75), degree = 2, boundary = c(0, 1), coeff_spline = c(0.1, 0.2, 0.3))
GENERATE_PHI <- function(x, knot, degree, boundary, coeff_spline) {
  phi_matrix <- bSpline(
    x = x,
    knots = knot,
    degree = degree,
    Boundary.knots = boundary,
    intercept = TRUE
  )
  phi_vector <- drop(phi_matrix %*% coeff_spline)
  return(list(
    phi_vector = phi_vector,
    phi_matrix = phi_matrix
  ))
}
