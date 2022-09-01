#' Estimate the variances at phase two for the regression coefficients by
#' using a bootstrap procedure
#'
#' @param brepsMax The number of bootstrap random samples.
#'
#' @param phase1_size A vector consisting of the sample sizes of all
#'     strata at the first phase.
#'
#' @param phase2_size A vector consisting of the sample sizes of all
#'     strata at the second phase.
#'
#' @param member_phase2 A vector consisting of the observations
#' which are selected to form those at the second phase.
#'
#' @param DELTA A matrix that each row indicates the interval-censoring
#'     type of the observation with the sample size by 3. For example,
#'     `(1, 0, 0)` indicates that `UVec` is observed but `VVec` reaches the
#'     right boundary, `(0, 1, 0)` indicates both `UVec` and `VVec` are
#'     observed and `(0, 0, 1)` indicates `UVec` reaches the left boundary
#'     but `VVec` is observed for the observation.
#'
#' @param OMEGA A vector that consisting of the weights for each observation
#'     after two-phase sampling. If the two-phase sampling is not used,
#'     the values of `OMEGA` are all 1's.
#'
#' @param Gama_initial An initial value of `Gama`.
#'
#' @param d The number of unknown coefficients of the regression parameters.
#'
#' @param qn The number of unknown coefficients of the B-spline basis,
#'     that is, `df` or `degree + length(knot) + 1` (because the intercept
#'     is included).
#'
#' @param ZMat ZMat The design matrix corresponding to the regression
#'     parameters with the sample size by `d`.
#'
#' @param UVec A vector consisting of the left observation visit given
#'     in the interval-censoring design.
#'
#' @param VVec A vector consisting of the right observation visit given
#'     in the interval-censoring design.
#'
#' @param knot The internal breakpoints (i.e., the breakpoints without two
#'     ending points) that define the spline. The default is `NULL`, which
#'     results in a basis for the ordinary polynomial regression. Typical
#'     values are the mean or median for one knot, quantiles for more knots.
#'     See also `Boundary.knots`.
#'
#' @param degree A non-negative integer degree of the piece-wise polynomial.
#'     The default value is 3 for cubic splines. Zero degree is allowed for
#'     this function, which is the only difference compared with `bs`
#'     function in package `splines`. Note that the order of a polynomial
#'     equals to the degree plus one.
#'
#' @param boundary A 2-dimensional boundary points at which to anchor the
#'     B-spline basis. By default, they are the range of the non-`NA` data.
#'     If both `knot` and `boundary` are supplied, the basis parameters do
#'     not depend on `x`. Data can extend beyond `boundary`.
#'
#' @param GENERATE_GRADIENT_HESSIAN A function numerically calculates the
#'    gradient vector and the Hessian matrix of the given weighted
#'    log-likelihood function under two-phase sampling design.
#'
#' @param GENERATE_PHI The function generates the B-spline basis matrix and
#'     the values of a linear combination of B-spline basis function at a
#'     given vector.
#'
#' @return A `brepsMax` by `length(Gama)` matrix has the bootstrap result.
#'
#' @export
#'
#' @references \url{https://onlinelibrary.wiley.com/doi/10.1111/sjos.12152}
#'
#' @importFrom splines2 bSpline
#' @importFrom extraDistr rmvhyper
BOOTSTRAP_PROCEDURE <- function(brepsMax, phase1_size, phase2_size,
    member_phase2, DELTA, OMEGA, Gama_initial, d, qn, ZMat, UVec, VVec,
    knot, degree, boundary, GENERATE_GRADIENT_HESSIAN, GENERATE_PHI) {

	breps <- 1
	ResMat_B <- matrix(0, nrow = brepsMax, ncol = d + qn)
	multiple_K <- phase1_size %/% phase2_size
	remainder_R <- phase1_size %% phase2_size
	selection_S <- (1 - remainder_R / phase2_size) *
	    (1 - remainder_R / (phase1_size - 1))

	while(breps <= brepsMax){
		OMEGA_B <- rep(0, times = nrow(DELTA))
		for(stra in seq(length(multiple_K))){
			if(runif(1) <= selection_S[stra]){
				boot_vec <- drop(rmvhyper(nn = 1,
										   n = rep(multiple_K[stra], times = phase2_size[stra]),
									       k = phase2_size[stra]))
				### it records how many multiple of the elements are
				###	sampled in bootstrap
				OMEGA_B[member_phase2[[stra]]] <- OMEGA[member_phase2[[stra]]] * boot_vec
			} else{
				boot_vec <- drop(rmvhyper(nn = 1,
										   n = rep(multiple_K[stra] + 1, times = phase2_size[stra]),
									       k = phase2_size[stra]))
				### it records how many multiple of the elements are
				### sampled in bootstrap
				OMEGA_B[member_phase2[[stra]]] <- OMEGA[member_phase2[[stra]]] * boot_vec
			}
		}

		RESULT_B <- GR(gama.ini = Gama_initial,
							  d = d,
							 qn = qn,
						   ZMat = ZMat,
						   UVec = UVec,
					       VVec = VVec,
				           knot = knot,
				         degree = degree,
				       boundary = boundary,
					      DELTA = DELTA,
					      OMEGA = OMEGA_B,
	  GENERATE_GRADIENT_HESSIAN = GENERATE_GRADIENT_HESSIAN,
			       GENERATE_PHI = GENERATE_PHI)

		ResMat_B[breps, ] <- RESULT_B$gama
		breps <- breps + 1
	}
	return(ResMat_B)
}

