#' Estimate both the coefficients of the regression parameters and the
#' coefficients of the B-spline basis functions and estimate the variance
#' of the estimated coefficents of the regression parameters by using
#' Generalized Rosen (GR) algorithm
#'
#' @param gama.ini An initial value of `Gama`.
#'
#' @param active.set A vector includes an active set corresponding to
#'     `gama.ini`.
#'
#' @param d The number of unknown coefficients of the regression parameters.
#'
#' @param qn The number of unknown coefficients of the B-spline basis,
#'     that is, `df` or `degree + length(knot) + 1` (because the intercept
#'     is included).
#'
#' @param ZMat The design matrix corresponding to the regression
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
#' @param GENERATE_GRADIENT_HESSIAN A function numerically calculates the
#'    gradient vector and the Hessian matrix of the given weighted
#'    log-likelihood function under two-phase sampling design.
#'
#' @param GENERATE_PHI The function generates the B-spline basis matrix and
#'     the values of a linear combination of B-spline basis function at a
#'     given vector.
#'
#' @return A list includes: a vector consists the estimated coefficients of
#'     the regression parameters, an error message corresponding to the
#'     Hessian matrix, a status message corresponding to the number of
#'     replications of updating `lamda` or `gama` during the calculation,
#'     and an estimated covariance matrix of the regression coefficients.
#'
#' @export
#'
#'@importFrom splines2 bSpline
GR <- function(gama.ini, active.set, d, qn, ZMat, UVec, VVec,
	knot, degree, boundary, DELTA, OMEGA, GENERATE_GRADIENT_HESSIAN,
	GENERATE_PHI){
	### 7/01: adding drop()
	alpha.dim <- qn
	beta.dim <- d
	converge.status <- 1
	error <- 0
	A1 <- cbind(rep(0, alpha.dim - 1), diag(1, nrow = alpha.dim - 1))
	A2 <- cbind(diag(- 1, nrow = alpha.dim - 1), rep(0, alpha.dim - 1))
	A.ori <- cbind(matrix(rep(0, (alpha.dim - 1) * beta.dim), ncol = beta.dim), A1 + A2)

	gama <- gama.ini
	gama.dim <- length(gama)
	active.set <- numeric(length = 0)
	active <- numeric(length = 0)
	lamda <- rep(1, alpha.dim - 1)
	A <- A.ori[active.set, , drop = FALSE]

	count2 <- 0
	### the number of updating times of 'lamda'
	while(max(lamda) > 0){

		delta <- 1
		count1 <- 0
		### the number of updating times of 'gama'
		while(max(abs(delta)) >= 1e-5){
		### if 'max(abs(delta)) < 1e-5' then meeting the stop condition
		### in Step 5

			active.set <- unique(append(active.set, active)[order(append(active.set, active))])

			GGH <- GENERATE_GRADIENT_HESSIAN(Gama = gama,
												d = beta.dim,
											   qn = alpha.dim,
											 ZMat = ZMat,
											 UVec = UVec,
										     VVec = VVec,
											 knot = knot,
										   degree = degree,
										 boundary = boundary,
											DELTA = DELTA,
											OMEGA = OMEGA,
									 GENERATE_PHI = GENERATE_PHI)
			U <- GGH$gradient
			W <- - GGH$hessian
			### obtain U = gradient vector
			### obtain W = negative Hessian = - hessian matrix
			if(is.infinite(max(W))) {error <- 1; break;}
			if(missing(W)){error <- 1; break;}
			if(sum(is.na(W)) > 0) {error <- 1; break;}
			if(min(abs(eigen(W)$values)) < 1e-5 | max(abs(eigen(W)$values)) > 1e20) {error <- 1; break;}
			### save error code 1 while something wrong with 'W', then break
			### in this situation: 'W' is not invertible
			W.inv <- solve(W)
			if(length(active.set) == 0){
			### if 'active.set' is empty, then the direction is just the
			### Newton-Raphson's direction
				direction <- drop(W.inv %*% U)
			} else{
			### if 'active.set' is not empty, then use the formula in Step 0
				A <- A.ori[active.set, , drop = FALSE]
				### always use the new matrix A corrsponding to 'active.set'
				direction <- drop((diag(1, gama.dim, gama.dim) - W.inv %*% t(A) %*%
					  solve(A %*% W.inv %*% t(A)) %*% A) %*% W.inv %*% U)
			### 7/01: adding drop()
			}

			ratio <- - drop(A.ori %*% gama) / drop(A.ori %*% direction)
			step <- ifelse(max(ratio, na.rm = TRUE) <= 0, 1,
						   min(ratio[ratio > 0], na.rm = TRUE))

			ksi <- min(step, 1)
			### 'ksi' = (1 / 2) ^ 'KAPPA' in the draft, 'KAPPA' starts with 0
			gama.update <- gama + ksi * direction
			### obtain the gradient vector U.update at 'gama.update'
			### and the gradient vector U at 'gama'
			while(crossprod(GENERATE_GRADIENT_HESSIAN(Gama = gama.update,
				d = beta.dim, qn = alpha.dim, ZMat = ZMat, UVec = UVec,
				VVec = VVec, knot = knot, degree = degree, boundary = boundary,
				DELTA = DELTA, OMEGA = OMEGA, GENERATE_PHI = GENERATE_PHI)$gradient) >
				  crossprod(GENERATE_GRADIENT_HESSIAN(Gama = gama,
				d = beta.dim, qn = alpha.dim, ZMat = ZMat, UVec = UVec,
				VVec = VVec, knot = knot, degree = degree, boundary = boundary,
				DELTA = DELTA, OMEGA = OMEGA, GENERATE_PHI = GENERATE_PHI)$gradient)){

				ksi <- ksi / 2
				### equivalently 'KAPPA <- KAPPA + 1'
				gama.update <- gama + ksi * direction
				if(ksi < 1e-5) break
				### break the loop if 'KAPPA' is large
			}
			### it breaks if '|| U.update || <= || U ||' and 'KAPPA' is smallest

			if(step > ksi){
				delta <- ksi * direction
				### 'gama.update - gama = ksi * direction'
			} else{
				delta <- step * direction
				### 'gama.update - gama = step * direction'
				active <- which(ratio == step)
				### pick up the possible index
			}
			gama <- gama + delta
			count1 <- count1 + 1
			### the number of updating times of 'gama'
			if(count1 > 20) {converge.status <- 0; break;}
			if(missing(direction)) break
			### if 'gama' cannot be updated to the final one within a few iterations,
			### it means the case will be discarded
		}
		### ending while 'direction' is small enough

		if(length(active.set) == 0){
		### if active set is empty, then 'break', meaning 'active.set' must
		### non-empty during the calculation
			break
		} else{
		### if active set has some elements
			lamda <- solve(A %*% W.inv %*% t(A)) %*% A %*% W.inv %*% U
			inactive <- which.max(lamda)
			active.set <- active.set[- inactive]
		}
		count2 <- count2 + 1
		### the number of updating times of 'lamda'
		if(count2 > 5){
			converge.status <- 0
			break
		}
	}
	### ending of adjusting 'lamda'
	GGH.update <- GENERATE_GRADIENT_HESSIAN(Gama = gama.update,
											   d = beta.dim,
											  qn = alpha.dim,
											ZMat = ZMat,
								   		    UVec = UVec,
								   		    VVec = VVec,
								    		knot = knot,
										  degree = degree,
										boundary = boundary,
										   DELTA = DELTA,
										   OMEGA = OMEGA,
									GENERATE_PHI = GENERATE_PHI)
	variance_observed <- GGH.update$variance_observed
	### update U, W by 'gama.update'
	return(list(gama = gama.update,
			   error = error,
     converge.status = converge.status,
   variance_observed = variance_observed))
}

