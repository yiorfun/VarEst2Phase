#' Calculate the gradient vector and the Hessian matrix of the weighted
#' log-likelihood function derived from a Cox proportional-hazards model
#' with interval censored data under two-phase sampling design
#'
#' @param Gama A vector whose the first `d` elements are identical with
#'     the `d`-dimensional regrssion coefficients and the last `qn`
#'     elements are identical with the `qn`-dimensional spline
#'     coefficients.
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
#' @param GENERATE_PHI The function generates the B-spline basis matrix and
#'     the values of a linear combination of B-spline basis function at a
#'     given vector.
#'
#' @return A list includes: a `d + qn` dimensional gradient vector evaluated
#'     at given `Gama`, a `(d + qn)` by `(d + qn)` Hessian matrix evaluated
#'     at given `Gama`, and a `d` by `d` estimated covariance matrix of the
#'     regression parameters at given `Gama`.
#'
#' @export
#'
#'@importFrom splines2 bSpline
GENERATE_COX_GRADIENT_HESSIAN <- function(Gama, d, qn, ZMat,
	UVec, VVec, knot, degree, boundary, DELTA, OMEGA, GENERATE_PHI){

	sample_size <- nrow(ZMat)
	Alpha <- Gama[1 : d]
	Beta <- Gama[(d + 1) : (d + qn)]

	PhiU_res <- GENERATE_PHI(x = UVec,
						  knot = knot,
					    degree = degree,
				      boundary = boundary,
			      coeff_spline = Beta)
	PhiUVec <- PhiU_res$phi_vector
	PhiUMat <- PhiU_res$phi_matrix
	PhiV_res <- GENERATE_PHI(x = VVec,
						  knot = knot,
					    degree = degree,
				      boundary = boundary,
			      coeff_spline = Beta)
	PhiVVec <- PhiV_res$phi_vector
	PhiVMat <- PhiV_res$phi_matrix

	#ZPT <- drop(ZMat %*% Alpha)
	#AU <- 1 + as.vector(exp(PhiUVec + ZPT))
	#BV <- 1 + as.vector(exp(PhiVVec + ZPT))
	#invAU <- 1 / AU
	#invBV <- 1 / BV
	#invAU2 <- invAU ^ 2 - invAU
	#invBV2 <- invBV ^ 2 - invBV

	ZPT <- drop(ZMat %*% Alpha)
	AUVec <- as.vector(exp(PhiUVec + ZPT))
	BUVec <- exp(- AUVec)
	AVVec <- as.vector(exp(PhiVVec + ZPT))
	BVVec <- exp(- AVVec)

	AUBU <- AUVec * BUVec
	AVBV <- AVVec * BVVec
	AUBU_OMBU <- AUBU / (1 - BUVec)
	AU_OMBU <- AUVec / (1 - BUVec)
	AUBU_BUMBV <- AUBU / (BUVec - BVVec)
	AVBV_BUMBV <- AVBV / (BUVec - BVVec)
	AUBUMAVBV_BUMBV <- AUBU_BUMBV - AVBV_BUMBV
	BUBVAUMAV_BUMBVSQ <- BUVec * BVVec * (AUVec - AVVec) / (BUVec - BVVec) ^ 2

	GradVec <- rep(0, d + qn)
	HessMat <- matrix(0, nrow = d + qn, ncol = d + qn)
	LDotOneMat <- matrix(0, nrow = sample_size, ncol = d)
	### n * d matrix to compute observed information,
	### the i-th row saves the values of LDotOne at each parameter given Xi
	LDotTwoMat <- matrix(0, nrow = sample_size, ncol = qn)
	### n * q_n matrix to compute observed information,
	### the i-th row saves the values of LDotTwo at each parameter given Xi

	Index_1 <- which(DELTA[ , 1] == 1)
	Index_2 <- which(DELTA[ , 2] == 1)
	Index_3 <- which(DELTA[ , 3] == 1)

	for(k in seq(d)){
		#Grad_1 <- as.vector(invAU[Index_1] * ZMat[Index_1, k])
		#Grad_2 <- as.vector((invAU[Index_2] + invBV[Index_2] - 1) * ZMat[Index_2, k])
		#Grad_3 <- as.vector((invBV[Index_3] - 1) * ZMat[Index_3, k])
		#GradVec[k] <- sum(OMEGA[Index_1] * Grad_1) + sum(OMEGA[Index_2] * Grad_2) + sum(OMEGA[Index_3] * Grad_3)
		#LDotOneMat[ , k] <- c(Grad_1, Grad_2, Grad_3)

		Grad_1 <- as.vector(ZMat[Index_1, k] * AUBU_OMBU[Index_1])
		Grad_2 <- as.vector((- ZMat[Index_2, k]) * AUBUMAVBV_BUMBV[Index_2])
		Grad_3 <- as.vector(ZMat[Index_3, k] * AVVec[Index_3])
		GradVec[k] <- sum(OMEGA[Index_1] * Grad_1) + sum(OMEGA[Index_2] * Grad_2) - sum(OMEGA[Index_3] * Grad_3)
		LDotOneMat[ , k] <- c(Grad_1, Grad_2, - Grad_3)

		for(kp in seq(d)){
			#Hess_1 <- sum(OMEGA[Index_1] * invAU2[Index_1] * ZMat[Index_1, k] * ZMat[Index_1, kp])
			#Hess_2 <- sum(OMEGA[Index_2] * (invAU2[Index_2] + invBV2[Index_2]) * ZMat[Index_2, k] * ZMat[Index_2, kp])
			#Hess_3 <- sum(OMEGA[Index_3] * invBV2[Index_3] * ZMat[Index_3, k] * ZMat[Index_3, kp])
			#HessMat[k, kp] <- Hess_1 + Hess_2 + Hess_3

			Hess_1 <- sum(ZMat[Index_1, k] * ZMat[Index_1, kp] * AUBU_OMBU[Index_1] * (1 - AU_OMBU[Index_1]))
			Hess_2 <- sum(ZMat[Index_2, k] * ZMat[Index_2, kp] * (AUBUMAVBV_BUMBV[Index_2] + BUBVAUMAV_BUMBVSQ[Index_2] * (AUVec[Index_2] - AVVec[Index_2])))
			Hess_3 <- sum(ZMat[Index_3, k] * ZMat[Index_3, kp] * AVVec[Index_3])
			HessMat[k, kp] <- Hess_1 - Hess_2 - Hess_3

		}
		for(jp in seq(qn)){
			#Hess_1 <- sum(OMEGA[Index_1] * invAU2[Index_1] * ZMat[Index_1, k] * PhiUMat[Index_1, jp])
			#Hess_2 <- sum(OMEGA[Index_2] * (invAU2[Index_2] * ZMat[Index_2, k] * PhiUMat[Index_2, jp] + invBV2[Index_2] * ZMat[Index_2, k] * PhiVMat[Index_2, jp]))
			#Hess_3 <- sum(OMEGA[Index_3] * invBV2[Index_3] * ZMat[Index_3, k] * PhiVMat[Index_3, jp])
			#HessMat[k, d + jp] <- Hess_1 + Hess_2 + Hess_3

			Hess_1 <- sum(ZMat[Index_1, k] * PhiUMat[Index_1, jp] * AUBU_OMBU[Index_1] * (1 - AU_OMBU[Index_1]))
			Hess_2 <- sum(ZMat[Index_2, k] * (AUBU_BUMBV[Index_2] * PhiUMat[Index_2, jp] - AVBV_BUMBV[Index_2] * PhiVMat[Index_2, jp] + BUBVAUMAV_BUMBVSQ[Index_2] * (AUVec[Index_2] * PhiUMat[Index_2, jp] - AVVec[Index_2] * PhiVMat[Index_2, jp])))
			Hess_3 <- sum(ZMat[Index_3, k] * PhiVMat[Index_3, jp] * AVVec[Index_3])
			HessMat[k, d + jp] <- Hess_1 - Hess_2 - Hess_3
		}
	}
	for(j in seq(qn)){
		#Grad_1 <- as.vector(invAU[Index_1] * PhiUMat[Index_1, j])
		#Grad_2 <- as.vector((((AU[Index_2] - 1) / (BV[Index_2] - AU[Index_2]) + invBV[Index_2]) * PhiVMat[Index_2, j] + ((1 - BV[Index_2]) / (BV[Index_2] - AU[Index_2]) + invAU[Index_2]) * PhiUMat[Index_2, j]))
		#Grad_3 <- as.vector((invBV[Index_3] - 1) * PhiVMat[Index_3, j])
		#GradVec[d + j] <- sum(OMEGA[Index_1] * Grad_1) + sum(OMEGA[Index_2] * Grad_2) + sum(OMEGA[Index_3] * Grad_3)
		#LDotTwoMat[ , j] <- c(Grad_1, Grad_2, Grad_3)

		Grad_1 <- as.vector(PhiUMat[Index_1, j] * AUBU_OMBU[Index_1])
		Grad_2 <- as.vector((- PhiUMat[Index_2, j] * AUBU_BUMBV[Index_2] + PhiVMat[Index_2, j] * AVBV_BUMBV[Index_2]))
		Grad_3 <- as.vector(PhiVMat[Index_3, j] * AVVec[Index_3])
		GradVec[d + j] <- sum(OMEGA[Index_1] * Grad_1) + sum(OMEGA[Index_2] * Grad_2) - sum(OMEGA[Index_3] * Grad_3)
		LDotTwoMat[ , j] <- c(Grad_1, Grad_2, - Grad_3)
		for(kp in seq(d)){
			#Hess_1 <- sum(OMEGA[Index_1] * invAU2[Index_1] * PhiUMat[Index_1, j] * ZMat[Index_1, kp])
			#Hess_3 <- sum(OMEGA[Index_3] * invBV2[Index_3] * PhiVMat[Index_3, j] * ZMat[Index_3, kp])
			#Hess_2 <- sum(OMEGA[Index_2] * (invBV2[Index_2] * PhiVMat[Index_2, j] * ZMat[Index_2, kp] + invAU2[Index_2] * PhiUMat[Index_2, j] * ZMat[Index_2, kp]))
			#HessMat[d + j, kp] <- Hess_1 + Hess_2 + Hess_3

			Hess_1 <- sum(PhiUMat[Index_1, j] * ZMat[Index_1, kp] * AUBU_OMBU[Index_1] * (1 - AU_OMBU[Index_1]))
			Hess_2 <- sum(ZMat[Index_2, kp] * (AUBU_BUMBV[Index_2] * PhiUMat[Index_2, j] - AVBV_BUMBV[Index_2] * PhiVMat[Index_2, j] + BUBVAUMAV_BUMBVSQ[Index_2] * (AUVec[Index_2] * PhiUMat[Index_2, j] - AVVec[Index_2] * PhiVMat[Index_2, j])))
			Hess_3 <- sum(PhiVMat[Index_3, j] * ZMat[Index_3, kp] * AVVec[Index_3])
			HessMat[d + j, kp] <- Hess_1 - Hess_2 - Hess_3
		}
		for(jp in seq(qn)){
			#Hess_1 <- sum(OMEGA[Index_1] * invAU2[Index_1] * PhiUMat[Index_1, j] * PhiUMat[Index_1, jp])
			#Hess_3 <- sum(OMEGA[Index_3] * invBV2[Index_3] * PhiVMat[Index_3, j] * PhiVMat[Index_3, jp])
			#Hess_2 <- sum(OMEGA[Index_2] * (((AU[Index_2] - 1) * (BV[Index_2] - 1) / (BV[Index_2] - AU[Index_2]) ^ 2) * (PhiUMat[Index_2, jp] * PhiVMat[Index_2, j] + PhiVMat[Index_2, jp] * PhiUMat[Index_2, j] - PhiVMat[Index_2, jp] * PhiVMat[Index_2, j] - PhiUMat[Index_2, jp] * PhiUMat[Index_2, j]) + invAU2[Index_2] * PhiUMat[Index_2, jp] * PhiUMat[Index_2, j] + invBV2[Index_2] * PhiVMat[Index_2, jp] * PhiVMat[Index_2, j]))
			#HessMat[d + j, d + jp] <- Hess_1 + Hess_2 + Hess_3

			Hess_1 <- sum(PhiUMat[Index_1, j] * PhiUMat[Index_1, jp] * AUBU_OMBU[Index_1] * (1 - AU_OMBU[Index_1]))
			Hess_2 <- sum((AUBU_BUMBV[Index_2] * PhiUMat[Index_2, j] * PhiUMat[Index_2, jp] - AVBV_BUMBV[Index_2] * PhiVMat[Index_2, j] * PhiVMat[Index_2, jp] + BUVec[Index_2] * BVVec[Index_2] * (AUVec[Index_2] * PhiUMat[Index_2, j] - AVVec[Index_2] * PhiVMat[Index_2, j]) * (AUVec[Index_2] * PhiUMat[Index_2, jp] - AVVec[Index_2] * PhiVMat[Index_2, jp]) / (BUVec[Index_2] - BVVec[Index_2]) ^ 2))
			Hess_3 <- sum(PhiVMat[Index_3, j] * PhiVMat[Index_3, jp] * AVVec[Index_3])
			HessMat[d + j, d + jp] <- Hess_1 - Hess_2 - Hess_3
		}
	}

	A11 <- matrix(0, nrow = d, ncol = d)
	A12 <- matrix(0, nrow = d, ncol = qn)
	A22 <- matrix(0, nrow = qn, ncol = qn)
	WEIGHTS <- c(OMEGA[Index_1], OMEGA[Index_2], OMEGA[Index_3])
	for(i in seq(sample_size)){
		A11 <- A11 + WEIGHTS[i] * tcrossprod(LDotOneMat[i, ])
		A12 <- A12 + WEIGHTS[i] * tcrossprod(LDotOneMat[i, ], LDotTwoMat[i, ])
		A22 <- A22 + WEIGHTS[i] * tcrossprod(LDotTwoMat[i, ])
	}
	### For a matrix A, A[1, ] is a column vector
	A11 <- A11 / sample_size
	A12 <- A12 / sample_size
	A21 <- t(A12)
	A22 <- A22 / sample_size
	ObsInf <- A11 - A12 %*% ginv(A22) %*% A21
	InverseInformation <- ginv(ObsInf)
	variance_observed <- InverseInformation / sample_size
	### usually, I^{-1}(theta_0) / n is the variance,

	return(list(gradient = GradVec,
				 hessian = HessMat,
       variance_observed = variance_observed))
}

