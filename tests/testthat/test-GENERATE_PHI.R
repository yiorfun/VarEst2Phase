test_that("GENERATE_PHI() function returns the Spline basis matrix and the vector of a linear combination of B-spline basis functions evaluated at a given vector", {
  TEST_RES <- GENERATE_PHI(x = c(0.1, 0.2, 0.3),
						            knot = c(0.25, 0.75),
					            degree = 2,
					          boundary = c(0, 1),
				        coeff_spline = c(0.1, 0.1, 0.1, 0.1, 0.1))
  expect_equal(TEST_RES$phi_vector, c(0.1, 0.1, 0.1))
})
