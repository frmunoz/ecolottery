context("Test tcor()")

test_that("tcor() function works", {
  
  # Invalid input
  expect_error(tcor())
  expect_error(tcor(-4))
  expect_error(tcor(NA))
  
  # Good dimensions and structure
  expect_equal(dim(tcor(10)), c(10, 2))
  expect_is(tcor(10), "data.frame")
  expect_equal(colnames(tcor(10)), c("t1", "t2"))
  
  # Make sure the correlation is enough
  traits = tcor(10000)
  expect_equal(cor(traits$t1, traits$t2), 0.5, tolerance = 0.02)
  
  # Extreme correlation values
  expect_silent(tcor(10, rho = 0))
  expect_silent(tcor(10, rho = 1))
  expect_error(tcor(10, rho = 4), "rho must belong to [0; 1] interval",
               fixed = TRUE)
  expect_error(tcor(10, rho = -4), "rho must belong to [0; 1] interval",
               fixed = TRUE)
  expect_error(tcor(10, rho = NA), "rho must belong to [0; 1] interval",
               fixed = TRUE)
  
  # When providing trait vector
  given_trait = rnorm(10000)
  corr_traits = tcor(10000, rho = 0.7, x = given_trait)
  expect_equal(colnames(corr_traits), c("trait", "t2"))
  expect_equal(cor(corr_traits$trait, corr_traits$t2), 0.7, tolerance = 0.02)
  
  expect_warning(tcor(100, x = given_trait),
                 "Provided trait vector x does not have length n!")
})