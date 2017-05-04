context("Test coalesc() structure and extreme cases")

test_that("coalesc() result structure", {
  
  # Default results when defining nothing
  suppressWarnings({
    res <- coalesc(100, 0.1, 50, filt = function(x) 1 - x) 
  })
  
  expect_is(res, "list", info = "Result must be a list")
  expect_is(res$com, "data.frame")
  expect_is(res$pool, "data.frame")
  expect_equal(ncol(res$com), 3)
  expect_equal(nrow(res$com), 100)
  expect_equal(ncol(res$pool), 3)
  expect_equal(nrow(res$pool), 5000)
  expect_named(res$com, c("ind", "sp", "tra1"))
  expect_named(res$pool, c("ind", "sp", "tra1"))
  expect_is(res$com$tra1, "numeric")
  expect_is(res$com$tra1, "numeric")
  
  # Result structure when defining pool
  pool <- cbind(1:10000, rep(1:500, each = 20), rep(runif(500), each = 20))
  
  res <- coalesc(50, 0.1, pool = pool)
  
  expect_is(res, "list", info = "Result must be a list")
  expect_is(res$com, "data.frame")
  expect_is(res$pool, "data.frame")
  expect_equal(ncol(res$com), 3)
  expect_equal(nrow(res$com), 50)
  expect_equal(ncol(res$pool), 3)
  expect_equal(nrow(res$pool), 10000)
  expect_named(res$com, c("ind", "sp", "tra1"))
  expect_named(res$pool, c("ind", "sp", "tra1"))
  expect_is(res$com$tra1, "numeric")
  expect_is(res$com$tra1, "numeric")
})


test_that("coalesc() behaves well with extreme cases", {
  
  pool <- cbind(1:10000, rep(1:500, each = 20), rep(runif(500), each = 20))
  
  ## Invalid inputs (negative parameters)
  # Bad J values
  expect_error(coalesc(J = -4, theta = 50), "J must be positive.",
               fixed = TRUE)
  expect_error(coalesc(0, theta = 50), "J must be positive.",
               fixed = TRUE)
  
  # Problematic m values
  # m too high
  expect_error(coalesc(10, m = 10, 40), fixed = TRUE,
               "The migration parameter takes values between 0 and 1")
  
  # negative m
  expect_error(coalesc(10, m = -1, 40), fixed = TRUE,
               "The migration parameter takes values between 0 and 1")
  
  # null m
  expect_silent(coalesc(10, m = 0, 40))
  
  # Bad theta values
  expect_error(coalesc(10, m = 0.1, -1), fixed = TRUE,
               "The theta parameter must be positive")
  # Provide both pool and theta
  expect_warning(coalesc(10, 0.1, 40, pool = pool, verbose = TRUE),
                 fixed = TRUE,
                 "Both a theta value and a regional pool provided, discarding theta")
  
  # Extreme pool of species (one individual one species, individual
  # from a single species)
  
  
  # Community bigger than pool of species
  
  # Community with a single individual
})
