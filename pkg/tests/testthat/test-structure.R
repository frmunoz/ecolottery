context("Tests Structure")

test_that("Result structure of the forward function with limiting similarity", {
  res <- forward(initial = data.frame(sp = sort(rep(as.character(1:10), 10)), trait = runif(100)), prob = 0.1, gens = 100, pool = sort(rep(as.character(1:10), 10)), filt = function(x) 0.5 - abs(0.5 - x), limit.sim = TRUE, keep = TRUE)
  expect_is(res, "list", info = "The result of the forward function must be a list")
  expect_is(res[[1]], "list", info = "The aa result (first item of the result given by the forward function) must be a list")
  expect_is(res[[1]][[1]], "data.frame", info = "Each element of the aa result (forward function) must be data.frame")
  
  expect_is(res[[2]], "numeric", info = "The limit.sim.t result (second item of the result given by the forward function) must be numeric")
  expect_is(res[[3]], "data.frame", info = "The pool result (third item of the result given by the forward function) must be a data.frame")
  
  expect_true(sum(names(res) == c("aa", "limit.sim.t", "pool")) == 3, info = "The result of the forward function must be named aa, limit.sim.t and pool")
})

test_that("Result structure of the coalesc function ", {
  res <- coalesc(100, 50, 0.1, filt = function(x) 1 - x) 
  
  expect_is(res, "matrix", info = "The result of the coalesc function must be a matrix")
  expect_true(ncol(res) == 3, info = "The result of the coalesc function must be a matrix with 3 columns")
})