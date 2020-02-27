context("Test forward() structure and extreme cases")

test_that("forward() works and returns desired result structure", {
  
  # Case with speciation
  set.seed(1)
  initial <- data.frame(ind = paste("init",1:100, sep="."),
                        sp = rep(1:10, each = 10),
                        trait = rep(runif(10), 10))
  
  sink(tempfile())
  res <- forward(initial, theta = 50, gens = 5, pool = NULL, limit.sim = TRUE)
  sink()
  
  expect_is(res, "list")
  expect_named(res, c("com", "sp_t", "dist.t", "pool", "call"))
  expect_is(res$com, "data.frame")
  expect_named(res$com, c("ind", "sp", "trait"))
  expect_is(res$com$trait, "numeric")
  expect_equal(dim(res$com), c(100, 3))

  # Cases with migration (an external pool is defined)
  pool <- data.frame(ind = paste("pool",1:10000, sep="."),
                     sp = rep(1:500, each = 20),
                     trait = rep(runif(500), 20))
  
  set.seed(1)
  initial <- data.frame(ind = paste("init",1:100, sep="."),
                        sp = rep(1:10, each = 10),
                        trait = rep(runif(10), 10))
  
  sink(tempfile())
  res <- forward(initial, m = 0.1, gens = 5, pool = pool, limit.sim = TRUE)
  sink()
  
  expect_is(res, "list")
  expect_named(res, c("com", "sp_t", "dist.t", "pool", "call"))
  expect_is(res$com, "data.frame")
  expect_is(res$pool, "data.frame")
  expect_named(res$com, c("ind", "sp", "trait"))
  expect_named(res$pool, c("ind", "sp", "trait"))
  expect_is(res$com$trait, "numeric")
  expect_is(res$pool$trait, "numeric")
  expect_equal(dim(res$com), c(100, 3))
  expect_equal(dim(res$pool), dim(pool))
  
  # Simulations when keeping all generations
  sink(tempfile())
  res_keep <- forward(initial, m = 0.1, gens = 5, pool = pool,
                      limit.sim = FALSE, keep = TRUE)
  sink()
  
  expect_is(res_keep, "list")
  expect_named(res_keep, c("com_t", "sp_t", "pool", "call"))
  expect_is(res_keep$com_t, "list")
  expect_is(res_keep$com_t[[1]], "data.frame")
  expect_is(res_keep$pool, "data.frame")
  expect_named(res_keep$com_t[[1]], c("ind", "sp", "trait"))
  expect_named(res_keep$pool, c("ind", "sp", "trait"))
  expect_is(res_keep$com_t[[1]]$trait, "numeric")
  expect_is(res_keep$pool$trait, "numeric")
  expect_equal(dim(res_keep$com_t[[1]]), c(100, 3))
  expect_length(res_keep$com_t, 5)
  expect_equal(dim(res_keep$pool), dim(pool))
  
  # With limiting similarity
  sink(tempfile())
  res_limit <- forward(initial, m = 0.1, gens = 5, pool = pool,
                      limit.sim = TRUE, keep = TRUE)
  sink()
  
  expect_is(res_limit, "list")
  expect_named(res_limit, c("com_t", "sp_t", "dist.t", "pool", "call"))
  expect_is(res_limit$com_t, "list")
  expect_is(res_limit$com_t[[1]], "data.frame")
  expect_is(res_limit$pool, "data.frame")
  expect_is(res_limit$dist.t, "numeric")
  expect_named(res_limit$com_t[[1]], c("ind", "sp", "trait"))
  expect_named(res_limit$pool, c("ind", "sp", "trait"))
  expect_is(res_limit$com_t[[1]]$trait, "numeric")
  expect_is(res_limit$pool$trait, "numeric")
  expect_equal(dim(res_limit$com_t[[1]]), c(100, 3))
  expect_length(res_limit$com_t, 5)
  expect_equal(dim(res_limit$pool), dim(pool))
  expect_length(res_limit$dist.t, 5)
})