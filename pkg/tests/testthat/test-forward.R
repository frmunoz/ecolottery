context("Test forward() structure and extreme cases")

test_that("forward() works and returns desired result structure", {
  pool <- data.frame(id = 1:10000,
                     sp = rep(1:500, each = 20),
                     trait = rep(runif(500), 20))
  
  set.seed(1)
  initial <- pool[sample(nrow(pool), size = 50, replace = FALSE),]
  
  sink(tempfile())
  res <- forward(initial, prob = 0.1, gens = 10, pool = pool, limit.sim = TRUE)
  sink()
  
  expect_is(res, "list")
  expect_named(res, c("com", "sp_t", "dist.t", "pool"))
  expect_is(res$com, "data.frame")
  expect_is(res$pool, "data.frame")
  expect_named(res$com, c("id", "sp", "trait"))
  expect_named(res$pool, c("id", "sp", "trait"))
  expect_is(res$com$trait, "numeric")
  expect_is(res$pool$trait, "numeric")
  expect_equal(dim(res$com), c(50, 3))
  expect_equal(dim(res$pool), dim(pool))
  
  
  # Simulations when keeping all generations
  sink(tempfile())
  res_keep <- forward(initial[, -1], prob = 0.1, gens = 4, pool = pool,
                      limit.sim = FALSE, keep = TRUE)
  sink()
  
  expect_is(res_keep, "list")
  expect_named(res_keep, c("com_t", "sp_t", "pool"))
  expect_is(res_keep$com_t, "list")
  expect_is(res_keep$com_t[[1]], "data.frame")
  expect_is(res_keep$pool, "data.frame")
  expect_named(res_keep$com_t[[1]], c("id", "sp", "trait"))
  expect_named(res_keep$pool, c("id", "sp", "trait"))
  expect_is(res_keep$com_t[[1]]$trait, "numeric")
  expect_is(res_keep$pool$trait, "numeric")
  expect_equal(dim(res_keep$com_t[[1]]), c(50, 3))
  expect_length(res_keep$com_t, 4)
  expect_equal(dim(res_keep$pool), dim(pool))
  
  # With limiting similarity
  sink(tempfile())
  res_limit <- forward(initial[, -1], prob = 0.1, gens = 4, pool = pool,
                      limit.sim = TRUE, keep = TRUE)
  sink()
  
  expect_is(res_limit, "list")
  expect_named(res_limit, c("com_t", "sp_t", "dist.t", "pool"))
  expect_is(res_limit$com_t, "list")
  expect_is(res_limit$com_t[[1]], "data.frame")
  expect_is(res_limit$pool, "data.frame")
  expect_is(res_limit$dist.t, "numeric")
  expect_named(res_limit$com_t[[1]], c("id", "sp", "trait"))
  expect_named(res_limit$pool, c("id", "sp", "trait"))
  expect_is(res_limit$com_t[[1]]$trait, "numeric")
  expect_is(res_limit$pool$trait, "numeric")
  expect_equal(dim(res_limit$com_t[[1]]), c(50, 3))
  expect_length(res_limit$com_t, 4)
  expect_equal(dim(res_limit$pool), dim(pool))
  expect_length(res_limit$dist.t, 4)
})