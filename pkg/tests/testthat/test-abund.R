context("Test abund()")

test_that("get_rel_abund() works", {
  given_com <- data.frame(ind = 1:4, sp = c(1, 1, 2, 3))
  
  expected_rel_abund <- data.frame(sp = as.character(1:3),
                                   ab = as.integer(c(2, 1, 1)),
                                   relab = c(0.5, 0.25, 0.25),
                                   stringsAsFactors = FALSE)
  
  comp_rel_abund <- get_rel_abund(given_com)
  
  expect_is(comp_rel_abund, "data.frame")
  expect_named(comp_rel_abund, c("sp", "ab", "relab"))
  expect_equal(dim(comp_rel_abund), c(3, 3))
  expect_equal(get_rel_abund(given_com),
               expected_rel_abund)
  
  # Extreme case single individual of a single species in community
  expect_equal(get_rel_abund(data.frame(ind = 1, sp = 1)), 
               data.frame(sp = "1", ab = 1, relab = 1,
                          stringsAsFactors = FALSE))
  
})

test_that("abund() works", {
  
  fwd_initial <- sort(rep(as.character(1:10), 10))
  fwd_final   <- forward(initial = fwd_initial, prob = 0.1, gens = 1000)
  fwd_keep    <- forward(initial = fwd_initial, prob = 0.1, gens = 1000,
                         keep = TRUE)
  suppressWarnings({
    coal_final  <- coalesc(100, 0.9, 50)
    
    ab_fwd  <- abund(fwd_final)
    ab_keep <- abund(fwd_keep)
    ab_coal <- abund(coal_final)
  })
  
  expect_length(ab_coal, 2)
  expect_named(ab_coal, c("com", "pool"))
  
  expect_length(ab_fwd, 2)
  expect_warning(abund(fwd_final), fixed = TRUE,
                 "Some communities were undefined; returning NA abundances")
  
  expect_length(ab_keep, 3)
  expect_length(ab_keep$com_t, 1000)
  expect_named(ab_keep, names(fwd_keep))
})