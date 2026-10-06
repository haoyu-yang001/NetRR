test_that("exposure counts treated neighbours in each group", {
  set.seed(1)
  A <- er_graph(30, 0.2)
  z <- sample(rep(0:1, 15))
  g <- sample(c("b", "a", "c"), 30, TRUE)
  expect_equal(drop(exposure(A, z)), colSums(A[z == 1, ]))
  E <- exposure(A, z, g)
  expect_equal(colnames(E), c("a", "b", "c"))
  for (k in colnames(E)) {
    expect_equal(E[, k], drop(A %*% (z * (g == k))), ignore_attr = TRUE)
  }
})

test_that("imbalance matches its definition and is scale invariant", {
  set.seed(2)
  A <- er_graph(40, 0.2)
  z <- sample(rep(0:1, 20))
  x <- rnorm(40)
  E <- cbind(exposure(A, z), x)
  d <- colMeans(E[z == 1, ]) - colMeans(E[z == 0, ])
  M <- 40 * 0.25 * drop(d %*% solve(cov(E)) %*% d)
  expect_equal(imbalance(A, z, X = x), M)
  expect_equal(imbalance(A, z, X = 10 * x + 3), M)
  expect_equal(imbalance(matrix(0, 40, 40), z), 0)
})

test_that("n_matched equals the overlap of interference profiles", {
  A <- matrix(c(0, 1, 0, 0,
                1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0), 4, 4)
  z <- c(1, 0, 0, 1)
  # exposures: 0, 1, 1, 0 -> profile 0 has one treated and no control
  expect_equal(n_matched(A, z), 0)
  expect_equal(n_matched(A, c(1, 1, 0, 0)), 1)
})

test_that("invalid input is rejected", {
  A <- er_graph(10, 0.3)
  expect_error(exposure(A[, -1], rep(0:1, 5)), "square")
  expect_error(exposure(A, rep(0:1, 4)), "0/1")
  expect_error(exposure(A, rep(1, 10)), "Both arms")
  expect_error(exposure(A, rep(0:1, 5), group = 1:3), "group")
  expect_error(netrr(A, n1 = 0), "n1")
})
