test_that("estimators recover the direct effect without noise", {
  set.seed(7)
  A <- er_graph(100, 0.1)
  g <- rep(1:2, 50)
  z <- sample(rep(0:1, 50))
  E <- exposure(A, z, g)
  x <- rnorm(100)
  Y <- 0.5 + 2 * z + E %*% c(1, 3) + 0.7 * x
  expect_equal(dte(Y, z, A, g, X = x, method = "lm"), 2)
  Y2 <- drop(1 + 2 * z + exp(E[, 1] / 5) + E[, 2]^2)
  est <- dte(Y2, z, A, g, method = "match")
  expect_equal(as.vector(est), 2)
  expect_equal(attr(est, "pairs"), n_matched(A, z, g))
})

test_that("difference in means matches its definition", {
  z <- c(1, 0, 1, 0)
  expect_equal(dte(c(3, 1, 5, 2), z), 4 - 1.5)
  expect_error(dte(1:4, z, method = "lm"), "`A`")
})
