test_that("netrr returns a balanced assignment of the requested size", {
  set.seed(3)
  A <- er_graph(80, 0.2)
  g <- rep(1:2, 40)
  x <- rnorm(80)
  z <- netrr(A, group = g, X = x, n1 = 30, pa = 0.05, nref = 200)
  expect_equal(sum(z), 30)
  expect_lte(attr(z, "imbalance"), attr(z, "threshold"))
  expect_equal(attr(z, "imbalance"), imbalance(A, as.vector(z), g, x))
})

test_that("netrr warns and returns the best draw when maxit is too small", {
  set.seed(4)
  A <- er_graph(60, 0.2)
  expect_warning(z <- netrr(A, pa = 0.001, nref = 100, maxit = 1), "best")
  expect_equal(sum(z), 30)
})

test_that("netmm reports the correct number of pairs and beats complete randomization", {
  set.seed(5)
  A <- er_graph(100, 0.05)
  g <- rep(1:2, 50)
  z <- netmm(A, group = g, n1 = 40, maxit = 300)
  expect_equal(sum(z), 40)
  expect_equal(attr(z, "pairs"), n_matched(A, as.vector(z), g))
  cr <- replicate(50, n_matched(A, sample(rep(1:0, c(40, 60))), g))
  expect_gt(attr(z, "pairs"), max(cr))
})

test_that("netmm handles weighted adjacency matrices", {
  set.seed(6)
  A <- er_graph(50, 0.1) * 0.1
  z <- netmm(A, maxit = 200)
  expect_equal(attr(z, "pairs"), n_matched(A, as.vector(z)))
})
