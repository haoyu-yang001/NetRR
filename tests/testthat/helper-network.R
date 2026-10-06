er_graph <- function(n, p) {
  A <- matrix(0, n, n)
  A[upper.tri(A)] <- stats::rbinom(n * (n - 1) / 2, 1, p)
  A + t(A)
}
