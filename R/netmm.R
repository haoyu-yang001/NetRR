#' Network maximized matching randomization (NetMM)
#'
#' Searches for an assignment with `n1` treated units that maximizes
#' \eqn{|P_0|}, the number of treated-control pairs with identical
#' interference ([n_matched()]).
#'
#' @details
#' The search starts from a complete randomization and keeps the best
#' assignment found. At each iteration it swaps \eqn{k} random treated units
#' with \eqn{k} random control units, where \eqn{k = 1 + \mathrm{Geom}(1/2)}
#' is capped at \eqn{\min(n_1, n - n_1)}. The candidate replaces the current
#' assignment when its \eqn{|P_0|} is at least as large. Any assignment can
#' be reached in one step with positive probability, so the result converges
#' to a global maximizer as `maxit` grows. The interference is updated
#' incrementally, so each iteration costs \eqn{O(nkK)}.
#'
#' @inheritParams netrr
#' @param maxit Number of search iterations.
#'
#' @return A 0/1 assignment vector with attribute `pairs`, its
#'   \eqn{|P_0|}.
#'
#' @inherit netrr references
#' @seealso [netrr()], [n_matched()], [dte()]
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' A <- matrix(rbinom(n * n, 1, 0.05), n, n)
#' A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
#' z <- netmm(A, maxit = 500)
#' attr(z, "pairs")
#' n_matched(A, sample(rep(0:1, n / 2)))  # complete randomization
#' @export
netmm <- function(A, group = NULL, n1 = floor(nrow(A) / 2), maxit = 1000) {
  n <- .check_A(A)
  .check_n1(n1, n)
  G <- .onehot(group, n)
  z <- .draw(n, n1)
  E <- as.matrix(A %*% (z * G))
  P <- .npairs(E, z)
  kmax <- min(n1, n - n1)
  for (t in seq_len(maxit)) {
    k <- min(1 + stats::rgeom(1, 0.5), kmax)
    idx <- c(.pick(which(z == 1), k), .pick(which(z == 0), k))
    z1 <- z
    z1[idx] <- 1 - z1[idx]
    E1 <- E + as.matrix(A[, idx, drop = FALSE] %*% ((z1[idx] - z[idx]) * G[idx, , drop = FALSE]))
    P1 <- .npairs(E1, z1)
    if (P1 >= P) {
      z <- z1
      E <- E1
      P <- P1
    }
  }
  structure(z, pairs = P)
}
