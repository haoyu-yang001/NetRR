#' Estimate the direct treatment effect
#'
#' @details
#' \itemize{
#'   \item `"dm"`: difference in means, \eqn{\bar Y_T - \bar Y_C}.
#'   \item `"lm"`: coefficient of `z` in the least-squares regression of `Y`
#'     on `z`, the interference [exposure()] and `X`.
#'   \item `"match"`: the matched estimator \eqn{\hat\eta_m = |P_0|^{-1}
#'     \sum_{(i,j) \in P_0} (Y_i - Y_j)}. Treated and control units with the
#'     same interference profile are paired at random, as many pairs as
#'     possible ([n_matched()]). Under an additive model its variance is
#'     \eqn{2\sigma^2 / |P_0|}.
#' }
#'
#' @inheritParams imbalance
#' @param Y Numeric outcome vector.
#' @param A Adjacency matrix, as in [exposure()]. Not needed for `"dm"`.
#' @param method One of `"dm"`, `"lm"` or `"match"`.
#'
#' @return The estimate. For `"match"` it has attribute `pairs`, the number
#'   of matched pairs.
#'
#' @inherit netrr references
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' A <- matrix(rbinom(n * n, 1, 0.05), n, n)
#' A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
#' z <- netrr(A, nref = 500)
#' Y <- 0.5 + 2 * z + 1.5 * exposure(A, z)[, 1] + rnorm(n)
#' dte(Y, z)
#' dte(Y, z, A, method = "lm")
#' dte(Y, z, A, method = "match")
#' @export
dte <- function(Y, z, A = NULL, group = NULL, X = NULL,
                method = c("dm", "lm", "match")) {
  method <- match.arg(method)
  if (length(Y) != length(z)) stop("`Y` and `z` must have the same length.", call. = FALSE)
  if (method == "dm") {
    .check_z(z, length(z))
    return(mean(Y[z == 1]) - mean(Y[z == 0]))
  }
  if (is.null(A)) stop("`A` is needed for method \"", method, "\".", call. = FALSE)
  E <- exposure(A, z, group)
  if (method == "lm") {
    D <- cbind(1, z, E, if (!is.null(X)) as.matrix(X))
    return(unname(stats::lm.fit(D, Y)$coefficients[2]))
  }
  d <- unlist(lapply(split(seq_along(z), .key(E)), function(i) {
    t1 <- i[z[i] == 1]
    t0 <- i[z[i] == 0]
    m <- min(length(t1), length(t0))
    Y[.pick(t1, m)] - Y[.pick(t0, m)]
  }), use.names = FALSE)
  if (length(d) == 0) warning("No treated-control pair has identical interference.", call. = FALSE)
  structure(mean(d), pairs = length(d))
}
