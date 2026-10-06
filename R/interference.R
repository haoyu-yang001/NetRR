#' Interference experienced by each unit
#'
#' Counts the treated neighbours of every unit, separately for each group of
#' neighbours: \eqn{a_{ik} = \sum_j A_{ij} z_j 1\{g_j = k\}}.
#'
#' @param A Adjacency matrix (`n x n`). `A[i, j]` is the weight of the link
#'   through which unit `j` affects unit `i`; use a symmetric 0/1 matrix for
#'   an undirected network. Sparse matrices from the \pkg{Matrix} package
#'   also work.
#' @param z 0/1 treatment assignment of length `n`.
#' @param group Optional group label of each unit, used when peer effects
#'   differ across groups of neighbours (heterogeneous interference). `NULL`
#'   puts all units in one group (homogeneous interference).
#'
#' @return An `n x K` matrix with one column per group.
#'
#' @examples
#' A <- matrix(c(0, 1, 1,
#'               1, 0, 0,
#'               1, 0, 0), 3, 3)
#' exposure(A, c(1, 0, 1))
#' exposure(A, c(1, 0, 1), group = c("a", "a", "b"))
#' @export
exposure <- function(A, z, group = NULL) {
  n <- .check_A(A)
  .check_z(z, n)
  as.matrix(A %*% (z * .onehot(group, n)))
}

#' Interference imbalance between the treatment and control groups
#'
#' The Mahalanobis distance used by [netrr()]:
#' \deqn{M = n p (1 - p) (\bar a_T - \bar a_C)^\top \mathrm{cov}(a)^{-1}
#'   (\bar a_T - \bar a_C),}
#' where \eqn{a} stacks the interference from [exposure()] and any extra
#' covariates `X`, \eqn{\bar a_T} and \eqn{\bar a_C} are its group means, and
#' \eqn{p} is the proportion treated. A Moore-Penrose inverse is used when
#' \eqn{\mathrm{cov}(a)} is singular. Under complete randomization \eqn{M} is
#' approximately \eqn{\chi^2} with as many degrees of freedom as columns of
#' \eqn{a}.
#'
#' @inheritParams exposure
#' @param X Optional numeric matrix or data frame of covariates to balance
#'   together with the interference.
#'
#' @return A non-negative number.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' A <- matrix(rbinom(n * n, 1, 0.1), n, n)
#' A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
#' imbalance(A, sample(rep(0:1, n / 2)))
#' @export
imbalance <- function(A, z, group = NULL, X = NULL) {
  E <- cbind(exposure(A, z, group), if (!is.null(X)) as.matrix(X))
  d <- colMeans(E[z == 1, , drop = FALSE]) - colMeans(E[z == 0, , drop = FALSE])
  p <- mean(z)
  length(z) * p * (1 - p) * drop(d %*% .pinv(stats::cov(E)) %*% d)
}

#' Number of matched pairs with identical interference
#'
#' Computes \eqn{|P_0|}, the largest number of one-to-one pairs of a treated
#' and a control unit whose interference vectors from [exposure()] are
#' identical. It equals \eqn{\sum_s \min(n_{T,s}, n_{C,s})} over the
#' distinct interference profiles \eqn{s}. [netmm()] maximizes it, and the
#' matched estimator in [dte()] averages over these pairs.
#'
#' @inheritParams exposure
#'
#' @return An integer.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' A <- matrix(rbinom(n * n, 1, 0.05), n, n)
#' A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
#' n_matched(A, sample(rep(0:1, n / 2)))
#' @export
n_matched <- function(A, z, group = NULL) {
  .npairs(exposure(A, z, group), z)
}
