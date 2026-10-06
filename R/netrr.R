#' Network interference balancing rerandomization (NetRR)
#'
#' Draws complete randomizations with `n1` treated units until the
#' interference imbalance [imbalance()] falls below a threshold.
#'
#' @details
#' The threshold is the `pa` quantile of the imbalance over `nref` complete
#' randomizations, so roughly a fraction `pa` of assignments is acceptable.
#' The first acceptable draw is returned. If none is found within `maxit`
#' draws, the most balanced draw seen is returned with a warning.
#'
#' @inheritParams imbalance
#' @param n1 Number of treated units.
#' @param pa Acceptance probability, between 0 and 1. Smaller values give
#'   better balance but need more draws.
#' @param nref Number of complete randomizations used to set the threshold.
#' @param maxit Maximum number of draws.
#'
#' @return A 0/1 assignment vector with attributes `imbalance` (its
#'   imbalance) and `threshold`.
#'
#' @references Yang, H., Qin, Y., Yu, Y., and Li, Y. (2026). Design
#'   strategies for networked experiments via interference balancing.
#'   *Journal of Business & Economic Statistics*.
#'   \doi{10.1080/07350015.2026.2658067}
#'
#' @seealso [netmm()], [dte()]
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' A <- matrix(rbinom(n * n, 1, 0.1), n, n)
#' A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
#' z <- netrr(A, pa = 0.05, nref = 500)
#' attr(z, "imbalance")
#'
#' # heterogeneous interference over two groups, plus a covariate
#' g <- rep(1:2, n / 2)
#' x <- rnorm(n)
#' z <- netrr(A, group = g, X = x, pa = 0.05, nref = 500)
#' @export
netrr <- function(A, group = NULL, X = NULL, n1 = floor(nrow(A) / 2),
                  pa = 0.1, nref = 1000, maxit = 10000) {
  n <- .check_A(A)
  .check_n1(n1, n)
  if (pa <= 0 || pa >= 1) stop("`pa` must be strictly between 0 and 1.", call. = FALSE)
  M <- function(z) imbalance(A, z, group, X)
  thr <- stats::quantile(replicate(nref, M(.draw(n, n1))), pa, type = 1, names = FALSE)
  best <- NULL
  bestM <- Inf
  for (i in seq_len(maxit)) {
    z <- .draw(n, n1)
    m <- M(z)
    if (m < bestM) {
      best <- z
      bestM <- m
    }
    if (m <= thr) break
  }
  if (bestM > thr) warning("No assignment met the threshold; returning the best one found.", call. = FALSE)
  structure(best, imbalance = bestM, threshold = thr)
}
