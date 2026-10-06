# Moore-Penrose pseudo-inverse of a symmetric matrix
.pinv <- function(S, tol = 1e-8) {
  s <- svd(S)
  k <- s$d > tol * s$d[1]
  s$v[, k, drop = FALSE] %*% (t(s$u[, k, drop = FALSE]) / s$d[k])
}

# n x K indicator matrix of the neighbour groups
.onehot <- function(group, n) {
  if (is.null(group)) group <- rep(1L, n)
  if (length(group) != n) stop("`group` must have one entry per node.", call. = FALSE)
  lev <- sort(unique(group))
  G <- outer(group, lev, "==") + 0
  colnames(G) <- as.character(lev)
  G
}

.check_A <- function(A) {
  if (length(dim(A)) != 2 || nrow(A) != ncol(A)) {
    stop("`A` must be a square adjacency matrix.", call. = FALSE)
  }
  nrow(A)
}

.check_z <- function(z, n) {
  if (length(z) != n || !all(z %in% 0:1)) {
    stop("`z` must be a 0/1 vector with one entry per node.", call. = FALSE)
  }
  if (sum(z) == 0 || sum(z) == n) {
    stop("Both arms need at least one unit.", call. = FALSE)
  }
}

.check_n1 <- function(n1, n) {
  if (length(n1) != 1 || n1 != round(n1) || n1 < 1 || n1 > n - 1) {
    stop("`n1` must be an integer between 1 and nrow(A) - 1.", call. = FALSE)
  }
}

# complete randomization with n1 treated units
.draw <- function(n, n1) sample(rep(1:0, c(n1, n - n1)))

# k random elements of x (safe when length(x) == 1)
.pick <- function(x, k) x[sample.int(length(x), k)]

# one key per unit identifying its interference profile
.key <- function(E) apply(round(E, 8), 1, paste, collapse = ",")

# |P0|: largest number of one-to-one treated-control pairs with equal interference
.npairs <- function(E, z) {
  tab <- table(.key(E), factor(z, levels = 0:1))
  sum(pmin(tab[, 1], tab[, 2]))
}
