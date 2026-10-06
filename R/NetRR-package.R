#' NetRR: Design of Networked Experiments via Interference Balancing
#'
#' Implements the randomization designs of Yang, Qin, Yu and Li (2026) for
#' estimating the direct treatment effect when units interfere through a
#' network.
#'
#' Main functions:
#' \itemize{
#'   \item [netrr()]: network interference balancing rerandomization.
#'   \item [netmm()]: network maximized matching randomization.
#'   \item [dte()]: difference-in-means, regression and matched estimators.
#'   \item [exposure()], [imbalance()], [n_matched()]: interference, its
#'     imbalance, and the number of matched pairs.
#' }
#'
#' @inherit netrr references
#' @keywords internal
"_PACKAGE"
