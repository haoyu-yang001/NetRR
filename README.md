# NetRR

<!-- badges: start -->
[![R-CMD-check](https://github.com/haoyu-yang001/NetRR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haoyu-yang001/NetRR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Randomization designs for experiments on a network, from Yang, Qin, Yu and Li (2026), "Design Strategies for Networked Experiments via Interference Balancing", *Journal of Business & Economic Statistics*, [doi:10.1080/07350015.2026.2658067](https://doi.org/10.1080/07350015.2026.2658067).

When units are connected, a unit's outcome depends on its own treatment and also on how many of its neighbours are treated (its *interference*). If the treatment and control groups differ in interference, the difference in means for the direct treatment effect is biased and noisy. NetRR offers two designs that balance interference:

* **NetRR** (`netrr()`): network interference balancing rerandomization. It keeps drawing complete randomizations until the Mahalanobis distance between the groups' mean interference is small.
* **NetMM** (`netmm()`): network maximized matching randomization. It searches for an assignment with as many treated-control pairs of identical interference as possible.

Neighbours can be split into groups with different peer effects (`group`), and NetRR can balance extra covariates at the same time (`X`).

## Installation

Requires R >= 4.0.0 and has no dependencies beyond base R.

```r
# install.packages("remotes")
remotes::install_github("haoyu-yang001/NetRR")
```

## Example

```r
library(NetRR)
set.seed(1)
n <- 200
A <- matrix(rbinom(n * n, 1, 0.1), n, n)
A <- 1 * ((A + t(A)) > 0); diag(A) <- 0      # undirected network

# NetRR: accept the most balanced 5% of complete randomizations
z <- netrr(A, pa = 0.05)
Y <- 0.5 + 2 * z + 1.5 * exposure(A, z)[, 1] + rnorm(n)
dte(Y, z)                     # difference in means
dte(Y, z, A, method = "lm")   # regression on the interference

# NetMM: maximize the number of matched pairs, then use the matched estimator
z <- netmm(A, maxit = 2000)
Y <- 0.5 + 2 * z + 1.5 * exposure(A, z)[, 1] + rnorm(n)
dte(Y, z, A, method = "match")
```

| Function | Purpose |
| :--- | :--- |
| `netrr()` | Network interference balancing rerandomization |
| `netmm()` | Network maximized matching randomization |
| `dte()` | Direct effect: difference in means, regression, or matched pairs |
| `exposure()` | Treated neighbours of each unit, per neighbour group |
| `imbalance()` | Mahalanobis distance of interference (and covariates) between arms |
| `n_matched()` | Number of treated-control pairs with identical interference |

## Citation

```r
citation("NetRR")
```

The code and data used in the paper are in its online supplementary material.
