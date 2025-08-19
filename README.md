
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ppwdeming

<!-- badges: start -->

[![R-CMD-check](https://github.com/jjkraker/ppwdeming/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jjkraker/ppwdeming/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ppwdeming is to provide functions for weighted Deming
regression, using weights modeled via precision profile (used commonly
in the realm of clinical chemistry). Functions are included for
implementing weights in situations of known and unknown precision
profile settings.

## Installation

You can install the development version of `ppwdeming` like so:

``` r
install.packages("ppwdeming")  # once available on CRAN
```

## Example

This is a basic example which shows you how to run the main functions:

``` r
# library
library(ppwdeming)

# parameter specifications
sigma <- 1
kappa <- 0.08
alpha <- 1
beta  <- 1.1
true  <- 8*10^((0:99)/99)
truey <- alpha+beta*true
# simulate single sample - set seed for reproducibility
set.seed(1039)
# specifications for predicate method
X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(100))
# specifications for test method
Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(100))

# fit RL with given sigma and kappa
RL_results <- PWD_RL(X,Y,sigma,kappa)
cat("\nWith given sigma and kappa, the estimated intercept is",
    signif(RL_results$alpha,4), "and the estimated slope is",
    signif(RL_results$beta,4), "\n")

# fit with RL precision profile to estimate parameters
RL_gh_fit  <- PWD_get_gh(X,Y,printem=TRUE)
# RL precision profile estimated parameters
cat("\nsigmahat=", signif(RL_gh_fit$sigma,6),
    "and kappahat=", signif(RL_gh_fit$kappa,6))

# run the residual analysis from the model output
post  <- PWD_resi(X, RL_gh_fit$resi, printem=TRUE)

# fit with RL precision profile to estimate parameters and variability
RL_inf <- PWD_inference(X,Y,MDL=12,printem=TRUE)
```

along with the outlier review:

``` r
# add some outliers
Y[c(1,2,100)] <- Y[c(1,2,100)] + c(7,4,-45)

# check for outliers, re-fit, and store output
outliers_assess <- PWD_outlier(X,Y,K=5)
```

An alternative example in which the precision profiles are known:

``` r
# parameter specifications
alpha <- 1
beta  <- 1.1
true  <- 8*10^((0:99)/99)
truey <- alpha+beta*true
# forms of precision profiles
gfun    <- function(true, gparms) {
  gvals = gparms[1]+gparms[2]*true^gparms[3]
  gvals
}
hfun    <- function(true, hparms) {
  hvals = hparms[1]+hparms[2]*true^hparms[3]
  hvals
}

# Loosely motivated by Vitamin D data set
g     <- 4e-16+0.07*true^1.27
h     <- 6e-2+7e-5*truey^2.2
# simulate single sample - set seed for reproducibility
set.seed(1039)
# specifications for predicate method
X     <- true +sqrt(g)*rnorm(100)
# specifications for test method
Y     <- truey+sqrt(h)*rnorm(100)

# fit with to estimate linear parameters
pwd_known_fit <- PWD_known(X, Y, gfun, hfun,
                           c(4e-16, 0.07, 1.27), c(6e-2, 7e-5, 2.2))
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.
`devtools::build_readme()` is handy for this. -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
