#' Estimate of Variance Profile Functions (proportional)
#' @name PWD_get_gh
#'
#' @description
#' This code estimates the variance profiles, assumed proportional, of
#' the Rocke-Lorenzato form;
#' also provides the resulting weighted Deming fit and residuals.
#'
#' @usage
#' PWD_get_gh(X, Y, lambda=1, epsilon=1.e-8, printem=TRUE)
#'
#' @param X		the vector of predicate readings,
#' @param Y		the vector of test readings,
#' @param lambda		*optional* (default of 1) - the ratio of the X to the Y precision profile (defaults to 1),
#' @param epsilon		*optional* (default of 1.e-8) - convergence tolerance limit,
#' @param printem	  *optional* - if TRUE, routine will print out results.
#'
#' @details
#' This workhorse routine optimizes the likelihood in the **unknown** *g*, *h*
#' setting over its *n*+4 parameters
#' (the two Rocke-Lorenzato precision profile parameters \eqn{\sigma}
#' and \eqn{\kappa}, the intercept \eqn{\alpha} and slope \eqn{\beta},
#' and the *n* latent true concentrations \eqn{\mu_i}).
#'
#' That is, the assumed forms are:
#'    * predicate precision profile model: \eqn{g_i = var(X_i) = \lambda\left(\sigma^2 + \left[\kappa\cdot \mu_i\right]^2\right)} and
#'    * test precision profile model: \eqn{h_i = var(Y_i) = \sigma^2 + \left[\kappa\cdot (\alpha + \beta\mu_i)\right]^2}.
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha }{the fitted intercept}
#'   \item{beta }{the fitted slope}
#'   \item{fity }{the vector of predicted Y}
#'   \item{mu }{the vector of estimated latent true values}
#'   \item{resi }{the vector of residuals}
#'   \item{sigma }{the estimate of the Rocke-Lorenzato \eqn{\sigma}}
#'   \item{kappa }{the estimate of the Rocke-Lorenzato \eqn{\kappa}}
#'   \item{like }{the -2 log likelihood L}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @example /Examples/PWD_get_gh_man_example.R
#'
#' @references Hawkins DM and Kraker JJ. Precision Profile Weighted Deming
#' Regression for Methods Comparison, on *Arxiv* (2025, <arxiv.org/abs/2508.02888>).
#'
#' @references Rocke DM, Lorenzato S (2012). A Two-Component Model for Measurement
#' Error in Analytical Chemistry.  *Technometrics*, **37:2**:176-184.
#'
#' @importFrom stats optimize
#'
#' @export

PWD_get_gh <- function(X, Y, lambda=1, epsilon=1.e-8, printem=TRUE) {

  key  <- order(X)
  sX   <- X[key]
  sY   <- Y[key]

  # Establish search ranges
  n    <- length(X)
  lowr <- 1:round(n/3)
  hir  <- round(2*n/3):n
  maxs <- max(abs(sY[lowr]-sX[lowr]))              # sigma
  maxk <- max(abs(log(sY[hir]/sX[hir])), na.rm=T)  # kappa
  #  print(paste("Entering PWD_get_GH, lambda=", lambda, maxs, maxk))
  gr   <- 1.618
  a    <- 0
  b    <- maxs
  fity <- X

  innr <- function(kappa, sigma, lambda) {
    do <- PWD_RL(X, Y, sigma, kappa, lambda)
    return(like=do$like)
  }

  h    <- maxs
  while (h > epsilon) {
    h        <- b - a
    c        <- b - h / gr
    d        <- a + h / gr
    vc       <- optimize(innr, c(0,maxk), c, lambda)
    fc       <- vc$objective
    vd       <- optimize(innr, c(0,maxk), d, lambda)
    fd       <- vd$objective
    if (fc < fd) {
      b <- d
    } else {
      a <- c
    }
  }
  sigma <- c
  kappa <- vc$minimum
  like    <- fc
  if (fd < fc) {
    sigma <- d
    kappa <- vd$minimum
    like    <- fd
  }
  #  print("Done optimizing")
  do <- PWD_RL(X, Y, sigma, kappa, lambda)
  if (printem) {
    with(do, cat(sprintf("alpha %6.3f beta %5.3f like %8.5g \n", alpha, beta, like)))
  }
  return(list(alpha=do$alpha, beta=do$beta, fity=do$fity, mu=do$mu, resi=do$resi,
              sigma=sigma, kappa=kappa, like=do$like))
}
