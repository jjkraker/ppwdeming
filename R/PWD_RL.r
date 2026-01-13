#' Weighted Deming -- Rocke-Lorenzato - known sigma, kappa
#' @name PWD_RL
#'
#' @description
#' This code fits the weighted Deming regression on
#' predicate readings (X) and test readings (Y),
#' with user-supplied Rocke-Lorenzato ("RL") parameters
#' sigma (\eqn{\sigma}) and kappa (\eqn{\kappa}).
#'
#' @usage
#' PWD_RL(X, Y, sigma, kappa, lambda=1, alpha=NA, beta=NA, mu=NA, epsilon=1e-8)
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param sigma		the RL \eqn{\sigma} parameter.
#' @param kappa		the RL \eqn{\kappa} parameter.
#' @param lambda		*optional* (default of 1) - the ratio of the `X` to
#' the `Y` precision profile.
#' @param alpha     *optional* (default of NA) - numeric, single value, initial estimate of \eqn{\alpha}.
#' @param beta      *optional* (default of NA) - numeric, single value, initial estimate of \eqn{\beta}.
#' @param mu        *optional* (default of NA) - numeric, vector of length of `X`, initial estimate of \eqn{\mu}.
#' @param epsilon		*optional*  (default of 1e-8) - convergence tolerance limit.
#'
#' @details The Rocke-Lorenzato precision profile model assumes the following
#' forms for the variances, with proportionality constant \eqn{\lambda}:
#'    * predicate precision profile model: \eqn{g_i = var(X_i) = \lambda\left(\sigma^2 + \left[\kappa\cdot \mu_i\right]^2\right)} and
#'    * test precision profile model: \eqn{h_i = var(Y_i) = \sigma^2 + \left[\kappa\cdot (\alpha + \beta\mu_i)\right]^2}.
#'
#' The algorithm uses maximum likelihood estimation.  Proportionality constant
#' \eqn{\lambda} is assumed to be known or estimated externally.
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha }{the fitted intercept}
#'   \item{beta }{the fitted slope}
#'   \item{fity }{the vector of predicted Y}
#'   \item{mu }{the vector of estimated latent true values}
#'   \item{resi }{the vector of residuals}
#'   \item{L }{the -2 log likelihood L}
#'   \item{innr }{the number of inner refinement loops executed}
#'   \item{error }{an error code if the iteration fails}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @examples
#' # library
#' library(ppwdeming)
#'
#' # parameter specifications
#' sigma <- 1
#' kappa <- 0.08
#' alpha <- 1
#' beta  <- 1.1
#' true  <- 8*10^((0:99)/99)
#' truey <- alpha+beta*true
#' # simulate single sample - set seed for reproducibility
#' set.seed(1039)
#' # specifications for predicate method
#' X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(100))
#' # specifications for test method
#' Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(100))
#'
#' # fit RL with given sigma and kappa
#' RL_results <- PWD_RL(X,Y,sigma,kappa)
#' cat("\nWith given sigma and kappa, the estimated intercept is",
#'     signif(RL_results$alpha,4), "and the estimated slope is",
#'     signif(RL_results$beta,4), "\n")
#'
#' @references Hawkins DM and Kraker JJ (in press). Precision Profile Weighted
#' Deming Regression for Methods Comparison. *The Journal of Applied Laboratory Medicine*.
#' <doi:10.1093/jalm/jfaf183>
#'
#' @references Hawkins DM (2014). A Model for Assay Precision.
#' *Statistics in Biopharmaceutical Research*, **6**, 263-269.
#' <doi:10.1080/19466315.2014.899511>
#'
#' @importFrom stats optim
#' @importFrom stats complete.cases
#'
#' @export

PWD_RL <- function(X, Y, sigma, kappa, lambda=1, alpha=NA, beta=NA, mu=NA, epsilon=1e-8) {
  whichmissing <- (!complete.cases(X)) | (!complete.cases(Y))
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  allY <- Y
  X <- X[!whichmissing]
  Y <- Y[!whichmissing]
  if(sum(!is.na(mu)) > 0) mu <- mu[!whichmissing]

  getmu <- function(X, Y, alpha, beta, g, h, mu, epsilon=1e-8) {
    diffr <- 2*epsilon
    innr  <- 0
    while (diffr > epsilon & innr < 100) {
      innr <- innr + 1
      old  <- mu
      fity <- alpha + beta*mu
      mu   <- (h * X + g * beta * (Y - alpha))/(h + g * beta^2)
      diffr <- sum((mu - old)^2)/sum(mu^2)
    }
    return(list(mu=mu, innr=innr))
  }

  # calculates L from alpha, beta, rho
  qform   <- function(par) {
    alpha <- par[1]
    beta  <- par[2]
    diffr <- 2*epsilon
    refine <- 0
    while(diffr > epsilon & refine < 100) {
      refine <- refine+1
      old <- mu
      fity  <- alpha + beta*mu
      g     <- lambda * (rho^2 + mu^2)
      h     <-           rho^2 + fity^2
      mu    <- getmu(X, Y, alpha, beta, g, h, mu)$mu
      diffr <- sum((old-mu)^2)/sum(mu^2)
    }
    W     <- sum((X-mu)^2/g+(Y-fity)^2/h)
    slgh  <- sum(log(g*h))
    kappa <- sqrt(W/tun)
    sigma <- rho*kappa
    L     <- tun*log(W) + slgh + A
    return(list(L=L, W=W, sigma=sigma, kappa=kappa, alpha=alpha,
                beta=beta, mu=mu, fity=fity))
  }

  # Wrapper
  wrapqform <- function(par) {
    do <- qform(par)
    do$L
  }

  rho   <- sigma/(kappa+1e-6)
  if(is.na(alpha))    alpha <- 0
  if(is.na(beta ))    beta  <- 1
  if (is.na(sum(mu))) mu    <- X

  n     <- length(X)
  tun   <- 2*n
  A     <- tun * (1-log(tun))
  par   <- c(alpha, beta)

  doit  <- optim(par, wrapqform)

  pars  <- doit$par
  wrap  <- qform(pars)

  resi <- Y - wrap$fity

  allresi = rep(NA, length(allX))
  allresi[!whichmissing] = resi
  allfity = rep(NA, length(allX))
  allfity[!whichmissing] = wrap$fity
  allmu = rep(NA, length(allX))
  allmu[!whichmissing] = wrap$mu

  return(list(alpha = pars[1], beta = pars[2], fity = allfity,
              mu = allmu, resi = allresi, L = wrap$L))
}
