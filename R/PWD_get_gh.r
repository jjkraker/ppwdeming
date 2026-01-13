#' Estimate of Variance Profile Functions (proportional)
#' @name PWD_get_gh
#'
#' @description
#' This code estimates the variance profiles, assumed proportional, of
#' the Rocke-Lorenzato form;
#' also provides the resulting weighted Deming fit and residuals.
#'
#' @usage
#' PWD_get_gh(X, Y, lambda = 1, rho=NA, alpha=NA, beta=NA, mu=NA,
#'            epsilon = 1e-8, printem=FALSE)
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param lambda		*optional* (default of 1) - the ratio of the `X` to
#' the `Y` precision profile.
#' @param rho       *optional* (default of NA) - numeric, single value or vector, initial estimate(s) of \eqn{\rho = \frac{\sigma}{\kappa}}.
#' @param alpha     *optional* (default of NA) - numeric, single value, initial estimate of \eqn{\alpha}.
#' @param beta      *optional* (default of NA) - numeric, single value, initial estimate of \eqn{\beta}.
#' @param mu        *optional* (default of NA) - numeric, vector of length of `X`, initial estimate of \eqn{\mu}.
#' @param epsilon		*optional* (default of 1.e-8) - convergence tolerance limit.
#' @param printem	  *optional* (default of FALSE) - if TRUE, routine will print out results as a `message`.
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
#' The search algorithm implements an efficient approach via reparameterization
#' to the ratio \eqn{\rho = \frac{\sigma}{\kappa}}.
#'
#' If initial estimates are not provided, the parameters are initialized as:
#'    * `alpha` and `beta` are initially intercept and slope from simple linear regression;
#'    * `rho` is initialized as the vector c(0.01, 1, 100); and
#'    * `mu` is initialized as the values of `X`.
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
#'   \item{L }{the -2 log likelihood L}
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
#' # fit with RL precision profile to estimate parameters
#' RL_gh_fit  <- PWD_get_gh(X,Y,printem=TRUE)
#' # RL precision profile estimated parameters
#' cat("\nsigmahat=", signif(RL_gh_fit$sigma,6),
#'     "and kappahat=", signif(RL_gh_fit$kappa,6), "\n")
#' # with estimated linear coefficients
#' cat("\nalphahat=", signif(RL_gh_fit$alpha,6),
#'     "and betahat=", signif(RL_gh_fit$beta,6), "\n")
#'
#' @references Hawkins DM and Kraker JJ (in press). Precision Profile Weighted
#' Deming Regression for Methods Comparison. *The Journal of Applied Laboratory Medicine*.
#' <doi:10.1093/jalm/jfaf183>
#'
#' @references Rocke DM, Lorenzato S (2012). A Two-Component Model for Measurement
#' Error in Analytical Chemistry.  *Technometrics*, **37:2**:176-184.
#'
#' @importFrom stats optimize
#' @importFrom stats lm
#' @importFrom stats coef
#'
#' @export

PWD_get_gh <- function (X, Y, lambda = 1, rho=NA,
                        alpha=NA, beta=NA, mu=NA, epsilon = 1e-8, printem=FALSE) {
  whichmissing <- (!complete.cases(X)) | (!complete.cases(Y))
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  allY <- Y
  X <- X[!whichmissing]
  Y <- Y[!whichmissing]
  if(sum(!is.na(mu)) > 0) mu <- mu[!whichmissing]

  # inner function takes old mu, alpha, beta, g, h and gets new mu
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

  # calculates L from alpha, beta, sigma, kappa
  qform   <- function(par) {
    rho   <- par[1]
    alpha <- par[2]
    beta  <- par[3]
    diffr <- 2*epsilon
    refine <- 0
    while(diffr > epsilon & refine < 100) {
      refine <- refine+1
      old <- mu
      fity  <- alpha + beta*mu
      resi  <- Y - alpha - beta*X
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
    L     <- L - 1000 * min(c(0, rho))
    return(list(L=L, W=W, sigma=sigma, kappa=kappa, alpha=alpha,
                beta=beta, mu=mu, fity=fity, resi=resi))
  }

  # Wrapper
  wrapqform <- function(par) {
    do <- qform(par)
    do$L
  }

  # Preliminary ranging
  if (is.na(alpha + beta)) {
    #generate starting values
    fitlm <- lm(Y~X)
    alpha <- coef(fitlm)[1]
    beta  <- coef(fitlm)[2]
  }
  if (any(is.na(rho)))     rho   <- c(0.01, 1, 100)
  if (is.na(sum(mu))) mu    <- X

  n     <- length(X)
  tun   <- 2*n
  A     <- tun * (1-log(tun))

  nrho  <- length(rho)
  best  <- 1e10
  for (mm in 1:nrho) {
    par   <- c(rho[mm], alpha, beta)
    doit  <- optim(par, wrapqform)
    L     <- doit$value
    if (L < best) {
      best <- L
      bestpars <- doit$par
    }
  }
  wrap <- qform(bestpars)
  resi <- Y - bestpars[2] - bestpars[3]*X

  allresi = rep(NA, length(allX))
  allresi[!whichmissing] = resi
  allfity = rep(NA, length(allX))
  allfity[!whichmissing] = wrap$fity
  allmu = rep(NA, length(allX))
  allmu[!whichmissing] = wrap$mu

  return(list(alpha = bestpars[2], beta = bestpars[3], fity = allfity,
              mu = allmu, resi = allresi, rho=wrap$sigma/wrap$kappa, sigma = wrap$sigma,
              kappa = wrap$kappa, L = wrap$L))
}
