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
#' PWD_RL(X, Y, sigma, kappa, lambda=1, epsilon=1e-6)
#'
#' @param X		the vector of predicate readings,
#' @param Y		the vector of test readings,
#' @param sigma		the RL sigma parameter,
#' @param kappa		the RL kappa parameter,
#' @param lambda		*optional* (default of 1) - the ratio of the X to
#' the Y precision profile,
#' @param epsilon		*optional* - convergence tolerance limit.
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
#'   \item{like }{the -2 log likelihood L}
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
#' @references Hawkins DM and Kraker JJ. Precision Profile Weighted Deming
#' Regression for Methods Comparison, on *Arxiv* (2025) <doi:10.48550/arXiv.2508.02888>
#'
#' @references Hawkins DM (2014). A Model for Assay Precision.
#' *Statistics in Biopharmaceutical Research*, **6**, 263-269.
#' <doi:10.1080/19466315.2014.899511>
#'
#' @importFrom stats optim
#'
#' @export

PWD_RL <- function(X, Y, sigma, kappa, lambda=1, epsilon=1e-6) {
  mu     <- X
  old    <- mu
  fity   <- mu
  diffr  <- 2*epsilon
  innr   <- 0
  error  <- ""
  best   <- 1e20
  beta   <- 1
  error  <- ""
  while(diffr > epsilon & innr < 26) {# weight depends on beta, so refine.
    innr  <- innr + 1
    g     <- lambda*(sigma^2 + (kappa*mu)^2)
    h     <- sigma^2 + (kappa*fity)^2
    w     <- 1/(h + beta^2 * g)
    sumw  <- sum(w)
    wsq   <- w^2
    xbar  <- sum(w * X) / sumw
    ybar  <- sum(w * Y) / sumw
    devx  <- X - xbar
    devy  <- Y - ybar
    sxxx  <- sum(wsq * devx^2 * g)
    sxxy  <- sum(wsq * devx^2 * h)
    sxyx  <- sum(wsq * devx * devy * g)
    sxyy  <- sum(wsq * devx * devy * h)
    syyx  <- sum(wsq * devy^2 * g)
    surd  <- (sxxy - syyx)^2 + 4 * sxyx * sxyy
    if (surd < 0) {
      error <- "Negative surd in WD_RL"
      break
    }
    beta   <- (syyx - sxxy + sqrt(surd)) / (2 * sxyx)
    alpha  <- ybar - beta *  xbar
    mu     <- w * (h * X + g * beta * (Y - alpha))
    fity   <- alpha + beta * mu
    resi   <- Y - alpha - beta*X
    L      <- sum((X-mu)^2/g + (Y-fity)^2/h + log(g*h))
    diffr <- sum((mu - old)^2) / sum(mu^2)
    old   <- mu
  }

  easyOK <- innr < 25
  if (!easyOK) {          # The fast method failed.  Use optim.

    getmu <- function(X, Y, sigma, kappa, lambda, alpha, beta, epsilon=1e-5) {
      mu <- X
      diffr <- 2*epsilon
      looper <- 0
      while(diffr > epsilon & looper < 100) {
        looper <- looper + 1
        old   <- mu
        fity  <- alpha + beta*mu
        g     <- lambda * (sigma^2 + (kappa*mu  )^2)
        h     <-           sigma^2 + (kappa*fity)^2
        mu    <- (h*X + g*beta*(Y-alpha))/(h + g*beta^2)
        diffr <- sum((mu-old)^2)/sum(mu^2)
      }

      fity  <- alpha + beta*mu
      g     <- lambda * (sigma^2 + (kappa*mu  )^2)
      h     <-           sigma^2 + (kappa*fity)^2

      return(list(mu=mu, g=g, h=h))
    }

    inner <- function(par) {
      alpha <- par[1]
      beta  <- par[2]
      gmu   <- getmu(X, Y, sigma, kappa, lambda, alpha, beta)
      mu    <- gmu$mu
      g     <- gmu$g
      h     <- gmu$h
      fity  <- alpha + beta * mu
      L <- sum((X-mu)^2/g + (Y-fity)^2/h + log(g*h))
      return(L)
    }

    dd    <- optim(c(0,1), inner)
    alpha <- dd$par[1]
    beta  <- dd$par[2]
    L     <- dd$value
    mu    <- getmu(X, Y, sigma, kappa, lambda, alpha, beta)$mu
    fity  <- alpha + beta*mu
    resi  <- Y - fity
  }
  return(list(alpha=alpha, beta=beta, fity=fity, mu=mu, resi=resi,
              like=L))
}
