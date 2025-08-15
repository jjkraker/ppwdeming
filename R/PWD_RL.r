#' Weighted Deming Regression -- with input Rocke-Lorenzato parameters
#' @name PWD_RL
#'
#' @description
#' This code fits the weighted Deming regression on
#' predicate readings (X) and test readings (Y),
#' with user-supplied Rocke-Lorenzato ("RL") parameters
#' sigma (\eqn{\sigma}) and kappa (\eqn{\kappa}).
#'
#' @usage
#' PWD_RL(X, Y, sigma, kappa, lambda=1, epsilon=1e-10)
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
#' @example /Examples/PWD_RL_man_example.R
#'
#' @references Hawkins DM and Kraker JJ. Precision Profile Weighted Deming
#' Regression for Methods Comparison, on *Arxiv* (2025, <arxiv.org/abs/2508.02888>).
#'
#' @references Hawkins DM (2014). A Model for Assay Precision.
#' *Statistics in Biopharmaceutical Research*, **6**, 263-269.
#' http://dx.doi.org/10.1080/19466315.2014.899511
#'
#' @export

PWD_RL <- function(X, Y, sigma, kappa, lambda=1, epsilon=1e-10) {
  mu     <- X
  old    <- mu
  fity   <- mu
  diffr  <- 2*epsilon
  innr   <- 0
  error  <- ""
  best   <- 1e20
  beta   <- 1
  allbet <- NULL
  alllik <- NULL
  while((innr < 5) | (diffr > epsilon & innr < 100)) {# weight depends on beta, so refine.
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
    like   <- sum((X-mu)^2/g + (Y-fity)^2/h + log(g*h))
    allbet <- c(allbet, round(beta,4))
    alllik <- c(alllik, like)
    if (like < best) {
      best   <- like
      besbet <- beta
      besalp <- alpha
      besmu  <- mu
      besg   <- g
      besh   <- h
      besfit <- fity
      besres <- resi
      when   <- innr
    }
    diffr <- sum((mu - old)^2) / sum(mu^2)
    old   <- mu
  }
  return(list(alpha=besalp, beta=besbet, fity=besfit, mu=besmu, resi=besres,
              like=best, innr=innr, error=error))
}
