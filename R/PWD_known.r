#' Weighted Deming Regression -- general weights
#' @name PWD_known
#'
#' @description
#' This code is used for the setting of known precision profiles implemented
#' in user-provided R functions called `gfun` and `hfun`.
#'
#' @usage
#' PWD_known(X, Y, gfun, hfun, gparms, hparms, epsilon=1e-10,
#'           MDL=NA, getCI=TRUE, printem=TRUE)
#'
#' @param X		the vector of predicate readings,
#' @param Y		the vector of test readings,
#' @param gfun	 a function with two arguments, a vector of size *n* and a vector of parameters,
#' @param hfun	 a function with two arguments, a vector of size *n* and a vector of parameters,
#' @param gparms		 a numeric vector containing any parameters referenced by `gfun`,
#' @param hparms		 a numeric vector containing any parameters referenced by `hfun`,
#' @param MDL		*optional* medical decision level(s),
#' @param getCI		*optional* - allows for jackknifed standard errors on the regression and MDL,
#' @param epsilon		*optional* convergence tolerance limit,
#' @param printem	  *optional* - if TRUE, routine will print out results.
#'
#' @details The functions `gfun` and `hfun` are allowed as inputs,
#' to support flexibility in specification of the forms of these variance functions.
#' The known precision profiles specified by the functions `gfun` and `hfun`,
#' when provided with estimated vectors of \eqn{\mu} and \eqn{\alpha + \beta\mu}
#' respectively and with any required parameters, will produce
#' the vectors g and h.  These vectors are then integrated into the
#' iterative estimation of the slope and intercept of the linear relationship
#' between predicate and test readings.
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha }{the fitted intercept}
#'   \item{beta }{the fitted slope}
#'   \item{cor }{the Pearson correlation between X and Y}
#'   \item{fity }{the vector of predicted Y}
#'   \item{mu }{the vector of estimated latent true values}
#'   \item{resi }{the vector of residuals}
#'   \item{scalr }{the vector of scaled residuals using the specified g and h}
#'   \item{like }{the -2 log likelihood L}
#'   \item{sealpha }{the jackknife standard error of alpha}
#'   \item{sebeta }{the jackknife standard error of beta}
#'   \item{covar }{the jackknife covariance between alpha and beta}
#'   \item{preMDL }{the predictions at the MDL(s)}
#'   \item{preMDLl }{the lower confidence limit(s) of preMDL}
#'   \item{preMDLu }{the upper confidence limit(s) of preMDL}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @example /Example/PWD_known_man_example.R
#'
#' @export

PWD_known <- function(X, Y, gfun, hfun, gparms, hparms, epsilon=1e-10, MDL=NA, getCI=TRUE, printem=TRUE) {

  alpha      <- 0
  beta       <- 1
  fullmu     <- X
  pseudalpha <- NULL
  pseudbeta  <- NULL
  n          <- length(X)
  top        <- 0
  sealpha    <- NA
  sebeta     <- NA
  covar      <- NA
  if (getCI) top <- n
  for (dele in 0:top) {
    tx  <- X
    ty  <- Y
    mu <- fullmu
    if (dele > 0) {
      tx  <- X[-dele]
      ty  <- Y[-dele]
      mu <- fullmu[-dele]
    }
    old  <- mu
    diff <- 2*epsilon
    while(diff > epsilon) {# weight depends on beta, so refine.
      fity  <- alpha + beta * mu
      g     <- gfun(mu  , gparms)
      h     <- hfun(fity, hparms)
      w     <- 1/(h + beta^2 * g)
      sumw  <- sum(w)
      wsq   <- w^2
      xbar  <- sum(w * tx) / sumw
      ybar  <- sum(w * ty) / sumw
      devx  <- tx - xbar
      devy  <- ty - ybar
      sxxx  <- sum(wsq * devx^2 * g)
      sxxy  <- sum(wsq * devx^2 * h)
      sxyx  <- sum(wsq * devx * devy * g)
      sxyy  <- sum(wsq * devx * devy * h)
      syyx  <- sum(wsq * devy^2 * g)
      surd  <- (sxxy - syyx)^2 + 4 * sxyx * sxyy
      beta  <- (syyx - sxxy + sqrt(surd)) / (2 * sxyx)
      alpha <- ybar - beta *  xbar
      mu    <- w * (h * tx + g * beta * (ty - alpha))
      diff  <- sum((mu - old)^2) / sum(mu^2)
      old   <- mu
    }
    if (dele == 0) {
      fullalpha <- alpha
      fullbeta  <- beta
      fullmu    <- mu
      fullxbar  <- xbar
      fullybar  <- ybar
      fity      <- alpha + beta * mu
      resi     <- Y - alpha - beta*X
      profl     <- sqrt(g + beta^2*h)
      scalr     <- resi / profl
      like <- sum((X-mu)^2/g + (Y-alpha-beta*mu)^2/h + log(g*h))
    } else {
      pseudalpha <- c(pseudalpha, n*fullalpha - (n-1)*alpha)
      pseudbeta  <- c(pseudbeta , n*fullbeta  - (n-1)*beta )
    }
  }
  alpha  <- fullalpha
  beta   <- fullbeta
  if (getCI) {
    sealpha <- sd(pseudalpha) / sqrt(n)
    sebeta  <- sd(pseudbeta ) / sqrt(n)
    covar   <- sealpha * sebeta * cor(pseudalpha, pseudbeta)
  }
  nMDL    <- 0
  preMDL  <- NA
  preMDLl <- NA
  preMDLu <- NA
  if (!is.na(sum(MDL))) {
    nMDL    <- length(MDL)
    preMDL  <- alpha + beta*MDL
    if (getCI) {
      tcut    <- qt(0.975, n-1)
      MoEpre  <- tcut*sqrt(sealpha^2 + (sebeta*MDL)^2 + 2*covar*MDL)
      preMDLl <- preMDL - MoEpre
      preMDLu	<- preMDL + MoEpre
    }
    #	print(data.frame(preMDL, preMDLl, preMDLu))
  }

  if (printem) {
    cat(sprintf("%9s %8s %8s %9s\n", "Parameter", "estimate", "se", "CI"))
    tcut <- qt(0.975, n-1)
    CI   <- fullalpha + tcut * sealpha * c(-1,1)
    cat(sprintf("Intercept %8.3f %8.3f (%7.3f, %6.3f)\n", fullalpha, sealpha, CI[1], CI[2]))
    CI   <- fullbeta + tcut * sebeta * c(-1,1)
    cat(sprintf("slope     %8.3f %8.3f (%7.3f, %6.3f)\n", fullbeta , sebeta, CI[1], CI[2]))
    if (nMDL > 0) {
      for (kk in 1:nMDL) {
        cat(sprintf("MDL %7.3f prediction %7.3f CI %7.3f %7.3f\n",
                    MDL[kk], preMDL[kk], preMDLl[kk], preMDLu[kk]))
      }
    }
  }
  corXY = cor(X,Y)

  return(list(alpha=fullalpha, beta=fullbeta, cor=corXY, fity=fity, mu=fullmu, resi=resi,
              scalr=scalr, like=like, sealpha=sealpha, sebeta=sebeta, covar=covar,
              preMDL=preMDL, preMDLl=preMDLl, preMDLu=preMDLu))
}

