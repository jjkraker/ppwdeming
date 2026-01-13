#' Weighted Deming Regression -- general weights
#' @name PWD_known
#'
#' @description
#' This code is used for the setting of known precision profiles implemented
#' in user-provided R functions called `gfun` and `hfun`.
#'
#' @usage
#' PWD_known(X, Y, gfun, hfun, gparms, hparms, epsilon=1e-8,
#'           MDL=NA, getCI=TRUE, printem=FALSE)
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param gfun	 a function with two arguments, a vector of size *n* and a vector of parameters.
#' @param hfun	 a function with two arguments, a vector of size *n* and a vector of parameters.
#' @param gparms		 a numeric vector containing any parameters referenced by `gfun`.
#' @param hparms		 a numeric vector containing any parameters referenced by `hfun`.
#' @param MDL		*optional*  (default of NA) - medical decision level(s).
#' @param getCI		*optional*  (default of TRUE) - allows for jackknifed standard errors on the regression and MDL.
#' @param epsilon		*optional*  (default of 1.e-8) - convergence tolerance limit.
#' @param printem	  *optional*  (default of FALSE) - if TRUE, routine will print out results as a `message`.
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
#'   \item{L }{the -2 log likelihood L}
#'   \item{sealpha }{the jackknife standard error of alpha}
#'   \item{sebeta }{the jackknife standard error of beta}
#'   \item{covar }{the jackknife covariance between alpha and beta}
#'   \item{preMDL }{the predictions at the MDL(s)}
#'   \item{preMDLl }{the lower confidence limit(s) of preMDL}
#'   \item{preMDLu }{the upper confidence limit(s) of preMDL}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @examples
#' # library
#' library(ppwdeming)
#'
#' # parameter specifications
#' alpha <- 1
#' beta  <- 1.1
#' true  <- 8*10^((0:99)/99)
#' truey <- alpha+beta*true
#' # forms of precision profiles
#' gfun    <- function(true, gparms) {
#'   gvals = gparms[1]+gparms[2]*true^gparms[3]
#'   gvals
#' }
#' hfun    <- function(true, hparms) {
#'   hvals = hparms[1]+hparms[2]*true^hparms[3]
#'   hvals
#' }
#'
#' # Loosely motivated by Vitamin D data set
#' g     <- 4e-16+0.07*true^1.27
#' h     <- 6e-2+7e-5*truey^2.2
#' # simulate single sample - set seed for reproducibility
#' set.seed(1039)
#' # specifications for predicate method
#' X     <- true +sqrt(g)*rnorm(100)
#' # specifications for test method
#' Y     <- truey+sqrt(h)*rnorm(100)
#'
#' # fit with to estimate linear parameters
#' pwd_known_fit <- PWD_known(X, Y, gfun, hfun,
#'                            gparms=c(4e-16, 0.07, 1.27),
#'                            hparms=c(6e-2, 7e-5, 2.2), MDL=12,
#'                            printem=TRUE)
#'
#' @importFrom stats complete.cases
#'
#' @export

PWD_known <- function(X, Y, gfun, hfun, gparms, hparms, epsilon=1e-8, MDL=NA, getCI=TRUE, printem=FALSE) {
  whichmissing <- (!complete.cases(X)) | (!complete.cases(Y))
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  allY <- Y
  X <- X[!whichmissing]
  Y <- Y[!whichmissing]

  getmu <- function(tx, ty, alpha, beta, g, h, tmu, epsilon=1e-8) {
    diffr <- 2*epsilon
    innr  <- 0
    while (diffr > epsilon & innr < 100) {
      innr <- innr + 1
      old  <- tmu
      fity <- alpha + beta*tmu
      tmu  <- (h * tx + g * beta * (ty - alpha))/(h + g * beta^2)
      diffr <- sum((tmu - old)^2)/sum(tmu^2)
    }
    return(list(mu=tmu, innr=innr))
  }

  # calculates L from alpha, beta
  qform   <- function(par) {
    alpha <- par[1]
    beta  <- par[2]
    diffr <- 2*epsilon
    refine <- 0
    while(diffr > epsilon & refine < 100) {
      refine <- refine+1
      old   <- tmu
      fity  <- alpha + beta*tmu
      resi  <- ty - alpha - beta*tx
      g     <- gfun(tmu,  gparms)
      h     <- hfun(fity,hparms)
      tmu   <- getmu(tx, ty, alpha, beta, g, h, tmu)$mu
      scalr <- resi / sqrt((h+g*beta^2))
      diffr <- sum((old-tmu)^2)/sum(tmu^2)
    }

    W     <- sum((tx-tmu)^2/g+(ty-fity)^2/h)
    slgh  <- sum(log(g*h))
    L     <- W + slgh
    return(list(L=L, alpha=alpha, beta=beta, mu=tmu, g=g, h=h, fity=fity, resi=resi, scalr=scalr))
  }

  # Wrapper
  wrapqform <- function(par) {
    do <- qform(par)
    do$L
  }

  alpha    <- 0              # Initial values
  beta     <- 1
  fullmu   <- X
  n        <- length(X)
  allpars  <- NULL
  top      <- 1
  sealpha  <- NA
  sebeta   <- NA
  covar    <- NA
  if (getCI) top <- n
  for (dele in 0:top) {
    tx  <- X
    ty  <- Y
    tmu <- fullmu
    if (dele > 0) {
      tx  <- X[-dele]
      ty  <- Y[-dele]
      tmu <- fullmu[-dele]
    }
    par  <- c(alpha, beta)
    do   <- optim(par, wrapqform)
    pars <- do$par
    if (dele == 0) {
      full      <- qform(pars)
      alpha     <- do$par[1]
      beta      <- do$par[2]
      L         <- do$value
      fullmu    <- full$mu
      fullg     <- full$g
      fullh     <- full$h
      fullfity  <- full$fity
      fullresi  <- full$resi
      fullscalr <- full$scalr
    } else {
      allpars <- rbind(allpars, do$par)
    }
  }
  if (getCI) {
    sealpha <- (n-1)*sd(allpars[,1]) / sqrt(n)
    sebeta  <- (n-1)*sd(allpars[,2]) / sqrt(n)
    covar   <- sealpha * sebeta * cor(allpars[,1], allpars[,2])
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
  }

  if (printem) {
    message(sprintf("%9s %8s %8s %9s", "Parameter", "estimate", "se", "CI"))
    tcut <- qt(0.975, n-1)
    CI   <- alpha + tcut * sealpha * c(-1,1)
    message(sprintf("Intercept %8.3f %8.3f (%7.3f, %6.3f)", alpha, sealpha, CI[1], CI[2]))
    CI   <- beta + tcut * sebeta * c(-1,1)
    message(sprintf("slope     %8.3f %8.3f (%7.3f, %6.3f)\n", beta , sebeta, CI[1], CI[2]))
    if (nMDL > 0) {
      for (kk in 1:nMDL) {
        message(sprintf("MDL %7.3f prediction %7.3f CI %7.3f %7.3f",
                        MDL[kk], preMDL[kk], preMDLl[kk], preMDLu[kk]))
      }
    }
    if(sum(whichmissing) > 0) message(sprintf("\t\t Fit on n = %i complete readings", sum(!whichmissing)))
  }
  corXY <- cor(X,Y)

  allfity = rep(NA, length(allX))
  allfity[!whichmissing] = full$fity
  allmu = rep(NA, length(allX))
  allmu[!whichmissing] = full$mu
  allresi = rep(NA, length(allX))
  allresi[!whichmissing] = full$resi
  allscalr = rep(NA, length(allX))
  allscalr[!whichmissing] = full$scalr

  return(list(alpha=alpha, beta=beta, cor=corXY, fity=allfity, mu=allmu,
              resi=allresi, scalr=allscalr, L=L, sealpha=sealpha, sebeta=sebeta,
              covar=covar, preMDL=preMDL, preMDLl=preMDLl, preMDLu=preMDLu))
}
