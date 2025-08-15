#' Linnet proportional CV weighted Deming
#' @name WD_Linnet
#'
#' @description
#' This routine, provided for convenience, makes Linnet’s constant CV fit.
#'
#' @usage
#' WD_Linnet(X, Y, lambda=1, MDL=NA, getCI=TRUE, epsilon=1e-9, printem=FALSE)
#'
#' @param X		the vector of predicate readings,
#' @param Y		the vector of test readings,
#' @param lambda		ratio of g function to h function,
#' @param MDL		optional medical decision limit(s),
#' @param getCI		if TRUE, generates jackknife standard errors,
#' @param epsilon		optional tolerance limit,
#' @param printem	if TRUE, prints results.
#'
#' @details
#' *This could be added upon literature review.*
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha }{the fitted intercept}
#'   \item{beta }{the fitted slope}
#'   \item{cor }{the Pearson correlation between X and Y}
#'   \item{sealpha }{the jackknife standard error of alpha}
#'   \item{sebeta }{the jackknife standard error of beta}
#'   \item{covar }{the jackknife covariance between alpha and beta}
#'   \item{preMDL }{the predictions at the MDL(s)}
#'   \item{preMDLl }{the lower confidence limit(s) of preMDL}
#'   \item{preMDLu }{the upper confidence limit(s) of preMDL}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @example /Examples/WD_Linnet_man_example.R
#'
#' @references Linnet K (1993). Evaluation of regression procedures for methods
#' comparison studies. *Clinical Chemistry*, **39**, 424-432.
#'
#' @importFrom stats cor qt sd
#'
#' @export

WD_Linnet <- function(X, Y, lambda=1, MDL=NA, getCI=TRUE, epsilon=1e-9, printem=FALSE) {
  pseuda <- NULL
  pseudb <- NULL
  n      <- length(X)
  top    <- n
  if (!getCI) top <- 0
  for (dele in 0:top) {
    x <- X
    y <- Y
    if (dele > 0) {
      x <- X[-dele]
      y <- Y[-dele]
    }
    newX   <- x
    newY   <- y
    diff <- 2*epsilon
    while (diff > epsilon) {
      w    <- ((x+y)/2)^(-2)
      sumw <- sum(w)
      xbar <- sum(w*x)/sumw
      ybar <- sum(w*y)/sumw
      uw   <- sum(w*(x-xbar)^2)
      qw   <- sum(w*(y-ybar)^2)
      pw   <- sum(w*(x-xbar)*(y-ybar))
      surd <- (uw - lambda*qw)^2 + 4*lambda*pw^2
      b    <- (lambda*qw - uw + sqrt(surd))/(2*lambda*pw)
      a    <- ybar - b * xbar
      d    <- y -a - b*x
      deno <- 1 + lambda * b^2
      oldX <- newX
      oldY <- newY
      newX    <- x + lambda * b * d / deno
      newY    <- y - d / deno
      diff <- sum((oldX-newX)^2+(oldY-newY)^2) / sum(oldX^2+oldY^2)
    }
    if (dele == 0) {
      fulla <- a
      fullb <- b
    } else {
      pseuda <- c(pseuda, n*fulla - (n-1)*a)
      pseudb <- c(pseudb, n*fullb - (n-1)*b)
    }
  }
  alpha   <- fulla
  beta    <- fullb
  sealpha <- NA
  sebeta  <- NA
  covar   <- NA
  nMDL    <- 0
  preMDL  <- NA
  preMDLl <- NA
  preMDLu <- NA
  if (getCI) {
    sealpha <- sd(pseuda) / sqrt(n)
    sebeta  <- sd(pseudb) / sqrt(n)
    covar   <- sealpha * sebeta * cor(pseuda, pseudb)
  }
  if (!is.na(sum(MDL))) {
    nMDL    <- length(MDL)
    preMDL  <- alpha + beta*MDL
    MoEpre  <- tcut*sqrt(sealpha^2 + (sebeta*MDL)^2 + 2*covar*MDL)
    preMDLl <- preMDL - MoEpre
    preMDLu	<- preMDL + MoEpre
  }
  if (getCI & printem) {
    cat(sprintf("Linnet weighted Deming\n\testimate\tse\tCI\n"))
    tcut <- qt(0.975, n-1)
    CIa <- alpha + sealpha * c(-tcut, tcut)
    CIb <- beta  + sebeta  * c(-tcut, tcut)
    cat(sprintf("Intercept\t%3.3f\t%3.3f\t%3.3f\t%3.3f\n", alpha, sealpha, CIa[1], CIa[2]))
    cat(sprintf("Slope\t%3.3f\t%3.3f\t%3.3f\t%3.3f\n", beta, sebeta, CIb[1], CIb[2]))
  }
  corXY = cor(X,Y)
  return(list(alpha=alpha, beta=beta, cor=corXY, sealpha=sealpha, sebeta=sebeta,
              covar=covar, preMDL=preMDL, preMDLl=preMDLl, preMDLu=preMDLu))
}
