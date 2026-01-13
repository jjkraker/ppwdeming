#' Linnet proportional CV weighted Deming
#' @name WD_Linnet
#'
#' @description
#' This routine, provided for convenience, makes Linnet’s constant CV fit.
#'
#' @usage
#' WD_Linnet(X, Y, lambda=1, MDL=NA, getCI=TRUE, epsilon=1e-10, printem=FALSE)
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param lambda		ratio of g function to h function.
#' @param MDL		*optional*  (default of NA) - medical decision limit(s).
#' @param getCI	  *optional*  (default of TRUE) - 	if TRUE, generates jackknife standard errors.
#' @param epsilon		*optional*  (default of 1e-10) -  tolerance limit.
#' @param printem	*optional*  (default of FALSE) - if TRUE, prints out results as a `message`.
#'
#' @details Note that in cases where sigma happens to come out zero,
#' Linnet’s constant CV fit differs from the precision-profile fit
#' since the underlying precision profile models are not the same.
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
#' @examples
#' # library
#' library(ppwdeming)
#'
#' # parameter specifications
#' alpha <- 1
#' beta  <- 1.1
#' true  <- 8*10^((0:99)/99)
#' truey <- alpha+beta*true
#' kappa <- 0.1
#'
#' # simulate single sample - set seed for reproducibility
#' set.seed(1039)
#' # specifications for predicate method
#' X     <- true *(1+kappa*rnorm(100))
#' # specifications for test method
#' Y     <- truey *(1+kappa*rnorm(100))
#'
#' # fit with to estimate linear parameters
#' wd_fit <- WD_Linnet(X, Y, MDL=12, printem=TRUE)
#' cat("\nThe Linnet constant-CV estimated intercept is",
#'     signif(wd_fit$alpha,4), "and the estimated slope is",
#'     signif(wd_fit$beta,4), "\n")
#'
#' @references Linnet K (1993). Evaluation of regression procedures for methods
#' comparison studies. *Clinical Chemistry*, **39**: 424-432.
#'
#' @importFrom stats cor qt sd
#' @importFrom stats complete.cases
#'
#' @export

WD_Linnet <- function(X, Y, lambda=1, MDL=NA, getCI=TRUE, epsilon=1e-10, printem=FALSE) {
  whichmissing <- (!complete.cases(X)) | (!complete.cases(Y))
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  allY <- Y
  X <- X[!whichmissing]
  Y <- Y[!whichmissing]

  allalpha <- NULL
  allbeta  <- NULL
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
      w    <- ((newX+newY)/2)^(-2)
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
      alpha <- a
      beta  <- b
    } else {
      allalpha <- c(allalpha, a)
      allbeta  <- c(allbeta , b)
    }
  }
  sealpha <- NA
  sebeta  <- NA
  covar   <- NA
  nMDL    <- 0
  preMDL  <- NA
  preMDLl <- NA
  preMDLu <- NA

  tcut <- qt(0.975, n-1)

  if (getCI) {
    sealpha <- (n-1)*sd(allalpha) / sqrt(n)
    sebeta  <- (n-1)*sd(allbeta ) / sqrt(n)
    covar   <- sealpha * sebeta * cor(allalpha, allbeta)
  }
  if (!is.na(sum(MDL))) {
    nMDL    <- length(MDL)
    preMDL  <- alpha + beta*MDL
    MoEpre  <- tcut*sqrt(sealpha^2 + (sebeta*MDL)^2 + 2*covar*MDL)
    preMDLl <- preMDL - MoEpre
    preMDLu	<- preMDL + MoEpre
  }
  if (getCI & printem) {
    message(sprintf("Linnet weighted Deming\n\t\test\tse\tCI"))
    CIa <- alpha + sealpha * c(-tcut, tcut)
    CIb <- beta  + sebeta  * c(-tcut, tcut)
    message(sprintf("Intercept\t%3.3f\t%3.3f\t%3.3f\t%3.3f", alpha, sealpha, CIa[1], CIa[2]))
    message(sprintf("Slope    \t%3.3f\t%3.3f\t%3.3f\t%3.3f", beta, sebeta, CIb[1], CIb[2]))
    if(sum(whichmissing) > 0) message(sprintf("\t\t Fit on n = %i complete readings", sum(!whichmissing)))
  }
  corXY = cor(X,Y)
  return(list(alpha=alpha, beta=beta, cor=corXY, sealpha=sealpha, sebeta=sebeta,
              covar=covar, preMDL=preMDL, preMDLl=preMDLl, preMDLu=preMDLu))
}
