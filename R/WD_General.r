#' Weighted Deming Regression
#' @name WD_General
#'
#' @description
#' This code fits the weighted Deming regression on predicate readings (X)
#' and test readings (Y).
#'
#' @usage
#' WD_General(X, Y, g, h, epsilon=1e-10)
#'
#' @param X		the vector of predicate readings,
#' @param Y		the vector of test readings,
#' @param g		the vector of variances of the X,
#' @param h		the vector of variances of the Y,
#' @param epsilon		*optional* convergence tolerance limit.
#'
#' @details For input vectors `g` and `h` containing the variances of
#' predicate readings `X` and test readings `Y`, respectively, iteratively fits
#' weighted Deming regression.
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha }{the fitted intercept}
#'   \item{beta }{the fitted slope}
#'   \item{cor }{the Pearson correlation between X and Y}
#'   \item{fity }{the vector of predicted Y}
#'   \item{mu }{the vector of estimated latent true values}
#'   \item{resi }{the vector of residuals}
#'   \item{like }{the -2 log likelihood L}
#'   \item{innr }{the number of inner refinement loops executed}
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
#' # Loosely motivated by Vitamin D data set
#' g     <- 4e-16+0.07*true^1.27
#' h     <- 6e-2+7e-5*truey^2.2
#'
#' # simulate single sample - set seed for reproducibility
#' set.seed(1039)
#' # specifications for predicate method
#' X     <- true +sqrt(g)*rnorm(100)
#' # specifications for test method
#' Y     <- truey+sqrt(h)*rnorm(100)
#'
#' # fit with to estimate linear parameters
#' wd_fit <- WD_General(X,Y,g,h)
#' cat("\nWith given g and h, the estimated intercept is",
#'     signif(wd_fit$alpha,4), "and the estimated slope is",
#'     signif(wd_fit$beta,4), "\n")
#'
#' @references Ripley BD and Thompson M (1987). Regression techniques for the detection
#' of analytical bias.  *Analyst*, **112**, 377-383.
#'
#' @export

WD_General <- function(X, Y, g, h, epsilon=1e-10) {
  old    <- X
  diff   <- 2*epsilon
  flag   <- any(is.na(g+h))
  like   <- 1e10
  innr   <- 0
  n      <- length(X)
  mu     <- rep(NA, n)
  fity   <- mu
  resi   <- mu
  best   <- 1e10
  beta   <- 1
  alllik <- NULL
  if (flag) {
    print("g h error in WD")
    print(summary(g))
    print(summary(h))
    print("X")
    print(X)
    print("Y")
    print(Y)
  } else {
    while(diff > epsilon) {# weight depends on beta, so refine.
      innr   <- innr+1
      w      <- 1/(h + beta^2 * g)
      sumw   <- sum(w)
      wsq    <- w^2
      xbar   <- sum(w * X) / sumw
      ybar   <- sum(w * Y) / sumw
      devx   <- X - xbar
      devy   <- Y - ybar
      sxxx   <- sum(wsq * devx^2 * g)
      sxxy   <- sum(wsq * devx^2 * h)
      sxyx   <- sum(wsq * devx * devy * g)
      sxyy   <- sum(wsq * devx * devy * h)
      syyx   <- sum(wsq * devy^2 * g)
      surd   <- (sxxy - syyx)^2 + 4 * sxyx * sxyy
      beta   <- (syyx - sxxy + sqrt(surd)) / (2 * sxyx)
      alpha  <- ybar - beta *  xbar
      mu     <- w * (h * X + g * beta * (Y - alpha))
      fity   <- alpha + beta*mu
      resi   <- Y - alpha - beta*X
      like   <- sum((X-mu)^2/g + (Y-fity)^2/h + log(g*h))
      alllik <- c(alllik, like)
      if (like < best) {
        best     <- like
        besalpha <- alpha
        besbeta  <- beta
        besmu    <- mu
        besresi  <- resi
        besfity  <- fity
      }
      diff   <- sum((mu - old)^2) / sum(mu^2)
      old   <- mu
    }
  }
  corXY = cor(X,Y)
  return(list(alpha=besalpha, beta=besbeta, cor=corXY, fity=besfity, mu=besmu,
              resi=besresi, like=best, innr=innr))
}
