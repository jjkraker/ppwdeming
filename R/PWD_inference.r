#' Weighted Deming Regression -- Inference
#' @name PWD_inference
#'
#' @description
#' This routine fits the regression, uses the jackknife to get its precision,
#' and optionally prints it out.  Currently implements Rocke-Lorenzato as
#' the variance profile model.
#'
#' @usage
#' PWD_inference(X, Y, lambda=1, MDL=NA, epsilon=1e-8, printem=TRUE)
#'
#' @param X		the vector of predicate readings,
#' @param Y		the vector of test readings,
#' @param lambda		*optional* (default of 1) - the ratio of the X to the Y precision profile.
#' @param MDL		    *optional* (default to missing) - medical decision level(s),
#' @param epsilon		*optional* (default of 1.e-8) - convergence tolerance limit,
#' @param printem	  *optional* - if TRUE, routine will print out results.
#'
#' @details  For the linear model relating the predicate and test readings,
#' the standard errors of the estimators \eqn{\hat{\alpha}},
#' \eqn{\hat{\beta}},  and their covariance are estimated by
#' the jackknife.  The estimates of the intercept and slope are output,
#'  along with their standard errors and covariance.
#'
#'  These estimates are further used
#'  to estimate the predictions at the input `MDL`.
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha }{the fitted intercept}
#'   \item{beta }{the fitted slope}
#'   \item{cor }{the Pearson correlation between X and Y}
#'   \item{fity }{the vector of predicted Y}
#'   \item{mu }{the vector of estimated latent true values}
#'   \item{resi }{the vector of residuals}
#'   \item{preresi }{the vector of leave-one-out predicted residuals}
#'   \item{sigma }{the estimate of the Rocke-Lorenzato \eqn{\sigma}}
#'   \item{kappa }{the estimate of the Rocke-Lorenzato \eqn{\kappa}}
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
#' # fit with RL precision profile to estimate parameters and variability
#' \dontrun{
#' RL_inf <- PWD_inference(X,Y,MDL=12,printem=TRUE)
#' }
#'
#' @references Hawkins DM and Kraker JJ. Precision Profile Weighted Deming
#' Regression for Methods Comparison, on *Arxiv* (2025, [arxiv.org/abs/2508.02888](https://arxiv.org/abs/2508.02888)).
#'
#' @references  Efron, B (1982). The jackknife, the bootstrap and other resampling plans.
#' Society for Industrial and Applied Mathematics.
#'
#' @export

PWD_inference <- function(X, Y, lambda=1, MDL=NA, epsilon=1e-8, printem=TRUE) {

  n <- length(X)
  pseudalpha <- NULL
  pseudbeta  <- NULL
  preresi    <- NULL
  for (dele in 0:n) {
    x  <- X
    y  <- Y
    if (dele > 0) {
      x   <- X[-dele]
      y   <- Y[-dele]
    }

    do <- PWD_get_gh(x, y, lambda, epsilon)
    if (dele == 0) {
      fullalpha  <- do$alpha
      fullbeta   <- do$beta
      mu         <- do$mu
      resi       <- do$resi
      fity       <- do$fity
      sigma      <- do$sigma
      kappa      <- do$kappa
      like       <- do$like
    } else {
      subalpha   <- do$alpha
      subbeta    <- do$beta
      pseudalpha <- c(pseudalpha, n*fullalpha - (n-1)*subalpha)
      pseudbeta  <- c(pseudbeta , n*fullbeta  - (n-1)*subbeta )
      preresi    <- c(preresi, Y[dele]-subalpha-subbeta*X[dele])
    }
  }
  alpha   <- fullalpha
  beta    <- fullbeta
  sealpha <- sd(pseudalpha) / sqrt(n)
  sebeta  <- sd(pseudbeta ) / sqrt(n)
  covar   <- sealpha * sebeta * cor(pseudalpha,pseudbeta)
  tcut <- qt(0.975, n-1)
  nMDL    <- 0
  preMDL  <- NA
  preMDLl <- NA
  preMDLu <- NA
  if (!is.na(sum(MDL))) {
    nMDL    <- length(MDL)
    preMDL  <- alpha + beta*MDL
    MoEpre  <- tcut*sqrt(sealpha^2 + (sebeta*MDL)^2 + 2*covar*MDL)
    preMDLl <- preMDL - MoEpre
    preMDLu	<- preMDL + MoEpre
  }
  if (printem) {
    cat(sprintf("%9s %8s %8s %9s\n", "Parameter", "estimate", "se", "CI"))
    CI   <- fullalpha + tcut * sealpha * c(-1,1)
    cat(sprintf("Intercept %8.3f %8.3f (%7.3f, %6.3f)\n", fullalpha, sealpha, CI[1], CI[2]))
    CI   <- fullbeta + tcut * sebeta * c(-1,1)
    cat(sprintf("slope     %8.3f %8.3f (%7.3f, %6.3f)\n", fullbeta , sebeta, CI[1], CI[2]))
    cat(sprintf("\nsigma %6.4f kappa %6.4f -2 log likelihood %7.3f\n", sigma, kappa, like))
    if (nMDL > 0) {
      for (kk in 1:nMDL) {
        cat(sprintf("MDL %7.3f prediction %7.3f CI %7.3f %7.3f\n",
                    MDL[kk], preMDL[kk], preMDLl[kk], preMDLu[kk]))
      }
    }
  }
  corXY = cor(X,Y)
  return(list(alpha=fullalpha, beta=fullbeta, cor=corXY, fity=fity, mu=mu,
              resi=resi, preresi=preresi, sigma=sigma, kappa=kappa, like=like,
              sealpha=sealpha, sebeta=sebeta, covar=covar,
              preMDL=preMDL, preMDLl=preMDLl, preMDLu=preMDLu))
}
