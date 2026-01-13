#' Weighted Deming Regression -- Inference
#' @name PWD_inference
#'
#' @description
#' This routine fits the regression, uses the jackknife to get its precision,
#' and optionally prints it out.  Implements Rocke-Lorenzato as
#' the variance profile model.
#'
#' @usage
#' PWD_inference(X, Y, lambda=1, rho=NA, alpha=NA, beta=NA, mu=NA, MDL=NA,
#'               epsilon=1e-8, printem=FALSE)
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param lambda		*optional* (default of 1) - the ratio of the `X` to
#' the `Y` precision profile.
#' @param rho       *optional* (default of NA) - numeric, single value or vector, initial estimate(s) of \eqn{\rho = \frac{\sigma}{\kappa}}.
#' @param alpha     *optional* (default of NA) - numeric, single value, initial estimate of \eqn{\alpha}.
#' @param beta      *optional* (default of NA) - numeric, single value, initial estimate of \eqn{\beta}.
#' @param mu        *optional* (default of NA) - numeric, vector of length of `X`, initial estimate of \eqn{\mu}.
#' @param MDL		    *optional* (default to missing) - medical decision level(s).
#' @param epsilon		*optional* (default of 1.e-8) - convergence tolerance limit.
#' @param printem	  *optional* (default of FALSE) - if TRUE, routine will print out results as a `message`.
#'
#' @details  For the linear model relating the predicate and test readings,
#' the standard errors of the estimators \eqn{\hat{\alpha}},
#' \eqn{\hat{\beta}},  and their covariance are estimated by
#' the jackknife.  The point estimates of the intercept and slope are output,
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
#' \donttest{RL_inf <- PWD_inference(X,Y,MDL=12,printem=TRUE)}
#'
#' @references Hawkins DM and Kraker JJ (in press). Precision Profile Weighted
#' Deming Regression for Methods Comparison. *The Journal of Applied Laboratory Medicine*.
#' <doi:10.1093/jalm/jfaf183>
#'
#' @references  Efron, B (1982). The jackknife, the bootstrap and other resampling plans.
#' Society for Industrial and Applied Mathematics.
#'
#' @importFrom stats complete.cases
#'
#' @export

PWD_inference <- function(X, Y, lambda=1, rho=NA, alpha=NA, beta=NA, mu=NA, MDL=NA,  epsilon=1e-8, printem=FALSE) {
  whichmissing <- (!complete.cases(X)) | (!complete.cases(Y))
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  allY <- Y
  X <- X[!whichmissing]
  Y <- Y[!whichmissing]
  if(sum(!is.na(mu)) > 0) mu <- mu[!whichmissing]

  n <- length(X)
  allalpha <- NULL
  allbeta  <- NULL
  preresi  <- NULL
  for (dele in 0:n) {
    x   <- X
    y   <- Y
    mud <- mu
    if (dele > 0) {
      x   <- X [-dele]
      y   <- Y [-dele]
      mud <- mu[-dele]
    }
    do <- PWD_get_gh(x, y, lambda, rho, alpha, beta, mud, epsilon)

    if (dele == 0) {
      resi  <- do$resi
      fity  <- do$fity
      sigma <- do$sigma
      kappa <- do$kappa
      rho   <- sigma/kappa
      alpha <- do$alpha
      beta  <- do$beta
      mu    <- do$mu
      L     <- do$L
    } else {
      allalpha <- c(allalpha, do$alpha)
      allbeta  <- c(allbeta , do$beta )
      preresi    <- c(preresi, Y[dele]-do$alpha-do$beta*X[dele])
    }
  }
  sealpha <- (n-1) * sd(allalpha) / sqrt(n)
  sebeta  <- (n-1) * sd(allbeta )/ sqrt(n)
  covar   <- sealpha * sebeta * cor(allalpha, allbeta)
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
    message(sprintf("%9s %8s %8s %9s", "Parameter", "estimate", "se", "CI"))
    CI   <- alpha + tcut * sealpha * c(-1,1)
    message(sprintf("Intercept %8.3f %8.3f (%7.3f, %6.3f)", alpha, sealpha, CI[1], CI[2]))
    CI   <- beta + tcut * sebeta * c(-1,1)
    message(sprintf("slope     %8.3f %8.3f (%7.3f, %6.3f)", beta , sebeta, CI[1], CI[2]))
    message(sprintf("\nsigma %6.4f kappa %6.4f -2 log likelihood %7.3f", sigma, kappa, L))
    if (nMDL > 0) {
      for (kk in 1:nMDL) {
        message(sprintf("MDL %7.3f prediction %7.3f CI %7.3f %7.3f",
                        MDL[kk], preMDL[kk], preMDLl[kk], preMDLu[kk]))
      }
    }
    if(sum(whichmissing) > 0) message(sprintf("\t\t Fit on n = %i complete readings", sum(!whichmissing)))
  }
  corXY = cor(X,Y)

  allpreresi = rep(NA, length(allX))
  allpreresi[!whichmissing] = preresi
  allresi = rep(NA, length(allX))
  allresi[!whichmissing] = resi
  allfity = rep(NA, length(allX))
  allfity[!whichmissing] = fity
  allmu = rep(NA, length(allX))
  allmu[!whichmissing] = mu

  return(list(alpha=alpha, beta=beta, cor=corXY, fity=allfity, mu=allmu,
              resi=allresi, preresi=allpreresi, sigma=sigma, kappa=kappa, L=L,
              sealpha=sealpha, sebeta=sebeta, covar=covar,
              preMDL=preMDL, preMDLl=preMDLl, preMDLu=preMDLu))
}
