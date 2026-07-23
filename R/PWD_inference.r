#' Precision-Profile Weighted Deming Regression -- Inference
#' @name PWD_inference
#'
#' @description
#' This routine fits the regression and uses the jackknife to get its precision.
#' Implements Rocke-Lorenzato as the variance profile model.
#'
#' @usage
#' PWD_inference(X, Y, lambda=1,
#'               rho=lifecycle::deprecated(),
#'               alpha=lifecycle::deprecated(), beta=lifecycle::deprecated(),
#'               mu=lifecycle::deprecated(),
#'               MDL=NA,
#'               quad=FALSE, epsilon=1e-8,
#'               printem=lifecycle::deprecated())
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param lambda		*optional* (default of 1) - the ratio of the `X` to the `Y` precision profile.
#' @param rho       `r lifecycle::badge("deprecated")` `rho = numvalue` initialization is no longer implemented in algorithm.
#' @param alpha     `r lifecycle::badge("deprecated")` `alpha = numvalue` initialization is no longer implemented in algorithm.
#' @param beta      `r lifecycle::badge("deprecated")` `beta = numvalue` initialization is no longer implemented in algorithm.
#' @param mu        `r lifecycle::badge("deprecated")` `mu = numvector` initialization is no longer implemented in algorithm.
#' @param quad      *optional* (default of FALSE) - logical, selects fitting a linear or a quadratic regression.
#' @param MDL		    *optional* (default to missing) - medical decision level(s).
#' @param epsilon		*optional* (default of 1.e-8) - convergence tolerance limit.
#' @param printem `r lifecycle::badge("deprecated")` `printem = TRUE` is no
#' longer supported; separate summary functions are provided for printing output.
#'
#' @details  For the linear model relating the predicate and test readings,
#' the standard errors of the estimators \eqn{\hat{\alpha}},
#' \eqn{\hat{\beta}}, (and \eqn{\hat{\gamma}}, if quadratic) and their covariance are estimated by
#' the jackknife.  The point estimates of the intercept and slope are output,
#'  along with their standard errors and covariance.
#'
#'  These estimates are further used
#'  to estimate the predictions at the input `MDL`, if appropriate.
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha }{the fitted intercept}
#'   \item{beta }{the fitted slope}
#'   \item{gamma}{the fitted coefficient of the square: if quad is FALSE, the value is zero}
#'   \item{quad}{logical, FALSE for linear or TRUE for quadratic regression}
#'   \item{cor }{the Pearson correlation between X and Y}
#'   \item{fity }{the vector of predicted Y}
#'   \item{mu }{the vector of estimated latent true values}
#'   \item{resi }{the vector of residuals}
#'   \item{preresi }{the vector of leave-one-out predicted residuals}
#'   \item{sigma }{the estimate of the Rocke-Lorenzato \eqn{\sigma}}
#'   \item{kappa }{the estimate of the Rocke-Lorenzato \eqn{\kappa}}
#'   \item{L }{the -2 log likelihood L}
#'   \item{sealpha }{the jackknife standard error of  \eqn{alpha}}
#'   \item{sebeta }{the jackknife standard error of  \eqn{beta}}
#'   \item{segamma }{the jackknife standard error of  \eqn{gamma}}
#'   \item{covar }{the jackknife covariance between `alpha` and `beta` (and `gamma`, if quadratic)}
#'   \item{CIalpha}{the level 95% confidence-interval estimate of \eqn{\alpha}}
#'   \item{CIbeta}{the level 95% confidence-interval estimate of \eqn{\beta}}
#'   \item{CIgamma}{the level 95% confidence-interval estimate of \eqn{\gamma}}
#'   \item{preMDL }{the predictions at the MDL(s)}
#'   \item{sepreMDL }{the jackknife standard error of predictions at the MDL}
#'   \item{preMDLl }{the lower confidence limit(s) of preMDL}
#'   \item{preMDLu }{the upper confidence limit(s) of preMDL}
#'   \item{whichmissing}{a logical vector identifying locations of missing X or Y values}
#'   \item{MDL}{vector of medical decision level(s)}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @examples
#' # library
#' library(ppwdeming)
#'
#' # parameter specifications
#' n <- 100
#'
#' sigma <- 1
#' kappa <- 0.08
#' alpha <- 1
#' beta  <- 1.1
#' quad  <- FALSE # or
#' quad  <- TRUE
#' gamma <- 0
#' if (quad) gamma <- -0.005
#' true  <- 8*10^((0:(n-1))/(n-1))
#' truey <- alpha+beta*true+gamma*true^2
#' # simulate single sample - set seed for reproducibility
#' set.seed(1039)
#' # specifications for predicate method
#' X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(n))
#' # specifications for test method
#' Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(n))
#'
#' # fit with RL precision profile to estimate parameters and variability
#' RL_inf <- PWD_inference(X,Y,MDL=12, quad=quad)
#'
#' # summary of results
#' summary(RL_inf)
#' summary(RL_inf, conflevel=.95)
#'
#' @references Hawkins DM and Kraker JJ (2026). Precision Profile Weighted
#' Deming Regression for Methods Comparison.
#' *The Journal of Applied Laboratory Medicine*, **11**(2), 379-392.
#' <doi:10.1093/jalm/jfaf183>
#'
#' @references  Efron, B (1982). The jackknife, the bootstrap and other resampling plans.
#' Society for Industrial and Applied Mathematics.
#'
#' @importFrom stats complete.cases
#'
#' @export

PWD_inference <- function(X, Y, lambda=1,
                          rho=lifecycle::deprecated(),
                          alpha=lifecycle::deprecated(), beta=lifecycle::deprecated(),
                          mu=lifecycle::deprecated(),
                          MDL=NA,
                          quad=FALSE, epsilon=1e-8,
                          printem=lifecycle::deprecated()) {
  if (lifecycle::is_present(printem)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_inference(printem)",
      details = "Argument printem no longer prints results-message.  \n See summary functions instead for printing output. \n Argument will be dropped in next release."
    )
  }
  if (lifecycle::is_present(rho)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_inference(rho)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(alpha)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_inference(alpha)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(beta)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_inference(beta)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(mu)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_inference(mu)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  rho=NA; alpha=NA; beta=NA; gamma=NA; mu=NA

  whichmissing <- (!complete.cases(X)) | (!complete.cases(Y))
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  allY <- Y
  X <- X[!whichmissing]
  Y <- Y[!whichmissing]
  if(sum(!is.na(mu)) > 0) mu <- mu[!whichmissing]

  n       <- length(X)
  alljack <- NULL
  preresi <- NULL
  for (dele in 0:n) {
    x   <- X
    y   <- Y
    mud <- mu
    if (dele > 0) {
      x   <- X [-dele]
      y   <- Y [-dele]
    }
    do <- PWD_get_gh(x, y, lambda, quad=quad, epsilon=epsilon)

    if (dele == 0) {
      resi  <- do$resi
      fity  <- do$fity
      sigma <- do$sigma
      kappa <- do$kappa
      gamma <- do$gamma
      rho   <- sigma/kappa
      alpha <- do$alpha
      beta  <- do$beta
      gamma <- do$gamma
      mu    <- do$mu
      L     <- do$L
    } else {
      jack    <- c(do$alpha, do$beta, do$gamma)
      alljack <- rbind(alljack, jack)
      rhs     <- c(1, X[dele], X[dele]^2)
      delfit  <- t(jack) %*% rhs
      preresi <- c(preresi, Y[dele] - delfit)
    }
  }
  covx    <- (n-1)^2 / n * cov(alljack)
  tcut <- qt(0.975, n-1)
  nMDL    <- 0
  preMDL  <- NA
  sepreMDL <- NA
  preMDLl <- NA
  preMDLu <- NA
  if (!is.na(sum(MDL))) {
    nMDL    <- length(MDL)
    terms   <- c(1, MDL, MDL^2)
    preMDL  <- alpha + beta*MDL + gamma*MDL^2
    varof   <- c(terms) %*% covx %*% c(terms)
    sepreMDL  <- sqrt(varof)
    MoEpre  <- tcut*sqrt(varof)
    preMDLl <- preMDL - MoEpre
    preMDLu	<- preMDL + MoEpre
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

  if (!quad) covar <- covx[1,2]
  if ( quad) covar <- c(covx[1, 2:3], covx[2,3])
  sealpha  <- sqrt(covx[1,1])
  sebeta   <- sqrt(covx[2,2])
  segamma  <- ifelse( quad, sqrt(covx[3,3]), NA)

  clevelused <- .95
  tcut <- qt(0.975, n-1)

  CIalpha   <- alpha + tcut * sealpha * c(-1,1)
  CIbeta  <- beta + tcut * sebeta * c(-1,1)
  CIgamma  <- gamma+c(-tcut,tcut)*segamma

  fullout <- list(alpha=alpha, beta=beta, gamma=gamma, quad=quad,
                  cor=corXY, fity=allfity, mu=allmu,
                  resi=allresi, preresi=allpreresi,
                  sigma=sigma, kappa=kappa, L=L,
                  sealpha=sqrt(covx[1,1]), sebeta=sqrt(covx[2,2]),
                  segamma=sqrt(covx[3,3]), covar=covar,
                  preMDL=preMDL, preMDLl=preMDLl, preMDLu=preMDLu,
                  clevelused=clevelused, CIalpha=CIalpha, CIbeta=CIbeta, CIgamma=CIgamma,
                  preMDL=preMDL, sepreMDL=sepreMDL, preMDLl=preMDLl, preMDLu=preMDLu,
                  whichmissing=whichmissing, MDL=MDL)

  class(fullout) <- c("pwdinf")

  return(fullout)

}
