#' Weighted Deming Regression -- Outlier scanning
#' @name PWD_outlier
#'
#' @description
#' This function tests for outliers from the fitted regression, and refits on
#' a sanitized data set (with outliers removed).
#'
#' @usage
#' PWD_outlier(X, Y, K, lambda=1, Pcut=0.01,
#'             rho=lifecycle::deprecated(),
#'             alpha=lifecycle::deprecated(), beta=lifecycle::deprecated(),
#'             mu=lifecycle::deprecated(),
#'             quad=FALSE)
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param K		the maximum number of outliers to seek.
#' @param lambda		*optional* (default of 1) - the ratio of the `X` to
#' the `Y` precision profile.
#' @param Pcut		  *optional*, default 0.01 (1%), cutoff for statistical significance of Bonferroni P.
#' @param rho       `r lifecycle::badge("deprecated")` `rho = numvalue` initialization is no longer implemented in algorithm.
#' @param alpha     `r lifecycle::badge("deprecated")` `alpha = numvalue` initialization is no longer implemented in algorithm.
#' @param beta      `r lifecycle::badge("deprecated")` `beta = numvalue` initialization is no longer implemented in algorithm.
#' @param mu        `r lifecycle::badge("deprecated")` `mu = numvector` initialization is no longer implemented in algorithm.
#' @param quad      *optional* (default of FALSE) - logical, selects fitting a linear or a quadratic regression.
#'
#' @details
#' The method is modeled on the Rosner sequential ESD outlier procedure and
#' assumes the sample is large enough to ignore the effect of random variability
#' in the parameter estimates on the distribution of the residuals.
#'
#' @returns A list containing the following components:
#'
#'   \item{ndrop}{the number of significant outliers}
#'   \item{drop}{a vector of the indices of the outliers}
#'   \item{cor }{the Pearson correlation between X and Y}
#'   \item{cleancor }{the Pearson correlation between cleaned X and Y (after outliers removed)}
#'   \item{scalr}{the scaled residuals of all cases from the sanitized fit and whose normal tail areas provide the basis for the outlier P values}
#'   \item{basepar}{the sigma, kappa, alpha, beta, gamma of the full data set}
#'   \item{lastpar}{the sigma, kappa, alpha, beta, gamma of the sanitized data set}
#'   \item{forward }{dataframe summarizing the forward identification of possible outliers}
#'   \item{backward }{dataframe summarizing the backward reinclusion of cases}
#'   \item{tee }{the t statistics of the final identified outliers}
#'   \item{BonP }{the Bonferroni P-value of the final identified outliers}
#'   \item{outlis }{dataframe containing the outlier cases, test statistics, and P-values}
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
#' true  <- 8*10^((0:(n-1))/(n-1))
#' truey <- alpha+beta*true
#' # simulate single sample - set seed for reproducibility
#' set.seed(1069)
#' # specifications for predicate method
#' X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(n))
#' # specifications for test method
#' Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(n))
#' # add some outliers
#' Y[c(1,2,100)] <- Y[c(1,2,100)] + c(-10,9,-50)
#'
#' # check for outliers and store output
#' outliers_assess <- PWD_outlier(X, Y, K=5)
#'
#' # summary of process
#' summary(outliers_assess)
#'
#' @references Hawkins DM and Kraker JJ (2026). Precision Profile Weighted
#' Deming Regression for Methods Comparison.
#' *The Journal of Applied Laboratory Medicine*, **11**(2), 379-392.
#' <doi:10.1093/jalm/jfaf183>
#'
#' @references  Hawkins DM (2008). *Outliers* in Wiley Encyclopedia of Clinical Trials,
#' eds R. D’Agostino, L. Sullivan, and J. Massaro. Wiley, New York.
#'
#' @importFrom stats pnorm
#' @importFrom stats complete.cases
#' @importFrom stats na.omit
#' @importFrom stats pt
#' @importFrom utils capture.output
#'
#' @export

PWD_outlier <- function(X, Y, K, lambda=1, Pcut=0.01,
                        rho=lifecycle::deprecated(),
                        alpha=lifecycle::deprecated(), beta=lifecycle::deprecated(),
                        mu=lifecycle::deprecated(),
                        quad=FALSE) {
  if (lifecycle::is_present(rho)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_outlier(rho)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(alpha)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_outlier(alpha)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(beta)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_outlier(beta)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(mu)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_outlier(mu)",
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

  N       <- length(X)
  keep    <- rep(FALSE, N)
  clean   <- 1:N
  nclean  <- N
  ndrop   <- 0
  initlis <- NULL
  initZ   <- NULL

  do_orig <- PWD_get_gh(X, Y, lambda, quad=quad)

  for (m in 1:K) {
    x        <- X[clean]
    y        <- Y[clean]
    do       <- PWD_get_gh(x, y, lambda, quad=quad)
    if (m == 1) do <- do_orig
    printres <- FALSE
    sigma    <- unname(do$sigma)
    kappa    <- unname(do$kappa)
    alpha    <- unname(do$alpha)
    beta     <- unname(do$beta )
    gamma    <- unname(do$gamma)
    L        <- unname(do$L    )
    resi     <- y - alpha - beta*x
    fitres   <- PWD_resi(x, resi)
    sigr     <- fitres$sigmar
    kapr     <- fitres$kappar
    scalr    <- fitres$scalr
    if (m == 1) {
      basepar  <- c(sigma, kappa, alpha, beta, gamma, L   )
      lastpar  <- basepar                 # If there are no outliers
      cor      <- cor(X,Y)
      allscalr <- scalr
    }
    meanr    <- mean(scalr)
    sdr      <- sd(scalr)
    Z        <- (scalr-meanr)/sdr
    whereis  <- (1:nclean)[abs(Z) == max(abs(Z))][1]
    idof     <- clean[ whereis]
    clean    <- clean[-whereis]
    nclean   <- nclean-1
    drop     <- c(drop, idof)
    ndrop    <- ndrop+1
    maxZ     <- Z[whereis]
    initlis  <- c(initlis, idof)
    initZ    <- c(initZ, maxZ)
  }

  drop       <- initlis
  backlis    <- NULL
  backP      <- NULL
  for (m in 1:K) {
    x         <- X[clean]
    y         <- Y[clean]
    do        <- PWD_get_gh(x, y, lambda, quad=quad)
    sigma     <- unname(do$sigma)
    kappa     <- unname(do$kappa)
    alpha     <- unname(do$alpha)
    beta      <- unname(do$beta )
    gamma     <- unname(do$gamma)
    L         <- unname(do$L    )
    resi      <- y - alpha - beta*x - gamma*x^2
    fitres    <- PWD_resi(x, resi)
    sigr      <- fitres$sigmar
    kapr      <- fitres$kappar
    scalr     <- fitres$scalr
    meanr     <- mean(scalr)
    sdr       <- sd(scalr)
    fitsusp   <- alpha + beta*X[drop] + gamma*X[drop]^2
    resisusp  <- Y[drop] - fitsusp
    profl     <- sqrt(sigr^2 + (kapr*fitsusp)^2)
    scalrsusp <- resisusp/profl
    tee       <- sqrt((nclean-1)/nclean) * (scalrsusp-meanr)/sdr
    rawp      <- pt(-abs(tee), nclean-1)
    BonP      <- 2*(nclean+1)*rawp
    maxbon    <- max(BonP)
    whereis   <- (1:ndrop)[BonP == maxbon][1]
    teeof     <- tee [whereis]
    idof      <- drop[whereis]
    backlis <- c(backlis, idof)
    backP   <- c(backP, maxbon)
    if (maxbon < Pcut) {
      lastpar <- c(sigma, kappa, alpha, beta, gamma, L   )
      break
    }
    drop    <- drop[-whereis]
    clean   <- c(clean, idof)
    nclean  <- nclean+1
    ndrop   <- ndrop-1
  }

  cleancor    <- cor(X[clean], Y[clean])
  keep[clean] <- TRUE
  if (ndrop > 0) {
    allscalr[clean] <- scalr
    allscalr[drop ] <- scalrsusp
  }
  forward <- data.frame(initlis, round(initZ,3))
  colnames(forward) <- c("case","Z")
  backward <- data.frame(backlis, backP)
  colnames(backward) <- c("case", "Bonf P")
  outlis <- NA
  if (ndrop > 0) {
    outlis   <- data.frame(drop, X[drop], Y[drop], round(tee,3), signif(BonP, 5))
    colnames(outlis) <- c("case", "X", "Y", "t", "Bonf P")
  }

  fullout <- list(ndrop=ndrop, drop=drop, cor=cor,
                  cleancor=cleancor,
                  scalr=allscalr, basepar=basepar, lastpar=lastpar,
                  forward=forward, backward=backward, Pcut=Pcut,
                  tee=tee, BonP=BonP, outlis=outlis)
class(fullout) <- c("pwdout", "list")

return(fullout)
}
