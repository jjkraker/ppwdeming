#' Weighted Deming Regression -- Outlier scanning
#' @name PWD_outlier
#'
#' @description
#' This function tests for outliers from the fitted regression, and refits on
#' a sanitized data set (with outliers removed).
#'
#' @usage
#' PWD_outlier(X, Y, K, lambda=1, Pcut=0.01, rho=NA, alpha=NA, beta=NA, mu=NA)
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param K		the maximum number of outliers to seek.
#' @param lambda		*optional* (default of 1) - the ratio of the `X` to
#' the `Y` precision profile.
#' @param Pcut		  *optional*, default 0.01 (1%), cutoff for statistical significance of Bonferroni P.
#' @param rho       *optional* (default of NA) - numeric, single value or vector, initial estimate(s) of \eqn{\rho = \frac{\sigma}{\kappa}}.
#' @param alpha     *optional* (default of NA) - numeric, single value, initial estimate of \eqn{\alpha}.
#' @param beta      *optional* (default of NA) - numeric, single value, initial estimate of \eqn{\beta}.
#' @param mu        *optional* (default of NA) - numeric, vector of length of `X`, initial estimate of \eqn{\mu}.
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
#'   \item{basepar}{the sigma, kappa, alpha, beta of the full data set}
#'   \item{lastpar}{the sigma, kappa, alpha, beta of the sanitized data set}
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
#' sigma <- 1
#' kappa <- 0.08
#' alpha <- 1
#' beta  <- 1.1
#' true  <- 8*10^((0:99)/99)
#' truey <- alpha+beta*true
#' # simulate single sample - set seed for reproducibility
#' set.seed(1069)
#' # specifications for predicate method
#' X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(100))
#' # specifications for test method
#' Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(100))
#' # add some outliers
#' Y[c(1,2,100)] <- Y[c(1,2,100)] + c(-10,9,-50)
#'
#' # check for outliers, re-fit, and store output
#' \donttest{outliers_assess <- PWD_outlier(X, Y, K=5)}
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

PWD_outlier <- function(X, Y, K, lambda=1, Pcut=0.01, rho=NA, alpha=NA, beta=NA, mu=NA) {
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

  do_orig <- PWD_get_gh(X, Y, lambda, rho, alpha, beta, mu)

  for (m in 1:K) {
    x        <- X[clean]
    y        <- Y[clean]
    do       <- PWD_get_gh(x, y, lambda)
    if (m == 1) do <- do_orig
    printres <- FALSE
    alpha    <- unname(do$alpha)
    beta     <- unname(do$beta)
    resi     <- y - alpha - beta*x
    fitres   <- PWD_resi(x, resi)
    sigr     <- fitres$sigmar
    kapr     <- fitres$kappar
    scalr    <- fitres$scalr
    if (m == 1) {
      basepar  <- c(do$sigma, do$kappa, unname(do$alpha), unname(do$beta), do$like)
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
    do        <- PWD_get_gh(x, y, lambda)
    alpha     <- unname(do$alpha)
    beta      <- unname(do$beta)
    resi      <- y - alpha - beta*x
    fitres    <- PWD_resi(x, resi)
    sigr      <- fitres$sigmar
    kapr      <- fitres$kappar
    scalr     <- fitres$scalr
    meanr     <- mean(scalr)
    sdr       <- sd(scalr)
    fitsusp   <- alpha + beta*X[drop]
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
      lastpar <- c(do$sigma, do$kappa, unname(do$alpha), unname(do$beta), do$like)
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

  return(list(ndrop=ndrop, drop=drop, cor=cor,
              cleancor=cleancor,
              scalr=allscalr, basepar=basepar, lastpar=lastpar,
              forward=forward, backward=backward,
              tee=tee, BonP=BonP, outlis=outlis))
}
