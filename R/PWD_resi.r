#' Fit Rocke-Lorenzato profile model to residuals
#' @name PWD_resi
#'
#' @description
#' This routine fits the Rocke-Lorenzato precision profile model to the
#' **residuals** from the fit.
#'
#' @usage
#' PWD_resi(true, resi, epsilon=1e-8, printem=FALSE)
#'
#' @param true  	the vector of values used to predict the precision – commonly X.
#' @param resi		the vector of residuals whose variance is thought to be a function of “true”.
#' @param epsilon		*optional* (default of 1e-8) - convergence tolerance limit.
#' @param printem	  *optional* (default of FALSE) - if TRUE, routine will print out results as a `message`.
#'
#' @details  The Rocke-Lorenzato precision profile model is
#' \deqn{SD^2 = \sigma_r^2 + (\kappa_r\cdot true)^2}
#' for the *residuals* from a precision-profile model fit.
#'
#' Under this model, the approach for reviewing residuals is to fit a
#' variance profile model to the residuals \eqn{r_i} themselves.  The output
#' of this function includes a maximum-likelihood estimate of the remaining
#' parameter in the special cases of:
#'    * constant variance (\eqn{\kappa_r} = 0); and
#'    * constant coefficient of variation (\eqn{\sigma_r} = 0).
#'
#' @returns A list containing the following components:
#'
#'   \item{sigmar}{the estimate of \eqn{\sigma_r}}
#'   \item{kappar}{the estimate of \eqn{\kappa_r}}
#'   \item{L}{the -2 log likelihood}
#'   \item{scalr}{the scaled residuals}
#'   \item{poolsig}{the maximum likelihood estimate of \eqn{\sigma_r} if \eqn{\kappa_r} = 0}
#'   \item{poolkap}{the maximum likelihood estimate of \eqn{\kappa_r} if \eqn{\sigma_r} = 0}
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
#' # fit the model and store output
#' RL_gh_fit  <- PWD_get_gh(X,Y)
#' # run the residual analysis from the model output
#' post  <- PWD_resi(X, RL_gh_fit$resi, printem=TRUE)
#'
#' @references Hawkins DM and Kraker JJ (in press). Precision Profile Weighted
#' Deming Regression for Methods Comparison. *The Journal of Applied Laboratory Medicine*.
#' <doi:10.1093/jalm/jfaf183>
#'
#' @references Hawkins DM (2014). A Model for Assay Precision.
#' *Statistics in Biopharmaceutical Research*, **6**, 263-269.
#' http://dx.doi.org/10.1080/19466315.2014.899511
#'
#' @importFrom stats optimize pchisq shapiro.test
#' @importFrom stats median
#' @importFrom stats complete.cases
#'
#' @export

PWD_resi    <- function(true, resi, epsilon=1e-8, printem=FALSE) {
  whichmissing <- (!complete.cases(true)) | (!complete.cases(resi))
  missingcases <- (1:length(true))[whichmissing]
  alltrue <- true
  true <- true[!whichmissing]
  resi <- resi[!whichmissing]

  n       <- length(true)
  absres  <- abs(resi)
  key     <- order(true)
  sortres <- absres[key]
  sorttru <- true  [key]
  ratio   <- sortres / (abs(sorttru)+epsilon)
  cvsq    <- ratio  ^2
  vars    <- sortres^2
  lowr    <- 1:round(n/3)
  hir     <- round(2*n/3):n
  maxsig  <- 5*median(sortres[lowr])
  maxkap  <- 5*median(ratio  [hir])

  innr     <- function(kap, sig) {
    modl   <- sig^2 + (kap*sorttru)^2
    sum((vars / modl + log(modl)))
  }

  a       <- 0
  b       <- maxsig
  gr      <- 1.618
  dif     <- b-a
  h       <- maxsig
  while (h > epsilon) {
    h  <- b - a
    c  <- b - h / gr
    d  <- a + h / gr
    vc <- optimize(innr, c(0,maxkap), c)
    fc <- vc$objective
    vd <- optimize(innr, c(0,maxkap), d)
    fd <- vd$objective
    if (fc < fd) {
      b <- d
    } else {
      a <- c
    }
  }
  sigma <- c
  kappa <- vc$minimum
  L  <- fc
  if (fd < fc) {
    sigma <- d
    kappa <- vd$minimum
    L  <- fd
  }

  profl   <- sqrt(sigma^2 + (kappa*true)^2)
  scalr   <- resi/profl

  proflf  <- sqrt(sigma^2 + (kappa*sorttru)^2)

  poolsig <- sqrt(mean(vars))
  poolkap <- sqrt(mean(cvsq))

  if (printem) {
    SW    <- shapiro.test(scalr)$p.value
    message(sprintf("Rocke-Lorenzato fit to residuals\nsigma %6.4f kappa %6.4f",
                sigma, kappa))
    if(sum(whichmissing) > 0) message(sprintf("\t Fit on n = %i complete residuals\n", sum(!whichmissing)))
    message(sprintf("P value for normality %6.4f",
                SW))
  }

  allscalr = rep(NA, length(alltrue))
  allscalr[!whichmissing] = scalr

  return(list(sigmar=sigma, kappar=kappa, L=L, scalr=scalr,
              poolsig=poolsig, poolkap=poolkap))
}

