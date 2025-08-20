#' Fit Rocke-Lorenzato profile model to residuals
#' @name PWD_resi
#'
#' @description
#' This routine fits the Rocke-Lorenzato precision profile model to the
#' **residuals** from the fit (via `PWD_inference`).
#'
#' @usage
#' PWD_resi(true, resi, epsilon=1e-5, printem=TRUE)
#'
#' @param true  	the vector of values used to predict the precision – commonly X,
#' @param resi		the vector of residuals whose variance is thought to be a function of “true”,
#' @param epsilon		*optional* (default of 1e-5) - convergence tolerance limit,
#' @param printem	  *optional* - if TRUE, routine will print out results.
#'
#' @details  The Rocke-Lorenzato precision profile model is
#' \deqn{SD^2 = \sigma_r^2 + (\kappa_r\cdot true)^2}
#' for the *residuals* from a precision-profile model fit.
#'
#'  Under this model, the approach for reviewing residuals is to fit a
#'  variance profile model to the residuals \eqn{r_i} themselves.
#'  This function includes a check for the special cases of
#'    * constant variance (\eqn{\kappa_r=0}) - in this case,
#'    one could switch to the simpler unweighted Deming model;
#'    * and of constant coefficient of variation (\eqn{\sigma_r=0}) - in this case,
#'    one could switch to the constant CV weighted Deming model.
#'
#'  using chi-squared tests.
#'
#' @returns A list containing the following components:
#'
#'   \item{sigmar}{the estimate of \eqn{\sigma_r}}
#'   \item{kappar}{the estimate of \eqn{\kappa_r}}
#'   \item{like}{the likelihood}
#'   \item{scalr}{the scaled residuals}
#'   \item{poolsig}{the maximum likelihood estimate of \eqn{\sigma_r} if \eqn{\kappa_r} =0}
#'   \item{poolkap}{the maximum likelihood estimate of \eqn{\kappa_r} if \eqn{\sigma_r} =0}
#'   \item{tests}{the chi-squared test statistics for \eqn{\kappa_r}=0 and for \eqn{\sigma_r}=0}
#'   \item{Pvals}{the P values for the two chi-squared tests}
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
#' RL_gh_fit  <- PWD_get_gh(X,Y,printem=FALSE)
#' # run the residual analysis from the model output
#' post  <- PWD_resi(X, RL_gh_fit$resi, printem=TRUE)
#'
#' @references Hawkins DM and Kraker JJ. Precision Profile Weighted Deming
#' Regression for Methods Comparison, on *Arxiv* (2025, [arxiv.org/abs/2508.02888](https://arxiv.org/abs/2508.02888)).
#'
#' @references Hawkins DM (2014). A Model for Assay Precision.
#' *Statistics in Biopharmaceutical Research*, **6**, 263-269.
#' http://dx.doi.org/10.1080/19466315.2014.899511
#'
#' @importFrom stats optimize pchisq shapiro.test
#'
#' @export

PWD_resi    <- function(true, resi, epsilon=1e-5, printem=TRUE) {
  n       <- length(true)
  absres  <- abs(resi)
  key     <- order(true)
  sortres <- round(absres[key],4)
  sorttru <- true  [key]
  vars    <- sortres^2
  cvsq    <- (sortres/sorttru)^2
  logcv   <- round(0.5 * log(cvsq),4)
  lowr    <- 1:round(n/3)
  hir     <- round(2*n/3):n
  maxsig  <- max(sortres[lowr])
  maxkap  <- 0.5*max(abs(logcv[hir]), na.rm=TRUE)

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
  like  <- fc
  if (fd < fc) {
    sigma <- d
    kappa <- vd$minimum
    like  <- fd
  }

  profl   <- sqrt(sigma^2 + (kappa*true)^2)
  scalr   <- resi/profl

  proflf  <- sqrt(sigma^2 + (kappa*sorttru)^2)

  poolsig <- sqrt(mean(vars))
  poolkap <- sqrt(mean(cvsq))
  consdl  <- innr(0, poolsig)
  concvl  <- innr(poolkap, 0)
  check   <- innr(vc$minimum, c)
  tests   <- c(consdl, concvl) - fc
  Pvals   <- pchisq(tests, 1, lower.tail=FALSE)
  if (printem) {
    SW    <- shapiro.test(scalr)$p.value
    cat(sprintf("Rocke-Lorenzato fit to residuals\nsigma %6.4f kappa %6.4f\n",
                sigma, kappa))
    cat(sprintf("P values for constant sd %6.4f cv %6.4f normality %6.4f\n",
                Pvals[1], Pvals[2], SW))
  }

  return(list(sigmar=sigma, kappar=kappa, like=like, scalr=scalr,
              poolsig=poolsig, poolkap=poolkap, tests=tests, Pvals=Pvals))
}

