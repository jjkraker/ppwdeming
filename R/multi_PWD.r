#' Multiple Instrument Weighted Deming Regression
#' @name multi_PWD
#'
#' @description
#' This code estimates the Rocke-Lorenzato precision profile parameters and the
#' regression of each instrument’s values on the underlying latent true
#' concentration.
#'
#' @usage
#' multi_PWD(X, names, refine=1)
#'
#' @param X		the *n* by *I* matrix of values reported by the *I* instruments on the *n* samples.
#' @param names      *optional* (default of NA) - labels of the *I* instruments.
#' @param refine       *optional* (default of 1) - the number of resubstitutions used to refine lambda estimates.
#'
#' @details  The function optimizes the likelihood over the parameter `rho`,
#' the ratio of the Rocke-Lorenzato parameters sigma to kappa: \eqn{\rho = \frac{\sigma}{\kappa}}.
#' Calls `multi_PWD_inner` to provide the estimates of the other parameters.
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha}{the vector of *I* estimated intercepts.}
#'   \item{beta}{the vector of *I* estimated slopes.}
#'   \item{fity}{the *n* by *I* matrix of fitted values.}
#'   \item{mu}{the vector of *n* estimated true concentrations.}
#'   \item{resi}{the *n* by *I* matrix of residuals.}
#'   \item{rho}{the estimate of \eqn{\rho}.}
#'   \item{sigma}{the estimate of \eqn{\sigma}.}
#'   \item{kappa}{the estimate of \eqn{\kappa}.}
#'   \item{lambda}{the vector of *I* estimated \eqn{\lambda_i}.}
#'   \item{L}{the -2 log likelihood}
#'   \item{scalr}{the *n* by *I* matrix of scaled residuals.}
#'   \item{whichmissing}{a logical vector identifying locations of row of X with missing values}
#'   \item{results}{a printer-ready data frame of the fitted regression}
#'   \item{names}{labels of the *I* instruments.}
#'
#' @references Hawkins DM and Kraker JJ. Multiple Instrument Methods Comparison
#' by Precision weighted Deming Regression,
#' on *Arxiv* (2026) <doi:10.48550/arXiv.2607.11776>
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @examples
#' # library
#' library(ppwdeming)
#'
#' # parameter specifications
#' n  <- 100
#' I  <- 4
#' alpha  <- c(0  , 0  , 2  , -2 )
#' beta   <- c(0.9, 1.1, 1.2, 0.8)
#' names   <- c("Inst_1", "Inst_2", "Inst_3", "Inst_4")
#' lambda <- c(1.8, 0.45, 0.9, 0.9)
#' rutlam <- sqrt(lambda)
#' sigma  <- 2
#' kappa  <- 0.08
#' true   <- 8*10^((0:(n-1))/(n-1))
#'
#' # simulated data
#' set.seed(1039)
#' X     <- NULL
#' for (i in 1:I)  {
#'   truey <- alpha[i] + beta[i]*true
#'   X <- cbind(X, truey+rutlam[i]*(sigma*rnorm(n) + truey*kappa*rnorm(n)))
#' }
#'
#' # fit with RL precision profile to estimate parameters
#' fit <- multi_PWD(X, names=names)
#'
#' # results
#' summary(fit)
#'
#' @references .
#'
#' @export

multi_PWD <- function(X, names=NA, refine=1) {

  whichmissing <- !complete.cases(X)
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  X <- X[!whichmissing,]

  mwrap <- function(logrho, X) {
    multi_PWD_inner(X, exp(logrho), refine=refine)$L
  }
  opper   <- optimize(mwrap, log(c(1e-5, 1e5)), X)
  rho     <- exp(opper$minimum)
  do      <- multi_PWD_inner(X, rho, refine=refine)
  results <- rbind(signif(do$alpha,4), round(do$beta,4), round(do$lambda,3))
  if(any(is.na(names))) names <- 1:dim(X)[2]
  colnames(results) <- names
  rownames(results) <- c("Intercept", "slope", "lambda")

  fullout <- list(alpha=unname(do$alpha), beta=unname(do$beta), fity=do$fity,
                  mu=do$mu, resi=do$resi,
                  rho=rho, sigma=rho*do$kappa, kappa=do$kappa,
                  lambda=do$lambda,
                  L=do$L, scalr=do$scalr,
                  whichmissing=whichmissing, results=results, names=names)

  class(fullout) <- c("pwdgetgh")

  return(fullout)

}
