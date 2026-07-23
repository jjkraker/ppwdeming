#' Multiple Instrument Weighted Deming Regression
#' @name PWD_multi
#'
#' @description
#' This code estimates the Rocke-Lorenzato precision profile parameters and the
#' regression of each instrument’s values on the underlying latent true
#' concentration.
#'
#' @usage
#' PWD_multi(X, names, refine=1)
#'
#' @param X		the *n* by *I* matrix of values reported by the *I* instruments on the *n* samples.
#' @param names      *optional* (default of NA) - labels of the *I* instruments.
#' @param refine       *optional* (default of 1) - the number of resubstitutions used to refine lambda estimates.
#'
#' @details  The function optimizes the likelihood over the parameter `rho`, calling `PWD_innermult` to provide the estimates of the other parameters.
#'
#' @returns A list containing the following components:
#'
#'   \item{L}{the -2 log likelihood}
#'   \item{rho}{the estimate of \eqn{\tho}.}
#'   \item{kappa}{the estimate of \eqn{\kappa}.}
#'   \item{lambda}{the vector of *I* estimated \eqn{\lambda_i}.}
#'   \item{alpha}{the vector of *I* estimated intercepts.}
#'   \item{beta}{the vector of *I* estimated slopes.}
#'   \item{mu}{the vector of *n* estimated true concentrations.}
#'   \item{fity}{the *n* by *I* matrix of fitted values.}
#'   \item{resi}{the *n* by *I* matrix of residuals.}
#'   \item{scalr}{the *n* by *I* matrix of scaled residuals.}
#'   \item{results}{a printer-ready data frame of the fitted regression}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @examples
#' # library
#' library(ppwdeming)
#'
#' # parameter specifications
#' n  <- 120
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
#' # application
#' fit <- PWD_multi(X)
#' with(fit, cat(sprintf("L %6.3f rho %6.4g kappa %6.4f\n", L, rho, kappa)))
#' print(fit$results)
#'
#' @references .
#'
#' @export

PWD_multi <- function(X, names=NA, refine=1) {

  mwrap <- function(logrho, X) {
    PWD_innermult(X, exp(logrho), refine=refine)$L
  }
  opper   <- optimize(mwrap, log(c(1e-5, 1e5)), X)
  rho     <- exp(opper$minimum)
  do      <- PWD_innermult(X, rho, refine=refine)
  results <- rbind(signif(do$alpha,4), round(do$beta,4), round(do$lambda,3))
  if(any(is.na(names))) names <- 1:dim(X)[2]
  colnames(results) <- names
  rownames(results) <- c("Intercept", "slope", "lambda")
  return(list(L=do$L, rho=rho, kappa=do$kappa,
              lambda=do$lambda, alpha=do$alpha, beta=do$beta,
              mu=do$mu, fity=do$fity, resi=do$resi, scalr=do$scalr, results=results))
}
