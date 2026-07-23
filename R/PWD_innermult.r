#' Parameter Estimation from rho
#' @name PWD_innermult
#'
#' @description
#' This is primarily a workhorse routine that PWD_multi calls and, given a
#' value of rho, estimates the other model parameters.
#'
#' @usage
#' PWD_innermult(X, rho, refine=1, eps=1e-9)
#'
#' @param X		the *n* by *I* matrix of values reported by the *I* instruments on the *n* samples.
#' @param rho the input value of \eqn{\rho = \frac{\sigma}{\kappa}}.
#' @param refine       *optional* (default of 1) - the number of resubstitutions used to refine lambda estimates.
#' @param eps      *optional* (default of 1e-9) - a convergence criterion.
#'
#' @details Using the specified \eqn{\rho}, the routine calculates the estimates of
#' the regression parameters and \eqn{\lambda}.
#'
#' @returns A list containing the following components:
#'
#'   \item{L}{the -2 log likelihood}
#'   \item{kappa}{the estimate of \eqn{\kappa}.}
#'   \item{alpha}{the vector of *I* estimated intercepts.}
#'   \item{beta}{the vector of *I* estimated slopes.}
#'   \item{lambda}{the vector of *I* estimated \eqn{\lambda_i}.}
#'   \item{mu}{the vector of *n* estimated true concentrations.}
#'   \item{fity}{the *n* by *I* matrix of fitted values.}
#'   \item{resi}{the *n* by *I* matrix of residuals.}
#'   \item{scalr}{the *n* by *I* matrix of scaled residuals.}
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
#' # application (to be added)
#'
#' @references .
#'
#' @export

PWD_innermult <- function(X, rho, refine=1, eps=1e-9)  {
  size   <- dim(X)
  n      <- size[1]
  I      <- size[2]
  nI     <- n*I

  # Initialize parameter estimates

  alpha  <- rep(0,I)
  beta   <- rep(1,I)
  lambda <- rep(1, I)
  mu     <- apply(X, 1, FUN=mean)
  ones   <- rep(1,n)
  for (lamloop in 0:refine) {
    g     <- NULL
    for (i in 1:I) {
      g  <- cbind(g , rho^2 + (alpha[i] + beta[i]*mu)^2)
    }
    differ <- 2*eps
    loopr  <- 0

    #   Loop to get alpha, beta, mu and kappa

    while(differ > eps & loopr < 100) {
      loopr  <- loopr + 1
      oldmu  <- mu
      V      <- NULL
      for (i in 1:I) {

        #       Fit WLS to get alpha beta

        ruter    <- 1/sqrt(g[,i])
        Z        <- ruter * cbind(X[,i]-mu, ones, mu)
        ZTZ      <- t(Z) %*% Z
        ZTZI     <- solve(ZTZ)
        R        <- ZTZI[1,1]   # 1/R is residual SS
        V        <- c(V, 1/R)
        coef     <- -ZTZI[2:3,1]/R
        alpha[i] <- coef[1]
        beta [i] <- coef[2]+1
      }

      #     Normalize alpha, beta

      alpha  <- alpha - mean(alpha)
      beta   <- beta  - mean(beta ) + 1

      #     Recalculate g

      for (i in 1:I) {
        fity   <- alpha[i]+beta[i]*mu
        g[,i]  <- rho^2 + fity^2
      }

      #     Recalculate mu

      for (j in 1:n) {
        denom <- lambda*g[j,]
        mu[j] <- sum((beta*(X[j,]-alpha)/denom)) / sum(beta^2/denom)
      }
      differ <- sum((oldmu-mu)^2)/sum(mu^2)
    }
    lambda <- V / mean(V)
    kappa  <- sqrt(mean(V/(lambda*n)))
    B <- sum(log(c(g))) + n*sum(log(lambda))
    L <- nI*(1-log(nI) + log(sum(V/lambda))) + B
  }

  # Wrapup

  fity  <- NULL
  resi  <- NULL
  scalr <- NULL
  for (i in 1:I) {
    fit   <- alpha[i] + beta[i]*mu
    fity  <- cbind(fity , fit)
    res   <- X[,i] - fit
    gi    <- lambda[i]*(rho^2+fit^2)
    scal1 <- res/(kappa*sqrt(gi))
    resi  <- cbind(resi, res)
    scalr <- cbind(scalr, scal1)
  }

  return(list(L=L, kappa = kappa, alpha=alpha, beta=beta,
              lambda=lambda, mu=mu, fity=fity, resi=resi, scalr=scalr))
}
