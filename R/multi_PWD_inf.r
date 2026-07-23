#' Jack-knife inference on the parameters of the Multiple Instrument fit
#' @name multi_PWD_inf
#'
#' @description
#' This code uses jackknifing of the Multiple-Instrument fit to get
#' standard errors of the regression parameters and lambda values.
#'
#' @usage
#' multi_PWD_inf(X, names=NA, refine=1)
#'
#' @param X		the *n* by *I* matrix of values reported by the *I* instruments on the *n* samples.
#' @param names      *optional* (default of NA) - labels of the *I* instruments.
#' @param refine       *optional* (default of 1) - the number of re-substitutions used to refine lambda estimates.
#'
#' @details  Uses jackknifing to get standard errors for the slopes `beta`,
#' intercepts `alpha`, and `lambda` values of the *I* instruments.
#'
#' @returns A list containing the following components:
#'
#'   \item{L}{the -2 log likelihood of the full-sample fit.}
#'   \item{rho}{the full-sample estimate of \eqn{\rho}.}
#'   \item{kappa}{the full-sample estimate of \eqn{\kappa}.}
#'   \item{lambda}{the vector of *I* estimated \eqn{\lambda_i}.}
#'   \item{lambdase}{the vector of *I* standard errors of the estimated \eqn{\lambda_i}.}
#'   \item{alpha}{the vector of *I* estimated intercepts.}
#'   \item{alphase}{the vector of *I* standard errors of the estimated intercepts.}
#'   \item{beta}{the vector of *I* estimated slopes.}
#'   \item{betase}{the vector of *I* standard errors of the estimated slopes.}
#'   \item{mu}{the vector of *n* estimated true concentrations.}
#'   \item{fity}{the *n* by *I* matrix of fitted values.}
#'   \item{resi}{the *n* by *I* matrix of residuals.}
#'   \item{whichmissing}{a logical vector identifying locations of row of X with missing values}
#'   \item{scalr}{the *n* by *I* matrix of scaled residuals.}
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
#' inf <- multi_PWD_inf(X, names=names)
#' inf$alpha; inf$beta
#'
#' # summary
#' summary(inf)
#'
#' @export

multi_PWD_inf <- function(X, names=NA, refine=1) {
  whichmissing <- !complete.cases(X)
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  X <- X[!whichmissing,]

  allalpha  <- NULL
  allbeta   <- NULL
  alllam    <- NULL
  n         <- dim(X)[1]
  I         <- dim(X)[2]
  if(any(is.na(names))) names <- 1:I
  for (dele in 0:n) {
    x       <- X
    if (dele > 0) x <- X[-dele,]
    do        <- multi_PWD(x, refine=refine)
    if (dele == 0) {
      baseL     <- do$L
      rho       <- do$rho
      kappa     <- do$kappa
      baselam   <- do$lambda
      basealpha <- do$alpha
      basebeta  <- do$beta
      mu        <- do$mu
      fity      <- do$fity
      resi      <- do$resi
      scalr     <- do$scalr
    } else {
      allalpha  <- rbind(allalpha, do$alpha )
      allbeta   <- rbind(allbeta , do$beta  )
      alllam    <- rbind(alllam  , do$lambda)
    }
  }

  tcut        <- qt(0.975, n-1)
  mult        <- sqrt(n-1)
  alphase     <- mult * apply(allalpha, 2, FUN=sd)
  betase      <- mult * apply(allbeta , 2, FUN=sd)
  lamse       <- mult * apply(alllam  , 2, FUN=sd)

  fullout <- list(alpha=basealpha, beta=basebeta,
                  fity=fity, mu=mu,
                  resi=resi,
                  rho=rho, sigma=rho*kappa, kappa=kappa, L=baseL,
                  lambda=baselam, selambda=lamse,
                  sealpha=alphase, sebeta=betase,
                  whichmissing=whichmissing,
                  scalr=scalr, names=names)

  class(fullout) <- c("pwdinf")

  return(fullout)
}
