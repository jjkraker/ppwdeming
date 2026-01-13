#' Weighted Deming Regression -- Outlier scanning
#' @name PWD_outlier
#'
#' @description
#' This function tests for outliers from the fitted regression, and refits on
#' a sanitized data set (with outliers removed).
#'
#' @usage
#' PWD_outlier(X, Y, K, lambda=1, Pcut=0.01, rho=NA, alpha=NA, beta=NA, mu=NA,
#'             printem=FALSE)
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
#' @param printem	  *optional* (default of FALSE) - if TRUE, routine will print out results as a `message`.
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
#'   \item{keep}{logical vector identifying which cases retained in sanitized data set}
#'   \item{basepar}{the sigma, kappa, alpha, beta of the full data set}
#'   \item{lastpar}{the sigma, kappa, alpha, beta of the sanitized data set}
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
#' # add some outliers
#' Y[c(1,2,100)] <- Y[c(1,2,100)] + c(7,4,-45)
#'
#' # check for outliers, re-fit, and store output
#' \donttest{outliers_assess <- PWD_outlier(X, Y, K=5, printem=TRUE)}
#'
#' @references Hawkins DM and Kraker JJ (in press). Precision Profile Weighted
#' Deming Regression for Methods Comparison. *The Journal of Applied Laboratory Medicine*.
#' <doi:10.1093/jalm/jfaf183>
#'
#' @references  Hawkins DM (2008). *Outliers* in Wiley Encyclopedia of Clinical Trials,
#' eds R. D’Agostino, L. Sullivan, and J. Massaro. Wiley, New York.
#'
#' @importFrom stats pnorm
#' @importFrom stats complete.cases
#' @importFrom stats na.omit
#' @importFrom utils capture.output
#'
#' @export

PWD_outlier <- function(X, Y, K, lambda=1, Pcut=0.01, rho=NA, alpha=NA, beta=NA, mu=NA, printem=FALSE) {
  whichmissing <- (!complete.cases(X)) | (!complete.cases(Y))
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  allY <- Y
  X <- X[!whichmissing]
  Y <- Y[!whichmissing]
  if(sum(!is.na(mu)) > 0) mu <- mu[!whichmissing]

  outlis   <- NULL
  N        <- length(X)
  keep     <- rep(TRUE, N)
  mapper   <- 1:N

  # Get initial fullsample fit
  full  <- PWD_get_gh(X, Y, lambda, rho, alpha, beta, mu)
  alpha <- full$alpha
  beta  <- full$beta
  rho   <- full$rho
  mu    <- full$mu
  if (printem) {
    message(sprintf("Outlier identification\nFull sample fit"))
    do <- PWD_inference(X, Y, lambda, rho, alpha, beta, mu, printem=FALSE)
  }
  # Forward identification of suspects
  for (m in 1:K) {
    x          <- X [keep]
    y          <- Y [keep]
    mud        <- mu[keep]
    do         <- PWD_get_gh(x, y, lambda, rho, alpha, beta, mud)
    printres   <- FALSE
    if (m == 1) {
      basepar  <- c(full$sigma, full$kappa, full$alpha, full$beta, full$L)
      lastpar  <- basepar                 # If there are no outliers
      printres <- TRUE & printem
    }

    rho        <- do$rho
    alpha      <- do$alpha
    beta       <- do$beta
    mu         <- do$mu
    resi       <- y - alpha - beta*x

    fitres     <- PWD_resi(x, resi, printem=printres)
    sigr       <- fitres$sigmar
    kapr       <- fitres$kappar
    scalr      <- fitres$scalr

    profl      <- resi/scalr
    ascal      <- abs(scalr)
    maxa       <- max(ascal)
    inx        <- (1:N)[ascal == max(ascal)][1]
    susp       <- mapper[inx]

    Bonmin     <- 2*(N-m+1)*pnorm(-maxa)
    if (printem) message(sprintf("Suspect %3.0f outlier Z %5.2f ", susp, scalr[inx]))
    outlis     <- c(outlis, susp)
    keep[susp] <- FALSE
    mapper     <- mapper[-inx]
  }
  x          <- X[keep]
  y          <- Y[keep]
  do         <- PWD_get_gh(x, y, lambda, rho, alpha, beta, mud)

  # Now do reinclusion
  if (printem) message(sprintf("\nOutlier reinclusion"))
  for (m in 1:K) {
    x          <- X[keep]
    y          <- Y[keep]
    do         <- PWD_get_gh(x, y, lambda, rho, alpha, beta, mud)
    sigr       <- fitres$sigmar
    kapr       <- fitres$kappar
    alpha      <- do$alpha
    beta       <- do$beta

    resis      <- Y - alpha - beta*X
    profl      <- sqrt(sigr^2 + (beta*kapr*X)^2)
    scalr      <- resis/profl
    BonP       <- 2*(N-m+1)*pnorm(-abs(scalr))
    outlis     <- (1:N)[!keep]
    minout     <- min(abs(scalr[!keep]))
    susp       <- (1:N)[!keep & abs(scalr) == minout]
    Bonmax     <- max(BonP[!keep])

    if(printem) {
      echoit  <- data.frame(outlis, X[!keep], Y[!keep],
                            round(scalr[!keep],3), round(BonP[!keep],5))
      colnames(echoit) <- c("case", "X", "Y", "outlier Z", "Bonferroni P")
      message(paste0(capture.output(echoit), collapse = "\n"))
      message(sprintf("Least suspect %3.0f Z %5.3f BonP %6.4f", susp,minout,Bonmax))
    }

    if (Bonmax < Pcut) {
      if (printem) message(sprintf("Any remaining suspects significant"))
      lastpar <- c(do$sigma, do$kappa, do$alpha, do$beta, do$L)

      break
    }
    if (printem) message(sprintf("Reinclude %1.0f\n", susp))
    keep[susp] <- TRUE      # Reinclude----
  }

  # End reinclusions

  ndrop  <- sum(!keep)
  listem <- (1:N)[!keep]
  if (ndrop > 0 & printem) {
    message(sprintf("\nFit to retained clean cases"))
    dofir  <- PWD_inference(x, y, lambda=1, rho, alpha, beta, mu[keep], printem=TRUE)
    if(sum(whichmissing) > 0) message(sprintf("\t\t From among n = %i complete readings", sum(!whichmissing)))
  }

  corXY = cor(X,Y)
  cleancorxy = cor(x,y)

  allscalr = rep(NA, length(allX))
  allscalr[!whichmissing] = scalr
  allkeep = rep(NA, length(allX))
  allkeep[!whichmissing] = keep
  alllistem <- as.numeric(na.omit((1:length(allX))[!allkeep]))

  return(list(ndrop=ndrop, drop=alllistem, cor=corXY, cleancor=cleancorxy, scalr=allscalr,
              keep=allkeep, basepar=as.numeric(basepar[1:4]), lastpar=as.numeric(lastpar[1:4])))
}
