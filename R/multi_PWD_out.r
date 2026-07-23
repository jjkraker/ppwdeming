#' Locate and test for outliers of the fit
#' @name multi_PWD_out
#'
#' @description
#' The routine uses a conventional multivariate outlier test statistic set,
#' in the framework of Rosner’s multiple outlier search methodology,
#' for the Multiple-Instrument fit. The vectors tested as multivariate outliers
#' are the vectors of scaled residuals of each case.
#'
#' @usage
#' multi_PWD_out(X, K=NA, names=NA, refine=1, Pcut=0.01)
#'
#' @param X		the *n* by *I* matrix of values reported by the *I* instruments on the *n* samples.
#' @param K    *optional* (default of n/20) the number of cases to seek in the forward outlier screen.
#' @param names      *optional* (default of NA) - labels of the *I* instruments.
#' @param refine       *optional* (default of 1) - the number of re-substitutions used to refine lambda estimates.
#' @param Pcut   *optional* (default of 0.01) the Bonferroni significance level for the outlier test.
#'
#' @details  A forward selection step identifies the *K* sequentially
#' most extreme cases using scaled residuals from a sequence of fits.
#' This is done by finding the Mahalanobis distance of each case’s vector
#' of scaled residuals from the mean vector of all cases currently thought
#' to be "clean". The case with the largest such distance is moved from
#' the "clean" list to the "suspect" list and the model refitted using the
#' new "clean" list.  This is repeated until `K` cases are on the suspect list.
#'
#' Backward re-inclusion follows.  The model is fitted to the clean cases, and
#' the Hotelling \eqn{T^2} statistic between each suspect and the set of clean
#' cases is computed.  The case with the minimum \eqn{T^2} is found and its
#' Bonferroni P value computed.  If it is not significant, the case is moved
#' from the suspect back to the clean list.
#' The process is repeated until either all suspects have been exonerated or
#' a set of statistically significant outliers remains.
#'
#' @returns A list containing the following components:
#'
#'   \item{ndrop}{the number of statistically significant outliers. }
#'   \item{drop}{the identities of the significant outliers.}
#'   \item{initlis}{the *K* vector of cases identified as potential outliers in forward selection. }
#'   \item{Dilis}{the *K* vector of Mahalanobis distances of each initial suspect.}
#'   \item{Tilis}{the *K* vector of transforms of the Dilis to Hotelling \eqn{T^2}.}
#'   \item{backlis}{the vector of cases re-included in the backward step.}
#'   \item{backP}{their outlier *P* values.}
#'   \item{bonfP}{The Bonferroni *P* values of the remaining suspect cases.}
#'   \item{forward}{a print-ready listing of the forward selection steps.}
#'   \item{backward}{a print-ready listing of the backward re-inclusion steps.}
#'   \item{outlis}{a print-ready list of the significant outliers’ scaled residuals and Bonferroni *P* values.}
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
#' # add outliers
#' X[1,4] <- X[1,4] + 15
#' X[n,2] <- X[n,2] / 1.5
#'
#' # application
#' out <- multi_PWD_out(X, K=4, names=names)
#'
#' # full summary of process
#' summary(out)
#'
#' @importFrom stats cov
#' @importFrom stats pf
#'
#' @export

multi_PWD_out <- function(X, K=NA, names=NA, refine=1, Pcut=0.01) {
  n        <- dim(X)[1]
  I        <- dim(X)[2]

  if(any(is.na(names))) names <- 1:I
  if(is.na(K)) K = trunc(n/20)

  # Forward

  clean  <- 1:n
  drop   <- NULL
  Dilis  <- NULL
  Tilis  <- NULL
  nclean <- n
  ndrop  <- 0
  for (k in 1:K) {
    do     <- multi_PWD(X[clean,], refine=refine)
    alpha  <- do$alpha
    beta   <- do$beta
    lambda <- do$lambda
    if(k == 1) {
      mu        <- do$mu
      basealpha <- alpha
      basebeta  <- beta
      baselam   <- lambda
    }
    kappa  <- do$kappa
    rho    <- do$rho
    sigma  <- rho*kappa
    scalr  <- do$scalr
    x      <- scalr
    xbar   <- apply(x, 2, FUN=mean)
    devi   <- t(t(x)-xbar)[,-1]
    covx   <- cov(x[,-1])
    covi   <- solve(covx)
    M      <- diag(devi %*% covi %*% t(devi))
    maxM   <- max(M)
    whic   <- (1:n)[M == maxM][1]
    orig   <- clean[whic]
    drop   <- c(drop, orig)
    Dilis  <- c(Dilis, maxM)
    Ti     <- (nclean-2)*maxM/((nclean-1)^2/nclean-maxM)
    Tilis  <- c(Tilis, Ti)
    clean  <- clean[-whic]
    nclean <- nclean - 1
    ndrop  <- ndrop + 1
  }

  initlis <- drop
  forward <- cbind(initlis, round(Dilis,3), round(Tilis,3))
  colnames(forward) <- c("case", "Mahal D", "Hotel T")

  # Reinclude

  backlis <- NULL
  backP   <- NULL
  for (k in 1:K) {
    do     <- multi_PWD(X[clean,], refine=refine)
    alpha  <- do$alpha
    beta   <- do$beta
    kappa  <- do$kappa
    rho    <- do$rho
    sigma  <- rho*kappa
    lambda <- do$lambda
    scalr  <- NULL
    resi   <- NULL
    for (i in 1:I) {
      fit   <- alpha[i] + beta[i]*mu
      g     <- sqrt(lambda[i]*(sigma^2 + (kappa*fit)^2))
      resi  <- cbind(resi, X[,i]-fit)
      scalr <- cbind(scalr, (X[,i]-fit)/g)
    }
    xbar    <- apply(scalr[clean,], 2, FUN=mean)[-1]
    covmat  <- cov(scalr[clean,-1])
    covi    <- solve(covmat)
    x       <- scalr[drop,-1]
    devi    <- t(t(x)-xbar)
    if (ndrop > 1) {
      M     <- diag(devi %*% covi %*% t(devi))
    } else {
      M     <- t(devi) %*% covi %*% devi
    }
    minM    <- min(M)
    Tsq     <- (nclean/(nclean+1))*M
    F       <- (nclean-1-I+1)/((nclean-1)*(I-1))*Tsq
    Pval    <- pf(F, I-1, nclean-1-I+1+1, lower.tail=FALSE)
    bonfP   <- (nclean + 1) * Pval
    whic    <- (1:ndrop)[M == minM][1]
    orig    <- drop[whic]
    backlis <- c(backlis, orig)
    backP   <- c(backP, max(bonfP))
    if (max(bonfP) < Pcut) break
    drop    <- drop[-whic]
    clean   <- c(clean, orig)
    ndrop   <- ndrop - 1
    nclean  <- nclean + 1
  }
  backward <- cbind(backlis, round(backP,4))
  colnames(backward) <- c("case", "Bonf P")

  outscal <- scalr[drop,]
  if (ndrop == 1) outscal <- matrix(outscal, 1)
  if (ndrop == 0) {
    outlis <- NA
  } else {
    cn     <- c("case", "Bonf P", names)
    rn     <- drop
    outlis <- cbind(drop, signif(bonfP,5), round(outscal,3))
    colnames(outlis) <- cn
    rownames(outlis) <- rn
  }

  fullout <- list(ndrop=ndrop, drop=drop, initlis=initlis,
              Dilis=Dilis, Tilis=Tilis, backlis=backlis, backP=backP,
              forward=forward, backward=backward, Pcut=Pcut,
              BonP=bonfP, outlis=outlis)

  class(fullout) <- c("pwdout", "list")

  return(fullout)
}
