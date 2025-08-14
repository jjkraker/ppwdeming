#' Weighted Deming Regression -- Outlier scanning
#' @name PWD_outlier
#'
#' @description
#' This function tests for outliers from the fitted regression, and refits on
#' a sanitized data set (with outliers removed).
#'
#' @usage
#' PWD_outlier(X, Y, K, lambda=1, Pcut=0.01, printem=TRUE)
#'
#' @param X		the vector of predicate readings,
#' @param Y		the vector of test readings,
#' @param K		the maximum number of outliers to seek,
#' @param lambda		*optional* the ratio of the X to the Y precision profile (defaults to 1),
#' @param Pcut		*optional*, default 0.01 (1%), cutoff for statistical significance of Bonferroni P,
#' @param printem	  *optional* - if TRUE, routine will print out results.
#'
#' @details
#' *To be added from paper*
#'
#' @returns A list containing the following components:
#'
#'   \item{ndrop}{the number of significant outliers}
#'   \item{drop}{a vector of the indices of the outliers}
#'   \item{cor }{the Pearson correlation between X and Y}
#'   \item{cleancor }{the Pearson correlation between cleaned X and Y (after outliers removed)}
#'   \item{scalr}{the scaled residuals of all cases from the sanitized fit}
#'   \item{keep}{logical vector identifying which cases retained in sanitized data set}
#'   \item{basepar}{the sigma, kappa, alpha, beta of the full data set}
#'   \item{lastpar}{the sigma, kappa, alpha, beta of the sanitized data set}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @example /Examples/PWD_outlier_man_example.R
#'
#' @references  Hawkins DM. *Outliers*. in Wiley Encyclopedia of Clinical Trials,
#' eds R. D’Agostino, L. Sullivan, and J. Massaro. 2008. Wiley, New York.
#'
#' @importFrom stats pnorm
#'
#' @export

PWD_outlier <- function(X, Y, K, lambda=1, Pcut=0.01, printem=TRUE) {
  outlis   <- NULL
  N        <- length(X)
  keep     <- rep(TRUE, N)
  mapper   <- 1:N
  # Forward identification of suspects

  if (printem) {
    cat(sprintf("Outlier identification\nFull sample fit\n"))
    do <- PWD_inference(X, Y, lambda, printem=TRUE)
  }
  for (m in 1:K) {
    x          <- X[keep]
    y          <- Y[keep]
    do         <- PWD_get_gh(x, y, lambda)
    printres   <- FALSE
    if (m == 1) {
      basepar  <- c(do$sigma, do$kappa, do$alpha, do$beta, do$like)
      lastpar  <- basepar                 # If there are no outliers
      printres <- TRUE & printem
    }

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
    if (printem) cat(sprintf("Suspect %3.0f outlier Z %5.2f \n", susp, scalr[inx]))
    outlis     <- c(outlis, susp)
    keep[susp] <- FALSE
    mapper     <- mapper[-inx]
  }
  x          <- X[keep]
  y          <- Y[keep]
  do         <- PWD_get_gh(x, y, lambda)

  # Now do reinclusion
  if (printem) cat(sprintf("Outlier reinclusion\n"))
  for (m in 1:K) {
    x          <- X[keep]
    y          <- Y[keep]
    do         <- PWD_get_gh(x, y, lambda)
    sigr         <- fitres$sigmar
    kapr         <- fitres$kappar
    alpha        <- do$alpha
    beta         <- do$beta

    resis        <- Y - alpha - beta*X
    profl        <- sqrt(sigr^2 + (beta*kapr*X)^2)
    scalr        <- resis/profl
    BonP         <- 2*(N-m+1)*pnorm(-abs(scalr))
    outlis       <- (1:N)[!keep]
    minout       <- min(abs(scalr[!keep]))
    susp         <- (1:N)[!keep & abs(scalr) == minout]
    Bonmax       <- max(BonP[!keep])

    if(printem) {
      echoit    <- data.frame(outlis, round(scalr[!keep],3), round(BonP[!keep],5))
      colnames(echoit) <- c("case", "outlier Z", "Bonferroni P")
      print(echoit)
      cat(sprintf("Least suspect %3.0f Z %5.3f BonP %6.4f\n", susp,minout,Bonmax))
    }

    if (Bonmax < Pcut) {
      if (printem) cat(sprintf("Any remaining suspects significant\n"))
      lastpar <- c(do$sigma, do$kappa, do$alpha, do$beta, do$like)

      break
    }
    cat(sprintf("Reinclude %1.0f\n", susp))
    keep[susp] <- TRUE      # Reinclude---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  }

  # End reinclusions

  ndrop  <- sum(!keep)
  listem <- (1:N)[!keep]
  if (ndrop > 0 & printem) {
    cat(sprintf("\nFit to retained clean cases\n"))
    dofir  <- PWD_inference(x, y, lambda=1, printem=TRUE)
  }

  corXY = cor(X,Y)
  cleancorxy = cor(x,y)
  return(list(ndrop=ndrop, drop=listem, cor=corXY, cleancor=cleancorxy, scalr=scalr,
              keep=keep, basepar=basepar, lastpar=lastpar))
}
