#' Summary information of output from `PWD_inference` or `multiPWD_inf`
#'
#' @description
#' The function, invisible to the user, is called by `summary()` to display the
#' formatted results of a model fitted by `PWD_inference` or by `multiPWD_inf`.
#'
#' @method summary pwdinf
#'
#' @param object    an object produced by a call to function  `PWD_inference` or `multiPWD_inf`
#' @param conflevel    *optional* (default of 0.95) - confidence level, between 0 to 1;
#'   value of `NA` will remove the confidence bounds.
#' @param digits    *optional* (default of 4) - number of decimal places to use for rounding displayed results.
#' @param ... Further arguments passed to or from other methods (required for
#'   compatibility with the generic \code{summary} function).
#'
#' @return `object` is returned invisibly
#'   per generic `summary` method conventions.
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
#' # fit with RL precision profile to estimate parameters and variability
#' RL_inf <- PWD_inference(X,Y,MDL=12)
#' summary(RL_inf)
#'
#' @importFrom stats qt
#'
#' @export

summary.pwdinf  <- function(object, conflevel=0.95, digits=4, ...) {

  n <- sum(!object$whichmissing)
  I <- length(object$alpha)
  if (is.numeric(conflevel)) {
    if (conflevel > 0 & conflevel < 1) {
      tcut <- qt(1-(1-conflevel)/2, n-1)
      CIalpha <- object$alpha[1] + c(-1,1)*tcut*object$sealpha[1]
      CIbeta <- object$beta[1] + c(-1,1)*tcut*object$sebeta[1]
      CIgamma <- c(ifelse(object$segamma == 0, NA,  object$gamma[1] - tcut*object$segamma[1]),
                   ifelse(object$segamma == 0, NA,  object$gamma[1] + tcut*object$segamma[1])
      )
    } else {
      stop("ERROR: conflevel must be number between 0 and 1.")
    }
  }

  ### add two-instrument estimates ###
  if (I == 1) {
    matcoef <- matrix(c(object$alpha, object$beta,
                        object$sealpha, object$sebeta),
                      ncol=2)

    if (object$segamma == 0) {
      rownames(matcoef) <- c("Intercept", "Slope")
      ests <- c(CIalpha[1], CIbeta[1])
      ses <- c(CIalpha[2], CIbeta[2])
    } else {
      matcoef <- rbind(matcoef, c(object$gamma, object$segamma))
      rownames(matcoef) <- c("Intercept", "Slope", "gamma")
      ests <- c(CIalpha[1], CIbeta[1], CIgamma[1])
      ses <- c(CIalpha[2], CIbeta[2], CIgamma[2])
    }
    colnames(matcoef) <- c("PWD Estimates", "SE")

    if (!is.na(conflevel)) {
      matcoef <- cbind(matcoef, ests, ses)
      colnames(matcoef) <- c("PWD Estimates", "SE", paste(" ", 100*conflevel, c("% lower bound","% upper bound"), sep=""))
    }
    cat("Two-Instrument Parameter Estimates:\n")
    print(round(matcoef,digits), na.print = "---")
  }

  ### print multiple-instrument estimates ###
  if (I > 1) {
    matcoef <- NULL
    for (i in 1:I) {
      alpha <- object$alpha[i]
      beta <- object$beta[i]
      lambda <- object$lambda[i]
      sealpha <- object$sealpha[i]
      sebeta <- object$sebeta[i]
      selambda <- object$selambda[i]

      veccoef <- c(alpha, sealpha,
                   beta, sebeta,
                   lambda, selambda)
      matcoef <- cbind(matcoef,veccoef)
    }
    rownames(matcoef) <- c("Intercept", " se alpha",
                           "Slope", " se beta",
                           "lambda", " se lambda")
    colnames(matcoef) <- object$names
    cat("\nMultiple Instrument Parameter Estimates:\n")
    print(round(matcoef,3), na.print = "---")
  }

  cat("\nVariance-parameter estimates are sigmahat=", round(object$sigma,digits),
      " and kappahat=", round(object$kappa,digits), ",\n", sep="")
  cat("\twith rhohat=", round(object$sigma/object$kappa,digits-1), "; ", sep="")
  cat(sprintf("\tat minimum -2 log likelihood = %7.3f.\n", object$L), sep="")
  if(sum(object$whichmissing) > 0) print(sprintf("\t Fit on n = %i complete readings", sum(!object$whichmissing)))

  class(object) <- "summary.pwdinf"
  invisible(object)
}
