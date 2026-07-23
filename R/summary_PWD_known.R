#' Summary information of output from `PWD_known`
#'
#' @description
#' The function, invisible to the user, is called by `summary()` to display the
#' formatted results of a model fitted by `PWD_known`.
#'
#' @method summary pwdknown
#'
#' @param object    an object produced by a call to function  `PWD_known`
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
#' alpha <- 1
#' beta  <- 1.1
#' true  <- 8*10^((0:99)/99)
#' truey <- alpha+beta*true
#' # forms of precision profiles
#' gfun    <- function(true, gparms) {
#'   gvals = gparms[1]+gparms[2]*true^gparms[3]
#'   gvals
#' }
#' hfun    <- function(true, hparms) {
#'   hvals = hparms[1]+hparms[2]*true^hparms[3]
#'   hvals
#' }
#'
#' # Loosely motivated by Vitamin D data set
#' g     <- 4e-16+0.07*true^1.27
#' h     <- 6e-2+7e-5*truey^2.2
#' # simulate single sample - set seed for reproducibility
#' set.seed(1039)
#' # specifications for predicate method
#' X     <- true +sqrt(g)*rnorm(100)
#' # specifications for test method
#' Y     <- truey+sqrt(h)*rnorm(100)
#'
#' # fit with to estimate linear parameters
#' pwd_known_fit <- PWD_known(X, Y, gfun, hfun,
#'                            gparms=c(4e-16, 0.07, 1.27),
#'                            hparms=c(6e-2, 7e-5, 2.2), MDL=12)
#' summary(pwd_known_fit)
#'
#' @importFrom stats qt
#'
#' @export

summary.pwdknown  <- function(object, digits=4, ...) {
  n <- length(object$whichmissing)

  matcoef <- matrix(c(object$alpha, object$beta,
                      object$sealpha, object$sebeta,
                      object$CIalpha[1], object$CIbeta[1],
                      object$CIalpha[2], object$CIbeta[2]),
                    ncol=4)
  rownames(matcoef) <- c("Intercept", "Slope")
  colnames(matcoef) <- c(" Estimate", " Std. Error", paste(" ", object$clevelused,"% lower bound", sep=""),
                         paste(" ", object$clevelused,"% upper bound", sep=""))
  print(round(matcoef,3), na.print = "---")

  cat(sprintf("\t with minimum -2 log likelihood =  %7.3f \n", object$L))
  if(sum(object$whichmissing) > 0) print(sprintf("\t Fit on n = %i complete readings", sum(!object$whichmissing)))

  if (!is.na(sum(object$MDL))) {
    predcoef <- matrix(c(object$preMDL,
                         object$sepreMDL,
                         object$preMDLl,
                         object$preMDLu),ncol=4)

    rownames(predcoef) <- c(paste("        MDL=", object$MDL, sep=""))
    colnames(predcoef) <- c(" Estimate", " Std. Error", paste(" ", object$clevelused,"% lower bound", sep=""),
                            paste(" ", object$clevelused,"% upper bound", sep=""))
    cat("\nPredictions at: \n")
    print(round(predcoef, digits), na.print = "---")
  }

  class(object) <- "summary.pwdknown"
  invisible(object)
}
