#' Summary information of output from `PWD_outlier` or `multi_PWD_out`
#'
#' @description
#' The function, invisible to the user, is called by `summary()` to display the
#' formatted results of a model fitted by `PWD_outlier` or by `multi_PWD_out`.
#'
#' @method summary pwdout
#'
#' @param object    an object produced by a call to function  `PWD_outlier` or `multi_PWD_out`
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
#' set.seed(1069)
#' # specifications for predicate method
#' X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(100))
#' # specifications for test method
#' Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(100))
#' # add some outliers
#' Y[c(1,2,100)] <- Y[c(1,2,100)] + c(-10,9,-50)
#'
#' # check for outliers, re-fit, and store output
#' outliers_assess <- PWD_outlier(X, Y, K=5)
#' summary(outliers_assess)
#'
#' @importFrom stats qt
#'
#' @export

summary.pwdout  <- function(object, ...) {
  cat("In order of identification, potential outliers are:\n")
  print(object$forward)
  cat("\nThe cases considered for re-inclusion are:\n")
  print(object$backward)
  cat("where the final case is determined to be an outlier,\n at level of significance", object$Pcut, ":\n")
  print(object$backward[nrow(object$backward),])
  cat("\nThus, the final", object$ndrop, "cases are identified as outliers:\n")
  print(object$outlis)

  class(object) <- "summary.pwdout"
  invisible(object)
}
