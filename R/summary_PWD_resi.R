#' Summary information of output from `PWD_resi`
#'
#' @description
#' The function, invisible to the user, is called by `summary()` to display the
#' formatted results of output from `PWD_resi`.
#'
#' @method summary pwdresi
#'
#' @param object    an object produced by a call to function  `PWD_resi`
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
#' # fit the model and store output
#' RL_gh_fit  <- PWD_get_gh(X,Y)
#' # run the residual analysis from the model output
#' post  <- PWD_resi(X, RL_gh_fit$resi)
#' summary(post)
#'
#' @importFrom stats qt shapiro.test
#'
#' @export

summary.pwdresi  <- function(object, digits=4, ...) {
  SW    <- shapiro.test(object$scalr)$p.value
  cat("Rocke-Lorenzato fit to residuals\n\t Estimates: sigma", round(object$sigma, digits),
      "kappa", round(object$kappa, digits), "\n")
  if(sum(object$whichmissing) > 0) cat("\t Fit on n =", sum(!object$whichmissing), "complete residuals\n")
  cat("P value for normality", round(SW, digits))

  class(object) <- "summary.pwdresi"
  invisible(object)
}
