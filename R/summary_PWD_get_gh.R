#' Summary information of output from `PWD_get_gh` or `multi_PWD`
#'
#' @description
#' The function, invisible to the user, is called by `summary()` to display the
#' formatted results of a model fitted by `PWD_inference` or by `multi_PWD`.
#'
#' @method summary pwdgetgh
#'
#' @param object    an object produced by a call to function  `PWD_inference` or `multi_PWD`
#' @param digits    *optional* (default of 6) - number of decimal places to use for rounding displayed results.
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
#' quad  <- FALSE   # for a linear fit
#' quad  <- TRUE    # for a quadratic fit
#' sigma <- 1
#' kappa <- 0.08
#' alpha <- 1
#' beta  <- 1.1
#' gamma <- 0
#' if (quad) gamma <- -0.005
#' true  <- 8*10^((0:99)/99)
#' truey <- alpha+beta*true+gamma*true^2
#' # simulate single sample - set seed for reproducibility
#' set.seed(1039)
#' # specifications for predicate method
#' X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(100))
#' # specifications for test method
#' Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(100))
#'
#' # fit with RL precision profile to estimate parameters
#' RL_gh_fit  <- PWD_get_gh(X,Y, quad=quad)
#'
#' # results
#' summary(RL_gh_fit)
#'
#' @importFrom stats qt
#'
#' @export

summary.pwdgetgh <- function(object, digits=6, ...) {

  cat("\nVariance-parameter estimates are sigmahat=", round(object$sigma,digits),
      " and kappahat=", round(object$kappa,digits), ",\n", sep="")
  cat("\twith rhohat=", round(object$sigma/object$kappa,digits-1), "; ", sep="")
  cat(sprintf("\tat minimum -2 log likelihood = %7.3f.\n", object$L), sep="")

  matcoef <- rbind(object$alpha, object$beta)

  if (!is.null(object$gamma)) {
    matcoef <- rbind(matcoef,ifelse(object$gamma==0, NA, object$gamma))
    rownames(matcoef) <- c("Intercept", "Slope", "gamma")} else {
      rownames(matcoef) <- c("Intercept", "Slope")
    }

  if (length(object$alpha) == 1) {
    colnames(matcoef) <- "estimates" } else {
      colnames(matcoef) <- object$names
    }

  cat("With point estimates of regression coefficients:\n")
  print(round(matcoef,digits), na.print = "---")

  class(object) <- "summary.pwdgetgh"
  invisible(object)
}
