#' Estimate of Variance Profile Functions (proportional)
#' @name PWD_get_gh
#'
#' @description
#' This code estimates the variance profiles, assumed proportional, of
#' the Rocke-Lorenzato form;
#' also provides the resulting weighted Deming fit and residuals.
#'
#' @usage
#' PWD_get_gh(X, Y, lambda = 1,
#'            rho=lifecycle::deprecated(),
#'            alpha=lifecycle::deprecated(), beta=lifecycle::deprecated(),
#'            mu=lifecycle::deprecated(),
#'            quad=FALSE, epsilon = 1e-8,
#'            printem=lifecycle::deprecated())
#'
#' @param X		the vector of predicate readings.
#' @param Y		the vector of test readings.
#' @param lambda		*optional* (default of 1) - the ratio of the `X` to
#' the `Y` precision profile.
#' @param rho       `r lifecycle::badge("deprecated")` `rho = numvalue` initialization is no longer implemented in algorithm.
#' @param alpha     `r lifecycle::badge("deprecated")` `alpha = numvalue` initialization is no longer implemented in algorithm.
#' @param beta      `r lifecycle::badge("deprecated")` `beta = numvalue` initialization is no longer implemented in algorithm.
#' @param mu        `r lifecycle::badge("deprecated")` `mu = numvector` initialization is no longer implemented in algorithm.
#' @param quad      *optional* (default of FALSE) - logical, selects fitting a linear or a quadratic regression.
#' @param epsilon		*optional* (default of 1.e-8) - convergence tolerance limit.
#' @param printem `r lifecycle::badge("deprecated")` `printem = TRUE` is no
#' longer supported; separate summary functions are provided for printing output.
#'
#' @details
#' This workhorse routine optimizes the likelihood in the **unknown** *g*, *h*
#' setting over its *n*+4 or *n*+5 parameters:
#' the two Rocke-Lorenzato precision profile parameters \eqn{\sigma}
#' and \eqn{\kappa}, the intercept \eqn{\alpha} and slope \eqn{\beta}
#' (and coefficien of the squared term \eqn{\gamma}, if quadratic),
#' and the *n* latent true concentrations \eqn{\mu_i}.
#'
#' That is, the assumed forms are:
#'    * predicate precision profile model: \eqn{g_i = var(X_i) = \lambda\left(\sigma^2 + \left[\kappa\cdot \mu_i\right]^2\right)} and
#'    * test precision profile model: \eqn{h_i = var(Y_i) = \sigma^2 + \left[\kappa\cdot (\alpha + \beta\mu_i)\right]^2}.
#'
#' The search algorithm implements an efficient approach via reparameterization
#' to the ratio \eqn{\rho = \frac{\sigma}{\kappa}}.
#'
#' @returns A list containing the following components:
#'
#'   \item{alpha }{the fitted intercept}
#'   \item{beta }{the fitted slope}
#'   \item{gamma}{the fitted coefficient of the square: if quad is FALSE, the value is zero}
#'   \item{quad}{logical, FALSE for linear or TRUE for quadratic regression}
#'   \item{fity }{the vector of predicted Y}
#'   \item{mu }{the vector of estimated latent true values}
#'   \item{resi }{the vector of residuals}
#'   \item{sigma }{the estimate of the Rocke-Lorenzato \eqn{\sigma}}
#'   \item{kappa }{the estimate of the Rocke-Lorenzato \eqn{\kappa}}
#'   \item{L }{the -2 log likelihood L}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @examples
#' # library
#' library(ppwdeming)
#'
#' # parameter specifications
#' n <- 100
#'
#' quad  <- FALSE   # for a linear fit
#' quad  <- TRUE    # for a quadratic fit
#' sigma <- 1
#' kappa <- 0.08
#' alpha <- 1
#' beta  <- 1.1
#' gamma <- 0
#' if (quad) gamma <- -0.005
#' true  <- 8*10^((0:(n-1))/(n-1))
#' truey <- alpha+beta*true+gamma*true^2
#' # simulate single sample - set seed for reproducibility
#' set.seed(1039)
#' # specifications for predicate method
#' X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(n))
#' # specifications for test method
#' Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(n))
#'
#' # fit with RL precision profile to estimate parameters
#' RL_gh_fit  <- PWD_get_gh(X,Y, quad=quad)
#'
#' # results
#' summary(RL_gh_fit)
#'
#' @references Hawkins DM and Kraker JJ (2026). Precision Profile Weighted
#' Deming Regression for Methods Comparison.
#' *The Journal of Applied Laboratory Medicine*, **11**(2), 379-392.
#' <doi:10.1093/jalm/jfaf183>
#'
#' @references Rocke DM, Lorenzato S (2012). A Two-Component Model for Measurement
#' Error in Analytical Chemistry.  *Technometrics*, **37:2**:176-184.
#'
#' @importFrom stats optimize
#' @importFrom stats lm
#' @importFrom stats coef
#'
#' @export

PWD_get_gh <- function(X, Y, lambda = 1,
                       rho=lifecycle::deprecated(),
                       alpha=lifecycle::deprecated(), beta=lifecycle::deprecated(),
                       mu=lifecycle::deprecated(),
                       quad=FALSE, epsilon = 1e-8,
                       printem=lifecycle::deprecated()) {

# updated as per code in June 29 email

  if (lifecycle::is_present(printem)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_get_gh(printem)",
      details = "Argument printem no longer prints results-message.  \n See summary functions instead for printing output. \n Argument will be dropped in next release."
    )
  }
  if (lifecycle::is_present(rho)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_get_gh(rho)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(alpha)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_get_gh(alpha)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(beta)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_get_gh(beta)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  if (lifecycle::is_present(mu)) {
    lifecycle::deprecate_warn(
      when = "3.0.0",
      what = "PWD_get_gh(mu)",
      details = "Initialization arguments rho, alpha, beta, and mu are no longer \n implemented in algorithm. \n Arguments will be dropped in next release."
    )
  }
  rho=NA; alpha=NA; beta=NA; gamma=NA; mu=NA

  whichmissing <- (!complete.cases(X)) | (!complete.cases(Y))
  missingcases <- (1:length(X))[whichmissing]
  allX <- X
  allY <- Y
  X <- X[!whichmissing]
  Y <- Y[!whichmissing]
  if(sum(!is.na(mu)) > 0) mu <- mu[!whichmissing]

  n     <- length(X)
  tun   <- 2*n
  A     <- tun * (1-log(tun))
  ones  <- rep(1, n)

  innerfit <- function(rho, X, Y, lambda=1, alpha, beta, gamma, mu, quad, epsilon=1e-7) {
    if (is.na(sum(mu))) mu    <- X
    diffr <- 2*epsilon
    innr  <- 0
    while (diffr > epsilon & innr < 100) {
      innr <- innr + 1
      old  <- mu
      fity <- alpha + beta*mu + gamma*mu^2
      g    <- rho^2 + mu^2
      h    <- rho^2 + fity^2
      if (!quad) {
        mu   <- (h * X + g * beta * (Y - alpha))/(h + g * beta^2)
      } else {
        mu   <- NULL
        for (i in 1:n) {
          x    <- X[i]
          y    <- Y[i]-alpha
          gi   <- g[i]
          hi   <- h[i]
          ccof <- c(-(x/gi+beta*y/hi),1/gi+(beta^2-2*gamma*y)/hi,
                    3*gamma*beta/hi, 2*gamma^2/hi)
          cub  <- polyroot(ccof)
          sele <- abs(Im(cub))<1e-5
          soln <- Re(cub[sele])[1]    # Smallest real root
          mu <- c(mu, soln)
        }
      }
      #     Check for possible mu degeneracy
      lenrat <- sum(mu^2) / sum(X^2)
      alpha  <- 0
      beta   <- 1
      gamma  <- 0
      if (lenrat > 0.5)  {   # Avoid mu degenerating to a point.
        musq   <- mu^2
        ruter  <- 1/sqrt(h)
        Z      <- ruter * cbind(Y-mu, ones, mu)
        if (quad) Z <- cbind(Z, ruter*musq)
        ZTZ    <- t(Z) %*% Z
        ZTZI   <- solve(ZTZ)
        R      <- ZTZI[1,1]
        coff   <- -ZTZI[-1,1]/R
        alpha <- coff[1]
        beta  <- coff[2]+1
        if (quad) gamma <- coff[3]
      }
      diffr <- sum((mu - old)^2)/sum(mu^2)
    }
    fity <- alpha + beta*mu + gamma*mu^2
    W     <- sum((X-mu)^2/g+(Y-fity)^2/h)
    kappa <- sqrt(W/tun)
    sigma <- rho*kappa
    slgh  <- sum(log(g*h))
    L     <- tun*log(W) + slgh + A
    #	print(paste(L, alpha, beta, gamma))
    return(list(L=L, W=W, alpha=alpha, beta=beta, gamma=gamma,
                sigma=sigma, kappa=kappa,mu=mu, fity=fity, resi=Y-fity))
  }    #  end of innerfit

  wrapper <- function(logrho, X, Y, lambda, alpha, beta, gamma, mu, quad) {
    innerfit(exp(logrho), X, Y, lambda, alpha, beta, gamma, mu, quad)$L
  }


  # Preliminary ranging
  if (is.na(alpha + beta)) {
    #generate starting values
    fitlm <- lm(Y~X)
    alpha <- coef(fitlm)[1]
    beta  <- coef(fitlm)[2]
    gamma <- 0
  }

  if (is.na(sum(mu))) mu    <- X

  # Ranging on rho
  nmesh  <- 20    # possible
  rmesh  <- seq(-5, 5, length.out=nmesh)
  profil <- NULL
  for (rho in rmesh) {
    LL <- wrapper(rho,X, Y, lambda, alpha, beta, gamma, mu, quad)
    profil <- c(profil, LL)
  }
  #  plot(rmesh, profil, xlab="ln(rho)", ylab="l", type="l")
  minof <- (1:nmesh)[profil == min(profil)]
  lowval <- rmesh[max(1,minof-1)]
  hival  <- rmesh[min(nmesh, minof+1)]
  #  cat(sprintf("Rho space %6.4f %6.4f\n", exp(lowval), exp(hival)))
  # End ranging

  fitit <- optimize(wrapper, c(-5,5), X=X, Y=Y, lambda=lambda,
                    alpha=alpha, beta=beta, gamma=gamma, mu=mu, quad=quad)
  bestrho <- exp(fitit$minimum)
  do      <- innerfit(bestrho, X, Y, lambda, alpha, beta, gamma, mu, quad)

  #  print(do)
  L     <- do$L
  alpha <- do$alpha
  beta  <- do$beta
  gamma <- do$gamma
  sigma <- do$sigma
  kappa <- do$kappa
  mu    <- do$mu
  fity  <- do$fity
  resi  <- do$resi

  allresi = rep(NA, length(allX))
  allresi[!whichmissing] = resi
  allfity = rep(NA, length(allX))
  allfity[!whichmissing] = fity
  allmu = rep(NA, length(allX))
  allmu[!whichmissing] = mu

  fullout <- list(alpha=unname(alpha), beta=unname(beta), gamma=unname(gamma), quad=quad,
                  fity=fity,mu=mu, resi=resi,
                  rho=sigma/kappa, sigma=sigma, kappa=kappa, L=L)

  class(fullout) <- c("pwdgetgh")

  return(fullout)

}
