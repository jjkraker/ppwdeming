# library
library(ppwdeming)

# parameter specifications
sigma <- 1
kappa <- 0.08
alpha <- 1
beta  <- 1.1
true  <- 8*10^((0:99)/99)
truey <- alpha+beta*true
# simulate single sample - set seed for reproducibility
set.seed(1039)
# specifications for predicate method
X     <- sigma*rnorm(100)+true *(1+kappa*rnorm(100))
# specifications for test method
Y     <- sigma*rnorm(100)+truey*(1+kappa*rnorm(100))

# fit RL with given sigma and kappa
RL_results <- PWD_RL(X,Y,sigma,kappa)
cat("\nWith given sigma and kappa, the estimated intercept is",
    signif(RL_results$alpha,4), "and the estimated slope is",
    signif(RL_results$beta,4), "\n")
