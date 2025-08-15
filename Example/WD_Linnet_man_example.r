# library
library(ppwdeming)

# parameter specifications
alpha <- 1
beta  <- 1.1
true  <- 8*10^((0:99)/99)
truey <- alpha+beta*true
kappa <- 0.1

# simulate single sample - set seed for reproducibility
set.seed(1039)
# specifications for predicate method
X     <- true *(1+kappa*rnorm(100))
# specifications for test method
Y     <- truey *(1+kappa*rnorm(100))

# fit with to estimate linear parameters
wd_fit <- WD_Linnet(X,Y,MDL=12)
cat("\nThe Linnet constant-CV estimated intercept is",
    signif(wd_fit$alpha,4), "and the estimated slope is",
    signif(wd_fit$beta,4), "\n")
