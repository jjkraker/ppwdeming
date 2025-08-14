# library
library(ppwdeming)

# parameter specifications
alpha <- 1
beta  <- 1.1
true  <- 8*10^((0:99)/99)
truey <- alpha+beta*true
# Loosely motivated by Vitamin D data set
g     <- 4e-16+0.07*true^1.27
h     <- 6e-2+7e-5*truey^2.2

# simulate single sample - set seed for reproducibility
set.seed(1039)
# specifications for predicate method
X     <- true +sqrt(g)*rnorm(100)
# specifications for test method
Y     <- truey+sqrt(h)*rnorm(100)

# fit with to estimate linear parameters
wd_fit <- WD_General(X,Y,g,h)
cat("\nWith given g and h, the estimated intercept is",
    signif(wd_fit$alpha,4), "and the estimated slope is",
    signif(wd_fit$beta,4), "\n")
