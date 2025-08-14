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
# add some outliers
Y[c(1,2,100)] <- Y[c(1,2,100)] + c(7,4,-45)

# check for outliers, re-fit, and store output
outliers_assess <- PWD_outlier(X,Y,K=5)
