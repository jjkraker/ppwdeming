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

# fit with RL precision profile to estimate parameters and variability
RL_inf <- PWD_inference(X,Y,MDL=12,printem=TRUE)
