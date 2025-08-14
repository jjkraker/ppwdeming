# library
library(ppwdeming)

# parameter specifications
alpha <- 1
beta  <- 1.1
true  <- 8*10^((0:99)/99)
truey <- alpha+beta*true
# forms of precision profiles
gfun    <- function(true, gparms) {
  gvals = gparms[1]+gparms[2]*true^gparms[3]
  gvals
}
hfun    <- function(true, hparms) {
  hvals = hparms[1]+hparms[2]*true^hparms[3]
  hvals
}

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
pwd_known_fit <- PWD_known(X, Y, gfun, hfun,
                           c(4e-16, 0.07, 1.27), c(6e-2, 7e-5, 2.2))
