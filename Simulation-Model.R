#### DOC STRING ===========================================
# This script proposes a model of eDNA degradation over time
# in a marine environment subject to UV exposure. The model
# will be used to simulate stochastic time series data (eDNA 
# concentration over time). 
#### END DOC STRING =======================================

#### Required packages ####
# install.packages("deSolve")
library(deSolve)

#### Parameters ####
# From Strickler et al., 2015, degradation rate under UV-B exposure ranges between 0.05 and 0.25 when UV-B is between 2 and 50 kJ/m2/day

#### Exponential decay models ####
#### ... Base model ####
ExponentialDecay <- function(t, y, p){
  # Parms
  lambda <- p[1]
  # State variables
  N <- y[1]
  # Eqn
  dN = -lambda*N
  list(c(dN))
}
# Parm values
lambda = 0.4
parms <- c(lambda)
# Initial values
N0 <- 200
# Time point solutions
t <- seq(0.1, 40, 0.1)
# Simulation
results <- lsoda(N0, t, ExponentialDecay, parms, verbose = FALSE)
plot(results)

#### ... When lambda ranges between 0.05 and 0.25 ####
lambda = seq(0.05, 0.25, 0.01)
# Initial values
N0 <- 1
# Time point solutions
t <- seq(0.1, 40, 0.1)
# Simulation
df <- data.frame(t)
for (i in 1:length(lambda)){
  results <- lsoda(N0, t, ExponentialDecay, lambda[i], verbose = FALSE)
  df <- cbind(df, results[ ,2])
}
plot(1, 1,pch = "", xlim = c(0, t[length(t)]), ylim = c(0, 1), xlab = "Time", ylab = "Concentration")
for (i in 2:ncol(df)){
  lines(t, df[ ,i])
}

#### ... When lambda is pulled from a normal distribution ####
set.seed(123)
lambdaDraws <- rnorm(n = length(t), mean = mean(lambda), sd = 0.04); max(lambdaDraws); min(lambdaDraws)

hist(lambdaDraws)
abline(v = mean(lambdaDraws), col = "red")

# Constant lambda
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/0.15)
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, type = "l", xlim = c(0, 20))

# Varying lambda
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/lambdaDraws[i])
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, type = "l", xlim = c(0, 20))

#### Inconsistent sampling ####








