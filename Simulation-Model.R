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

#### Base exponential decay model ####
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
t <- seq(0.1, 20, 0.1)
# Simulation
results <- lsoda(N0, t, ExponentialDecay, parms, verbose = FALSE)
plot(results)

## OR ##
# N(t) = N0*exp^(-t/tau)

## When lambda ranges between 0.05 and 0.25
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

#### Model ####
