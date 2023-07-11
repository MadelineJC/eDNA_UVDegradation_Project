#### DOC STRING ===========================================
# This script proposes a model of eDNA degradation over time
# in a marine environment subject to UV exposure. The model
# will be used to simulate stochastic time series data (eDNA 
# concentration over time). 
#### END DOC STRING =======================================

#### Required packages ####
# install.packages("deSolve")
library(tidyverse) # For like, 1 function, very sorry, lol
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
lambda = 0.25
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
t <- seq(0, 10, 0.1)
set.seed(123)
lambdaDraws <- rnorm(n = length(t), mean = mean(lambda), sd = 0.04); max(lambdaDraws); min(lambdaDraws)

hist(lambdaDraws)
abline(v = mean(lambdaDraws), col = "firebrick2")

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
##  Where "1" is a day
# Sampling every hour for 10 days
ConstSamp <- seq(0, 10, 1/24); length(ConstSamp) 

# Sampling every hour for 2 days, every 6 hours for 2 days, every 12 hours for 2 days, and once/day for 2 days 
Scenario1 <- c(seq(0, 2, 1/24), seq(2 + 1/4, 4, 6/24), seq(4 + 1/2, 7, 12/24), seq(8, 10)); length(Scenario1)
plot(1, 1, type = "l", xlim = c(0, 10), yaxt = "n", main = "Scenario 1", xlab = "Day", ylab = " ")
abline(v = Scenario1, col = "firebrick2")

# Sampling every hour for 12 hours, every 4 hours for 2 days, every 12 hours for 5 days, and once/day for 2 days
Scenario2 <- c(seq(0, 1/2, 1/24), seq(0.5 + 4/24, 2 + 4/24, 4/24), seq(2.5, 7.5, 12/24), seq(8, 10, 1)); length(Scenario2)
plot(1, 1, type = "l", xlim = c(0, 10), yaxt = "n", main = "Scenario 2", xlab = "Day", ylab = " ")
abline(v = Scenario2, col = "firebrick2")

# Sampling every 4 hours for two days, every 12 hours for 4 days, and once/day for 2 days
Scenario3 <- c(seq(0, 2, 6/24), seq(2.5, 7, 12/24), seq(8, 10, 1)); length(Scenario3)
plot(1, 1, type = "l", xlim = c(0, 10), yaxt = "n", main = "Scenario 3", xlab = "Day", ylab = " ")
abline(v = Scenario3, col = "firebrick2")

# Sampling every 2 hours for half a day, every 6 hours for 2 days, every 12 hours for 5 days, and once/day for 2 days
Scenario4 <- c(seq(0, 0.5, 2/24), seq(0.5 + 6/24, 2, 6/24), seq(2.5, 7.5, 12/24), seq(8, 10, 1)); length(Scenario3)
plot(1, 1, type = "l", xlim = c(0, 10), yaxt = "n", main = "Scenario 3", xlab = "Day", ylab = " ")
abline(v = Scenario4, col = "firebrick2")

#### ... Constant lambda ####
## Constant sampling; 241 sampling events over 10 days
t <- ConstSamp
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/0.15)
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Constant sampling; 241 sampling events over 10 days", xlab = "Day")

## Scenario 1; 66 sampling events over 10 days
t <- Scenario1
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/0.15)
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Scenario 1; 66 sampling events over 10 days", xlab = "Day")

## Scenario 2; 37 sampling events over 10 days
t <- Scenario2
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/0.15)
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Scenario 2; 37 sampling events over 10 days", xlab = "Day")

## Scenario 3; 22 sampling events over 10 days
t <- Scenario3
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/0.15)
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Scenario 3; 22 sampling events over 10 days", xlab = "Day")

## Scenario 4; 22 sampling events over 10 days
t <- Scenario4
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/0.15)
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Scenario 4; 22 sampling events over 10 days", xlab = "Day")

#### ... Varying lambda ####
lambdaDraws <- rep(lambdaDraws, 3)
## Constant sampling; 241 sampling events over 10 days
t <- ConstSamp
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/lambdaDraws[i])
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Constant sampling; 241 sampling events over 10 days", xlab = "Day")

## Scenario 1; 66 sampling events over 10 days
t <- Scenario1
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/lambdaDraws[i])
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Scenario 1; 66 sampling events over 10 days", xlab = "Day")

## Scenario 2; 37 sampling events over 10 days
t <- Scenario2
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/lambdaDraws[i])
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Scenario 2; 37 sampling events over 10 days", xlab = "Day")

## Scenario 3; 22 sampling events over 10 days
t <- Scenario3
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/lambdaDraws[i])
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Scenario 3; 22 sampling events over 10 days", xlab = "Day")

## Scenario 4; 22 sampling events over 10 days
t <- Scenario4
N_prev <- 1; N_t <- c()
for (i in 1:length(t)){
  N_next <- N_prev*exp(-t[i]/lambdaDraws[i])
  N_t[i] <- N_next
  N_next <- N_prev
}
plot(N_t, xlim = c(0, length(t)), main = "Scenario 4; 22 sampling events over 10 days", xlab = "Day")

#### Model fitting ####
## NOTE: We should think about this because I could be missing something...
# Jaime used mle2 to fit a Laplace model, and people on the internet seem to think it's a good bet for an exponential model as well
# This seems promising: https://stats.stackexchange.com/questions/240455/fitting-exponential-regression-model-by-mle 
# I'm using the nls () function (nonlinear least squares), which determines the nonlinear least-squares estimates of the parameters of a nonlinear model.
#### ... lambda = 0.25 ####
#### ... ... Constant sampling ####
t <- ConstSamp # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1 ####
t <- Scenario1 # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4 ####
t <- Scenario4 # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... lambda = 0.15 ####
#### ... ... Constant sampling ####
t <- ConstSamp # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1 ####
t <- Scenario1 # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4 ####
t <- Scenario4 # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... lambda = 0.05 ####
#### ... ... Constant sampling ####
t <- ConstSamp # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1 ####
t <- Scenario1 # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4 ####
t <- Scenario4 # Sampling events
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### Missing data ####
#### ... lambda = 0.25 ####
#### ... ... Constant sampling: 80% of data ####
t <- ConstSamp # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Constant sampling: 50% of data ####
t <- ConstSamp # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1: 80% of data ####
t <- Scenario1 # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1: 50% of data ####
t <- Scenario1 # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4: 80% of data ####
t <- Scenario4 # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4: 50% of data ####
t <- Scenario4 # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.25 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... lambda = 0.15 ####
#### ... ... Constant sampling: 80% of data ####
t <- ConstSamp # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Constant sampling: 50% of data ####
t <- ConstSamp # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1: 80% of data ####
t <- Scenario1 # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1: 50% of data ####
t <- Scenario1 # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4: 80% of data ####
t <- Scenario4 # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4: 50% of data ####
t <- Scenario4 # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.15 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.05)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... lambda = 0.05 ####
#### ... ... Constant sampling: 80% of data ####
t <- ConstSamp # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Constant sampling: 50% of data ####
t <- ConstSamp # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1: 80% of data ####
t <- Scenario1 # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 1: 50% of data ####
t <- Scenario1 # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4: 80% of data ####
t <- Scenario4 # Sampling events
t <- sort(sample(t, length(t)*0.8)) # Remove 20% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### ... ... Scenario 4: 50% of data ####
t <- Scenario4 # Sampling events
t <- sort(sample(t, length(t)*0.5)) # Remove 50% of observations
y0_obs <- 1 # Starting value
r_obs <- 0.05 # Constant rate of decay
y_obs <- y0_obs*exp(-t/r_obs) # Generate observations
plot(t, y_obs, pch = 16, xlab = "Day", ylab = "Observed Data") # Viz

# Additive error
y_obs_err <- y0_obs*exp(-t/r_obs) + rnorm(length(t), sd = 0.1)
plot(t, y_obs_err, main = "Additive error", pch = 16)
lines(t, y0_obs*exp(-t/r_obs), lwd = 2) 

# Fitting with nls()
fit <- nls(y_obs_err ~ y0_est*exp(-t/r_est), 
           start = list(y0_est = 1, r_est = 0.25)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; r_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(t, resid(fit), pch = 16)
abline(h = 0, lty = 2)

y_est <- y0_est*exp(-t/r_est)

plot(t, y_obs, type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(t, y_obs_err, col = "cornflowerblue", pch = 16)
points(t, y_est, col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)

#### Including process/sampling error for barrel selection/sampling ####
#### ... Getting constant degradation rates for each of 9 barrels ####
# Required objects
min <- 0.05; max <- 0.25; mean <- (min + max)/2; print(mean)
num_barrels <- 9
# Taking a look at distribution
hist(rnorm(1000, mean = mean, sd = 0.1))
# Sampling 9 lambdas
lambdas <- c()
set.seed(123)
for (i in 1:num_barrels) {
  lambdas[i] <- round(rnorm(1, mean = mean, sd = 0.01), digits = 2)
  while(lambdas[i] < min | lambdas[i] > max)
  {lambdas[i] <- round(rnorm(1, mean = mean, sd = 0.01), digits = 2)}
}; print(lambdas); print(mean(lambdas))

#### ... Setting up sampling regimes for each of 9 barrels ####
all_times <- seq(0 + 2/24, 10, 2/24); num_samples <- round(length(all_times)/num_barrels) # Full time series divided by number of barrels; we can change this to make it more realistic
ts <- matrix(nrow = num_samples + 1, ncol = num_barrels) # Matrix of sampling regimes wherein each column contains the sampling regime for a different barrel
ts[1, ] <- rep(0, ncol(ts))
set.seed(123)
for (i in 1:num_barrels){
  t <- sort(sample(all_times, num_samples, replace = FALSE))
  ts[2:nrow(ts), i] <- t
  all_times <- all_times[!all_times %in% t]
}; colnames(ts) <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9")
t1 <- ts[ ,1]; t2 <- ts[ ,2]; t3 <- ts[ ,3]; t4 <- ts[ ,4]; t5 <- ts[ ,5]; t6 <- ts[ ,6]; t7 <- ts[ ,7]; t8 <- ts[ ,8]; t9 <- ts[ ,9] # Come up with an automated way of doing this
all_times <- sort(c(0, ts[2:nrow(ts), ]))

# Observation without error
y0_obs <- 1e6
y_obs <- y0_obs*exp(-all_times/mean(lambdas))
plot(y_obs)

# Setting up barrel volumes to implement sampling error
samp_vol <- 5; barrel_vol <- 100
vol <- c(barrel_vol)
for (i in 2:length(t1)){
  vol[i] <- vol[i - 1] - samp_vol
}

## Barrel 1
# Setting up initial abundance, and empty vector to fill with observations
y_prev <- 1e6; y1_obs_err <- c(); y1_obs_err[1] <- y_prev; j <- 1
# Loop to generate data with (1) process error, (2) environmental stochasticity, and (3) sampling error
for (i in 2:length(t1)){
  y_d <- rbinom(1, y_prev, exp(-(t1[i] - t1[i-1])/lambdas[1]*exp(-rnorm(1, 0, 0.1)))) # Determine eDNAs that "die" via exponential distribution (process error) and environmental stochasticity
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) # Sample eDNAs from remaining with sampling error via sampled volume vs total remaining volume
  y_next <- y_prev - y_d - y_select # eDNAs in barrel at next time step is equal to those available at previous time step, minus those that died, and those that were selected
  y_prev <- y_next # Re-assigning y_prev for next iteration
  y1_obs_err[i] <- y_select # Add sampled eDNAs to observation vector
}

## Barrel 2
y_prev <- 1e6; y2_obs_err <- c(); y2_obs_err[1] <- y_prev; j <- 1
for (i in 2:length(t2)){
  y_d <- rbinom(1, y_prev, exp(-(t2[i] - t2[i-1])/lambdas[2]*exp(-rnorm(1, 0, 0.1)))) 
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) 
  y_next <- y_prev - y_d - y_select 
  y_prev <- y_next 
  y2_obs_err[i] <- y_select 
}

## Barrel 3
y_prev <- 1e6; y3_obs_err <- c(); y3_obs_err[1] <- y_prev; j <- 1
for (i in 2:length(t3)){
  y_d <- rbinom(1, y_prev, exp(-(t3[i] - t3[i-1])/lambdas[3]*exp(-rnorm(1, 0, 0.1)))) 
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) 
  y_next <- y_prev - y_d - y_select 
  y_prev <- y_next 
  y3_obs_err[i] <- y_select 
}

## Barrel 4
y_prev <- 1e6; y4_obs_err <- c(); y4_obs_err[1] <- y_prev; j <- 1
for (i in 2:length(t4)){
  y_d <- rbinom(1, y_prev, exp(-(t4[i] - t4[i-1])/lambdas[4]*exp(-rnorm(1, 0, 0.1)))) 
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) 
  y_next <- y_prev - y_d - y_select 
  y_prev <- y_next 
  y4_obs_err[i] <- y_select 
}

## Barrel 5
y_prev <- 1e6; y5_obs_err <- c(); y5_obs_err[1] <- y_prev; j <- 1
for (i in 2:length(t5)){
  y_d <- rbinom(1, y_prev, exp(-(t5[i] - t5[i-1])/lambdas[5]*exp(-rnorm(1, 0, 0.1)))) 
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) 
  y_next <- y_prev - y_d - y_select 
  y_prev <- y_next 
  y5_obs_err[i] <- y_select 
}

## Barrel 6
y_prev <- 1e6; y6_obs_err <- c(); y6_obs_err[1] <- y_prev; j <- 1
for (i in 2:length(t6)){
  y_d <- rbinom(1, y_prev, exp(-(t6[i] - t6[i-1])/lambdas[6]*exp(-rnorm(1, 0, 0.1)))) 
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) 
  y_next <- y_prev - y_d - y_select 
  y_prev <- y_next 
  y6_obs_err[i] <- y_select 
}

## Barrel 7
y_prev <- 1e6; y7_obs_err <- c(); y7_obs_err[1] <- y_prev; j <- 1
for (i in 2:length(t7)){
  y_d <- rbinom(1, y_prev, exp(-(t7[i] - t7[i-1])/lambdas[7]*exp(-rnorm(1, 0, 0.1)))) 
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) 
  y_next <- y_prev - y_d - y_select 
  y_prev <- y_next 
  y7_obs_err[i] <- y_select 
}

## Barrel 8
y_prev <- 1e6; y8_obs_err <- c(); y8_obs_err[1] <- y_prev; j <- 1
for (i in 2:length(t8)){
  y_d <- rbinom(1, y_prev, exp(-(t8[i] - t8[i-1])/lambdas[8]*exp(-rnorm(1, 0, 0.1)))) 
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) 
  y_next <- y_prev - y_d - y_select 
  y_prev <- y_next 
  y8_obs_err[i] <- y_select 
}

## Barrel 9
y_prev <- 1e6; y9_obs_err <- c(); y9_obs_err[1] <- y_prev; j <- 1
for (i in 2:length(t9)){
  y_d <- rbinom(1, y_prev, exp(-(t9[i] - t9[i-1])/lambdas[9]*exp(-rnorm(1, 0, 0.1)))) 
  y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i]) 
  y_next <- y_prev - y_d - y_select 
  y_prev <- y_next 
  y9_obs_err[i] <- y_select 
}

# Combining generated data by barrel
b1 <- cbind(t1, y1_obs_err); b2 <- cbind(t2, y2_obs_err); b3 <- cbind(t3, y3_obs_err); b4 <- cbind(t4, y4_obs_err); b5 <- cbind(t5, y5_obs_err); b6 <- cbind(t6, y6_obs_err); b7 <- cbind(t7, y7_obs_err); b8 <- cbind(t8, y8_obs_err); b9 <- cbind(t9, y9_obs_err)
full_ts <- rbind(b1, b2[2:nrow(b2), ], b3[2:nrow(b3), ], b4[2:nrow(b4), ], b5[2:nrow(b5), ], b6[2:nrow(b6), ], b7[2:nrow(b7), ], b8[2:nrow(b8), ], b9[2:nrow(b9), ])
full_ts <- as.data.frame(full_ts); colnames(full_ts) <- c("t", "y_obs_err")
# Ordering generated data by sampling time
full_ts_ordered <- full_ts[with(full_ts, order(t)), ]
# Plotting generated data
plot(full_ts_ordered$t, full_ts_ordered$y_obs_err)

# Fit generated data to exponential model
y_obs_err <- full_ts_ordered$y_obs_err
fit <- nls(y_obs_err ~ y0_est*exp(-all_times/lambda_est), 
           start = list(y0_est = 1e6, lambda_est = 0.15)) # Making starting value for r_est bad ~on purpose~
coef(fit); y0_est <- coef(fit)[1]; lambda_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(all_times, resid(fit), pch = 16); abline(h = 0, lty = 2)

y0_est <- as.numeric(y0_est); lambda_est <- as.numeric(lambda_est)
y_est <- y0_est*exp(-all_times/lambda_est)
plot(y_est)

plot(all_times[1:20], y_obs[1:20], type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(all_times[1:20], y_obs_err[1:20], col = "cornflowerblue", pch = 16)
points(all_times[1:20], y_est[1:20], col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex=0.8)


