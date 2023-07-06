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

#### Including process error for barrel selection ####
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
  lambdas[i] <- round(rnorm(1, mean = mean, sd = 0.1), digits = 2)
  while(lambdas[i] < min | lambdas[i] > max)
  {lambdas[i] <- round(rnorm(1, mean = mean, sd = 0.1), digits = 2)}
}; print(lambdas); print(mean(lambdas))

#### ... Setting up sampling regimes for each of 9 barrels ####
all_times <- seq(0 + 2/24, 10, 2/24); num_samples <- round(length(all_times)/num_barrels) # Full time series divided by number of barrels; we can change this to make it more realistic
ts <- matrix(nrow = num_samples, ncol = num_barrels) # Matrix of sampling regimes wherein each column contains the sampling regime for a different barrel
set.seed(123)
for (i in 1:num_barrels){
  t <- sort(sample(all_times, num_samples, replace = FALSE))
  ts[ ,i] <- t
  all_times <- all_times[!all_times %in% t]
}; colnames(ts) <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9")
t1 <- ts[ ,1]; t2 <- ts[ ,2]; t3 <- ts[ ,3]; t4 <- ts[ ,4]; t5 <- ts[ ,5]; t6 <- ts[ ,6]; t7 <- ts[ ,7]; t8 <- ts[ ,8]; t9 <- ts[ ,9] # Come up with an automated way of doing this



#### ... Setting up loop to generate eDNA concentration data for 9 barrels ####
# Let's try to do this for one barrel first

# General structure
N_prev <- 1; N_t <- c()
for (i in 1:length(ts[ ,1])){
  N_next <- N_prev*exp(-t1[i]/lambdas[1]) # Were t1 is the sampling regime for barrel 1, and lambda[1] is the constant rate of decay for barrel 1
  N_t[i] <- N_next
  N_next <- N_prev
}; cbind(t1, N_t)


# ?rbinom says that the arguments needed are n (number of observations), size (number of trials), and prob (probability of success on each trial) 
# What's the difference between n and size?
# n is how many runs of the simulation you do; i.e., the line below reads, Do 1 run, of N_prev coin flips, which exp(-t1[i]/lambdas[1]) probability of success.
N_d <- rbinom(1, N_prev, exp(-t1[i]/lambdas[1]))
N_next <- N_prev - N_d
N_prev <- N_next

N_prev <- 1e6; N_t <- c(); N_t[1] <- N_prev
for (i in 1:length(t1)){
  N_d <- rbinom(1, N_prev, exp(-t1[i]*lambdas[1]))
  N_next <- N_prev - N_d
  N_prev <- N_next
  N_t[i + 1] <- N_prev
}

## NOTE: t1[i] should be time elapsed

all_times <- sort(c(t1, t2, t3, t4, t5, t6, t7, t8, t9)); t_elap <- 2/24
N_prev <- 1e6; N_t <- c(); N_t[1] <- N_prev; j <- 1
for (i in all_times){
  if (i %in% t1){
    N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[1]))
  } else{
    if (i %in% t2){
      N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[2]))
    } else {
      if (i %in% t3){
        N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[3]))
      } else {
        if (i %in% t4){
          N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[4]))
        } else {
          if (i %in% t5){
            N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[5]))
          } else {
            if (i %in% t6){
              N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[6]))
            } else {
              if (i %in% t7){
                N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[7]))
              } else {
                if (i %in% t8){
                  N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[8]))
                } else {
                  if (i %in% t9){
                    N_d <- rbinom(1, N_prev, exp(-t_elap/lambdas[9]))
                  } else {print(":(")}
                }
              }
            }
          }
        }
      }
    }
  }
  N_next <- N_prev - N_d
  N_prev <- N_next
  N_t[j + 1] <- N_prev
  j <- j + 1
}; t <- c(0, all_times); data <- cbind(t, N_t)

data <- as.data.frame(data)
plot(data$t[1:25], data$N_t[1:25], type = "l")









