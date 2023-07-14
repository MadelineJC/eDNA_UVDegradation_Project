#### [INSERT A DOC STRING HERE] ####

library(tidyverse)

#### Getting constant degradation rates for each of 9 barrels ####
# Required objects
min <- 0.05; max <- 0.34; mean <- (min + max)/2; print(mean)
num_barrels <- 9
# Taking a look at distribution
hist(rnorm(1000, mean = mean, sd = 0.1))
# Sampling 9 decay rates
r_vec <- c()
set.seed(123)
for (i in 1:num_barrels) {
  r_vec[i] <- round(rnorm(1, mean = mean, sd = 0.01), digits = 2)
  while(r_vec[i] < min | r_vec[i] > max)
  {r_vec[i] <- round(rnorm(1, mean = mean, sd = 0.01), digits = 2)}
}; print(r_vec); print(mean(r_vec))

#### Setting up sampling regimes for each of 9 barrels ####
all_times <- seq(0 + 12/24, 60, 12/24); num_samples <- floor(length(all_times)/num_barrels) # Full time series (in days) divided by number of barrels; we can change this to make it more realistic; used floor() instead of round() to always round down
# all_times <- seq(1, 60, 1); num_samples <- floor(length(all_times)/num_barrels) # DAYS
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

#### Observation without error ####
y0_obs <- 1e6
y_obs <- y0_obs*exp(-all_times*mean(r_vec))
plot(x = all_times, y = y_obs)
abline(h = y0_obs*0.05)

#### Setting up barrel volumes to implement sampling error ####
samp_vol <- 5; barrel_vol <- 150 # Under these conditions, you can sample each barrel a maximum of 40 times
vol <- c(barrel_vol)
for (i in 2:length(t1)){
  vol[i] <- vol[i - 1] - samp_vol
}

#### Generating data ####
###### PUT SOMETHING HERE ####

# Combining generated data by barrel
b1 <- cbind(t1, y1_obs_err); b2 <- cbind(t2, y2_obs_err); b3 <- cbind(t3, y3_obs_err); b4 <- cbind(t4, y4_obs_err); b5 <- cbind(t5, y5_obs_err); b6 <- cbind(t6, y6_obs_err); b7 <- cbind(t7, y7_obs_err); b8 <- cbind(t8, y8_obs_err); b9 <- cbind(t9, y9_obs_err)
full_ts <- rbind(b1, b2[2:nrow(b2), ], b3[2:nrow(b3), ], b4[2:nrow(b4), ], b5[2:nrow(b5), ], b6[2:nrow(b6), ], b7[2:nrow(b7), ], b8[2:nrow(b8), ], b9[2:nrow(b9), ])
full_ts <- as.data.frame(full_ts); colnames(full_ts) <- c("t", "y_obs_err")
# Ordering generated data by sampling time
full_ts_ordered <- full_ts[with(full_ts, order(t)), ]
# Plotting generated data
plot(full_ts_ordered$t, full_ts_ordered$y_obs_err)

#### Fitting generated data to exponential model ####
y_obs_err <- full_ts_ordered$y_obs_err
fit <- nls(y_obs_err ~ y0_est*exp(-all_times/lambda_est), 
           start = list(y0_est = 1e6, lambda_est = 0.15))
coef(fit); y0_est <- coef(fit)[1]; lambda_est <- coef(fit)[2]
# Plotting residuals to see variance in error
plot(all_times, resid(fit), pch = 16); abline(h = 0, lty = 2)

#### Estimated abundances ####
y0_est <- as.numeric(y0_est); lambda_est <- as.numeric(lambda_est)
y_est <- y0_est*exp(-all_times/lambda_est)
plot(y_est)

#### Plotting (1) observations without error, (2) generated data, and (3) estimated data for comparison ####
plot(all_times[1:20], y_obs[1:20], type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(all_times[1:20], y_obs_err[1:20], col = "cornflowerblue", pch = 16)
points(all_times[1:20], y_est[1:20], col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex = 0.8)


#### Marty edits ####

# (1) Separate the time structure of the process model from the time structure of the sampling model
# (2) The process model should have a standardized and constant time step, maybe equal to 1 hr if the degradation rate is in units per hr. 
  # So the process model should go through full time series, and sampling model should only be applied on select time steps?
# (3) The random variable for environmental noise should occur at each time step of the process model rather than the sampling model 

## Barrel 1
# Setting up initial abundance, and empty vector to fill with observations
t_elap <- 0.5 # Units are in hours so 10 day * 24hrs, for example; CHANGE ALL_TIMES
y_prev <- 1e6; y1_obs_err <- c(); y1_obs_err[1] <- y_prev
# Loop to generate data with (1) process error, (2) environmental stochasticity, and (3) sampling error
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[1]*exp(-rnorm(1, 0, 0.1)))) # Determine eDNAs that "die" via exponential distribution (process error) and environmental stochasticity
  if (all_times[i] %in% t1){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[i - 1]) # Sample eDNAs from remaining with sampling error via sampled volume vs total remaining volume
    y_next <- y_prev - y_d - y_select # eDNAs in barrel at next time step is equal to those available at previous time step, minus those that died, and those that were selected
    y_prev <- y_next # Re-assigning y_prev for next iteration
    y1_obs_err[i] <- y_select # Add sampled eDNAs to observation vector
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y1_obs_err[i] <- y_prev
  }
}; plot(y1_obs_err)
First0 <- which(y1_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y1_obs_err[1:First0])

sampled <- which(all_times %in% t1)
y1_obs_err <- y1_obs_err[sampled]

# Initially, sample every barrel, every hour, until they're empty





