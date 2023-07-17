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
samp_vol <- 5; barrel_vol <- 100 # Under these conditions, you can sample each barrel a maximum of 40 times
vol <- c(barrel_vol)
for (i in 2:length(t1 + 1)){
  vol[i] <- vol[i - 1] - samp_vol
}; vol <- vol[2:length(vol)]

#### Sampling initial abundance ####
set.seed(123); y_init <- 1e6
y_samp_init <- y_select <- rbinom(1, y_init, samp_vol/(vol[1] + samp_vol))
y_samp_init <- c(0.0, y_samp_init)

#### Generating data ####
## Barrel 1
# Set seed because reproducibility
set.seed(123)
# Setting up initial abundance, and empty vector to fill with observations
t_elap <- 0.5 # Units are in hours so 10 day * 24hrs, for example; CHANGE ALL_TIMES
y_prev <- 1e6; y1_obs_err <- c(); y1_obs_err[1] <- y_prev; y1_samp <- c(); j <- 2
# Loop to generate data with (1) process error, (2) environmental stochasticity, and (3) sampling error
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[1]*exp(-rnorm(1, 0, 0.1)))) # Determine eDNAs that "die" via exponential distribution (process error) and environmental stochasticity
  if (all_times[i] %in% t1){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) # Sample eDNAs from remaining with sampling error via sampled volume vs total remaining volume
    y1_samp[j - 1] <- y_select; j <- j + 1 # Add sampled eDNAs to observation vector
    y_next <- y_prev - y_d - y_select # eDNAs in barrel at next time step is equal to those available at previous time step, minus those that died, and those that were selected
    y_prev <- y_next # Re-assigning y_prev for next iteration
    y1_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y1_obs_err[i] <- y_prev
  }
}; plot(all_times, y1_obs_err)
First0 <- which(y1_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y1_obs_err[1:First0])
y1_samp <- cbind(t1[2:length(t1)], y1_samp)

## Barrel 2
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y2_obs_err <- c(); y2_obs_err[1] <- y_prev; y2_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[2]*exp(-rnorm(1, 0, 0.1)))) 
  if (all_times[i] %in% t2){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) 
    y2_samp[j - 1] <- y_select; j <- j + 1 
    y_next <- y_prev - y_d - y_select 
    y_prev <- y_next 
    y2_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y2_obs_err[i] <- y_prev
  }
}; plot(all_times, y2_obs_err)
First0 <- which(y2_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y2_obs_err[1:First0])
y2_samp <- cbind(t2[2:length(t2)], y2_samp)

## Barrel 3
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y3_obs_err <- c(); y3_obs_err[1] <- y_prev; y3_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[3]*exp(-rnorm(1, 0, 0.1)))) 
  if (all_times[i] %in% t3){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) 
    y3_samp[j - 1] <- y_select; j <- j + 1 
    y_next <- y_prev - y_d - y_select 
    y_prev <- y_next 
    y3_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y3_obs_err[i] <- y_prev
  }
}; plot(all_times, y3_obs_err)
First0 <- which(y3_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y3_obs_err[1:First0])
y3_samp <- cbind(t3[2:length(t3)], y3_samp)

## Barrel 4
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y4_obs_err <- c(); y4_obs_err[1] <- y_prev; y4_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[4]*exp(-rnorm(1, 0, 0.1)))) 
  if (all_times[i] %in% t4){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) 
    y4_samp[j - 1] <- y_select; j <- j + 1 
    y_next <- y_prev - y_d - y_select 
    y_prev <- y_next 
    y4_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y4_obs_err[i] <- y_prev
  }
}; plot(all_times, y4_obs_err)
First0 <- which(y4_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y4_obs_err[1:First0])
y4_samp <- cbind(t4[2:length(t4)], y4_samp)

## Barrel 5
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y5_obs_err <- c(); y5_obs_err[1] <- y_prev; y5_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[5]*exp(-rnorm(1, 0, 0.1)))) 
  if (all_times[i] %in% t5){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) 
    y5_samp[j - 1] <- y_select; j <- j + 1 
    y_next <- y_prev - y_d - y_select 
    y_prev <- y_next 
    y5_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y5_obs_err[i] <- y_prev
  }
}; plot(all_times, y5_obs_err)
First0 <- which(y5_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y5_obs_err[1:First0])
y5_samp <- cbind(t5[2:length(t5)], y5_samp)

## Barrel 6
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y6_obs_err <- c(); y6_obs_err[1] <- y_prev; y6_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[6]*exp(-rnorm(1, 0, 0.1)))) 
  if (all_times[i] %in% t6){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) 
    y6_samp[j - 1] <- y_select; j <- j + 1 
    y_next <- y_prev - y_d - y_select 
    y_prev <- y_next 
    y6_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y6_obs_err[i] <- y_prev
  }
}; plot(all_times, y6_obs_err)
First0 <- which(y6_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y6_obs_err[1:First0])
y6_samp <- cbind(t6[2:length(t6)], y6_samp)

## Barrel 7
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y7_obs_err <- c(); y7_obs_err[1] <- y_prev; y7_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[7]*exp(-rnorm(1, 0, 0.1)))) 
  if (all_times[i] %in% t7){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) 
    y7_samp[j - 1] <- y_select; j <- j + 1 
    y_next <- y_prev - y_d - y_select 
    y_prev <- y_next 
    y7_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y7_obs_err[i] <- y_prev
  }
}; plot(all_times, y7_obs_err)
First0 <- which(y7_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y7_obs_err[1:First0])
y7_samp <- cbind(t7[2:length(t7)], y7_samp)

## Barrel 8
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y8_obs_err <- c(); y8_obs_err[1] <- y_prev; y8_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[8]*exp(-rnorm(1, 0, 0.1)))) 
  if (all_times[i] %in% t8){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) 
    y8_samp[j - 1] <- y_select; j <- j + 1 
    y_next <- y_prev - y_d - y_select 
    y_prev <- y_next 
    y8_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y8_obs_err[i] <- y_prev
  }
}; plot(all_times, y8_obs_err)
First0 <- which(y8_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y8_obs_err[1:First0])
y8_samp <- cbind(t8[2:length(t8)], y8_samp)

## Barrel 9
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y9_obs_err <- c(); y9_obs_err[1] <- y_prev; y9_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, t_elap*exp(-1*r_vec[9]*exp(-rnorm(1, 0, 0.1)))) 
  if (all_times[i] %in% t9){
    y_select <- rbinom(1, y_prev - y_d, samp_vol/vol[j - 1]) 
    y9_samp[j - 1] <- y_select; j <- j + 1 
    y_next <- y_prev - y_d - y_select 
    y_prev <- y_next 
    y9_obs_err[i] <- y_prev
  } else {
    y_next <- y_prev - y_d
    y_prev <- y_next
    y9_obs_err[i] <- y_prev
  }
}; plot(all_times, y9_obs_err)
First0 <- which(y9_obs_err == 0)[1]
plot(x = all_times[1:First0], y = y9_obs_err[1:First0])
y9_samp <- cbind(t9[2:length(t9)], y9_samp)

full_ts <- rbind(y_samp_init, y1_samp, y2_samp, y3_samp, y4_samp, y5_samp, y6_samp, y7_samp, y8_samp, y9_samp)
full_ts <- as.data.frame(full_ts); colnames(full_ts) <- c("t", "y_samp")
# Ordering generated data by sampling time
full_ts_ordered <- full_ts[with(full_ts, order(t)), ]
# Plotting generated data
plot(full_ts_ordered$t, full_ts_ordered$y_samp)

#### Fitting generated data to exponential model ####
#### ... TO DO! ####

#### Estimated abundances ####
y0_est <- as.numeric(y0_est); r_est <- as.numeric(r_est)
y_est <- y0_est*exp(-all_times*r_est)
plot(y_est)

#### Plotting (1) observations without error, (2) generated data, and (3) estimated data for comparison ####
plot(all_times[1:20], y_obs[1:20], type = "l", xlab = "Day", ylab = "eDNA Concentration")
points(all_times[1:20], y_obs_err[1:20], col = "cornflowerblue", pch = 16)
points(all_times[1:20], y_est[1:20], col = "firebrick2", pch = 18)
legend("topright", legend=c("Deterministic data", "Observed data", "Estimated data"),
       col=c("black", "cornflowerblue", "firebrick2"), lty = c(1, NA, NA), pch = c(NA, 16, 18), cex = 0.8)














# Initially, sample every barrel, every hour, until they're empty





