#### [INSERT A DOC STRING HERE] ####

library(tidyverse)
library(bayesplot)

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
y_samp_init <- c(NA, 0.0, y_samp_init)

#### Generating data ####
## Barrel 1
# Set seed because reproducibility
set.seed(123)
# Setting up initial abundance, and empty vector to fill with observations
t_elap <- 0.5 # Units are in days, so 0.5 is half a day
y_prev <- 1e6; y1_obs_err <- c(); y1_obs_err[1] <- y_prev; y1_samp <- c(); j <- 2
# Loop to generate data with (1) process error, (2) environmental stochasticity, and (3) sampling error
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[1]*exp(-rnorm(1, 0, 0.1)))) # Determine eDNAs that "die" via exponential distribution (process error) and environmental stochasticity
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
plot(x = all_times, y = y1_obs_err)
y1_samp <- cbind(t1[2:length(t1)], y1_samp)
Barrel <- c(rep(1, nrow(y1_samp))); y1_samp <- cbind(Barrel, y1_samp)

## Barrel 2
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y2_obs_err <- c(); y2_obs_err[1] <- y_prev; y2_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[2]*exp(-rnorm(1, 0, 0.1)))) 
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
plot(x = all_times, y = y2_obs_err)
y2_samp <- cbind(t2[2:length(t2)], y2_samp)
Barrel <- c(rep(2, nrow(y2_samp))); y2_samp <- cbind(Barrel, y2_samp)

## Barrel 3
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y3_obs_err <- c(); y3_obs_err[1] <- y_prev; y3_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[3]*exp(-rnorm(1, 0, 0.1)))) 
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
plot(x = all_times, y = y3_obs_err)
y3_samp <- cbind(t3[2:length(t3)], y3_samp)
Barrel <- c(rep(3, nrow(y3_samp))); y3_samp <- cbind(Barrel, y3_samp)

## Barrel 4
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y4_obs_err <- c(); y4_obs_err[1] <- y_prev; y4_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[4]*exp(-rnorm(1, 0, 0.1))))
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
plot(x = all_times, y = y4_obs_err)
y4_samp <- cbind(t4[2:length(t4)], y4_samp)
Barrel <- c(rep(4, nrow(y4_samp))); y4_samp <- cbind(Barrel, y4_samp)

## Barrel 5
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y5_obs_err <- c(); y5_obs_err[1] <- y_prev; y5_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[5]*exp(-rnorm(1, 0, 0.1))))
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
plot(x = all_times, y = y5_obs_err)
y5_samp <- cbind(t5[2:length(t5)], y5_samp)
Barrel <- c(rep(5, nrow(y5_samp))); y5_samp <- cbind(Barrel, y5_samp)

## Barrel 6
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y6_obs_err <- c(); y6_obs_err[1] <- y_prev; y6_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[6]*exp(-rnorm(1, 0, 0.1))))
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
plot(x = all_times, y = y6_obs_err)
y6_samp <- cbind(t6[2:length(t6)], y6_samp)
Barrel <- c(rep(6, nrow(y6_samp))); y6_samp <- cbind(Barrel, y6_samp)

## Barrel 7
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y7_obs_err <- c(); y7_obs_err[1] <- y_prev; y7_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[7]*exp(-rnorm(1, 0, 0.1))))
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
plot(x = all_times, y = y7_obs_err)
y7_samp <- cbind(t7[2:length(t7)], y7_samp)
Barrel <- c(rep(7, nrow(y7_samp))); y7_samp <- cbind(Barrel, y7_samp)

## Barrel 8
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y8_obs_err <- c(); y8_obs_err[1] <- y_prev; y8_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[8]*exp(-rnorm(1, 0, 0.1))))
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
plot(x = all_times, y = y8_obs_err)
y8_samp <- cbind(t8[2:length(t8)], y8_samp)
Barrel <- c(rep(8, nrow(y8_samp))); y8_samp <- cbind(Barrel, y8_samp)

## Barrel 9
set.seed(123)
t_elap <- 0.5 
y_prev <- 1e6; y9_obs_err <- c(); y9_obs_err[1] <- y_prev; y9_samp <- c(); j <- 2
for (i in 2:length(all_times)){
  y_d <- rbinom(1, y_prev, 1 - exp(-t_elap*r_vec[9]*exp(-rnorm(1, 0, 0.1))))
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
plot(x = all_times, y = y9_obs_err)
y9_samp <- cbind(t9[2:length(t9)], y9_samp)
Barrel <- c(rep(9, nrow(y9_samp))); y9_samp <- cbind(Barrel, y9_samp)

full_ts <- rbind(y_samp_init, y1_samp, y2_samp, y3_samp, y4_samp, y5_samp, y6_samp, y7_samp, y8_samp, y9_samp)
full_ts <- as.data.frame(full_ts); colnames(full_ts) <- c("barrel", "t", "y_samp")
# Ordering generated data by sampling time
full_ts_ordered <- full_ts[with(full_ts, order(t)), ]
# Plotting generated data
colours <- c("black", "firebrick", "darkorange", "goldenrod2", "forestgreen", "cornflowerblue", "purple3", "orchid", "salmon")
colidx <- c("white"); barrel <- as.numeric(full_ts_ordered$barrel)
for (i in 2:length(full_ts_ordered$t)){
  colidx[i] <- colours[barrel[i]]
}
plot(full_ts_ordered$t, full_ts_ordered$y_samp, col = colidx,
     xlab = "Day", ylab = "eDNA Abundance (Sampled)")

plot(x = all_times, y = y1_obs_err, col = "black", type = "l", xlab = "Day", ylab = "eDNA Abundance (Unsampled)")
lines(x = all_times, y = y2_obs_err, col = "firebrick")
lines(x = all_times, y = y3_obs_err, col = "darkorange")
lines(x = all_times, y = y4_obs_err, col = "goldenrod2")
lines(x = all_times, y = y5_obs_err, col = "forestgreen")
lines(x = all_times, y = y6_obs_err, col = "cornflowerblue")
lines(x = all_times, y = y7_obs_err, col = "purple3")
lines(x = all_times, y = y8_obs_err, col = "orchid")
lines(x = all_times, y = y9_obs_err, col = "salmon")

#### Fitting generated data to exponential model ####
#### ... From Steph's script: ####
#### ... ... Simulate data: ####
a<-1.2			# growth rate
b<-1.2/1200	# density dependence
sig<-0.8		# process error
N<-100			# timesteps

y<-numeric(N+1)	# vector of population sizes
y1<-600			# initial condition

# Simulate data:
set.seed(9823)
y[1]<-y1
for(i in 2:(N+1)) y[i]<-y[i-1]*exp(rnorm(n=1, a-b*y[i-1], sig))

S<-y[1:N];		# spawners
Y<-log(y[2:(N+1)]/S)	# productivity = log(recruits/spawner)

plot(y, type="l", xlab="time", ylab="spawners"); abline(h=a/b, lty=3)
plot(S, Y, xlab="spawners", ylab="productivity"); abline(a=a, b=-b, col=2); abline(h=0)

#### ... ... ... Fit with JAGS: Single population ####
# Need to write model function IN JAGS LANGUAGE (see user manual)
ricker.model<-function(){	
  # define priors on parameter (mean and PRECISION)
  a ~ dnorm(0, 0.0001)
  b ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100)
  tau <- pow(sigma, -2)	
  
  # simulate model and likelihood
  for(i in 1:N){ 					# for each data point
    Y.hat[i] <- a-b*S[i] 		# expected value based on parameters
    Y[i] ~ dnorm(Y.hat[i], tau) # probability of observation Y given expectation Y.hat
  }
}

fit3<-jags.fit(data=list('Y'=Y, 'S'=S, 'N'=N), params=c("a", "b", "sigma"), model=ricker.model, n.chains=3, n.adapt=100)

plot(fit3)
summary(fit3) 

#### ... ... Simulate hierarchical data: ####
a<-1.2								# growth rate
b<-1.2/1200							# density dependence
sig<-0.3							# process error
sigP<-0.8							# sd in growth rate among populations

Nt<-50								# number of timesteps
Np<-20								# number of populations
thetaP<-rnorm(Np, mean=0, sd=sigP)	# population random-effect on growth rate

t<-rep(1:(Nt+1), Np)				# vector of timesteps
P<-rep(1:Np, each=(Nt+1))			# vector of population numbers
y<-numeric((Nt+1)*Np)				# vector of population sizes (to be filled)
y[which(t==1)]<-runif(Np, 100, 1000)# initialize with random number between 100 and 1000

# Simulate data:
for(i in 1:Np){ #for each population
  for(j in 2:(Nt+1)){ #for each timestep
    y[which(P==i&t==j)]<-y[which(P==i&t==(j-1))]*exp(rnorm(n=1, (a+thetaP[i])-b*y[which(P==i&t==(j-1))], sig))
  }
}

# Make dataframe of spawners and productivity
data<-data.frame(S=numeric(Nt*Np), Y=numeric(Nt*Np), Pop=numeric(Nt*Np))
N<-Nt*Np

k<-0
for(i in 1:Np){
  for(j in 1:Nt){
    k<-k+1
    data$S[k]<- y[which(P==i&t==j)]
    data$Y[k]<- log(y[which(P==i&t==(j+1))]/y[which(P==i&t==j)])
    data$Pop[k]<-i
  }
}

xyplot(Y~S|Pop, data=data)

#### ... ... ... Fit with JAGS: Hierarchical ####
# Transform data to matrices to make model syntax easier
Np <- 20 # Because 20 populations
Nt <- 50 # Because number of time points
S<-matrix(data$S, nrow=Np, ncol=Nt, byrow=TRUE) # Where S is number of spawners
Y<-matrix(data$Y, nrow=Np, ncol=Nt, byrow=TRUE) # Where Y is productivity (= log(recruits/spawner))

H.ricker.model<-function(){
  # Priors on parameters
  a ~ dnorm(0, 0.0001)
  b ~ dnorm(0, 0.0001)
  sig ~ dunif(0, 100)
  sigP ~ dunif(0, 100)
  tauP <- pow(sigP, -2)
  tau <- pow(sig, -2)
  
  # Hierarchical part - a different thetaP for each population:
  for(i in 1:Np){
    thetaP[i] ~ dnorm(0, tauP)
  }
  
  # Model simulation and likelihood
  for(j in 1:Np){ #for each population
    for(i in 1:Nt){ #for each time step
      Y.hat[j,i] <- a + thetaP[j] - b*S[j,i] # calculate model prediction
      Y[j,i] ~ dnorm(Y.hat[j,i], tau) #likelihood
    } #end i
  } # end j	
}


hfit2<-jags.fit(data=list('S'=S, 'Y'=Y, 'Np'=Np, 'Nt'=Nt), params=c("a", "b", "sig", "sigP"), model=H.ricker.model)

S<-summary(hfit2)
par(mfrow=c(1,1))
for(i in 1:7){
  hist(hfit2[[1]][,i], main=colnames(hfit2[[1]])[i], col="#00000030", border=NA, xlab="", freq=FALSE)
  hist(hfit2[[2]][,i], main=colnames(hfit2[[2]])[i], col="#FF000030", add=TRUE, freq=FALSE, border=NA)
  hist(hfit2[[3]][,i], main=colnames(hfit2[[3]])[i], col="#00FF0030", add=TRUE, freq=FALSE, border=NA)
  abline(v=S[[1]][i,1], lty=2)
  abline(v=c(a,b,sig,sigP)[i])
  abline(v=hfit1.par[i], lty=3)
}
legend("topright", lty=c(1,2,3), c("Real", "JAGS", "lmer"))


#### For Ning's stuff ####
#### ... Data manip ####
## Steph made matrices wherein rows were populations, and columns were time points. We have to switch that up a bit potentially because populations aren't sampled at same time points?
y <- matrix(nrow = num_barrels, ncol = length(all_times)); colnames(y) <- all_times
t1 <- t1[2:length(t1)]; t2 <- t2[2:length(t2)]; t3 <- t3[2:length(t3)]; t4 <- t4[2:length(t4)]; t5 <- t5[2:length(t5)]; t6 <- t6[2:length(t6)]; t7 <- t7[2:length(t7)]; t8 <- t8[2:length(t8)]; t9 <- t9[2:length(t9)]
b1 <- 1; b2 <- 1; b3 <- 1; b4 <- 1; b5 <- 1; b6 <- 1; b7 <- 1; b8 <- 1; b9 <- 1
for (i in 1:length(all_times)){
  if (all_times[i] %in% t1){
    y[1, i] <- y1_samp[b1, 3]
    b1 <- b1 + 1
  } else {
    if (all_times[i] %in% t2){
      y[2, i] <- y2_samp[b2, 3]
      b2 <- b2 + 1
    } else {
      if (all_times[i] %in% t3){
        y[3, i] <- y3_samp[b3, 3]
        b3 <- b3 + 1
      } else {
        if (all_times[i] %in% t4){
          y[4, i] <- y4_samp[b4, 3]
          b4 <- b4 + 1
        } else {
          if (all_times[i] %in% t5){
            y[5, i] <- y5_samp[b5, 3]
            b5 <- b5 + 1
          } else {
            if (all_times[i] %in% t6){
              y[6, i] <- y6_samp[b6, 3]
              b6 <- b6 + 1
            } else {
              if (all_times[i] %in% t7){
                y[7, i] <- y7_samp[b7, 3]
                b7 <- b7 + 1
              } else {
                if (all_times[i] %in% t8){
                  y[8, i] <- y8_samp[b8, 3]
                  b8 <- b8 + 1
                } else {
                  if (all_times[i] %in% t9){
                    y[9, i] <- y9_samp[b9, 3]
                    b9 <- b9 + 1
                  } else { print(":(") }
                }
              }
            }
          }
        }
      }
    }
  }
}
y[ ,1] <- rep(full_ts_ordered$y_samp[1], nrow(y)) # Might have to change this such that there's a different initial sample from each bucket but for now, keeping as is because shouldn't make a huge difference

#### ... JAGS ####
exp_decay_model <- function(){
  
  # Priors on parameters
  alpha ~ dunif(0, 1)
  alpha_tau ~ dunif(0, 100)
  alpha_tau_pop ~ dunif(0, 100)
  
  # Hierarchical part - a different theta_pop for each population:
  for(i in 1:num_barrels){
    theta_pop[i] ~ dnorm(0, alpha_tau_pop)
  }
  
  # Model simulation
  for (i in 1:num_barrels){ # For each population
    for (j in 1:ncol(y)){ # For each time step
      if (!is.na(y[i, j])){
        yhat[i, j] <- yhat[i, j]*exp(-t_elap*alpha_pop*theta_pop[i]) # Calculate model prediction
      }
    }
  }
  
  # Likelihood
  for (i in 1:num_barrels){ 
    for (j in 1:length(all_times)){ 
      if(!is.na(y[i, j])){
        y[i, j] ~ dpois(yhat[i, j]) # Likelihood
      }
    } 
  }
}






























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















