###Laplace model dispersal simulations###
###September 2022###
###Jaime Grimm###

####Housekeeping####
setwd("~/Documents/PhD/Code")
require(tidyverse)
require(bbmle)
require(ggpubr)


####Laplace function####
#First write a function that defines the laplace model
#This is the deterministic component of the model

laplace <- function(bk,a1,a2,x,xk){ #where bk is a scaling factor, a1 and a2 are exponential coefficients, x is a sample and 
  #xk is the location of the farm
  if (x<xk){
    y <- (bk*exp(-a1*abs(xk-x))) #this is the upstream function
  }
  else if (x>xk){
    y <- (bk*exp(-a2*abs(x-xk))) #and this is the downstream function
  }
  return (y)
}

####NLL function####
#Here, we are writing a function that calculates the negative log likehihood for a set of parameters
#Note, we are keeping the likelihood function as a poission distribution
#First, we define our parameters into a numeric vector because dpois works on vectors.

params = list()

#The function we want to minimize is:
laplaceNLL = function(params, y, x) { #inside the parentheses here we pass info from the outside environment into the function
  bk=params[[1]]  #First we define our parameters within the function
  a1=params[[2]]
  a2=params[[3]]
  xk=params[[4]]
  c=params[[5]] #c is the wild-borne pathogen value
  
  y.hat <- c() #If we're calculating something within the function, we should define it in the function
  
  for (i in 1:length(x)) {
    if (x[i]<xk){
      y.hat[i] <- (bk*exp(-a1*abs(xk-x[i]))+c)
    }
    else if (x[i]>xk){
      y.hat[i] <- (bk*exp(-a2*abs(x[i]-xk))+c)
    }
  }
  
  #NLL<- -sum(dpois(x = y, lambda = y.hat, log=TRUE)) #poisson
  NLL <- -sum(dnbinom(x=y, size=5, prob =, mu = y.hat , log=TRUE))
  return(NLL)
}

parnames(laplaceNLL) = c("bk", "a1", "a2", "xk", "c")

####Simulations####
#Here, I'm writing a function to simulate datasets with some base parameters
#Then, I will call this function to simulate scenarios where I specify changes
simulate <- function(x=c(-5000,-3000, -2200, -1200, -700, -400, -250,-150, -100,-75, -50,-40, -35, 
                         -20, -10, -5, 5, 10, 20, 35, 40, 50, 75, 100, 150, 250, 400, 700, 1200, 2200, 3000, 5000), 
                     bk=100, a1=0.0019, a2=0.0009, xk=0, c=5){
  y <- numeric(length(x))
  for (i in 1:length(x))
    y[i] <- rnbinom(n=1, size=((laplace(bk,a1,a2,x[i],xk)+c)*5), mu=laplace(bk,a1,a2,x[i],xk)+c) #here, size is k - the dispersion parameter
  return(y)
} 

num.sims <- 100
mle.params <- matrix(nrow=num.sims, ncol=5) #create an empty matrix that will store the mle parameters for each simulation

x.sim <- c(-5000,-3000, -2200, -1200, -700, -400, -250,-150, -100,-75, -50,-40, -35, 
           -20, -10, -5, 5, 10, 20, 35, 40, 50, 75, 100, 150, 250, 400, 700, 1200, 2200, 3000, 5000)
params.start <- c(bk=100, a1=0.0019, a2=0.0009, c=5)

for (i in 1:num.sims){
  y.sim <- simulate(x.sim)
  mle.params[i,] <- coef(mle2(minuslogl = laplaceNLL, start=params.start, fixed= list(xk=0), data= list(x=x.sim, y=y.sim)))
} 
#

x.true <- seq(-5000, 5000, length.out=1000) #high resolution to get a nice smooth curve
y.true <- numeric(length(x.true))
for (l in 1:length(x.true)) {  
  y.true[l] <- laplace(100,0.0019,0.0009, x.true[l], 0) + 5 #hand the laplace function the parameter estimates used to simulate scenario
}

y.mle <- numeric(length(x.true))

#par(mfrow=c(3,4))
foo.plot <- plot(x.true, y.true, type="l", xlim = c(-5000, 5000), ylim = c(0, 200), xlab="Location relative to origin", 
                 ylab="Pathogen Density", lwd=2.5,cex.lab=1.2, cex.axis=1, cex.sub=4, col="black")
#legend(1000, 145, legend=c("Farm contribution", "Wild host contribution", "Sum of components"),
#       col=c("black", "black",rgb(70/255, 115/255, 120/255, alpha = 1)), lty=1:2, lwd=2:2:2, cex=0.8,
#       family="serif")


#title("bk = 200", line=-1, cex.main=1)

for (i in 1:num.sims){
  p <- mle.params[i,]
  for (j in 1:length(x.true)) {
    y.mle[j] <- (laplace(p[1],p[2],p[3],x.true[j],p[4]) + p[5])
  }
  lines(x.true, y.mle, col = rgb(70/255, 115/255, 120/255, alpha = 0.2))
}
#abline(v=)
#legend(100, 100, legend=c("Simulations (N=100)", "Data generating model"),
#       col=c(rgb(70/255, 115/255, 120/255, alpha = 1), "black"), lty=1:1, lwd=1:2, cex=0.8)

#Plot distribution of simulated parameter estimates against true estimate
#bk
bk.plot <- hist(mle.params[,1], xlab="bk", main = "", col="white")
abline(v=100, col=rgb(70/255, 115/255, 120/255, alpha = 1), lwd=3)
legend(105, 60, legend=("True bk value"), col=(rgb(70/255, 115/255, 120/255, alpha = 1)), lty=1, lwd=2) 

#a1
a1.plot <- hist(mle.params[,2], xlab="a1", main = "", col="white")
abline(v=0.0019, col=rgb(70/255, 115/255, 120/255, alpha = 1), lwd=3)
legend(0.0030, 20, legend=("True a1 value"), col=(rgb(70/255, 115/255, 120/255, alpha = 1)), lty=1, lwd=2) 

#a2
a2.plot <- hist(mle.params[,3], xlab="a2", main = "", col="white")
abline(v=0.0009, col=rgb(70/255, 115/255, 120/255, alpha = 1), lwd=3)
legend(-5, 25, legend=("True parameter value"), col=(rgb(70/255, 115/255, 120/255, alpha = 1)), lty=1, lwd=2) 

#xk
xk.plot <- hist(mle.params[,4], xlab="xk", main = "", col="white")
abline(v=0, col=rgb(70/255, 115/255, 120/255, alpha = 1), lwd=3)
legend(5, 50, legend=("True xk value"), col=(rgb(70/255, 115/255, 120/255, alpha = 1)), lty=1, lwd=2) 




#Let's start by making a series of scenarios with varying sample sizes
x1<- c(-2200, -20, 35, 4000)
x2 <- c(-2200, -700, -50,-20, -5, 35, 100, 400, 1200, 4000)
x3 <- c(-4000, -2200, -1200, -700, -400, -250, -100, -50, -35, -20, -10, -5, 5, 10, 20, 35, 50, 100, 250, 400, 700,
   1200, 2200, 4000)
x4 <- c(-4000,-3000, -2200, -1200, -700, -400, -250,-150, -100,-75, -50,-40, -35, -20, -10, -5, 5, 10, 20, 35, 40, 50, 75, 100, 150, 
        250, 400, 700, 1200, 2200, 3000, 4000)
x5 <- c(-2200, -700, -50,-20, -5, 35, 100, 400, 1200, 4000)

x.s1 <- simulate(x=x1) #n=4
x.s2 <- simulate(x=x2) #n=10
x.s3 <- simulate(x=x3) #n=24
x.s4 <- simulate(x=x4) #n=34
x.s5 <- simulate(x=x5) #n=10, clustered

####Maximum likelihood estimation####

parnames(laplaceNLL) = c("bk", "a1", "a2", "xk")
start <- c(bk=105, a1=0.019, a2=0.019, xk=3)
#Now, we put each of these simulations into the mle to get parameter estimates for each scenario
x.s1.est <- mle2(minuslogl = laplaceNLL, start=start, data= list(x= x1, y=x.s1)) 
x.s2.est <- mle2(minuslogl = laplaceNLL, start=start, data= list(x= x2, y=x.s2)) 
x.s3.est <- mle2(minuslogl = laplaceNLL, start=start, data= list(x= x3, y=x.s3)) 
x.s4.est <- mle2(minuslogl = laplaceNLL, start=start, data= list(x= x4, y=x.s4)) 
x.s5.est <- mle2(minuslogl = laplaceNLL, start=start, data= list(x= x5, y=x.s5)) 


####Produce model lines for plots####
##S1##
x.s1.true.x <- seq(-4000, 4000, length.out=1000) #high resolution to get a nice smooth curve
x.s1.true.y <- numeric(length(x.s1.true.x))

for (i in 1:100){
  for (l in 1:length(x.s1.true.x)) {  
    x.s1.true.y[l] <- laplace(100,0.0019,0.0009, x.s1.true.x[l], 0) #hand the laplace function the parameter estimates used to simulate scenario
  }
}

#Next, we get y estimates for the fitted/estimated model
x.s1.em.x <- seq(-4000, 4000, length.out=1000)
x.s1.em.y <- numeric(length(x.s1.em.x))

for (l in 1:length(x.s1.em.x)) {  
  x.s1.em.y[l] <- laplace(coef(x.s1.est)[1],coef(x.s1.est)[2],coef(x.s1.est)[3],x.s1.em.x[l],coef(x.s1.est)[4]) #hand the laplace function the parameter estimates used to simulate scenario
}

##S2##
x.s2.true.x <- seq(-4000, 4000, length.out=1000) #high resolution to get a nice smooth curve
x.s2.true.y <- numeric(length(x.s2.true.x))

for (l in 1:length(x.s2.true.x)) {  
  x.s2.true.y[l] <- laplace(100,0.0019,0.0009, x.s2.true.x[l], 0) #hand the laplace function the parameter estimates used to simulate scenario
}

#Next, we get y estimates for the fitted/estimated model
x.s2.em.x <- seq(-4000, 4000, length.out=1000)
x.s2.em.y <- numeric(length(x.s2.em.x))

for (l in 1:length(x.s2.em.x)) {  
  x.s2.em.y[l] <- laplace(coef(x.s2.est)[1],coef(x.s2.est)[2],coef(x.s2.est)[3],x.s2.em.x[l],coef(x.s2.est)[4]) #hand the laplace function the parameter estimates used to simulate scenario
}

##S3##
x.s3.true.x <- seq(-4000, 4000, length.out=1000) #high resolution to get a nice smooth curve
x.s3.true.y <- numeric(length(x.s3.true.x))

for (l in 1:length(x.s3.true.x)) {  
  x.s3.true.y[l] <- laplace(100,0.0019,0.0009, x.s3.true.x[l], 0) #hand the laplace function the parameter estimates used to simulate scenario
}

#Next, we get y estimates for the fitted/estimated model
x.s3.em.x <- seq(-4000, 4000, length.out=1000)
x.s3.em.y <- numeric(length(x.s3.em.x))

for (l in 1:length(x.s3.em.x)) {  
  x.s3.em.y[l] <- laplace(coef(x.s3.est)[1],coef(x.s3.est)[2],coef(x.s3.est)[3],x.s3.em.x[l],coef(x.s3.est)[4]) #hand the laplace function the parameter estimates used to simulate scenario
}

##S4##
x.s4.true.x <- seq(-4000, 4000, length.out=1000) #high resolution to get a nice smooth curve
x.s4.true.y <- numeric(length(x.s4.true.x))

for (l in 1:length(x.s4.true.x)) {  
  x.s4.true.y[l] <- laplace(100,0.0019,0.0009, x.s4.true.x[l], 0) #hand the laplace function the parameter estimates used to simulate scenario
}

#Next, we get y estimates for the fitted/estimated model
x.s4.em.x <- seq(-4000, 4000, length.out=1000)
x.s4.em.y <- numeric(length(x.s4.em.x))

for (l in 1:length(x.s4.em.x)) {  
  x.s4.em.y[l] <- laplace(coef(x.s4.est)[1],coef(x.s4.est)[2],coef(x.s4.est)[3],x.s4.em.x[l],coef(x.s4.est)[4]) #hand the laplace function the parameter estimates used to simulate scenario
}

##S5##
x.s5.true.x <- seq(-4000, 4000, length.out=1000) #high resolution to get a nice smooth curve
x.s5.true.y <- numeric(length(x.s5.true.x))

for (l in 1:length(x.s5.true.x)) {  
  x.s5.true.y[l] <- laplace(100,0.0019,0.0009, x.s5.true.x[l], 0) #hand the laplace function the parameter estimates used to simulate scenario
}

#Next, we get y estimates for the fitted/estimated model
x.s5.em.x <- seq(-4000, 4000, length.out=1000)
x.s5.em.y <- numeric(length(x.s5.em.x))

for (l in 1:length(x.s5.em.x)) {  
  x.s5.em.y[l] <- laplace(coef(x.s5.est)[1],coef(x.s5.est)[2],coef(x.s5.est)[3],x.s5.em.x[l],coef(x.s5.est)[4]) #hand the laplace function the parameter estimates used to simulate scenario
}

###Plot###
par(mfrow=c(1,1))
##S1##
x1.plot <- plot(x.s1~x1, xlim = c(-4000, 4000), ylim = c(0, 105), xlab="Sample location relative to farm", 
               ylab="Pathogen density") 
lines(x.s1.true.y~x.s1.true.x) 
lines(x.s1.em.y~x.s1.em.x, lty = 2, col="red")
legend(950, -95, legend=c("Data generating model", "Estimated model"), col=c("black", "red"), lty=1:2, cex=0.6)

##S2##
x2.plot <- plot(x.s2~x2, xlim = c(-4000, 4000), ylim = c(0, 105), xlab="Sample location relative to farm", 
               ylab="Pathogen density") 
lines(x.s2.true.y~x.s2.true.x) 
lines(x.s2.em.y~x.s2.em.x, lty = 2, col="red")
#legend(950, 95, legend=c("Data generating model", "Estimated model"), col=c("black", "red"), lty=1:2, cex=0.6)

##S3##
x3.plot <- plot(x.s3~x3, xlim = c(-4000, 4000), ylim = c(0, 105), xlab="Sample location relative to farm", 
               ylab="Pathogen density") 
lines(x.s3.true.y~x.s3.true.x) 
lines(x.s3.em.y~x.s3.em.x, lty = 2, col="red")
#legend(950, 95, legend=c("Data generating model", "Estimated model"), col=c("black", "red"), lty=1:2, cex=0.6)

##S4##
x4.plot <- plot(x.s4~x4, xlim = c(-4000, 4000), ylim = c(0, 105), xlab="Sample location relative to farm", 
                ylab="Pathogen density") 
lines(x.s4.true.y~x.s4.true.x) 
lines(x.s4.em.y~x.s4.em.x, lty = 2, col="red")
#legend(950, 95, legend=c("Data generating model", "Estimated model"), col=c("black", "red"), lty=1:2, cex=0.6)

##S5##
x5.plot <- plot(x.s5~x5, xlim = c(-4000, 4000), ylim = c(0, 105), xlab="Sample location relative to farm", 
                ylab="Pathogen density") 
lines(x.s5.true.y~x.s5.true.x) 
lines(x.s5.em.y~x.s5.em.x, lty = 2, col="red")
#legend(950, 95, legend=c("Data generating model", "Estimated model"), col=c("black", "red"), lty=1:2, cex=0.6)



