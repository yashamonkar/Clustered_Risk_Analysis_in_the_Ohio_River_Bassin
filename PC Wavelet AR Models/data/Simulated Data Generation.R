setwd("~/Correlated Risk Analysis/Decadal Influences/PC Wavelet AR Models")

#Simulating a streamflow field. 

#The objective of this field is to simulate a naturalized version of the 
#simulated streamflow. 
N_yrs <- 200 #This is a good number
N_stations <- 30 # Number of different stations. 
t <- seq(1,N_yrs,1)

sim_data <- matrix(NA,ncol = N_stations, nrow <- N_yrs)
regime <- rbinom(N_stations, 1, 0.75)
for(i in 1:N_stations) {
  if(regime[i]==1) {
    periods <- c(rnorm(2,12,0.35),rnorm(2,50,1.5))
  } else {
    periods <- c(rnorm(2,6,0.25),rnorm(2,32,1))
  } 
  weights <- c(runif(2,0.9,1),runif(2,0.75,0.85))
   temp <- rep(0,N_yrs) 
   for(j in 1:length(periods)) { temp <- temp + weights[j]*cos(2*pi*t/periods[j])}
   temp <- scale(temp) + rnorm(N_yrs,0,0.5)
   sim_data[,i] <- temp
   }


write.table(sim_data, "data/Simulated Series.txt", sep=" ")


