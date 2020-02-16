#### Code for HMM with Emission Probability as Normal. 

####################Intialization##############
.libPaths("/rigel/cwc/users/yva2000/rpackages/")
library(rstan)
library(dplyr)
options(mc.cores = parallel::detectCores()) #Local Multicore CPU excess 
print("Library Loading Complete")


#Fake Data Simulation
pi1 <- matrix(c(0.54,0.44,0.02), ncol = 3)
A <- matrix(c(0.54,0.44,0.02, 0.15,0.60,0.25, 0.40,0.60,0.0), ncol = 3, byrow = TRUE)
lambda <- c(0.77,4.5,12)

#2. Hidden Path
z <- vector("numeric", 100)

z[1] <- sample(1:3, size = 1, prob = pi1)
for(t in 2:100) {
  z[t] <- sample(1:3, size = 1, prob = A[z[t-1],])}

#3. Observations
y <- vector("numeric", 100)
for(t in 1:100){
  y[t] <- rpois(1, lambda[z[t]])}
plot(y, type='l', main = "Simulated Data")





# code available in hmm_example.R
stan_data <- list(N = length(y),
                  K = 3,
                  y = y)
hmm_fit <- stan("Poisson_HMM.stan", data = stan_data, iter = 4e3, chains = 6)
print(hmm_fit, pars = "z_star", include = FALSE, probs = c(0.05,0.95))

# extract samples
samples <- as.matrix(hmm_fit)
print(dim(samples))
theta <- samples[,grep("^theta",colnames(samples))]
print(dim(theta))
write.table(theta, "Transition_Parameters.txt", sep = " ")
mu <- samples[,grep("^mu",colnames(samples))]
print(dim(mu))
write.table(mu, "Emission_Parameters.txt", sep = " ")
z_star <- samples[,grep("^z_star",colnames(samples))]
print(dim(z_star))
write.table(z_star, "Hidden_States.txt", sep = " ")

# simulate observations for each iteration in the sample
y_hat <- list()
for (i in 1:nrow(samples)) {
  psi_seq <- sapply(z_star[i,], function(x){mu[i,x]})
  y_hat[[i]] <- rpois(length(psi_seq), psi_seq)
}

# plot
indxs <- sample(length(y_hat), 200, replace = FALSE)
plot(y, type = "n",
     main = "Observed vs Predicted Output",
     ylab = "Observation Value",
     xlab = "Time",
     ylim = c(-1,20))
for (i in indxs) {
  lines(y_hat[[i]], col = "#ff668890")
}
lines(y, lwd = 2)
legend("bottomright", c("Observed","Simulated"), col = c("#000000","#ff668890"), lty = c(1,1), lwd = c(2,1), cex = 0.8)

# plot
indxs <- sample(length(y_hat), 200, replace = FALSE)
og_dens <- density(y)
plot(og_dens$x,og_dens$y, type = "n",
     main = "Observed vs Predicted Output",
     ylab = "Observation Value",
     xlab = "Time",
     xlim = c(-1,30),
     ylim = c(0, 0.4))
for (i in indxs) {
  dens <- density(y_hat[[i]])
  lines(dens$x,dens$y, col = "#ff668890")
}
lines(og_dens$x,og_dens$y, lwd = 2)
legend("bottomright", c("Observed","Simulated"), col = c("#000000","#ff668890"), lty = c(1,1), lwd = c(2,1), cex = 0.8)
