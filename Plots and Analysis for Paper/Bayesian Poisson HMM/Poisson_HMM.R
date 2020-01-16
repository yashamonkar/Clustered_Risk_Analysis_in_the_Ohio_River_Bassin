#### Code for HMM with Emission Probability as Normal. 

####################Intialization##############
.libPaths("/rigel/cwc/users/yva2000/rpackages/")
library(rstan)
library(dplyr)
options(mc.cores = parallel::detectCores()) #Local Multicore CPU excess 
print("Library Loading Complete")


#Reading the Data. 
input_data <- read.table("data/Annual_Maximum_Streamflows.txt", sep=" ", header = TRUE)
Site_Info <- read.table('data/Site_Information.txt', sep=" ", header = TRUE)

#Read the site data.
Ann_Max <- input_data
Site_info <- Site_Info
start_year <- min(Ann_Max$Year);end_year <- max(Ann_Max$Year)

#Converting Site Specific to Regional
Ann_Max$Year <- NULL
thres_quantile <- apply( Ann_Max, 2, quantile, probs = c(0.9), na.rm = TRUE)
Ann_Max_Proxy <- matrix(NA, nrow = nrow(Ann_Max), ncol=ncol(Ann_Max))
for(i in 1:ncol(Ann_Max_Proxy)){
  Ann_Max_Proxy[,i] <- ifelse(Ann_Max[,i]>=thres_quantile[i], 1, 0)
}

###Hidden Markov Model - Poisson Process
ann_count <- rowSums(Ann_Max_Proxy)
Count_prox <- as.data.frame(ann_count)
colnames(Count_prox) <- c("Count")

# code available in hmm_example.R
stan_data <- list(N = length(ann_count),
                  K = 3,
                  y = Count_prox$Count)
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
plot(Count_prox$Count, type = "n",
     main = "Observed vs Predicted Output",
     ylab = "Observation Value",
     xlab = "Time",
     ylim = c(-1,20))
for (i in indxs) {
  lines(y_hat[[i]], col = "#ff668890")
}
lines(Count_prox$Count, lwd = 2)
legend("bottomright", c("Observed","Simulated"), col = c("#000000","#ff668890"), lty = c(1,1), lwd = c(2,1), cex = 0.8)

# plot
indxs <- sample(length(y_hat), 200, replace = FALSE)
og_dens <- density(Count_prox$Count)
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
