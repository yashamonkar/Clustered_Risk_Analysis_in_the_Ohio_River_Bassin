#### Code for HMM with Emission Probability as Normal. 

####################Intialization##############
.libPaths("/rigel/cwc/users/yva2000/rpackages/")
library(rstan)
library(dplyr)
options(mc.cores = parallel::detectCores()) #Local Multicore CPU excess 
print("Library Loading Complete")
source("Get_Simulation_Skill.R")

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
mu <- samples[,grep("^mu",colnames(samples))]
z_star <- samples[,grep("^z_star",colnames(samples))]

# simulate observations for each iteration in the sample
y_hat <- list()
for (i in 1:nrow(samples)) {
  psi_seq <- sapply(z_star[i,], function(x){mu[i,x]})
  y_hat[[i]] <- rpois(length(psi_seq), psi_seq)
}

y_hat <- matrix(unlist(y_hat), ncol = nrow(samples))

get_SimSkill(og_data = Count_prox$Count,
             Sim_Mat = y_hat,
             name = "Exceedance Counts", 
             moments = TRUE,
             prob_dens = TRUE, 
             cumm_dens = TRUE,
             auto = TRUE,
             wavelt = TRUE,
             lagged = TRUE,
             N_Sims = 2000)