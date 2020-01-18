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

#Data Manipulation
years <- input_data$Year;input_data$Year <- NULL 
input_data <- log(input_data) #Brings it closer to multivariate normal

#Principal Component Analysis
pcs <- prcomp(input_data, scale = TRUE)
var <- cumsum(pcs$sdev^2)
p <- pcs$x[,1]


# code available in hmm_example.R
stan_data <- list(N = length(p),
                  K = 3,
                  y = p)
hmm_fit <- stan("PCA_HMM.stan", data = stan_data, iter = 4e3, chains = 6)
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
  y_hat[[i]] <- rnorm(length(psi_seq), psi_seq, 2)
}

# plot
indxs <- sample(length(y_hat), 200, replace = FALSE)
plot(p, type = "n",
     main = "Observed vs Predicted Output",
     ylab = "Observation Value",
     xlab = "Time",
     ylim = c(-15,15))
for (i in indxs) {
  lines(y_hat[[i]], col = "#ff668890")
}
lines(p, lwd = 2)
legend("bottomright", c("Observed","Simulated"), col = c("#000000","#ff668890"), lty = c(1,1), lwd = c(2,1), cex = 0.8)

# plot
indxs <- sample(length(y_hat), 200, replace = FALSE)
og_dens <- density(p)
plot(og_dens$x,og_dens$y, type = "n",
     main = "Observed vs Predicted Output",
     ylab = "Observation Value",
     xlab = "Time",
     xlim = c(-15,15),
     ylim = c(0, 0.4))
for (i in indxs) {
  dens <- density(y_hat[[i]])
  lines(dens$x,dens$y, col = "#ff668890")
}
lines(og_dens$x,og_dens$y, lwd = 2)
legend("bottomright", c("Observed","Simulated"), col = c("#000000","#ff668890"), lty = c(1,1), lwd = c(2,1), cex = 0.8)
#################




######Complete Plots#######

###Emission Parameters
par(mfrow = c(1,3), mar = c(4,2.5,5,2))
for(i in 1:ncol(mu)){
  boxplot(mu[,i], main = colnames(mu)[i], cex.main=1)
}
mtext("Estimated Poisson Rate Parameters", side = 3, line = -1.5, outer = TRUE)

###Complete Transition Matrix
par(mfcol = c(sqrt(ncol(theta)),sqrt(ncol(theta))), mar = c(2,4,2,2))
for(i in 1:ncol(theta)){
  boxplot(theta[,i], 
          main = colnames(theta)[i] ,ylim = c(0,1), cex.main = 0.95, 
          outline = FALSE)
}
par(mfrow = c(1,1), mar = c(2,4,2,2))

###Mean, Sd and Maximum
mean_sim <- sd_sim <- rep(NA,length(y_hat))
for(i in 1:length(y_hat)){
  mean_sim[i] <- mean(y_hat[[i]])
  sd_sim[i] <- sd(y_hat[[i]])
}

par(mfrow=c(1,2), mar = c(2,2,4,2))
boxplot(mean_sim, main = 'Mean')
abline(h= mean(p), lwd =2, col='red')

boxplot(sd_sim, main = "Standard Deviation")
abline(h= sd(p), lwd =2, col='red')

mtext("Simulation Skill", side = 3, line = -1.5, outer = TRUE)
par(mfrow=c(1,1))

#####Probability Distribution Function
og_density <- density(p)
consolid_sims <- as.list(1)
consolid_points <- as.list(1)
for(j in 1:length(y_hat)) {
  sims_dens <- density(y_hat[[j]])
  consolid_sims[[j]] <- sims_dens$y
  consolid_points[[j]] <- sims_dens$x
}
consolid_points <- unlist(consolid_points)
consolid_sims <- unlist(consolid_sims)
consolid_df <- data.frame(x=consolid_points,y=consolid_sims)
consolid_df$x <- cut(consolid_df$x, breaks = seq(-15,25,.5))
mid_points <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", consolid_df$x) ),
                    upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", consolid_df$x) ))
consolid_df$x <- rowMeans(mid_points)
og_df <- data.frame(x1=og_density$x, y1= og_density$y)

plo <- ggplot(og_df, aes(x=x1,y=y1))+
  geom_line(size=1.25)+
  scale_x_continuous(limits=c(-15,25)) +
  geom_boxplot(consolid_df, mapping = aes(y = y, x = x,group = x), outlier.shape = NA,outlier.colour = NA)+ 
  ggtitle("Simulated Annual Count Exceedances") +
  xlab("Count Exceedances") + 
  ylab("Probability Density") +
  theme(plot.title = element_text(size = 10, face = "bold"))
print(plo)


####Uncertainty
output_sims <- matrix(unlist(y_hat), ncol = length(p), byrow = TRUE)

#Getting the percentiles
lower_percentile <- apply(output_sims, 2, function(x) quantile(x, probs=.05))
upper_percentile <- apply(output_sims, 2, function(x) quantile(x, probs=.95))
median_percentile <- apply(output_sims, 2, median)

par(mfrow=c(1,1), mar = c(4,4,3,1))
plot(1937:2017, p, type='l',col='red',
     lwd = 2,  main = "PC-1",
     xlab = "Year", ylab= "PC Score", ylim = c(-15,20))
polygon(c(1937:2017,rev(1937:2017)),c(lower_percentile,rev(upper_percentile)),col="gray")
lines(1937:2017, median_percentile, lwd = 2)
lines(1937:2017, p, col='red', lwd = 2)
legend('topright', legend = c("Observed Values","Median-Simulations","90% Uncertainty Interval"),
       lty = 1, col = c('red','black','grey'), lwd = 2, cex = 0.75)
