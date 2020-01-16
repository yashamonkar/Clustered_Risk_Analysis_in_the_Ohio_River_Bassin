#### Code for HMM Output Visualizations####
#We take the input from the generated Quantities in Stan and compare them with the 'depmixS4' plot too

####################Intialization##############
library(dplyr)
library(ggplot2)
print("Library Loading Complete")


####Reading the Data. 
pois_par <- read.table("Emission_Parameters.txt", sep=" ", header = TRUE)
transition_matrix <- read.table('Transition_Parameters.txt', sep=" ", header = TRUE)
hidden_states <- read.table("Hidden_States.txt", sep=" ", header = TRUE)
input_data <- read.table("~/Correlated Risk Analysis/Plots and Analysis for Paper/STAN MODELS TEST/stan_hmm/data/Annual_Maximum_Streamflows.txt", sep=" ", header = TRUE)


####Aggregating the Data. 
input_data$Year <- NULL
thres_quantile <- apply( input_data, 2, quantile, probs = c(0.9), na.rm = TRUE)
Ann_Max_Proxy <- matrix(NA, nrow = nrow(input_data), ncol=ncol(input_data))
for(i in 1:ncol(Ann_Max_Proxy)){
  Ann_Max_Proxy[,i] <- ifelse(input_data[,i]>=thres_quantile[i], 1, 0)
}
ann_count <- rowSums(Ann_Max_Proxy)
Count_prox <- as.data.frame(ann_count)
colnames(Count_prox) <- c("Count")


#####Emission Probability Distribution####
depmix_pois_par <- c(0.77, 4.56, 12.34)
bx_title <- c("Dry State", "Intermediate State", "Wet State")
par(mfrow = c(1,3), mar = c(4,2.5,5,2))
for(i in 1:ncol(pois_par)){
  boxplot(pois_par[i], main = bx_title[i], cex.main=1)
  abline(h=depmix_pois_par[i], col='red', lwd = 2)
}
mtext("Estimated Poisson Rate Parameters", side = 3, line = -1.5, outer = TRUE)

par(mfrow = c(1,1), mar = c(3,3,3,3))
for(i in 1:ncol(pois_par)){
  plot(density(pois_par[,i]), 
       main = paste0("Density of Posterior Predictive Distributions for ", bx_title[i]),
       cex.main = 0.75)
  abline(v=depmix_pois_par[i], col='red', lwd = 2)
  legend('topright', legend = c("DepmixS4 estimate"), col = 'red', lwd = 2, cex = 0.75)
}


#####Transition Probability Distribution####

#The Dry-to-Wet and Wet-to-Wet transitions
par(mfrow = c(1,2), mar = c(2,4,3,1))
boxplot(transition_matrix[,7], 
        main = "Dry to Wet",ylim = c(0,1), 
        ylab = "Probability of Transition", outline = FALSE)
abline(h=0, col='red', lwd = 2)

boxplot(transition_matrix[,9], 
        main = "Wet to Wet",ylim = c(0,1),
        outline = FALSE)
abline(h=0, col='red', lwd = 2)



#Complete Transition Matrix
par(mfrow = c(3,3), mar = c(2,4,2,2))
boxplot(transition_matrix[,1], 
        main = "Dry to Dry",ylim = c(0,1), cex.main = 0.95, 
        ylab = "Probability of Transition", outline = FALSE)
abline(h=0.54, col='red', lwd = 2)

boxplot(transition_matrix[,4], 
        main = "Dry to Intermediate",ylim = c(0,1), cex.main = 0.95, 
        outline = FALSE)
abline(h=0.46, col='red', lwd = 2)

boxplot(transition_matrix[,7], 
        main = "Dry to Wet",ylim = c(0,1), cex.main = 0.95, 
        outline = FALSE)
abline(h=0.0, col='red', lwd = 2)

boxplot(transition_matrix[,2], 
        main = "Intermediate to Dry",ylim = c(0,1), cex.main = 0.95, 
        ylab = "Probability of Transition", outline = FALSE)
abline(h=0.47, col='red', lwd = 2)

boxplot(transition_matrix[,5], 
        main = "Intermediate to Intermediate",ylim = c(0,1), cex.main = 0.95, 
        outline = FALSE)
abline(h=0.36, col='red', lwd = 2)

boxplot(transition_matrix[,8], 
        main = "Intermediate to Wet",ylim = c(0,1), cex.main = 0.95, 
        outline = FALSE)
abline(h=0.17, col='red', lwd = 2)

boxplot(transition_matrix[,3], 
        main = "Wet to Dry",ylim = c(0,1), cex.main = 0.95, 
        ylab = "Probability of Transition", outline = FALSE)
abline(h=0.69, col='red', lwd = 2)

boxplot(transition_matrix[,6], 
        main = "Wet to Intermediate",ylim = c(0,1), cex.main = 0.95, 
        outline = FALSE)
abline(h=0.32, col='red', lwd = 2)

boxplot(transition_matrix[,9], 
        main = "Wet to Wet",ylim = c(0,1), cex.main = 0.95, 
        outline = FALSE)
abline(h=0.0, col='red', lwd = 2)



###Hidden States
library(tidyr)
par(mfrow = c(1,1), mar = c(2,2,2,2))
dry <- hidden_states %>%
  gather(x, value,) %>%
  group_by(x)%>%
  tally(value == 1)
dry <- dry$n/(dim(hidden_states)[1])

intermediate <- hidden_states %>%
  gather(x, value,) %>%
  group_by(x)%>%
  tally(value == 2)
intermediate <- intermediate$n/(dim(hidden_states)[1])

wet <- hidden_states %>%
  gather(x, value,) %>%
  group_by(x)%>%
  tally(value == 3)
wet <- wet$n/(dim(hidden_states)[1])

combined <- as.data.frame(cbind(Count_prox$Count, dry, intermediate, wet))
cent <- rep(NA,nrow(combined))
for(i in 1:nrow(combined)){
  cent[i] = max(dry[i],wet[i],intermediate[i])
}

modefunc <- function(x){
  tabresult <- tabulate(x)
  themode <- which(tabresult == max(tabresult))
  if(sum(tabresult == max(tabresult))>1) themode <- NA
  return(themode)
}
est_state <- apply(hidden_states, 2, modefunc)

plot(1937:2017, combined$V1, type='l',
     main = "Count Exceedances")
points(1937:2017, combined$V1, pch=(14+3*est_state), col='red',
       cex = (1-cent)*4)
########################################



#####Simulate observations for each iteration in the sample########
y_hat <- list()
for (i in 1:nrow(hidden_states)) {
  psi_seq <- sapply(hidden_states[i,], function(x){pois_par[i,x]})
  y_hat[[i]] <- rpois(length(psi_seq), psi_seq)
}


#Mean, Sd and Maximum
mean_sim <- sd_sim <- max_sim <- rep(NA,length(y_hat))
for(i in 1:length(y_hat)){
  mean_sim[i] <- mean(y_hat[[i]])
  sd_sim[i] <- sd(y_hat[[i]])
  max_sim[i] <- max(y_hat[[i]])
}

par(mfrow=c(1,3), mar = c(2,2,4,2))
boxplot(mean_sim, main = 'Mean')
abline(h= mean(Count_prox$Count), lwd =2, col='red')

boxplot(sd_sim, main = "Standard Deviation")
abline(h= sd(Count_prox$Count), lwd =2, col='red')

boxplot(max_sim, main = 'Maximum')
abline(h= max(Count_prox$Count), lwd =2, col='red')

mtext("Simulation Skill", side = 3, line = -1.5, outer = TRUE)
par(mfrow=c(1,1))


#####Probability Distribution Function
og_density <- density(Count_prox$Count)
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
consolid_df$x <- cut(consolid_df$x, breaks = seq(-1,25,.5))
mid_points <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", consolid_df$x) ),
                    upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", consolid_df$x) ))
consolid_df$x <- rowMeans(mid_points)
og_df <- data.frame(x1=og_density$x, y1= og_density$y)

plo <- ggplot(og_df, aes(x=x1,y=y1))+
  geom_line(size=1.25)+
  scale_x_continuous(limits=c(-2,25)) +
  geom_boxplot(consolid_df, mapping = aes(y = y, x = x,group = x), outlier.shape = NA,outlier.colour = NA)+ 
  ggtitle("Simulated Annual Count Exceedances") +
  xlab("Count Exceedances") + 
  ylab("Probability Density") +
  theme(plot.title = element_text(size = 10, face = "bold"))
print(plo)


####Simulations across the board##
output_sims <- matrix(unlist(y_hat), ncol = 81, byrow = TRUE)

par(mfrow = c(1,1), mar = c(4,4,4,2))
plot(1937:2017, Count_prox$Count, type='l',
     main = "Count Exceedances", xlab = "Year",
     ylab= "Exccedances", ylim = c(-3,25))

#Getting the percentiles
lower_percentile <- apply(output_sims, 2, function(x) quantile(x, probs=.05))
upper_percentile <- apply(output_sims, 2, function(x) quantile(x, probs=.95))
median_percentile <- apply(output_sims, 2, median)

#Plotting the CDF
par(mfrow=c(1,1), mar = c(4,4,3,1))
plot(1937:2017, Count_prox$Count, type='l',col='red',
     lwd = 2,  main = "Observed and Simulated Exceedances",
     xlab = "Year", ylab= "Exccedances", ylim = c(-1,20))
polygon(c(1937:2017,rev(1937:2017)),c(lower_percentile,rev(upper_percentile)),col="gray")
lines(1937:2017, median_percentile, lwd = 2)
lines(1937:2017, Count_prox$Count, col='red', lwd = 2)
legend('topright', legend = c("Observed Values","Median-Simulations","90% Uncertainty Interval"),
       lty = 1, col = c('red','black','grey'), lwd = 2, cex = 0.75)

