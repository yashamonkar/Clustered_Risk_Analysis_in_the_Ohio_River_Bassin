#### Code for HMM Output Visualizations####
#We take the input from the generated Quantities in Stan and compare them with the 'depmixS4' plot too

####################Intialization##############
library(dplyr)
library(ggplot2)
library(gridExtra)
print("Library Loading Complete")


####Reading the Data. 
pois_par <- read.table("Emission_Parameters.txt", sep=" ", header = TRUE)
transition_matrix <- read.table('Transition_Parameters.txt', sep=" ", header = TRUE)
hidden_states <- read.table("Hidden_States.txt", sep=" ", header = TRUE)
input_data <- read.table("~/Correlated Risk Analysis/Plots and Analysis for Paper/data/Annual_Maximum_Streamflows.txt", sep=" ", header = TRUE)


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

pdf("HMM_Checks.pdf")
#####Emission Probability Distribution####
bx_title <- c("Dry State", "Intermediate State", "Wet State")
par(mfrow = c(1,3), mar = c(4,2.5,5,2))
for(i in 1:ncol(pois_par)){
  boxplot(pois_par[i], main = bx_title[i], cex.main=1.5)
}
mtext("Estimated Poisson Rate Parameters", side = 3, line = -1.5, outer = TRUE)

violin_dataset <- data.frame(mu1 = pois_par[,1], mu2 = pois_par[,2],mu3 = pois_par[,3], ind = 1)
p_dry <- ggplot(violin_dataset, aes(x=ind, y=mu1)) +
  geom_violin() + labs(y = " ", x = " ") +
  ggtitle(c("Dry")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 20, face = "bold"))

p_int <- ggplot(violin_dataset, aes(x=ind, y=mu2)) +
  geom_violin() + labs(y = " ", x = " ") +
  ggtitle(c("Intermediate")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 20, face = "bold"))

p_wet <- ggplot(violin_dataset, aes(x=ind, y=mu3)) +
  geom_violin(mapping = NULL) + labs(y = " ", x = " ") +
  ggtitle(c("Wet")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 20, face = "bold"))

grid.arrange(p_dry, p_int, p_wet,ncol=3)
#########################################################################

#####Transition Probability Distribution####

#Complete Transition Matrix
par(mfrow = c(3,3), mar = c(2,4,2,2))
boxplot(transition_matrix[,1], 
        main = "Dry to Dry",ylim = c(0,1), cex.main = 1.2, 
        ylab = "Probability of Transition", outline = FALSE)

boxplot(transition_matrix[,4], 
        main = "Dry to Intermediate",ylim = c(0,1), cex.main = 1.2, 
        outline = FALSE)

boxplot(transition_matrix[,7], 
        main = "Dry to Wet",ylim = c(0,1), cex.main = 1.2, 
        outline = FALSE)

boxplot(transition_matrix[,2], 
        main = "Intermediate to Dry",ylim = c(0,1), cex.main = 1.2, 
        ylab = "Probability of Transition", outline = FALSE)

boxplot(transition_matrix[,5], 
        main = "Intermediate to Intermediate",ylim = c(0,1), cex.main = 1.2, 
        outline = FALSE)

boxplot(transition_matrix[,8], 
        main = "Intermediate to Wet",ylim = c(0,1), cex.main = 1.2, 
        outline = FALSE)

boxplot(transition_matrix[,3], 
        main = "Wet to Dry",ylim = c(0,1), cex.main = 1.2, 
        ylab = "Probability of Transition", outline = FALSE)

boxplot(transition_matrix[,6], 
        main = "Wet to Intermediate",ylim = c(0,1), cex.main = 1.2, 
        outline = FALSE)

boxplot(transition_matrix[,9], 
        main = "Wet to Wet",ylim = c(0,1), cex.main = 1.2, 
        outline = FALSE)

#Violin Plots
violin_dataset <- transition_matrix
violin_dataset$ind <-1
d_to_d <- ggplot(violin_dataset, aes(x=ind, y=theta.1.1.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Dry|Dry)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 12, face = "bold"))

d_to_i <- ggplot(violin_dataset, aes(x=ind, y=theta.1.2.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Intermediate|Dry)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 12, face = "bold"))

d_to_w <- ggplot(violin_dataset, aes(x=ind, y=theta.1.3.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Wet|Dry)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 12, face = "bold"))

i_to_d <- ggplot(violin_dataset, aes(x=ind, y=theta.2.1.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Dry|Intermediate)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 12, face = "bold"))

i_to_i <- ggplot(violin_dataset, aes(x=ind, y=theta.2.2.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Intermediate|Intermediate)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 9, face = "bold"))

i_to_w <- ggplot(violin_dataset, aes(x=ind, y=theta.2.3.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Wet|Intermediate)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 12, face = "bold"))

w_to_d <- ggplot(violin_dataset, aes(x=ind, y=theta.3.1.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Dry|Wet)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 12, face = "bold"))

w_to_i <- ggplot(violin_dataset, aes(x=ind, y=theta.3.2.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Intermediate|Wet)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 12, face = "bold"))

w_to_w <- ggplot(violin_dataset, aes(x=ind, y=theta.3.3.)) +
  geom_violin() + labs(y = "P(Transition)", x = " ") + ylim(0, 1) +
  ggtitle(c("P(Wet|Wet)")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 12, face = "bold"))



grid.arrange(d_to_d, d_to_i, d_to_w,
             i_to_d, i_to_i, i_to_w,
             w_to_d, w_to_i, w_to_w,
             ncol=3, nrow = 3)

#####Simulate observations for each iteration in the sample########
y_hat <- list()
for (i in 1:nrow(hidden_states)) {
  psi_seq <- sapply(hidden_states[i,], function(x){pois_par[i,x]})
  y_hat[[i]] <- rpois(length(psi_seq), psi_seq)
}


####Simulations across the board##
output_sims <- matrix(unlist(y_hat), ncol = nrow(input_data), byrow = TRUE)

par(mfrow = c(1,1), mar = c(4,4,4,2))
plot(1935:2019, Count_prox$Count, type='l',
     main = "Count Exceedances", xlab = "Year",
     ylab= "Exccedances", ylim = c(-1,20))

#Getting the percentiles
lower_percentile <- apply(output_sims, 2, function(x) quantile(x, probs=.05))
upper_percentile <- apply(output_sims, 2, function(x) quantile(x, probs=.95))
median_percentile <- apply(output_sims, 2, median)

#Plotting the Density
par(mfrow=c(1,1), mar = c(4,4,3,1))
plot(1935:2019, Count_prox$Count, type='l',col='red',
     lwd = 2,  main = "Observed and Simulated Exceedances",
     xlab = "Year", ylab= "Exccedances", ylim = c(-1,20))
polygon(c(1935:2019,rev(1935:2019)),c(lower_percentile,rev(upper_percentile)),col="gray")
lines(1935:2019, median_percentile, lwd = 2)
lines(1935:2019, Count_prox$Count, col='red', lwd = 3)
legend('topright', legend = c("Observed Values","Median-Simulations","90% Uncertainty Interval"),
       lty = 1, col = c('red','black','grey'), lwd = 2, cex = 0.75)

dev.off()