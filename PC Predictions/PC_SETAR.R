#####
#1. The objective of this file is to make predictions in PC space. 
#2. Use a dimension lowering algorithm - linear PCs
#3. Make predictions and simulations on those PCs. 
#4. Check if the relevant indices and the spectral structure has been replicated. 


#########Settting the Path
setwd("~/Correlated Risk Analysis/Decadal Influences/PC Predictions")


###Loading Libraries##########
library(forecast)
library(tsDyn)


###Reading the Data and Initial Manipulation###
input_data <- read.table("data/Max_Annual_Streamflow.txt", sep="", header = TRUE)
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)
######################################2


####Principal Component Analysis#############
data_pca <- prcomp(input_data, scale = TRUE)
var <- cumsum(data_pca$sdev^2)
pdf(file = "plots/SETAR/Variance explained.pdf")
plot(var/max(var),pch=19, main = "Variance explained by PCs",xlab = "PC's",ylab="Fraction Variance explained")
npcs <- 3 #This is the selected number of PC'S
abline(h = var[npcs]/max(var), lty = 2, col='red')
dev.off()
pcs_sel <- data_pca$x[,1:npcs]
#################################################2


###############Diagnostics##########
pdf(file = "plots/SETAR/PC Diagnostics.pdf")
for(i in 1:ncol(pcs_sel)) {

x <- pcs_sel[,i]

#Plotting the Data
  plot(pcs_sel[,i], main = paste0("PC - ", i), type = 'l', 
       xlab = "Years", ylab = " ")

#Time Reversibility
par(mfrow=c(2,1))
par(mar = c(1,1,1,1))
   x <- pcs_sel[,i]
  plot(x, main = paste0("PC - ", i), type = 'l', 
       xlab = "Years", ylab = " ")
  plot(x[length(x):1], type="l", ax=F)
  box()

#Lag Plots
par(mfrow=c(2,2), mar=c(4,1,1,1))
autopairs(x, lag=1, type="regression")
autopairs(x, lag=2, type="regression")
autopairs(x, lag=3, type="regression")
autopairs(x, lag=4, type="regression")
autopairs(x, lag=5, type="regression")
autopairs(x, lag=6, type="regression")
autopairs(x, lag=7, type="regression")
autopairs(x, lag=8, type="regression")


#Histograms
par(mfrow=c(1,1), mar=c(4,2,4,1))
  hist(pcs_sel[,i], main = paste0("PC - ", i), 
       xlab = " ")


#PACF and ACF's
par(mfrow=c(2,1), mar=c(2,4,4,1))
acf(x, main = paste0("PC -",i ))
pacf(x, main = paste0("PC -",i )) 

#Lag Plots
lag.plot(x, lags=3, layout=c(1,3))
par(mfrow=c(1,1))
}
dev.off()
######################SETAR Model################
#SETAR is Self Exciting Threshold AutoRegressive Model 
#Reading is Howell Tong - Non-Linear Time Series. 

for(i in 1:ncol(pcs_sel)){
x <- pcs_sel[,i]
if( i< 3) {x <- scale(x)}
file_nam <- paste0("plots/SETAR/Simulated PC_",i,".pdf")
pdf(file = file_nam)
for(k in 1:5) {
#SETAR
embd <- k
selections <- selectSETAR(x, m=embd, thDelay=0:(embd-1))
sel_model <- selections$firstBests
mod.setar <-  setar(x, m = embd, thDelay = sel_model[1],
                    th = sel_model[2], mL = sel_model[4],
                    mH = sel_model[5])
summary(mod.setar)
#plot(mod.setar)

#SETAR Simulations
N_Sims <- 1000
mean_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
sd_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
min_sim <- max_sim <- matrix(NA, ncol = 1, nrow = N_Sims)

ml <- as.vector(head(mod.setar$coefficients,sel_model[4]+1)) #Getting the lower regime coefficients
mh <- as.vector(head(mod.setar$coefficients,-1)) 
mh <- tail(mh,-sel_model[4]-1) #Getting the upper regime coefficients
ML <- MH <- rep(0,embd+1)
ML[1:length(ml)] <- ml
MH[1:length(mh)] <- mh
B <- c(ML,MH)
for(j in 1:N_Sims) {
  sims.setar <- setar.sim(B = B,
                          lag = embd,
                          nthresh = 1, 
                          Thresh = tail(mod.setar$coefficients,1), 
                          type = "simul", 
                          n = 81)
  mean_sim[j,] <- mean(sims.setar$serie)
  sd_sim[j,] <- sd(sims.setar$serie)
  min_sim[j,] <- min(sims.setar$serie)
  max_sim[j,] <- max(sims.setar$serie)
}

par(mfrow=c(1,4))
par(mar = c(4, 2, 1.5, 1))

boxplot(mean_sim, 
        ylim = c(1,-1))
title(paste0("Mean-PC ", i))
abline(h=mean(x), col = 'red')

boxplot(sd_sim, 
        ylim = c(.5,3))
title(paste0("SD - PC ", i))
abline(h=sd(x), col = 'red')

boxplot(max_sim, 
        ylim = c(0,20))
title(paste0("Max - PC ", i))
abline(h=max(x), col = 'red')

boxplot(min_sim, 
        ylim = c(-20,2))
title(paste0("Min - PC ", i))
abline(h=min(x), col = 'red')
par(mfrow=c(1,1))

}
dev.off()

 }

##################Simulating the PCs Data. 
N_Sims <- 1000 


###First PC###
x <- pcs_sel[,1]
x_1st <- scale(x,scale = TRUE)
embd_1st <- 4
selections_1st <- selectSETAR(x_1st, m=embd_1st, thDelay=0:(embd_1st-1))
sel_model_1st <- selections_1st$firstBests
mod.setar_1st <-  setar(x_1st, m = embd_1st, thDelay = sel_model_1st[1],
                    th = sel_model_1st[2], mL = sel_model_1st[4],
                    mH = sel_model_1st[5])
ml <- as.vector(head(mod.setar_1st$coefficients,sel_model_1st[4]+1)) #Getting the lower regime coefficients
mh <- as.vector(head(mod.setar_1st$coefficients,-1)) 
mh <- tail(mh,-sel_model_1st[4]-1) #Getting the upper regime coefficients
ML <- MH <- rep(0,embd_1st+1)
ML[1:length(ml)] <- ml
MH[1:length(mh)] <- mh
B_1st <- c(ML,MH)
sims.setar <- setar.sim(B = B_1st,
                        lag = embd_1st,
                        nthresh = 1, 
                        Thresh = tail(mod.setar_1st$coefficients,1), 
                        type = "simul", 
                        n = 81)


###Second PC###
x <- pcs_sel[,2]
x_2nd <- scale(x,scale = TRUE)
embd_2nd <- 3
selections_2nd <- selectSETAR(x_2nd, m=embd_2nd, thDelay=0:(embd_2nd-1))
sel_model_2nd <- selections_2nd$firstBests
mod.setar_2nd <-  setar(x_2nd, m = embd_2nd, thDelay = sel_model_2nd[1],
                        th = sel_model_2nd[2], mL = sel_model_2nd[4],
                        mH = sel_model_2nd[5])
ml <- as.vector(head(mod.setar_2nd$coefficients,sel_model_2nd[4]+1)) #Getting the lower regime coefficients
mh <- as.vector(head(mod.setar_2nd$coefficients,-1)) 
mh <- tail(mh,-sel_model_2nd[4]-1) #Getting the upper regime coefficients
ML <- MH <- rep(0,embd_2nd+1)
ML[1:length(ml)] <- ml
MH[1:length(mh)] <- mh
B_2nd <- c(ML,MH)


###Third PC###
x_3rd <- pcs_sel[,3]
embd_3rd <- 4
selections_3rd <- selectSETAR(x_3rd, m=embd_3rd, thDelay=0:(embd_3rd-1))
sel_model_3rd <- selections_3rd$firstBests
mod.setar_3rd <-  setar(x_3rd, m = embd_3rd, thDelay = sel_model_3rd[1],
                        th = sel_model_3rd[2], mL = sel_model_3rd[4],
                        mH = sel_model_3rd[5])
ml <- as.vector(head(mod.setar_3rd$coefficients,sel_model_3rd[4]+1)) #Getting the lower regime coefficients
mh <- as.vector(head(mod.setar_3rd$coefficients,-1)) 
mh <- tail(mh,-sel_model_3rd[4]-1) #Getting the upper regime coefficients
ML <- MH <- rep(0,embd_3rd+1)
ML[1:length(ml)] <- ml
MH[1:length(mh)] <- mh
B_3rd <- c(ML,MH)



####Combined Simulations###
N_sims <- 1000
mean_sim <- matrix(NA, ncol = ncol(input_data), nrow = N_Sims)
sd_sim <- matrix(NA, ncol = ncol(input_data), nrow = N_Sims)
min_sim <- max_sim <- matrix(NA, ncol = ncol(input_data), nrow = N_Sims)


for(sim in 1:N_sims) {
  sim_pcs <- matrix(NA, ncol = npcs, nrow = length(x))
  
  sims.setar_1st <- setar.sim(B = B_1st,
                        lag = embd_1st,
                        nthresh = 1, 
                        Thresh = tail(mod.setar_1st$coefficients,1), 
                        type = "simul", 
                        n = 81)
  sim_pcs[,1] <- sims.setar_1st$serie
  sim_pcs[,1] <- sim_pcs[,1]*sd(pcs_sel[,1]) #Resaling back

  sims.setar_2nd <- setar.sim(B = B_2nd,
                        lag = embd_2nd,
                        nthresh = 1, 
                        Thresh = tail(mod.setar_2nd$coefficients,1), 
                        type = "simul", 
                        n = 81)
  sim_pcs[,2] <- sims.setar_2nd$serie
  sim_pcs[,2] <- sim_pcs[,2]*sd(pcs_sel[,2]) #Resaling back

  sims.setar_3rd <- setar.sim(B = B_3rd,
                        lag = embd_3rd,
                        nthresh = 1, 
                        Thresh = tail(mod.setar_3rd$coefficients,1), 
                        type = "simul", 
                        n = 81)
  sim_pcs[,3] <- sims.setar_3rd$serie
  
  #Converting to Actual Field Space
  PC_Simulations <- sim_pcs #These are the predictions
  nComp = ncol(PC_Simulations)
  Simulations_Scaled =  PC_Simulations %*% t(data_pca$rotation[,1:nComp])
  for(i in 1:ncol(Simulations_Scaled)) {   Simulations_Scaled[,i] <- scale(Simulations_Scaled[,i], center = FALSE , scale=1/data_pca$scale[i]) }
  for(i in 1:ncol(Simulations_Scaled)) {   Simulations_Scaled[,i] <- scale(Simulations_Scaled[,i], center = -1 * data_pca$center[i], scale=FALSE)}

  #Computing the Site Specific Error
  mean_sim[sim,] <- colMeans(Simulations_Scaled)
  sd_sim[sim,] <- apply(Simulations_Scaled, 2, sd)
  max_sim[sim,] <- apply(Simulations_Scaled, 2, max)
  min_sim[sim,] <- apply(Simulations_Scaled, 2, min)
  
}

#####Plotting the Simulations by Moments######
##Plotting the means##
pdf(file = "plots/SETAR/Site Mean Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))

for(j in 1:ncol(input_data)) {

  boxplot(mean_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(mean_sim[,j], .5)-.2,quantile(mean_sim[,j], .95)+.2))
  title(paste0("Mean-Site ", j), cex = 0.5)
  abline(h=mean(input_data[,j]), col = 'red')
  
}
dev.off()

##Plotting the Standard Deviation###
pdf(file = "plots/SETAR/Site Standard Deviation Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(sd_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(sd_sim[,j], .5)-.2,quantile(sd_sim[,j], .95)+.2))
  title(paste0("SD-Site ", j), cex = 0.5)
  abline(h=sd(input_data[,j]), col = 'red')
  
}
dev.off()

##Plotting the Max###
pdf(file = "plots/SETAR/Site Max Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(max_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(max_sim[,j], .5)-.2,quantile(max_sim[,j], .95)+.2))
  title(paste0("Max-Site ", j), cex = 0.5)
  abline(h=max(input_data[,j]), col = 'red')
  
}
dev.off()

##Plotting the Min###
pdf(file = "plots/SETAR/Site Min Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(min_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(min_sim[,j], .5)-.2,quantile(min_sim[,j], .95)+.2))
  title(paste0("Min-Site ", j), cex = 0.5)
  abline(h=min(input_data[,j]), col = 'red')
  
}
dev.off()

#####Plotting the Simulations by Sites######
pdf(file = "plots/SETAR/Site Specific Simulation Moments.pdf")
par(mfrow=c(1,4))
par(mar = c(4, 2, 1.5, 1))

for(j in 1:ncol(input_data)) {
  
  #Mean
  boxplot(mean_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(mean_sim[,j], .5)-.2,quantile(mean_sim[,j], .95)+.2))
  title(paste0("Mean-Site ", j), cex = 0.5)
  abline(h=mean(input_data[,j]), col = 'red')
  
  #Standard Deviation
  boxplot(sd_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(sd_sim[,j], .5)-.2,quantile(sd_sim[,j], .95)+.2))
  title(paste0("SD-Site ", j), cex = 0.5)
  abline(h=sd(input_data[,j]), col = 'red')
  
  #Max
  boxplot(max_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(max_sim[,j], .5)-.2,quantile(max_sim[,j], .95)+.2))
  title(paste0("Max-Site ", j), cex = 0.5)
  abline(h=max(input_data[,j]), col = 'red')
  
  #Min
  boxplot(min_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(min_sim[,j], .5)-.2,quantile(min_sim[,j], .95)+.2))
  title(paste0("Min-Site ", j), cex = 0.5)
  abline(h=min(input_data[,j]), col = 'red')
}

dev.off()





##########Notes##############
For the First PC
Warning message:
  Possible unit root in the low  regime. Roots are: 0.9761 1.2779
Scaling helps. 
1,4,5 do good, but all are good too.
thDelay          th  pooled-AIC          mL          mH 
0.0000000  -0.7570516 223.8116432   1.0000000   1.0000000 

thDelay          th  pooled-AIC          mL          mH 
2.0000000   0.4894034 211.4874546   4.0000000   4.0000000 

thDelay          th  pooled-AIC          mL          mH 
0.0000000  -0.7015116 193.8520008   5.0000000   1.0000000


For the Second PC
Scaling helps. helps
m = 1,2,3 are the best onses. 
thDelay         th pooled-AIC         mL         mH 
0.000000  -0.431934 213.574713   1.000000   1.000000 
thDelay         th pooled-AIC         mL         mH 
0.000000  -0.431934 210.112013   1.000000   2.000000
thDelay          th  pooled-AIC          mL          mH 
0.0000000  -0.2869277 205.1642087   3.0000000   2.0000000



For the third PC. 
m = 4 is great. 
thDelay         th pooled-AIC         mL         mH 
3.000000  -1.501299 266.082602   4.000000   1.000000




##################################################









