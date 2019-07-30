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
#pdf
plot(var/max(var),pch=19, main = "Variance explained by PCs",xlab = "PC's",ylab="Fraction Variance explained")
npcs <- 3 #This is the selected number of PC'S
abline(h = var[npcs]/max(var), lty = 2, col='red')
#dev.off()

pcs_sel <- data_pca$x[,1:npcs]
#write.table(pcs_sel, "results/Selected PCs.txt", sep = " ")
#################################################2


###############Diagnostics##########
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
######################SETAR Model################
#SETAR is Self Exciting Threshold AutoRegressive Model 
#Reading is Howell Tong - Non-Linear Time Series. 
for(i in 1:ncol(pcs_sel)){
x <- pcs_sel[,i]

#Simple AR Model
mod.ar <- ar(x, aic = TRUE, order.max = 10)
ar_order <- mod.ar$order
resd <- mod.ar$resid[!is.na(mod.ar$resid)]
#acf(resd)
#pacf(resd)
#plot(forecast(mod.ar, h = 20))

#SETAR
if(ar_order>-1){ 
  ar_order = 1}
mod.setar <- setar(x, m=ar_order, d = 1, steps = 1,nthresh = 1)
summary(mod.setar)
#plot(mod.setar)
n_ahead <- 20
predict.setar <-  predict(mod.setar, n.ahead=n_ahead, type="bootstrap", n.boot=200)
plot(1:81, x, type ='l', 
     xlim = c(0,n_ahead+85), 
     main = paste0("Embedding Order is ", ar_order))
lines(82:(81+n_ahead), predict.setar$pred, col ='red')
lines(82:(81+n_ahead), predict.setar$se[,2], col ='red',lty = 2)
lines(82:(81+n_ahead), predict.setar$se[,1], col ='red',lty = 2)


#SETAR Simulations
N_Sims <- 1000
mean_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
sd_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
min_sim <- max_sim <- matrix(NA, ncol = 1, nrow = N_Sims)

B <- as.vector(head(mod.setar$coefficients,-1))
for(j in 1:N_Sims) {
sims.setar <- setar.sim(B = B,lag = ar_order,nthresh = 1, Thresh = tail(mod.setar$coefficients,1), type = "simul", n = 81)
mean_sim[j,] <- mean(sims.setar$serie)
sd_sim[j,] <- sd(sims.setar$serie)
min_sim[j,] <- min(sims.setar$serie)
max_sim[j,] <- max(sims.setar$serie)
}

par(mfrow=c(1,4))
par(mar = c(4, 2, 1.5, 1))

boxplot(mean_sim)
title(paste0("Mean-PC ", i))
abline(h=mean(x), col = 'red')

boxplot(sd_sim)
title(paste0("SD - PC ", i))
abline(h=sd(x), col = 'red')

boxplot(max_sim)
title(paste0("Max - PC ", i))
abline(h=max(x), col = 'red')

boxplot(min_sim)
title(paste0("Max - PC ", i))
abline(h=min(x), col = 'red')
par(mfrow=c(1,1))

}

#############GARCH Modelling#############














#Multiple Models Comparision
mod <- list()
max_delay <- 7
max_embedding <- 6
for(del in 1:max_delay) {
  for(embd in 1:max_embedding) {
  mod_name <- paste0("Del_",del,"_Embed_Dim_",embd)
mod[[mod_name]] <- setar(x, m=embd, d = del, nthresh = 2)
  }
}  

AICs <- matrix(sapply(mod, AIC), ncol = max_delay, nrow = max_embedding)
par(mfrow=c(1,1), mar = c(4,4,3,1))
plot(1:max_delay, AICs[1,], lty = 1, type = 'l', 
     ylim = range(AICs)-5, ylab = "AIC",
     xlab = "Delay(d)", col = 1, 
     main = "AIC values for SETAR Models with m and d")
for(embd in 2:max_embedding){
  lines(1:max_delay, AICs[embd,],lty = embd, col = embd)
}
legend("bottomleft", legend = c(1:max_embedding), lty = c(1:max_embedding), col = c(1:max_embedding))
#egend("bottomright", lty=1:(length(frc.test)+1), col=1:(length(frc.test)+1),legend=c("observed",names(frc.test)), cex = 0.6) }



sapply(mod, AIC)
sapply(mod, MAPE)

#Out of Sample Forecasting
set.seed(10)
mod.test <- list()
x.train <- window(x, end=70)
x.test <- window(x, start=71)
mod.test[["linear"]] <- linear(x.train, m=ar_order)
mod.test[["setar"]] <- setar(x.train, m=ar_order, thDelay=1, d = 3)
mod.test[["lstar"]] <- lstar(x.train, m=ar_order, thDelay=1, trace=FALSE, control=list(maxit=1e5))
mod.test[["nnet"]] <- nnetTs(x.train, m=ar_order, size=3, control=list(maxit=1e5))
mod.test[["aar"]] <- aar(x.train, m=ar_order)
sapply(mod.test, AIC)
sapply(mod.test, MAPE)
frc.test <- lapply(mod.test, predict, n.ahead=10)
plot(x.test,ylim=range(x-2), type ='l')
for(i in 1:length(frc.test)) {
lines(1:10, frc.test[[i]], lty=i+1, col=i+1) }
legend("bottomright", lty=1:(length(frc.test)+1), col=1:(length(frc.test)+1),legend=c("observed",names(frc.test)), cex = 0.6) }


