#####
#1. The objective of this file is to make predictions in PC space. 
#2. Use a dimension lowering algorithm - linear PCs
#3. Make predictions on those PCs using locfit i.e local regression. 
#4. Check if the relevant indices and the spectral structure has been replicated. 


#########Settting the Path
setwd("~/Correlated Risk Analysis/Decadal Influences/PC Predictions")


###Loading Libraries##########
library(ggplot2)
library(locfit)
library(sp)


###Reading the Data and Initial Manipulation###
input_data <- read.table("data/Max_Annual_Streamflow.txt", sep="", header = TRUE)
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)
######################################2


####Principal Component Analysis#############
data_pca <- prcomp(input_data, scale = TRUE)
var <- cumsum(data_pca$sdev^2)
#pdf(file = "plots/SETAR/Variance explained.pdf")
plot(var/max(var),pch=19, main = "Variance explained by PCs",xlab = "PC's",ylab="Fraction Variance explained")
npcs <- 3 #This is the selected number of PC'S
abline(h = var[npcs]/max(var), lty = 2, col='red')
#dev.off()
pcs_sel <- data_pca$x[,1:npcs]
#################################################2


#####Locfit for each PC##########
pred_horizon <- 9 #Years This is a decent forecast horizon.

for(jk in 1:ncol(pcs_sel)) {
  pc <- scale(pcs_sel[,jk])
  lag_loc <- c(1, 2, 3) #These are the lag terms. 
  pc = as.matrix(pc)
  feature_vectors <- matrix(NA, ncol = length(lag_loc)+1, nrow = length(pc)-max(lag_loc))
  for(i in 1:(ncol(feature_vectors)-2)) {
    feature_vectors[,i+1] <-  tail(head(pc,-lag_loc[i]),-max(lag_loc)+lag_loc[i])
  }
  feature_vectors[,ncol(feature_vectors)] <- head(pc,-max(lag_loc))
  feature_vectors[,1] <- tail(pc,-max(lag_loc))
  feature_test <- tail(feature_vectors, pred_horizon)
  feature_train <- head(feature_vectors, -pred_horizon)

  x_l <- feature_train[,2:ncol(feature_train)]
  y_l <- feature_train[,1]
  
  #GCV Plots
  gcv_mat <- gcvplot(x_l,y_l,alpha = seq(0.5,0.9,by=0.05),maxk = 1000)
  plot(gcv_mat$alpha, gcv_mat$values, 
       main = paste0("GCV plot for PC_",jk),
       xlab = "alpha", ylab = 'GCV Score')
  best_alpha <- gcv_mat$alpha[which.min(gcv_mat$values)]
  fit <- locfit.raw(x_l,y_l,maxk = 1000, alpha = best_alpha)
  #plot(fit,get.data=TRUE)
  pred <- predict(fit, feature_test[,2:ncol(feature_test)], se = T, interval = "prediction")

#Plotting the Predictions 
st <- 50
plot(st:length(pc), pc[st:length(pc)], type='l',
     main = paste0("PC_",jk))
lines(tail(1:length(pc),pred_horizon), tail(pc,pred_horizon), lwd = 2)
lines(tail(1:length(pc),pred_horizon), pred$fit, lwd = 2,col='red')
lines(tail(1:length(pc),pred_horizon), pred$fit-2*pred$se.fit, lwd = 2,col='red',lty = 2)
lines(tail(1:length(pc),pred_horizon), pred$fit+2*pred$se.fit, lwd = 2,col='red',lty = 2)
legend('bottomleft', legend = c("True Data","Predicted Data"), 
      col = c('black','red'), lwd = 2, lty = 1)
}

