#######File Objective##########
1. Fitting a SETAR Model to individual PCs.
2. If there is a trend in the data detrending it using loess. loess
3. Trying to simulate the data and checking if the relevant metrics are recovered. 


###Getting to the directory####
setwd("~/Correlated Risk Analysis/Decadal Influences/PC Predictions")



###Loading the requiste directories####
library(dplyr)
library(corrplot)
library(stargazer)

####Reading the Data####
input_data <- read.table("data/Max_Annual_Streamflow.txt", sep="", header = TRUE)
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)


####Principal Component Analysis####
data_pca <- prcomp(input_data, scale = TRUE)
var <- cumsum(data_pca$sdev^2)
#pdf(file = "plots/Climate_Indices/Variance by PCs.pdf")
plot(var/max(var),pch=19, main = "Variance explained by PCs",xlab = "PC's",ylab="Fraction Variance explained")
npcs <- 3 #This is the selected number of PC'S
abline(h = var[npcs]/max(var), lty = 2, col='red')
#dev.off()
pcs_sel <- data_pca$x[,1:npcs]


####PC Visualization####
#pdf(file = "plots/Climate_Indices/PC Visualizaion.pdf")
for(i in 1:npcs) {
  x <- pcs_sel[,i]
  #Plotting the Data
  plot(pcs_sel[,i], main = paste0("PC - ", i), type = 'l', 
       xlab = "Years", ylab = " ")
  
  #PACF and ACF's
  par(mfrow=c(2,1), mar=c(2,4,4,1))
  acf(x, main = paste0("PC -",i ))
  par(mar = c(2,4,2,1))
  pacf(x, main = paste0("PC -",i )) 
  par(mfrow=c(1,1), mar=c(4,2,2,2))
}
#dev.off()

#######Climate Indices############
#Reading in the climate indices. 
nino_34 <- read.table("data/iersst_nino3.4a.dat.txt", header = TRUE, sep ="", dec=".")
pdo <- read.table("data/ipdo_a.txt", header = TRUE, sep ="", dec=".")
nao <- read.table("data/inao_a.txt", header = TRUE, sep ="", dec=".")
amo <- read.table("data/iamo_hadsst_ts_a.txt", header = TRUE, sep ="", dec=".")

#Renaming. 
colnames(pdo) <- c("Time","PDO_Index")
colnames(nao) <- c("Time","NAO_Index")
colnames(nino_34) <- c("Time","Nino_Index")
colnames(amo) <- c("Time","AMO_Index")
Months <- c("Jan","Feb","Mar","April","May","June", "Jul","Aug","Sept","Oct","Nov","Dec")
target <- c("Feb", "Mar", "April","May")
years <- c(1936, 2018)

#Getting the combined PDO-NAO Index.
nao_temp <- nao[23:dim(nao)[1],] #Starting from 1924. 
nao_temp$Month <- rep(Months,195)
nao_temp$Year <- rep(1824:2018,each= 12)
nao_temp <- nao_temp %>% filter(Month %in% target)
nao_temp <- nao_temp %>% group_by(Year) %>% summarise(NAO_Index = mean(NAO_Index))
nao_temp$NAO_detrended <- nao_temp$NAO_Index - mean(nao_temp$NAO_Index)
nao_temp <- nao_temp %>% filter(
  Year > years[1],
  Year < years[2]
)


pdo_temp <- pdo[12:dim(pdo)[1],]
pdo_temp$Month <- rep(Months,117)
pdo_temp$Year <- rep(1901:2017,each= 12)
pdo_temp <- pdo_temp %>% filter(Month %in% target)
pdo_temp <- pdo_temp %>% group_by(Year) %>% summarise(PDO_Index = mean(PDO_Index))
pdo_temp$PDO_detrended <- pdo_temp$PDO_Index - mean(pdo_temp$PDO_Index)
pdo_temp <- pdo_temp %>% filter(
  Year > years[1],
  Year < years[2]
)

nino_temp <- nino_34[12:1967,] #Start from 1855
nino_temp$Month <- rep(Months,163)
nino_temp$Year <- rep(1855:2017,each= 12)
nino_temp <- nino_temp %>% filter(Month %in% target)
nino_temp <- nino_temp %>% group_by(Year) %>% summarise(ENSO_Index = mean(Nino_Index))
nino_temp$ENSO_detrended <- nino_temp$ENSO_Index - mean(nino_temp$ENSO_Index)
nino_temp <- nino_temp %>% filter(
  Year > years[1],
  Year < years[2]
)



#Combining them together. 
climate_indices <- cbind(nino_temp,pdo_temp[,2:3],nao_temp[,2:3])
climate_indices$PDO_NAO <- climate_indices$NAO_detrended*climate_indices$PDO_Index
climate_indices$ENSO_NAO <- climate_indices$ENSO_detrended*climate_indices$NAO_detrended
climate_indices$ENSO_PDO <- climate_indices$ENSO_detrended*climate_indices$PDO_Index
climate_indices$PDO_detrended <- NULL
climate_indices$NAO_detrended <- NULL
climate_indices$Year <- NULL
climate_indices$ENSO_detrended <- NULL
colnames(climate_indices) <- c("ENSO","PDO","NAO","PDO_NAO","ENSO_NAO","ENSO_PDO")

climate_indices$ENSO_NAO <- scale(climate_indices$ENSO_NAO)
climate_indices$PDO_NAO <- scale(climate_indices$PDO_NAO)

#########Plotting the Correlations########
#pdf(file = "plots/Climate_Indices/PC Correlation Matrix.pdf")
for(i in 1:npcs) {
  temp <-  cbind(pcs_sel[,i],climate_indices)
  colnames(temp)[1] <- paste0("PC_",i)
  M<-cor(temp)
  corrplot(M, type="upper")
}
#dev.off()

########Regression Predictions#########
#pdf(file = "plots/Climate_Indices/PC Predictions.pdf")
for(i in 1:npcs) {
  x <- pcs_sel[,i]
  n_ahead <- 9
  clim_data <- climate_indices[1:(nrow(climate_indices)-n_ahead),]
  train <- head(x,-n_ahead)
  test <- tail(x,n_ahead)
  training_set <- as.data.frame(cbind(train,clim_data))
  yr1 <- 1937;yrn <- yr1+length(x)-1;yr2 <- yrn-n_ahead;yr3 <- yr2+1
  mod.lm <- lm(train~.,training_set)
  #print(summary(mod.lm))
  par(mfrow=c(2,2), mar = c(1,1,3,1))
  #plot(mod.lm)
  
  act_clim_data <- climate_indices[(nrow(climate_indices)-n_ahead+1):nrow(climate_indices),]
  pred <- predict(mod.lm,act_clim_data,n_ahead)
  
  par(mfrow=c(1,1), mar = c(3,3,3,2))
  yrs_start <- 50
  yrs <- (1936+yrs_start):2017
  test_yrs <- tail(yrs, n_ahead)
  plot(yrs, x[yrs_start:length(x)], type='l',
       main = paste0("Regression with Climate Indices for PC ",i))
  lines(test_yrs, tail(x,n_ahead), col = 'blue', lwd = 2.5)
  lines(test_yrs, pred$fit, col ='red', lwd = 2.5)
  legend("bottomleft", legend = c("Testing Data", "Training Data", "Predictions"),
         lwd = c(2.5,1,2.5), col = c("blue","black","red"), lty = 1)
  
  print(paste0("The model AIC is for PC_",i, " is ", AIC(mod.lm)))
}
#dev.off()

#Predicting just based on NAO and ENSO_PDO. 
#pdf(file = "plots/Climate_Indices/PC Prediction NAO ENSO_PDO.pdf")
for(i in 1:npcs) {
  x <- pcs_sel[,i]
  n_ahead <- 9
  clim_data <- climate_indices[1:(nrow(climate_indices)-n_ahead),]
  train <- head(x,-n_ahead)
  test <- tail(x,n_ahead)
  training_data <- data.frame(train = train, NAO = clim_data$NAO, ENSO_PDO = clim_data$ENSO_PDO)
  yr1 <- 1937;yrn <- yr1+length(x)-1;yr2 <- yrn-n_ahead;yr3 <- yr2+1
  mod.lm <- lm(train~., training_data)
  print(summary(mod.lm))
  par(mfrow=c(2,2), mar = c(1,1,3,1))
  #plot(mod.lm)
  
  clim_data <- climate_indices[(nrow(climate_indices)-n_ahead+1):nrow(climate_indices),]
  act_clim_data <- data.frame(NAO = clim_data$NAO, ENSO_PDO = clim_data$ENSO_PDO)
  pred <- predict(mod.lm,act_clim_data,n_ahead)
  
  par(mfrow=c(1,1), mar = c(3,3,3,2))
  yrs_start <- 50
  yrs <- (1936+yrs_start):2017
  test_yrs <- tail(yrs, n_ahead)
  plot(yrs, x[yrs_start:length(x)], type='l',
       main = paste0("Regression with Climate Indices for PC ",i))
  lines(test_yrs, tail(x,n_ahead), col = 'blue', lwd = 2.5)
  lines(test_yrs, pred$fit, col ='red', lwd = 2.5)
  legend("bottomleft", legend = c("Testing Data", "Training Data", "Predictions"),
         lwd = c(2.5,1,2.5), col = c("blue","black","red"), lty = 1)
  
  }
#dev.off()

#Predicting just based on just NAO. 
#pdf(file = "plots/Climate_Indices/PC Predictions NAO.pdf")
for(i in 1:npcs) {
  x <- pcs_sel[,i]
  n_ahead <- 9
  clim_data <- climate_indices[1:(nrow(climate_indices)-n_ahead),]
  train <- head(x,-n_ahead)
  test <- tail(x,n_ahead)
  training_data <- data.frame(train = train, NAO = clim_data$NAO)
  yr1 <- 1937;yrn <- yr1+length(x)-1;yr2 <- yrn-n_ahead;yr3 <- yr2+1
  mod.lm <- lm(train~., training_data)
  print(summary(mod.lm))
  par(mfrow=c(2,2), mar = c(1,1,3,1))
  plot(mod.lm)
  
  clim_data <- climate_indices[(nrow(climate_indices)-n_ahead+1):nrow(climate_indices),]
  act_clim_data <- data.frame(NAO = clim_data$NAO)
  pred <- predict(mod.lm,act_clim_data,n_ahead)
  
  par(mfrow=c(1,1), mar = c(3,3,3,2))
  yrs_start <- 50
  yrs <- (1936+yrs_start):2017
  test_yrs <- tail(yrs, n_ahead)
  plot(yrs, x[yrs_start:length(x)], type='l',
       main = paste0("Regression with Climate Indices for PC ",i))
  lines(test_yrs, tail(x,n_ahead), col = 'blue', lwd = 2.5)
  lines(test_yrs, pred$fit, col ='red', lwd = 2.5)
  legend("bottomleft", legend = c("Testing Data", "Training Data", "Predictions"),
         lwd = c(2.5,1,2.5), col = c("blue","black","red"), lty = 1)
}
#dev.off()

######Making Predictions on the real space based on the PCs#########

#Getting the regression table output. 
x <- pcs_sel[,1]
n_ahead <- 9
clim_data <- climate_indices[1:(nrow(climate_indices)-n_ahead),]
train <- head(x,-n_ahead)
reg_data <- data.frame(PC1 = train, clim_data) 
mod.lm1 <- lm(PC1~.,reg_data)
mod.lm2 <- lm(PC1~NAO+ENSO_PDO,reg_data)
mod.lm3 <- lm(PC1~NAO, reg_data)
stargazer(mod.lm1,mod.lm2,mod.lm3, title="Results", align=TRUE, type = "text", out='.txt')

#Selecting the 1st PC
x <- pcs_sel[,1]
n_ahead <- 9
clim_data <- climate_indices[1:(nrow(climate_indices)-n_ahead),]
train <- head(x,-n_ahead)
test <- tail(x,n_ahead)
training_data <- data.frame(train = train, NAO = clim_data$NAO, ENSO_PDO = clim_data$ENSO_PDO)
mod.lm <- lm(train~., training_data)
clim_data <- climate_indices[(nrow(climate_indices)-n_ahead+1):nrow(climate_indices),]
act_clim_data <- data.frame(NAO = clim_data$NAO, ENSO_PDO = clim_data$ENSO_PDO)
pred <- predict(mod.lm,act_clim_data,n_ahead)
yr1 <- 1937;yrn <- yr1+length(x)-1;yr2 <- yrn-n_ahead;yr3 <- yr2+1
plot(yr1:yrn, x, type='l',
     main = paste0("PC ",1))
lines(yr3:yrn, pred$fit, lwd = 2, col ='red')






#PC Loadings
library("biwavelet")
library("plotrix")
library("maps")
loadings <- data_pca$rotation 
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
ju <- abs(loadings[,1])
#pdf(file = "plots/Climate_Indices/Spatial Distribution of PCs.pdf")
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"), boundary = TRUE)
points(site_info$dec_long_va,site_info$dec_lat_va,pch=19,cex=1,col=color.scale(ju,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
title(paste0("Spatial Distribution of PC 1"))
legend("bottomright", legend = c("high","low"), col = c("blue","red"), cex =0.6, pch =19)
#dev.off()

#Reconstructing the streamflow field
PC_Predictions <- pred$fit #These are the predictions
nComp = 1
Predictions_Scaled =  PC_Predictions %*% t(data_pca$rotation[,1])
for(i in 1:ncol(Predictions_Scaled)) {   
  Predictions_Scaled[,i] <- scale(Predictions_Scaled[,i], center = FALSE , scale=1/data_pca$scale[i]) }

for(i in 1:ncol(Predictions_Scaled)) {   
  Predictions_Scaled[,i] <- scale(Predictions_Scaled[,i], center = -1 * data_pca$center[i], scale=FALSE)
}
Predictions_Scaled <- exp(Predictions_Scaled)

#Computing Site - Specific MSE
True_Values <- exp(tail(input_data,n_ahead))
site_MSE <- bias <- rep(NA,ncol(Predictions_Scaled))
for(i in 1:ncol(Predictions_Scaled)) {
  site_MSE[i] <- mean((Predictions_Scaled[,i]-True_Values[,i])^2)
  bias[i] <- sign(sum(Predictions_Scaled[,i]-True_Values[,i]))
}


#Plotting the Data
for(i in 1:ncol(Predictions_Scaled)) {
  if(bias[i] == -1) {bias[i] = 19 #Triangle
  } else {bias[i] = 17 } #Over Predictions 
}

#pdf(file = "plots/Climate_Indices/Bias in prediction.pdf")
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(site_info$dec_long_va,site_info$dec_lat_va,
       col=color.scale(log(site_info$drain_area_va),c(1,0.5,0),c(0,0,1),color.spec="rgb"),
       cex=1,
       pch=bias)
title("Error in Streamflow Prediction")
legend("bottomright", c("Color - Drainage Area","Size - Error in MSE adjusted","Shape - Bias"), cex = 0.6)
#dev.off()

########3###Comparision against base mark prediction Skill Testing################
#1. Long Term Mean. 
#2. Regression directly to the dataset. 




#########Long Term Mean####################
#pdf(file = "plots/Climate_Indices/Accuracy against LTM.pdf")
All_Predictions <- Predictions_Scaled
training_set <- exp(head(input_data,-n_ahead))
testing_set <- True_Values
ltm <- matrix(colMeans(training_set), nrow = 1, ncol = ncol(training_set))
ltm_skill <- matrix(NA, nrow = n_ahead, ncol = ncol(testing_set))
for(i in 1:ncol(ltm_skill)) {
  for(j in 1:nrow(ltm_skill)) {
    if(abs(ltm[i]-testing_set[j,i]) > abs(All_Predictions[j,i]-testing_set[j,i])) { ltm_skill[j,i] = 1
    } else { ltm_skill[j,i] = 0 
    } 
  }
}

par(mfrow=c(1,1))
for(i in 1:nrow(ltm_skill)) {
  skill <- ltm_skill[i,]
  for(j in 1:length(skill)) { 
    if(skill[j]==0) {skill[j] = c("red")
    } else { skill[j] = c("blue")
    }
  }
  
  map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
  points(site_info$dec_long_va,site_info$dec_lat_va,
         pch=19,
         cex=1,
         col=skill)
  title(paste0("Skill Testing vs Mean for Year ", i))
  legend("bottomright", c("Correct Prediction", "Wrong Prediction"), cex = 0.6, pch = 19, col = c("blue","red"))
}

skill <- colSums(ltm_skill)
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(site_info$dec_long_va,site_info$dec_lat_va,
       pch=19,
       cex=skill*1.5/n_ahead,
       col=color.scale(site_info$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"))
title("Combined Skill vs Long Term Mean")
legend("bottomright", c("Color - Drainage Area","Size - Skill"), cex = 0.6)
par(mar = c(4,4,4,1))
plot(1:n_ahead, rowMeans(ltm_skill),type='l',
     xlab = "Year Ahead",
     ylab = "Accurate Predictions", 
     main = "Prediction Accuracy against LTM")
legend("bottomleft", legend = c(paste0("Mean Acc is ", round(sum(ltm_skill)/(n_ahead*ncol(input_data)),2))),
       cex = 0.75)

#dev.off()

#########Regression to individual Stations#######
#pdf(file = "plots/Climate_Indices/Accuracy against Raw Regression.pdf")
reg_raw <- matrix(NA, nrow = n_ahead, ncol = ncol(testing_set))
for(i in 1:ncol(reg_raw)) {
  temp <- exp(input_data[,i])
  temp_test <- head(temp, - n_ahead)
  clim_data <- climate_indices[1:(nrow(climate_indices)-n_ahead),]
  temp_training <- data.frame(train = temp_test, NAO = clim_data$NAO, ENSO_PDO = clim_data$ENSO_PDO)
  mod.lm <- lm(train~., temp_training)
  
  clim_data <- climate_indices[(nrow(climate_indices)-n_ahead+1):nrow(climate_indices),]
  act_clim_data <- data.frame(NAO = clim_data$NAO, ENSO_PDO = clim_data$ENSO_PDO)
  pred <- predict(mod.lm,act_clim_data,n_ahead)
  reg_raw[,i] <- pred$fit 
  
}

reg_skill <- matrix(NA, nrow = n_ahead, ncol = ncol(reg_raw))
for(i in 1:ncol(reg_skill)) {
  for(j in 1:nrow(reg_skill)) {
    if(abs(reg_raw[j,i]-testing_set[j,i]) > abs(All_Predictions[j,i]-testing_set[j,i])) { reg_skill[j,i] = 1
    } else { reg_skill[j,i] = 0 
    } 
  }
}

par(mfrow=c(1,1))
for(i in 1:nrow(reg_raw)) {
  skill <- reg_skill[i,]
  for(j in 1:length(skill)) { 
    if(skill[j]==0) {skill[j] = c("red")
    } else { skill[j] = c("blue")
    }
  }
  
  map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
  points(site_info$dec_long_va,site_info$dec_lat_va,
         pch=19,
         cex=1,
         col=skill)
  title(paste0("Skill Testing vs Raw Regression for Year ", i))
  legend("bottomright", c("Correct Prediction", "Wrong Prediction"), cex = 0.6, pch = 19, col = c("blue","red"))
}

skill <- colSums(reg_skill)
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(site_info$dec_long_va,site_info$dec_lat_va,
       pch=19,
       cex=skill*1.5/n_ahead,
       col=color.scale(site_info$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"))
title("Combined Skill vs Raw Regression")
legend("bottomright", c("Color - Drainage Area","Size - Skill"), cex = 0.6)
par(mar = c(4,4,4,1))
plot(1:n_ahead, rowMeans(reg_skill),type='l',
     xlab = "Year Ahead",
     ylab = "Accurate Predictions", 
     main = "Prediction Accuracy against Raw Regressions")
legend("bottomleft", legend = c(paste0("Mean Acc is ", round(sum(reg_skill)/(n_ahead*ncol(input_data)),2))),
       cex = 0.75)
#dev.off()
############################################################


