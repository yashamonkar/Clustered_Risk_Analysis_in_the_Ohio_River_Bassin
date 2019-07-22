#Files uses PC's to reduce dimensionality
#The wavelets are used to get information on the non-stationary PCs. 
#The individual signals are reconstructed using the AR Models. 
#Simulate this to get the structure again. 
setwd("~/Correlated Risk Analysis/Decadal Influences/PC Wavelet AR Models")



####Loading the libraries#####
library(dataRetrieval)
library(dplyr)
library(dams)
library(NbClust)
library(factoextra)
library(maps)
library(corrplot)
library(extRemes)
library(xts)
library(cluster)
library(foreign)
library(nnet)
library(forecast)
##################################3


################Reading the Data##############
input_data <- read.table("data/Max_Annual_Streamflow.txt", sep="", header = TRUE)
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)
######################################


###########Wavelet Function######## 
wavelet=function(Y,dj=0.025){
  
  #Y is time series to be analyzed
  DT=1# is timestep for annual data, 1
  pad=1
  #dj=0.025
  param=6
  #pad data ? 0=F, 1=T
  #dj= spacing between discrete scales (.025)
  #param = wavenumber (6)
  
  s0=2*DT
  
  n1 = length(Y)
  J1=floor((log2(n1*DT/s0))/dj)
  
  
  #....construct time series to analyze, pad if necessary
  x = Y - mean(Y)
  
  
  if (pad == 1){
    base2 = trunc(log(n1)/log(2) + 0.4999)   # power of 2 nearest to N
    x = c(x, rep(0, 2^(base2 + 1) - n1))
  }
  n = length(x)
  
  #....construct wavenumber array used in transform [Eqn(5)]
  k = (1:trunc(n/2))
  k = k*((2*pi)/(n*DT))
  k = c(0, k, -rev(k[1:floor((n-1)/2)]))
  
  #....compute FFT of the (padded) time series
  f = fft(x)    # [Eqn(3)]
  
  #....construct SCALE array & empty PERIOD & WAVE arrays
  scale = s0*2^((0:J1)*dj)
  period = scale;
  wave = matrix(data=0, ncol=n, nrow=J1+1)  # define the wavelet array
  wave = as.complex(wave)  # make it complex
  wave=matrix(data=wave, ncol=n, nrow=J1+1)
  
  # loop through all scales and compute transform
  for(a1 in 1:(J1+1)){
    scl=scale[a1]		
    
    nn = length(k);
    k0 = param
    expnt = -(scl*k - k0)^2/(2*(k > 0))
    norm = sqrt(scl*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
    daughter = norm*exp(expnt)
    daughter = daughter*(k > 0)    # Heaviside step function
    fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
    coi = fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
    dofmin = 2                                   # Degrees of freedom
    
    out <- list(daughter=daughter, fourier_factor=fourier_factor,coi=coi,dofmin=dofmin)
    
    daughter=out$daughter
    fourier_factor=out$fourier_factor
    coi=out$coi
    dofmin=out$dofmin	
    wave[a1,] = fft((f*daughter), inverse = TRUE)/(length(f*daughter))  # wavelet transform[Eqn(4)]
  }
  
  period = fourier_factor*scale
  
  coi = coi*c(1:(floor(n1 + 1)/2), rev(1:floor(n1/2))) * DT
  
  wave = wave[,1:n1]  # get rid of padding before returning
  power=abs(wave)^2
  ncol=length(power[1,])
  nrow=length(scale)
  avg.power=apply(power,1,mean)
  result=list(wave=wave, period=period, scale=scale, power=power, coi=coi,nc=ncol,nr=nrow,p.avg=avg.power)
  return(result)
}
CI=function(conf, dat,type){
  
  #### enter confidence as decimal 0-1
  #### two types of tests available 1) red noise enter: "r" , white noise enter: "w"
  # requires the wavelet function
  
  na=length(dat)
  wlt=wavelet(dat)
  
  if(type=="r"){
    
    zz=arima(dat/10^10, order = c(1, 0, 0))
    alpha=zz$coef[1]
    print(alpha)
  } else{
    alpha=0
  }
  
  ps=wlt$period
  LP= length(ps)
  
  freq = 1/ps
  
  CI=1:LP    ## confidence interval
  
  for(i in 1:LP){
    
    P=(1-(alpha^2))/(1+(alpha^2)-(2*alpha*cos(2*pi*freq[i])))    # descrete fourier power spectrum page 69 [qn 16] ( torrence and compo)... null hypothesis test
    df=2*sqrt(1+((na/(2.32*ps[i]))^2))
    CI[i] =P*(qchisq(conf, df)/df)*var(dat)          #divide by 2 removes power of 2.....for mean no chi dist.[ eqn 17]
  }
  
  
  list(sig=CI)
  
}
#These include functions to compute the wavelets and the confidence intervals on the wavelets. 
#############END Wavelet Function###################




############################PCA######################
#We have 81 years and 30 stations
#This makes the dimensionality of the problem too large to handle. 
#Therefore we use PCA to reduce the dimensionality

#Keeping last 5 years for predictions. 
predict_ahead <- 0
yr1 <- 1;yr2 <- dim(input_data)[1]-predict_ahead;yr3 <- yr2+1;yr4 <- dim(input_data)[1]


input_data_pca <- prcomp(input_data, scale = TRUE)
var <- cumsum(input_data_pca$sdev^2)
pdf(file='plots/PC Variance Explained.pdf')
plot(var/max(var),pch=19, main = "Variance explained by PCs",xlab = "PC's",ylab="Fraction Variance explained")
dev.off()
pdf(file = 'plots/Selected PCs.pdf')
npcs <- 5 #NUmber of Selected PCs
par(mfrow=c(3,2))
par(mar = c(4, 1, 1, 1))
for(i in 1:npcs){plot(input_data_pca$x[,i],typ='l', ylab = NA, xlab = paste0("PC ",i))
}
dev.off()
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))

#Plotting the De-trended PC's
pdf("plots/Detrended PCs.pdf")
par(mfrow=c(2,2));par(mar = c(4, 3, 3, 1))
detrended_PC <- matrix(NA, ncol = npcs, nrow = nrow(input_data))
trend_coeff <- matrix(NA, ncol = 2, nrow = npcs)
for(i in 1:npcs) {
  p=input_data_pca$x[,i]
  fit <- lm(p~c(yr1:yr2))
  
  #Saving the data.
  trend <- c(yr1:yr2)*fit$coefficients[2]+fit$coefficients[1]
  detrended_PC[,i] <- p-trend
  trend_coeff[i,] <- fit$coefficients
  
  plot(yr1:yr2,p, main = paste0("PC ", i),type='l',xlab = "Year")
  abline(fit,col='blue')
  lines(lowess(yr1:yr2,p,f=1/9),lwd=2,col="red")
  detrended_PC[,i] <- p-trend
  plot(yr1:yr2,detrended_PC[,i], main = paste0("Detrended PC ", i),type='l',xlab = "Year")
  lines(lowess(yr1:yr2,detrended_PC[,i],f=1/9),lwd=2,col="red")
  
  
}
par(mfrow=c(1,1))
dev.off()



prin_comp  <- detrended_PC
write.table(prin_comp, "results/Selected PCs.txt", sep = " ")
#####################################

##################Wavelet Analysis###############
# We have reduced the dimensionality of the data using PCA. 
#Now we try to uncover further spatial aspects using Wavelet Analysis.
# We first fit a wavelet to each of the PC's 

#Step 1:- Getting the significant periods (over red noise). 
library("biwavelet")
library("plotrix")
library("maps")
library("stats") 
pdf(file = 'plots/PC Wavelet Spectrums.pdf')
sig_scales <- as.list(1)
par(mfrow=c(3,2))
par(mar = c(4, 1, 1.5, 1))
for(i in 1:ncol(prin_comp)) {
  p <- prin_comp[,i]
  wlt <- wavelet(p)
  C_r <- CI(0.9,p,'r')
  C_w <- CI(0.9,p,'w')
  plot(wlt$period,wlt$p.avg,xlim=c(0,75),main=paste0("Global Wavelet Spectrum for PC ", i),xlab="Period",ylab="Variance")
  lines(wlt$period,wlt$p.avg)
  lines(wlt$period,C_r$sig,col='red')
  lines(wlt$period,C_w$sig, col ='black')
  sig_scales[[i]] <- which(wlt$p.avg > C_r$sig)
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
dev.off()

#Step 2:- Reconstructing each PC. 
pdf(file = 'plots/PC Reconstructions.pdf')
par(mfrow=c(3,2))
par(mar = c(4, 1, 1.5, 1))
for(i in 1:ncol(prin_comp)) {
  p <- prin_comp[,i]
  wlt <- wavelet(p)
  Cd <- 0.776;psi0 <- pi^(-.025);dj=0.025 #From the Torrence and Compo
  reconst <- matrix(NA, ncol = length(p), nrow = length(wlt$scale))
  for(j in 1:ncol(reconst)) {reconst[,j] <- Re(wlt$wave[,j])/(wlt$scale[i]^0.5)}
  p_reconst <- colSums(reconst)*dj/(Cd*psi0) 
  plot(yr1:yr2, p_reconst, type='l',col='red', main = paste0("Reconstruced PC ", i),
       xlab = "Year", ylab = "Streamflow(Scaled)")
  lines(yr1:yr2,p)
  legend('topright', legend = c("Signal","Constructed"), cex = 0.6, col = c('black','red'),lty = 1)
  print(paste0("The correlation between signal and reconstructed TS for PC ", i, " is ", round(cor(p_reconst,p),2)))
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
dev.off()
#These look good now we can seperate the significant signals and the noise


#Step 3:- Decomposing each PC into signals and noise
#For this we will have to divide the seperate regions for clustering. 
pdf(file = 'plots/PCs - Signals and Noise.pdf')
par(mfrow=c(3,2))
par(mar = c(4, 1, 1.5, 1))
signals <- list(NA)
for(i in 1:ncol(prin_comp)) {
  
  #Fitting the Wavelet
  p <- prin_comp[,i]
  wlt <- wavelet(p)
  Cd <- 0.776;psi0 <- pi^(-.025);dj=0.025 #From the Torrence and Compo
  reconst <- matrix(NA, ncol = length(p), nrow = length(wlt$scale))
  for(j in 1:ncol(reconst)) {reconst[,j] <- Re(wlt$wave[,j])/(wlt$scale[i]^0.5)}
  
  
  #Getting the significant portions
  temp <- sig_scales[[i]]
  temp_list <- list(NA)
  n_signals <- length(which(diff(temp)>1))+1 #Number of seperate signals
  if (n_signals > 1) { st <- 1
  for(j in 1:n_signals) {
    en <- c(which(diff(temp)>1),length(temp))
    en <- en[j]
    temp_list[[j]] <- temp[st:en] 
    st <- en+1
  } 
  }
  if(n_signals < 2) { temp_list <- list(temp)}
  n_clust <- length(rapply(temp_list, length))
  
  if(rapply(temp_list, length) == 0) {
    plot(yr1:yr2, p, type='l', main = paste0("Signal and Noise PC ", i))
    p_reconst <- colSums(reconst)*dj/(Cd*psi0) 
    lines(yr1:yr2,p_reconst, col ='blue', lty = 2)
    legend('bottomright', legend = c("PC","Signal","Noise"), cex = 0.6, col = c('black','red','blue'),lty = c(1,1,2))
    
  } else {
    #Plotting the original data. 
    plot(yr1:yr2, p, type='l', main = paste0("Signal and Noise PC ", i), xlab = "Year")
    
    #Plotting the significant terms.
    for(jk in 1:length(rapply(temp_list, length))) {
      clust_members <- temp_list[[jk]]
      t <- reconst[clust_members,]
      t_reconst <- colSums(t)*dj/(Cd*psi0)
      lines(yr1:yr2,t_reconst, col ='red')
    }
    
    #Modelling the noise terms.
    noise <- unlist(temp_list)
    t_noise <- reconst[-noise,]
    t_reconst <- colSums(t_noise)*dj/(Cd*psi0)
    lines(yr1:yr2,t_reconst, col ='blue', lty = 2)
    legend('bottomright', legend = c("PC","Signals","Noise"), cex = 0.6, col = c('black','red','blue'),lty = c(1,1,2))
    print(paste0("The correlation between signal and noise TS for PC ", i, " is ", round(cor(t_reconst,p),2)))
 }
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
dev.off()
###############################################################3



######################Fitting ARMA Models###################
pdf(file='plots/Individual AR Models.pdf')
signals <- list(NA)
tot_pc_pred <- matrix(NA,nrow=predict_ahead,ncol=ncol(prin_comp))
par(mfrow=c(2,2))
par(mar = c(4, 2, 1.5, 1))
for(i in 1:ncol(prin_comp)) {
  
  #Fitting the Wavelet
  p <- prin_comp[,i]
  wlt <- wavelet(p)
  Cd <- 0.776;psi0 <- pi^(-.025);dj=0.025 #From the Torrence and Compo
  reconst <- matrix(NA, ncol = length(p), nrow = length(wlt$scale))
  for(j in 1:ncol(reconst)) {reconst[,j] <- Re(wlt$wave[,j])/(wlt$scale[i]^0.5)}
  
  
  #Getting the significant portions
  temp <- sig_scales[[i]]
  temp_list <- list(NA)
  n_signals <- length(which(diff(temp)>1))+1 #Number of seperate signals
  if (n_signals > 1) { st <- 1
  for(j in 1:n_signals) {
    en <- c(which(diff(temp)>1),length(temp))
    en <- en[j]
    temp_list[[j]] <- temp[st:en] 
    st <- en+1
  } 
  }
  if(n_signals < 2) { temp_list <- list(temp)}
  n_clust <- length(rapply(temp_list, length))
  
  #Case with no significant regions.
  if(rapply(temp_list, length) == 0) {
    
    #Creating the matrix to store the time series.
    p_reconst <- colSums(reconst)*dj/(Cd*psi0) 
    breakdown_ts <- matrix(data = p_reconst, ncol = 1, nrow = length(p))
    predict_ts <- matrix(data=NA, ncol = 1, nrow = predict_ahead)
    
    #Fitting an ARIMA Model. 
    fit <- ar(breakdown_ts[,1], aic = TRUE, order.max = 10)
    ord <- arimaorder(fit)
    plot(forecast(fit,h=20), xlab = paste0("Reconstructed PC ", i))
    predict_ts <- forecast(fit,h=5)$mean
    tot_pc_pred[,i] <- predict_ts
    
    resd <- fit$resid[!is.na(fit$resid)]
    #Diagnostics
    #par(mfrow = c(2,2))
    #plot(density(resd), main ="Resd");acf(resd, main = "ACF");pacf(resd, main = "PACF")
    #par(mfrow=c(1,1))
    
    #Getting the predictions
    #plot(2001:2017,prin_comp[65:81,i], type='l', main = paste0("Predictions for PC ", i)
    #     ,xlab = "Year", ylab = "PC")
    #lines((2017-predict_ahead+1):2017, predict_ts,col='red')
    #legend('bottomright', legend = c("Real","predicted"), lty = 1, col =c('black','red'), cex = 0.6)
    
  } else {#Creating the matrix to store the signals. 
    breakdown_ts <- matrix(NA, ncol = length(rapply(temp_list, length))+1, nrow = length(p))
    predict_ts <- matrix(NA, ncol = length(rapply(temp_list, length))+1, nrow = predict_ahead)
    for(jk in 1:length(rapply(temp_list, length))) {
      clust_members <- temp_list[[jk]]
      t <- reconst[clust_members,]
      t_reconst <- colSums(t)*dj/(Cd*psi0)
      breakdown_ts[,jk] <- t_reconst}
    #Modelling the noise terms.
    noise <- unlist(temp_list)
    t_noise <- reconst[-noise,]
    t_reconst <- colSums(t_noise)*dj/(Cd*psi0)
    breakdown_ts[,ncol(breakdown_ts)] <- t_reconst
    
    #Fitting ARIMA Models to each component. 
    for(jks in 1:ncol(breakdown_ts)) {
      fit <- ar(breakdown_ts[,jks], order.max = 10, aic = TRUE)
      ord <- arimaorder(fit)
      nota <- paste0("Reconstructed PC ", i, "'s signal ", jks)
      if(jks == ncol(breakdown_ts)){nota <- paste0("Reconstructed PC ", i, "'s noise ")}
      plot(forecast(fit,h=20), xlab =nota)
      predict_ts[,jks] <- forecast(fit,h=5)$mean
      #resd <- fit$resid[!is.na(fit$resid)]
      #Diagnostics
      #par(mfrow = c(2,2))
      #plot(density(resd), main ="Resd");acf(resd, main = "ACF");pacf(resd, main = "PACF")
      #par(mfrow=c(1,1))
      
    }
    predict_ts <- rowSums(predict_ts)
    tot_pc_pred[,i] <- predict_ts
    #Getting the predictions
    #plot(2001:2017,prin_comp[65:81,i], type='l', main = paste0("Predictions for PC ", i)
    #    ,xlab = "Year", ylab = "PC")
    #lines((2017-predict_ahead+1):2017, predict_ts,col='red')
    #legend('bottomright', legend = c("Real","predicted"), lty = 1, col =c('black','red'), cex = 0.6)
    
    
  }
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
dev.off()

#Saving the Simulations
N_Sims <- 1000 #This is the number of simulations. 
signals <- list(NA)
tot_pc_sims <- matrix(NA,nrow=N_Sims*dim(input_data)[1],ncol=ncol(prin_comp))
for(i in 1:ncol(prin_comp)){
  
  #Fitting the Wavelet
  p <- prin_comp[,i]
  wlt <- wavelet(p)
  Cd <- 0.776;psi0 <- pi^(-.025);dj=0.025 #From the Torrence and Compo
  reconst <- matrix(NA, ncol = length(p), nrow = length(wlt$scale))
  for(j in 1:ncol(reconst)) {reconst[,j] <- Re(wlt$wave[,j])/(wlt$scale[i]^0.5)}
  
  
  #Getting the significant portions
  temp <- sig_scales[[i]]
  temp_list <- list(NA)
  n_signals <- length(which(diff(temp)>1))+1 #Number of seperate signals
  if (n_signals > 1) { st <- 1
  for(j in 1:n_signals) {
    en <- c(which(diff(temp)>1),length(temp))
    en <- en[j]
    temp_list[[j]] <- temp[st:en] 
    st <- en+1
  } 
  }
  if(n_signals < 2) { temp_list <- list(temp)}
  n_clust <- length(rapply(temp_list, length))
  
  #Case with no significant regions.
  if(rapply(temp_list, length) == 0) {
    
    #Creating the matrix to store the time series.
    p_reconst <- colSums(reconst)*dj/(Cd*psi0) 
    breakdown_ts <- matrix(data = p_reconst, ncol = 1, nrow = length(p))
    sims_ts <- matrix(data=NA, ncol = 1, nrow = dim(input_data)[1])
    
    #Fitting an ARIMA Model. 
    fit <- ar(breakdown_ts[,1], aic = TRUE, order.max = 10)
    ord <- arimaorder(fit)
    #plot(forecast(fit,h=20), xlab = paste0("Reconstructed PC ", i))
    for(jt in 1:N_Sims) {
    sims_ts <- simulate(fit,dim(input_data)[1])
    
    #Adding the trend (which was subtracted)
    trend <- trend_coeff[i,1] + trend_coeff[i,2]*c(yr1:yr4)
    low <- jt*dim(input_data)[1] - dim(input_data)[1] + 1
    high <- jt*dim(input_data)[1]
    tot_pc_sims[low:high,i] <- sims_ts + trend }
    
  } else {#Creating the matrix to store the signals. 
    breakdown_ts <- matrix(NA, ncol = length(rapply(temp_list, length))+1, nrow = length(p))
        for(jk in 1:length(rapply(temp_list, length))) {
      clust_members <- temp_list[[jk]]
      t <- reconst[clust_members,]
      t_reconst <- colSums(t)*dj/(Cd*psi0)
      breakdown_ts[,jk] <- t_reconst}
    #Modelling the noise terms.
    noise <- unlist(temp_list)
    t_noise <- reconst[-noise,]
    t_reconst <- colSums(t_noise)*dj/(Cd*psi0)
    breakdown_ts[,ncol(breakdown_ts)] <- t_reconst
    
    for(jt in 1:N_Sims) {
      sims_ts <- matrix(NA, ncol = length(rapply(temp_list, length))+1, nrow = dim(input_data)[1])
    #Fitting ARIMA Models to each component. 
    for(jks in 1:ncol(breakdown_ts)) {
      fit <- ar(breakdown_ts[,jks], order.max = 10, aic = TRUE)
      sims_ts[,jks] <- simulate(fit, dim(input_data)[1])
      
    }
    sims_ts <- rowSums(sims_ts)
    
    #Adding the trend (which was subtracted)
    trend <- trend_coeff[i,1] + trend_coeff[i,2]*c(yr1:yr4)
    low <- jt*dim(input_data)[1] - dim(input_data)[1] + 1
    high <- jt*dim(input_data)[1]
    tot_pc_sims[low:high,i] <- sims_ts + trend
  }
  
  }

}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))



####Plotting the moments
mean_sim <- matrix(NA, ncol = npcs, nrow = N_Sims)
sd_sim <- matrix(NA, ncol = npcs, nrow = N_Sims)
min_sim <- max_sim <- matrix(NA, ncol = npcs, nrow = N_Sims)
for(i in 1:npcs) {
  for(jt in 1:N_Sims){
    low <- jt*dim(input_data)[1] - dim(input_data)[1] + 1
    high <- jt*dim(input_data)[1]
    temp <- tot_pc_sims[low:high,i]
    mean_sim[jt,i] <- mean(temp)
    sd_sim[jt,i] <- sd(temp)
    min_sim[jt,i] <- min(temp)
    max_sim[jt,i] <- max(temp)
    
  }
}

#Plotting the histograms
pdf(file = 'plots/Simulation Replications.pdf')
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(i in 1:npcs) {
  boxplot(mean_sim[,i])
  title(paste0("Mean-PC ", i))
  abline(h=mean(input_data_pca$x[,i]), col = 'red')
}

for(i in 1:npcs) {
  boxplot(sd_sim[,i])
  title(paste0("SD - PC ", i))
  abline(h=sd(input_data_pca$x[,i]), col = 'red')
}

for(i in 1:npcs) {
  boxplot(max_sim[,i])
  title(paste0("Max - PC ", i))
  abline(h=max(input_data_pca$x[,i]), col = 'red')
}

for(i in 1:npcs) {
  boxplot(min_sim[,i])
  title(paste0("min - PC ", i))
  abline(h=min(input_data_pca$x[,i]), col = 'red')
}
dev.off()












