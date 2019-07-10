#Files uses PC's to reduce dimensionality
#The wavelets are used to get information on the non-stationary PCs. 
#The individual signals are reconstructed using the AR Models. 



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


#########Getting the Streamflow Sites#########
#Reterieving Sites
sites_one <- whatNWISsites(bBox = c(-88.5, 38.5, -82.5, 41.5), parameterCd = c("00060"), hasDataTypeCd = "dv")
sites_two <- whatNWISsites(bBox = c(-82.5, 38.5, -79.5, 41.5), parameterCd = c("00060"), hasDataTypeCd = "dv")
sites_three <- whatNWISsites(bBox = c(-88.5, 37.5, -80.5, 38.5), parameterCd = c("00060"), hasDataTypeCd = "dv")
sites_four <- whatNWISsites(bBox = c(-88.5, 36.5, -81.5, 37.5), parameterCd = c("00060"), hasDataTypeCd = "dv")
sites_five <- whatNWISsites(bBox = c(-88.5, 35.5, -83.5, 36.5), parameterCd = c("00060"), hasDataTypeCd = "dv")


#Getting the site numbers
siteNumbers <- data.frame(sites = c(sites_one$site_no, sites_two$site_no, sites_three$site_no, sites_four$site_no, sites_five$site_no))
site_INFO <- readNWISsite(siteNumbers$sites)


#Getting the required sites
req_sites <- site_INFO %>% filter(drain_area_va > (5791*0.25)) #Criteria in Dave's Paper
req_sites$SUM_NA <- NA
req_sites$Record_Length <- NA
req_sites$Flow_Change <- NA 


#Plotting the Data
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(req_sites$dec_long_va,req_sites$dec_lat_va,pch=19,cex=1)
################################################################3


##############Subsetting the Streamflow Sites###########

#Getting the Record Length
for( i in 1:dim(req_sites)[1]) {
  parameterCd <- "00060"  # Discharge
  startDate <- "1937-01-01"
  endDate <- "2017-12-31"
  statCd = "00003"
  discharge <- as.matrix(readNWISdv(req_sites$site_no[i], parameterCd, startDate, endDate, statCd))
  req_sites$Record_Length[i] <- dim(discharge)[1]
  req_sites$Flow_Change[i] <- ncol(discharge)
  if(ncol(discharge) < 4) {  
    req_sites$SUM_NA[i] <- 100
    } else {
  req_sites$SUM_NA[i] <- sum(is.na(discharge[,4]))
    }
}

final_sites <- req_sites %>% filter(Record_Length == max(req_sites$Record_Length))

#Plotting the Data
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(final_sites$dec_long_va,final_sites$dec_lat_va,pch=19,cex=1)

print(paste0("The number of stream gauge stations for the 80 years analysis is - ", dim(final_sites)[1]))
#########################################################3

#######Getting the annual maximum#############
num_sites <- dim(final_sites)[1]
num_days <- max(req_sites$Record_Length) #Hard Code
streamflow <- as.data.frame(matrix(NA,nrow=num_days,ncol = num_sites))

for(i in 1:num_sites){
  parameterCd <- "00060"  # Discharge ft3/sec
  startDate <- "1937-01-01"
  endDate <- "2017-12-31"
  statCd = "00003"
  discharge <- readNWISdv(final_sites$site_no[i], parameterCd, startDate, endDate, statCd)
  streamflow[,i] <- discharge[,4]
  
}

#Getting the max of each year. 
dates_year <- format(seq(as.Date("1937-01-01"), as.Date("2017-12-31"), by="days"),"%Y")
streamflow$year <- dates_year
max_annual <- streamflow %>% group_by(year) %>% summarise_all(funs(max))
max_annual$year <- NULL
max_annual <- as.data.frame(max_annual)
###############################################################3



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
max_annual_pca <- prcomp(max_annual, scale = TRUE)
var <- cumsum(max_annual_pca$sdev^2)
pdf(file='Individual AR Models.pdf')
plot(var/max(var),pch=19, main = "Variance explained by PCs",xlab = "PC's")
dev.off()
pdf(file = 'Selected PCs.pdf')
par(mfrow=c(3,2))
par(mar = c(4, 1, 1, 1))
for(i in 1:5){plot(max_annual_pca$x[,i],typ='l', ylab = NA, xlab = paste0("PC ",i))
}
dev.off()
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
pc_loadings  <- max_annual_pca$x[,1:5]
#write.table(pc_loadings, "Selected PCs.txt", sep = " ")
#####################################

##################Wavelet Analysis###############
# We have reduced the dimensionality of the data using PCA. 
#Now we try to uncover further spatial aspects using Wavelet Analysis.
# We first fit a wavelet to each of the PC's 

#Step 1:- Getting the significant periods (over red noise). 
pdf(file = 'PC Wavelet Spectrums.pdf')
predict_ahead <- 5
pc_training <- head(pc_loadings,-predict_ahead);pc_testing <- tail(pc_loadings,predict_ahead)
yr1 <- 1937; yr2 <- 2017-predict_ahead; yr3 <- yr2+1; yr4 <- 2017
sig_scales <- as.list(1)
par(mfrow=c(3,2))
par(mar = c(4, 1, 1.5, 1))
for(i in 1:ncol(pc_loadings)) {
  p <- pc_training[,i]
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
pdf(file = 'PC Reconstructions.pdf')
par(mfrow=c(3,2))
par(mar = c(4, 1, 1.5, 1))
for(i in 1:ncol(pc_loadings)) {
  p <- pc_training[,i]
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
pdf(file = 'PCs - Signals and Noise.pdf')
par(mfrow=c(3,2))
par(mar = c(4, 1, 1.5, 1))
signals <- list(NA)
for(i in 1:ncol(pc_loadings)) {
  
  #Fitting the Wavelet
  p <- pc_training[,i]
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
pdf(file='Individual AR Models.pdf')
signals <- list(NA)
tot_pc_pred <- matrix(NA,nrow=predict_ahead,ncol=ncol(pc_loadings))
par(mfrow=c(2,2))
par(mar = c(4, 2, 1.5, 1))
for(i in 1:ncol(pc_loadings)) {
  
  #Fitting the Wavelet
  p <- pc_training[,i]
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
    #plot(2001:2017,pc_loadings[65:81,i], type='l', main = paste0("Predictions for PC ", i)
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
    #plot(2001:2017,pc_loadings[65:81,i], type='l', main = paste0("Predictions for PC ", i)
    #    ,xlab = "Year", ylab = "PC")
    #lines((2017-predict_ahead+1):2017, predict_ts,col='red')
    #legend('bottomright', legend = c("Real","predicted"), lty = 1, col =c('black','red'), cex = 0.6)
    
    
  }
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
dev.off()

#Plotting the Predictions.
pdf(file = 'PC Predictions.pdf')
signals <- list(NA)
tot_pc_pred <- matrix(NA,nrow=predict_ahead,ncol=ncol(pc_loadings))
par(mfrow=c(3,2))
par(mar = c(4, 2, 1.5, 1))
for(i in 1:ncol(pc_loadings)) {
  
  #Fitting the Wavelet
  p <- pc_training[,i]
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
    #plot(forecast(fit,h=20), xlab = paste0("Reconstructed PC ", i))
    predict_ts <- forecast(fit,h=5)$mean
    tot_pc_pred[,i] <- predict_ts
    
    resd <- fit$resid[!is.na(fit$resid)]
    #Diagnostics
    #par(mfrow = c(2,2))
    #plot(density(resd), main ="Resd");acf(resd, main = "ACF");pacf(resd, main = "PACF")
    #par(mfrow=c(1,1))
    
    #Getting the predictions
    plot(2001:2017,pc_loadings[65:81,i], type='l', main = paste0("Predictions for PC ", i)
         ,xlab = "Year", ylab = "PC")
    lines((2017-predict_ahead+1):2017, predict_ts,col='red')
    legend('bottomright', legend = c("Real","predicted"), lty = 1, col =c('black','red'), cex = 0.6)
    
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
      #plot(forecast(fit,h=20), xlab =nota)
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
    plot(2001:2017,pc_loadings[65:81,i], type='l', main = paste0("Predictions for PC ", i)
        ,xlab = "Year", ylab = "PC")
    lines((2017-predict_ahead+1):2017, predict_ts,col='red')
    legend('bottomright', legend = c("Real","predicted"), lty = 1, col =c('black','red'), cex = 0.6)
    
    
  }
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
dev.off()


###########Spatial Distribution of MSE##############

#PC Loadings
library("biwavelet")
library("plotrix")
library("maps")
loadings <- max_annual_pca$rotation 
pdf(file = 'Spatial Distribution of PCs.pdf')
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
for(i in 1:dim(pc_loadings)[2]) {
  ju <- loadings[,i]
  map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"), boundary = TRUE)
  points(final_sites$dec_long_va,final_sites$dec_lat_va,pch=19,cex=1,col=color.scale(ju,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
  title(paste0("Spatial Distribution of PC ", i))
  legend("bottomright", legend = c("high","low"), col = c("blue","red"), cex =0.6, pch =19)
}
dev.off()




PC_Predictions <- tot_pc_pred #These are the predictions
nComp = ncol(PC_Predictions) 
Predictions_Scaled =  PC_Predictions %*% t(max_annual_pca$rotation[,1:nComp])
for(i in 1:nComp) {   
  Predictions_Scaled[,i] <- scale(Predictions_Scaled[,i], center = FALSE , scale=1/max_annual_pca$scale[i]) }

for(i in 1:nComp) {   
  Predictions_Scaled[,i] <- scale(Predictions_Scaled[,i], center = -1 * max_annual_pca$center[i], scale=FALSE)
}

#Computing Site - Specific MSE
True_Values <- tail(max_annual,5)
site_MSE <- bias <- rep(NA,ncol(Predictions_Scaled))
for(i in 1:ncol(Predictions_Scaled)) {
  site_MSE[i] <- mean((Predictions_Scaled[,i]-True_Values[,i])^2)
  bias[i] <- sign(sum(Predictions_Scaled[,i]-True_Values[,i]))
}

#Plotting the Data
site_MSE_adj <- as.numeric(site_MSE/(colMeans(max_annual)^2)) + 0.3
for(i in 1:ncol(Predictions_Scaled)) {
  if(bias[i] == -1) {bias[i] = 19 #Triangle
  } else {bias[i] = 17 } #Over Predictions 
  }

pdf(file = 'Prediction Error.pdf')
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(final_sites$dec_long_va,final_sites$dec_lat_va,
       col=color.scale(final_sites$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"),
       cex=round(site_MSE_adj,4),
       pch=bias)
title("Error in Streamflow Prediction")
legend("bottomright", c("Color - Drainage Area","Size - Error in MSE adjusted","Shape - Bias"), cex = 0.6)
dev.off()





