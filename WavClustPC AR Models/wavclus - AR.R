#Files uses PC's to reduce dimensionality
#The wavelets are used to get information on the non-stationary PCs. 
#The individual signals are reconstructed using the AR Models. 
setwd("~/Correlated Risk Analysis/Decadal Influences/WavClustPC AR Models")



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


##############Hierarchical Clustering################
#Keeping last 5 years for predictions. 
predict_ahead <- 5
training_set <- head(max_annual,-predict_ahead)
testing_set <- tail(max_annual,predict_ahead)

nr <- nrow(training_set)
nc <- ncol(training_set)
yr1 <- 1937;yr2 <- 2017-predict_ahead;yr3 <- yr2+1;yr4 <- 2017 
np <- length(wavelet(training_set[,1])$p.avg) #Number of wavelet scales. 
w.arr=array(NA,dim=c(nc,np)) #Storing the wavelet scales
#Standardize to make sure the global wavelets are comparable. 
training_set_scaled=training_set
for (i in 1:nc)training_set_scaled[,i]=(training_set[,i]-mean(training_set[,i]))/sd(training_set[,i])
for ( i in 1:nc)w.arr[i,]=wavelet(training_set_scaled[,i])$p.avg

#Hierarchical Clustering

par(mfrow=c(1,1))
clus=hclust(dist(w.arr),method="ward.D2")
plot(clus)
(cls <- identify(clus))
nclus=length(cls)
ju=rep(NA,nc)
for (i in 1:nclus)ju[cls[[i]]]=rep(i,length(cls[[i]]))
par(mfrow=c(1,1))

#Plotting the Data
pdf("plots/Cluster Distribution in the Basin.pdf")
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(final_sites$dec_long_va,final_sites$dec_lat_va,pch=19,cex=1, col=ju)
title("Spatial Distribution of Clusters")
dev.off()


################Principal Component Analysis################
#Here we extract the 1st PC of each cluster. 

#Showing variance explained by PCs
pdf("plots/Variance Explained by PCs.pdf")
for (i in 1:nclus) {
  tpd=training_set[,cls[[i]]]
  pcw=prcomp(tpd, scale = TRUE)
  var <- cumsum(pcw$sdev^2)
  plot(var/max(var),pch=19, main = paste0("Variance explained by PCs of Cluster ",i), xlab = "PC's", ylab = "Fraction Total Variance")

}
dev.off()


##Plotting the Wavelet Spectrums.
library("biwavelet")
library("plotrix")
library("maps")
library("stats") 
pdf(file = 'plots/PC Wavelet Spectrums.pdf')
PCs <- matrix(NA,ncol=nclus,nrow=dim(training_set)[1])
for (i in 1:nclus) {
  tpd=training_set[,cls[[i]]]
  pcw=prcomp(tpd, scale = TRUE)
  p=pcw$x[,1]
  PCs[,i] <- p
  wlt=wavelet(p)
  Cw=CI(0.9,p,"w")
  C=CI(0.9,p,"r")
  par(mfrow=c(2,2));par(mar=c(4,2,3,0.5))
  plot(yr1:yr2,p,xlab="Year",main=paste("Cluster ",as.character(i), " PC 1"))
  lines(lowess(yr1:yr2,p,f=1/9),lwd=2,col="red")
  plot(wlt$period,wlt$p.avg,xlim=c(0,75),main="Global Wavelet Spectrum",xlab="Period",ylab="Variance") 
  lines(wlt$period,wlt$p.avg);
  lines(wlt$period,Cw$sig)
  lines(wlt$period,C$sig,col="red")
  wt1=wt(cbind(yr1:yr2,p))
  plot(wt1, type="power.corr.norm", xlab="Year",main=paste("Wavelet Cluster ",as.character(i), " PC 1"))
  jj=rep(1,nc);jj[ju==i]=19
  jc=rep(0,nc);jc[cls[[i]]]=abs(pcw$rotation[,1])
  map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
  points(final_sites$dec_long_va,final_sites$dec_lat_va,pch=jj,cex=1,col=color.scale(jc,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
  
}
dev.off()
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
PCs <- as.data.frame(PCs)
for(i in 1:ncol(PCs)) {
  colnames(PCs)[i] <- c(paste0("Cluster_",i,"_PC_1"))
}
write.table(PCs,"results/PCs.txt",sep=" ")


##Reconstructing the Wavelets - Sanity Check
pdf(file = 'plots/Cluster PC Reconstructions.pdf')
for(i in 1:nclus) {
  tpd <- training_set[,cls[[i]]]
  pcw <- prcomp(tpd, scale = TRUE)
  p=pcw$x[,1]
  wlt <- wavelet(p)
  Cd <- 0.776;psi0 <- pi^(-.025);dj=0.025 #From the Torrence and Compo
  reconst <- matrix(NA, ncol = length(p), nrow = length(wlt$scale))
  for(j in 1:ncol(reconst)) {reconst[,j] <- Re(wlt$wave[,j])/(wlt$scale[i]^0.5)}
  p_reconst <- colSums(reconst)*dj/(Cd*psi0) 
  plot(yr1:yr2, p_reconst, type='l',col='red', main = paste0("Reconstruced Cluster ", i, " PC"),
       xlab = "Year")
  lines(yr1:yr2,p)
  legend('topright', legend = c("Signal","Constructed"), cex = 0.6, col = c('black','red'),lty = 1)
  print(paste0("The correlation between signal and reconstructed TS for Cluster PC ", i, " is ", round(cor(p_reconst,p),2)))
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
dev.off()


################Auto Regression#################

#Getting the significant parts.
sig_scales <- as.list(1)
for(i in 1:nclus) {
  tpd <- training_set[,cls[[i]]]
  pcw <- prcomp(tpd, scale = TRUE)
  p=pcw$x[,1]
  wlt <- wavelet(p)
  C_r <- CI(0.9,p,'r')
  C_w <- CI(0.9,p,'w')
  plot(wlt$period,wlt$p.avg,xlim=c(0,75),main=paste0("Global Wavelet Spectrum for Cluster ", i),xlab="Period",ylab="Variance")
  lines(wlt$period,wlt$p.avg)
  lines(wlt$period,C_r$sig,col='red')
  lines(wlt$period,C_w$sig, col ='black')
  sig_scales[[i]] <- which(wlt$p.avg > C_r$sig)
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))


##Decomposing each PC into signals and noise
pdf("plots/Cluster PCs Reconstructions with Signal and Noise.pdf")
signals <- list(NA)
for(i in 1:nclus) {
  
  #Fitting the Wavelet
  tpd <- training_set[,cls[[i]]]
  pcw <- prcomp(tpd, scale = TRUE)
  p=pcw$x[,1]
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
    plot(yr1:yr2, p, type='l', main = paste0("Signal and Noise for Cluster ", i, " PC"))
    p_reconst <- colSums(reconst)*dj/(Cd*psi0) 
    lines(yr1:yr2,p_reconst, col ='blue', lty = 2)
    legend('bottomright', legend = c("PC","Signal","Noise"), cex = 0.6, col = c('black','red','blue'),lty = c(1,1,2))
    
  } else {
    #Plotting the original data. 
    plot(yr1:yr2, p, type='l', main = paste0("Signal and Noise for Cluster ", i, " PC 1"), xlab = "Year")
    
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


##Fitting ARMA Models
pdf(file='plots/ Individual AR Models.pdf')
signals <- list(NA)
tot_pc_pred <- matrix(NA,nrow=predict_ahead,ncol=nclus)
#par(mfrow=c(2,2))
#par(mar = c(4, 2, 1.5, 1))
for(i in 1:nclus) {
  
  #Fitting the Wavelet
  tpd <- training_set[,cls[[i]]]
  pcw <- prcomp(tpd, scale = TRUE)
  p=pcw$x[,1]
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
      nota <- paste0("Reconstructed Cluster ", i, "'s PC 1 signal ", jks)
      if(jks == ncol(breakdown_ts)){nota <- paste0("Reconstructed Cluster ", i, "'s PC 1 noise ")}
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


###Plotting the Predictions.
pdf(file = 'plots/Cluster PC Predictions.pdf')
signals <- list(NA)
tot_pc_pred <- matrix(NA,nrow=predict_ahead,ncol=nclus)
for(i in 1:nclus) {
  
  #Fitting the Wavelet
  tpd <- training_set[,cls[[i]]]
  pcw <- prcomp(tpd, scale = TRUE)
  p=pcw$x[,1]
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
    plot(yr1:yr2,p, type='l', main = paste0("Predictions for Cluster ", i, " PC 1")
         ,xlab = "Year", ylab = "PC", xlim = c(yr1-1,yr4+1))
    lines(yr3:yr4, predict_ts,col='red')
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
      nota <- paste0("Reconstructed Cluster ", i, "'s PC 1 signal ", jks)
      if(jks == ncol(breakdown_ts)){nota <- paste0("Reconstructed Cluster ", i, "'s PC 1 noise ")}
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
    #Getting the predictions
    plot(yr1:yr2,p, type='l', main = paste0("Predictions for Cluster ", i, " PC 1")
         ,xlab = "Year", ylab = "PC", xlim = c(yr1-1,yr4+1))
    lines(yr3:yr4, predict_ts,col='red')
    legend('bottomright', legend = c("Real","predicted"), lty = 1, col =c('black','red'), cex = 0.6)
    
    
  }
}
par(mfrow=c(1,1));par(mar = c(4, 3, 3, 1))
dev.off()


###########Computing the Predictions##############

PC_Predictions <- tot_pc_pred #These are the predictions
All_Predictions <- matrix(NA, nrow=5, ncol = dim(testing_set)[2])
test <- testing_set #True Data. 
sze=rep(1,nc) #Size parameter MSE
jj=rep(1,nc) #Bias.
for(i in 1:nclus){

nComp = 1 
tpd <- training_set[,cls[[i]]]
pcw <- prcomp(tpd, scale = TRUE)

Predictions_Scaled =  PC_Predictions[,i] %*% t(pcw$rotation[,1:nComp])  #Predicted PC values * PC1 loadings

for(j in 1:ncol(Predictions_Scaled)) {   
  Predictions_Scaled[,j] <- scale(Predictions_Scaled[,j], center = FALSE , scale=1/pcw$scale[j]) }

for(j in 1:ncol(Predictions_Scaled)) {   
  Predictions_Scaled[,j] <- scale(Predictions_Scaled[,j], center = -1 * pcw$center[j], scale=FALSE)
}

#Computing Site - Specific MSE
True_Values <- testing_set[,cls[[i]]]
site_MSE <- bias <- rep(NA,ncol(True_Values))
for(j in 1:ncol(True_Values)) {
  site_MSE[j] <- mean((Predictions_Scaled[,j]-True_Values[,j])^2)
  bias[j] <- sign(sum(Predictions_Scaled[,j]-True_Values[,j]))
}

#Saving the Predictions
for(j in 1:length(cls[[i]])) {
  All_Predictions[,cls[[i]][j]] <- Predictions_Scaled[,j]
}
#Plotting Help
site_MSE_adj <- as.numeric(site_MSE/(colMeans(True_Values)^2))*6
for(j in 1:ncol(Predictions_Scaled)) {
  if(bias[j] == -1) {bias[j] = 19 #Triangle
  } else {bias[j] = 17 } #Over Predictions 
}

#Cluster members
ju=rep(NA,nc)
for (j in 1:nclus)ju[cls[[j]]]=rep(j,length(cls[[j]]))

jj[cls[[i]]]=bias
sze[cls[[i]]]=site_MSE_adj
}

pdf("plots/Streamflow Predictions.pdf")
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(final_sites$dec_long_va,final_sites$dec_lat_va,
       pch=jj,
       cex=sze,
       col=color.scale(final_sites$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"))
title("Error in Streamflow Prediction")
legend("bottomright", c("Color - Drainage Area","Size - Error in MSE adjusted","Shape - Bias"), cex = 0.6)
dev.off()



########3###Comparision against base mark prediction Skill Testing################
#1. Long Term Mean. 
#2. AR applied to the raw test series. 


#########Long Term Mean####################
ltm <- matrix(colMeans(training_set), nrow = 1, ncol = ncol(testing_set))
ltm_skill <- matrix(NA, nrow = predict_ahead, ncol = ncol(testing_set))
for(i in 1:ncol(ltm_skill)) {
  for(j in 1:nrow(ltm_skill)) {
    if(abs(ltm[i]-testing_set[j,i]) > abs(All_Predictions[j,i]-testing_set[j,i])) { ltm_skill[j,i] = 1
    } else { ltm_skill[j,i] = 0 
      } 
  }
}

par(mfrow=c(1,1))
pdf("plots/Skill vs Long Term Mean.pdf")
for(i in 1:nrow(ltm_skill)) {
  skill <- ltm_skill[i,]
  for(j in 1:length(skill)) { 
    if(skill[j]==0) {skill[j] = c("red")
    } else { skill[j] = c("blue")
      }
  }

map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(final_sites$dec_long_va,final_sites$dec_lat_va,
       pch=19,
       cex=1,
       col=skill)
title(paste0("Skill Testing vs Mean for Year ", i))
legend("bottomright", c("Correct Prediction", "Wrong Prediction"), cex = 0.6, pch = 19, col = c("blue","red"))
}

skill <- colSums(ltm_skill)
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(final_sites$dec_long_va,final_sites$dec_lat_va,
       pch=19,
       cex=skill/2,
       col=color.scale(final_sites$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"))
title("Combined Skill vs Long Term Mean")
legend("bottomright", c("Color - Drainage Area","Size - Skill"), cex = 0.6)
dev.off()



################AR on raw time series####################
ar_raw <- matrix(NA, nrow = predict_ahead, ncol = ncol(testing_set))

pdf("plots/Raw Time Series AR Models.pdf")
par(mfrow=c(3,2))
par(mar = c(4, 2, 1.5, 1))
for(i in 1:ncol(ar_raw)) {
  temp_raw <- training_set[,i]
  fit <- ar(temp_raw, order.max = 10, aic = TRUE)
  ord <- arimaorder(fit)
  nota <- paste0("Station ", i)
  plot(forecast(fit,h=20), xlab =nota)
  ar_raw[,i] <- forecast(fit,h=5)$mean
}
dev.off()


ar_skill <- matrix(NA, nrow = predict_ahead, ncol = ncol(testing_set))
for(i in 1:ncol(ar_skill)) {
  for(j in 1:nrow(ar_skill)) {
    if(abs(ar_raw[j,i]-testing_set[j,i]) > abs(All_Predictions[j,i]-testing_set[j,i])) { ar_skill[j,i] = 1
    } else { ar_skill[j,i] = 0 
    } 
  }
}

par(mfrow=c(1,1))
pdf("plots/Skill vs Raw TS AR Models.pdf")
for(i in 1:nrow(ar_raw)) {
  skill <- ar_skill[i,]
  for(j in 1:length(skill)) { 
    if(skill[j]==0) {skill[j] = c("red")
    } else { skill[j] = c("blue")
    }
  }
  
  map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
  points(final_sites$dec_long_va,final_sites$dec_lat_va,
         pch=19,
         cex=1,
         col=skill)
  title(paste0("Skill Testing vs AR for Year ", i))
  legend("bottomright", c("Correct Prediction", "Wrong Prediction"), cex = 0.6, pch = 19, col = c("blue","red"))
}

skill <- colSums(ar_skill)
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(final_sites$dec_long_va,final_sites$dec_lat_va,
       pch=19,
       cex=skill/2,
       col=color.scale(final_sites$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"))
title("Combined Skill vs Raw AR")
legend("bottomright", c("Color - Drainage Area","Size - Skill"), cex = 0.6)
dev.off()
############################################################
