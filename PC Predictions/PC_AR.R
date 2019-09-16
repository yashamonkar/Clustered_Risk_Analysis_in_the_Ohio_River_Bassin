#####
#1. The objective of this file is to make predictions in PC space. 
#2. Use a dimension lowering algorithm - linear PCs
#3. Make predictions and simulations on those PCs. 
#4. Check if the relevant indices and the spectral structure has been replicated. 


#########Settting the Path
setwd("~/Correlated Risk Analysis/Decadal Influences/PC Predictions")


###Loading Libraries##########
library(forecast)
library(e1071)

###Reading the Data and Initial Manipulation###
input_data <- read.table("data/Max_Annual_Streamflow.txt", sep="", header = TRUE)
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)
######################################2


####Principal Component Analysis#############
data_pca <- prcomp(input_data, scale = TRUE)
var <- cumsum(data_pca$sdev^2)
#pdf(file = "plots/AR/Variance explained.pdf")
plot(var/max(var),pch=19, main = "Variance explained by PCs",xlab = "PC's",ylab="Fraction Variance explained")
npcs <- 3 #This is the selected number of PC'S
abline(h = var[npcs]/max(var), lty = 2, col='red')
##dev.off()
pcs_sel <- data_pca$x[,1:npcs]
#################################################2


###############Diagnostics##########
#pdf(file = "plots/AR/PC Diagnostics.pdf")
for(i in 1:ncol(pcs_sel)) {

x <- pcs_sel[,i]

#Plotting the Data
  plot(pcs_sel[,i], main = paste0("PC - ", i), type = 'l', 
       xlab = "Years", ylab = " ")



#Histograms
par(mfrow=c(1,1), mar=c(4,2,4,1))
  hist(pcs_sel[,i], main = paste0("PC - ", i), 
       xlab = " ")


#PACF and ACF's
par(mfrow=c(2,1), mar=c(2,4,4,1))
acf(x, main = paste0("PC -",i ))
pacf(x, main = paste0("PC -",i )) 
par(mfrow = c(1,1))

}
##dev.off()
######################AR Model################
#pdf(file = "plots/AR/Simulated PC.pdf")
for(i in 1:ncol(pcs_sel)){
x_s <- pcs_sel[,i]
x_s <- scale(x_s)
x <- x_s-min(x_s)*1.01
og_density <- density(x_s)

#Auto-Regressive Modelling
mod.ar <- ar(x, aic = TRUE)
x_p <- x
lambda <- BoxCox.lambda(x_p, method = "loglik")
x_t <- (-1+x_p^lambda)/lambda
mod_ar <- ar(x_t, aic = TRUE)

#AR Simulations
N_Sims <- 1000
mean_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
sd_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
min_sim <- max_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
consolid_sims <- as.list(1)
consolid_points <- as.list(1)

for(j in 1:N_Sims) {
  sims <- arima.sim(n = length(x),
                    list(order = arimaorder(mod_ar),
                         ar = mod_ar$ar))
  sims <- min(x_s)*1.01+(lambda*sims+1)^(1/lambda)
  sims <- na.omit(as.numeric(sims))
  mean_sim[j,1] <- mean(sims)
  sd_sim[j,1] <- sd(sims)
  max_sim[j,1] <- max(sims)
  min_sim[j,1] <- min(sims)
  sims_dens <- density(sims)
  consolid_sims[[j]] <- sims_dens$y
  consolid_points[[j]] <- sims_dens$x
   }
consolid_points <- unlist(consolid_points)
consolid_sims <- unlist(consolid_sims)
consolid_df <- data.frame(x=consolid_points,y=consolid_sims)
consolid_df$x <- cut(consolid_df$x, breaks = seq(-6,6,.2))
mid_points <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", consolid_df$x) ),
                    upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", consolid_df$x) ))
consolid_df$x <- rowMeans(mid_points)
og_df <- data.frame(x1=og_density$x, y1= og_density$y)
  

  
plo <- ggplot(og_df, aes(x=x1,y=y1))+
  geom_line(size=1.25)+
  scale_x_continuous(limits=c(-5,5)) +
  geom_boxplot(consolid_df, mapping = aes(y = y, x = x,group = x), outlier.shape = NA,outlier.colour = NA)+ 
  ggtitle(paste0("Simulated PDF \n for PC",i)) +
  xlab(" ") + 
  ylab("Probability Density")
print(plo)


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
##dev.off()


##################Simulating the PCs Data. 
N_Sims <- 1000

x <- scale(pcs_sel[,1])
mod_ar_1 <- ar(x, aic = TRUE)

x <- scale(pcs_sel[,2])
mod_ar_2 <- ar(x, aic = TRUE)

x <- scale(pcs_sel[,3])
mod_ar_3 <- ar(x, aic = TRUE)





####Combined Simulations###
N_sims <- 1000
sel_pcs <- npcs
mean_sim <- matrix(NA, ncol = ncol(input_data), nrow = N_Sims)
sd_sim <- matrix(NA, ncol = ncol(input_data), nrow = N_Sims)
min_sim <- max_sim <- matrix(NA, ncol = ncol(input_data), nrow = N_Sims)


for(sim in 1:N_sims) {
  sim_pcs <- matrix(NA, ncol = sel_pcs, nrow = dim(pcs_sel)[1])
  for(j in 1:npcs) {
  x <- scale(pcs_sel[,1])
  mod_ar <- ar(x, aic = TRUE)
  sims <- arima.sim(n = length(x),
                       list(order = arimaorder(mod_ar),
                            ar = mod_ar$ar))
  sims <- sims*sd(pcs_sel[,j]) #Resaling back
  sim_pcs[,j] <- sims }

  #Converting to Actual Field Space
  PC_Simulations <- sim_pcs #These are the predictions
  nComp = sel_pcs
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
#pdf(file = "plots/AR/Site Mean Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))

for(j in 1:ncol(input_data)) {

  boxplot(mean_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(mean_sim[,j], .5)-.2,quantile(mean_sim[,j], .95)+.2))
  title(paste0("Mean-Site ", j), cex = 0.5)
  abline(h=mean(input_data[,j]), col = 'red')
  
}
##dev.off()

##Plotting the Standard #deviation###
#pdf(file = "plots/AR/Site Standard deviation Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(sd_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(sd_sim[,j], .5)-.2,quantile(sd_sim[,j], .95)+.2))
  title(paste0("SD-Site ", j), cex = 0.5)
  abline(h=sd(input_data[,j]), col = 'red')
  
}
##dev.off()

##Plotting the Max###
#pdf(file = "plots/AR/Site Max Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(max_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(max_sim[,j], .5)-.2,quantile(max_sim[,j], .95)+.2))
  title(paste0("Max-Site ", j), cex = 0.5)
  abline(h=max(input_data[,j]), col = 'red')
  
}
##dev.off()

##Plotting the Min###
#pdf(file = "plots/AR/Site Min Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(min_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(min_sim[,j], .5)-.2,quantile(min_sim[,j], .95)+.2))
  title(paste0("Min-Site ", j), cex = 0.5)
  abline(h=min(input_data[,j]), col = 'red')
  
}
##dev.off()

#####Plotting the Simulations by Sites######
#pdf(file = "plots/AR/Site Specific Simulation Moments.pdf")
par(mfrow=c(1,4))
par(mar = c(4, 2, 1.5, 1))

for(j in 1:ncol(input_data)) {
  
  #Mean
  boxplot(mean_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(mean_sim[,j], .5)-.2,quantile(mean_sim[,j], .95)+.2))
  title(paste0("Mean-Site ", j), cex = 0.5)
  abline(h=mean(input_data[,j]), col = 'red')
  
  #Standard #deviation
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
##dev.off()
##########################################################

##################Spatial Distribution of Errors#########################
pdf(file = "plots/AR/Spatial Distribution of Errors.pdf")
par(mfrow=c(1,1))
#Plotting the Distribution of Data. 
library("biwavelet")
library("plotrix")
library("maps")

#Means
mean_diff <- apply(mean_sim, 2, median, na.rm = TRUE)-colMeans(input_data)
bias <- rep(NA,ncol(input_data))
for(i in 1:length(bias)) {
  if(mean_diff[i] < 0) {bias[i] = 1}
  else {bias[i] = 19}
}
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(site_info$dec_long_va,site_info$dec_lat_va,
       col=color.scale(site_info$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"),
       cex=abs(mean_diff)*200,
       pch=bias)
title("Error in Mean Streamflow Prediction")
legend("bottomright", c("Color - Drainage Area","Size - Error in MSE adjusted","Shape - Bias"), cex = 0.6)



#SD
sd_diff <- apply(sd_sim, 2, median, na.rm = TRUE)-apply(input_data, 2, sd, na.rm = TRUE)
bias <- rep(NA,ncol(input_data))
for(i in 1:length(bias)) {
  if(sd_diff[i] < 0) {bias[i] = 1}
  else {bias[i] = 19}
}
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(site_info$dec_long_va,site_info$dec_lat_va,
       col=color.scale(site_info$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"),
       cex=abs(sd_diff)*15,
       pch=bias)
title("Error in Standard Deviation Streamflow Prediction")
legend("bottomright", c("Color - Drainage Area","Size - Error in MSE adjusted","Shape - Bias"), cex = 0.6)


#Max
max_diff <- apply(max_sim, 2, median, na.rm = TRUE)-apply(input_data, 2, max, na.rm = TRUE)
bias <- rep(NA,ncol(input_data))
for(i in 1:length(bias)) {
  if(max_diff[i] < 0) {bias[i] = 1}
  else {bias[i] = 19}
}
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(site_info$dec_long_va,site_info$dec_lat_va,
       col=color.scale(site_info$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"),
       cex=abs(max_diff)*3,
       pch=bias)
title("Error in Maximum Deviation Streamflow Prediction")
legend("bottomright", c("Color - Drainage Area","Size - Error in MSE adjusted","Shape - Bias"), cex = 0.6)



#Min
min_diff <- apply(min_sim, 2, median, na.rm = TRUE)-apply(input_data, 2, min, na.rm = TRUE)
bias <- rep(NA,ncol(input_data))
for(i in 1:length(bias)) {
  if(min_diff[i] < 0) {bias[i] = 1}
  else {bias[i] = 19}
}
map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(site_info$dec_long_va,site_info$dec_lat_va,
       col=color.scale(site_info$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"),
       cex=abs(min_diff)*2,
       pch=bias)
title("Error in Minimum Deviation Streamflow Prediction")
legend("bottomright", c("Color - Drainage Area","Size - Error in MSE adjusted","Shape - Bias"), cex = 0.6)
#dev.off()


##############Plotting the Wavelet Spectra Distribution###########

#Functions for the Wavelet Spectrum and Confidence Intervals
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
### Confidence level function
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


N_Sims <- 1000 

#pdf(file = "plots/AR/Global Wavelet Simulations.pdf")
par(mfrow = c (1,1), mar = c(4,4,4,1))
for(i in 1:npcs) {
  x <- scale(pcs_sel[,i])
  wlt_og=wavelet(x)
  Cw=CI(0.9,x,"w");
  C=CI(0.9,x,"r");
  plot(wlt_og$period,wlt_og$p.avg,xlim=c(0,75),
       main=paste0("Global Wavelet Spectrum for PC",i),
       xlab="Period",ylab="Variance", col ='red',
       ylim = c((min(wlt_og$p.avg)-2),(max(wlt_og$p.avg)+2))) 
  lines(wlt_og$period,wlt_og$p.avg, col ='red');
  lines(wlt_og$period,Cw$sig, lty = 2);
  lines(wlt_og$period,C$sig,col="red", lty = 2);
  avg_pow_matrix <- matrix(NA, nrow = N_Sims, ncol = length(wlt_og$p.avg))
  
  for(j in 1:N_Sims) {
  
  mod_ar <- ar(x, aic = TRUE)
  sims <- arima.sim(n = length(x),
                    list(order = arimaorder(mod_ar),
                         ar = mod_ar$ar))
  wlt=wavelet(sims);
  avg_pow_matrix[j,] <- wlt$p.avg
  }
  lower_percentile <- apply(avg_pow_matrix, 2, function(x) quantile(x, probs=.05))
  upper_percentile <- apply(avg_pow_matrix, 2, function(x) quantile(x, probs=.95))
  median_percentile <- apply(avg_pow_matrix, 2, median)
  lines(wlt_og$period, lower_percentile)
  lines(wlt_og$period, upper_percentile)
  lines(wlt_og$period, median_percentile, lwd = 2)
}
##dev.off()