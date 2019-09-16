#####
#1. The objective of this file is to make predictions in PC space. 
#2. Use a dimension lowering algorithm - linear PCs
#3. Make predictions and simulations on those PCs. 
#4. Check if the relevant indices and the spectral structure has been replicated. 


#########Settting the Path
setwd("~/Correlated Risk Analysis/Decadal Influences/PC Predictions")


###Loading Libraries##########
library(ggplot2)


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

######################KNN Simulations################
#k=number of neighbors,ns= number of simulations
#otherwise a v similar structure
knnsim=function(y,x,xtest,k,ns,w=NULL){
  x=as.matrix(x)
  xtest=as.matrix(xtest)
  y=as.matrix(y)
  if(is.null(w))w=rep(1,ncol(x))
  if(nrow(y)!= nrow(x))
    print('error: lengths of y and x differ')
  if(ncol(x)!= ncol(xtest))
    print('error: col lengths of x and xtest differ')
  
  na=rep(NA,nrow(xtest))
  yknn=matrix(NA,nrow=ns,ncol(y))
  
  yk=na
  
  yr=seq(1,nrow(y))
  
  for(i in 1:nrow(xtest)){
    a=matrix(NA,nrow(x),ncol(x))
    for(n in 1:ncol(x))
      a[,n]=100*w[n]*(x[,n]- xtest[i,n])^2
    
    #Finding the cumulative distance matrix
    c=rowSums(a,na.rm=T)
    
    yk[rank(c,ties.method='first')]=yr
    
    j=rank(c)		#equal prob for equidistant neighbours
    sj=sum(j[1:k]^(-1))
    pj=(j[1:k]^(-1))/sj #Probability of sampling
    
    ynp=sample(yk[1:k],ns,replace=T,prob=pj) #yk[1:k] are the closet neighbours. 
    for(p in 1:ncol(y)) 
      yknn[,p]=y[ynp,p]
    
  }
  return(yknn)
}

##################Simulating the PCs Data################ 
for(i in 1:ncol(pcs_sel)) {
  y <- scale(pcs_sel[,i])
  x <- 1:length(y)

og_density <- density(y)
consolid_sims <- as.list(1)
consolid_points <- as.list(1)
sims_all <- as.list(1)
x_tests <- seq(1,length(x),by = 0.5)
for(j in 1:length(x_tests)) {
  xtest <- x_tests[j]
  sims_all[[j]] <-  sims <- knnsim(y,x,xtest,9,100)
  sims_dens <- density(sims)
  consolid_sims[[j]] <- sims_dens$y
  consolid_points[[j]] <- sims_dens$x
}
consolid_points <- unlist(consolid_points)
consolid_sims <- unlist(consolid_sims)
consolid_df <- data.frame(x=consolid_points,y=consolid_sims)
consolid_df$x <- cut(consolid_df$x, breaks = seq(-4,4,.2))
mid_points <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", consolid_df$x) ),
                    upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", consolid_df$x) ))
consolid_df$x <- rowMeans(mid_points)
og_df <- data.frame(x1=og_density$x, y1= og_density$y)

plo <- ggplot(og_df, aes(x=x1,y=y1))+
  geom_line(size=1.25)+
  scale_x_continuous(limits=c(-3,3)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_boxplot(consolid_df, mapping = aes(y = y, x = x,group = x), outlier.shape = NA,outlier.colour = NA)+ 
  ggtitle(paste0("Simulated PDF \n for PC ", i)) +
  xlab(" ") + 
  ylab("Probability Density")
print(plo)

}


#Simulation Accuracy
N_Sims <- 1000
for(i in 1:ncol(pcs_sel)) {
  y <- scale(pcs_sel[,i])
  x <- 1:length(y)
  mean_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
  sd_sim <- matrix(NA, ncol = 1, nrow = N_Sims)
  min_sim <- max_sim <- matrix(NA, 1, nrow = N_Sims)
  skew_sim <- matrix(NA, 1, nrow = N_Sims)
  for(j in 1:N_Sims) {
    sims_all <- as.list(1)
    x_tests <- seq(1,length(x),by = 1)
    for(k in 1:length(x_tests)) {
      xtest <- x_tests[k]
      sims_all[[k]] <- knnsim(y,x,xtest,9,1)}
    sims_all <- unlist(sims_all)
    mean_sim[j,1] <- mean(sims_all)
    sd_sim[j,1] <- sd(sims_all)
    max_sim[j,1] <- max(sims_all)
    min_sim[j,1] <- min(sims_all)}
  
  par(mfrow=c(1,4))
  
  boxplot(mean_sim, 
          ylim = c(quantile(mean_sim, .5)-.02,quantile(mean_sim, .95)+.02))
  title(paste0("mean"), cex = 0.5)
  abline(h=mean(y), col = 'red')
  
  boxplot(sd_sim, 
          ylim = c(quantile(sd_sim, .5)-.1,quantile(sd_sim, .95)+.1))
  title(paste0("SD"), cex = 0.5)
  abline(h=sd(y), col = 'red')
  
  boxplot(max_sim, 
          ylim = c(quantile(max_sim, .5)-.2,quantile(max_sim, .95)+.2))
  title(paste0("max"), cex = 0.5)
  abline(h=max(y), col = 'red')
  
  boxplot(min_sim, 
          ylim = c(quantile(min_sim, .5)-.2,quantile(min_sim, .95)+.2))
  title(paste0("Min"), cex = 0.5)
  abline(h=min(y), col = 'red')
  par(mfrow=c(1,1))

}









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
                        thDelay = thDelay_1st,
                        n = 81)
  sim_pcs[,1] <- sims.setar_1st$serie
  sim_pcs[,1] <- sim_pcs[,1]*sd(pcs_sel[,1]) #Resaling back

  sims.setar_2nd <- setar.sim(B = B_2nd,
                        lag = embd_2nd,
                        nthresh = 1, 
                        Thresh = tail(mod.setar_2nd$coefficients,1), 
                        type = "simul", 
                        thDelay = thDelay_2nd,
                        n = 81)
  sim_pcs[,2] <- sims.setar_2nd$serie
  sim_pcs[,2] <- sim_pcs[,2]*sd(pcs_sel[,2]) #Resaling back

  sims.setar_3rd <- setar.sim(B = B_3rd,
                        lag = embd_3rd,
                        nthresh = 1, 
                        Thresh = tail(mod.setar_3rd$coefficients,1), 
                        type = "simul", 
                        thDelay = thDelay_3rd,
                        n = 81)
  sim_pcs[,3] <- sims.setar_3rd$serie
  sim_pcs[,3] <- sim_pcs[,3]*sd(pcs_sel[,3])
  
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
#pdf(file = "plots/SETAR/Site Mean Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))

for(j in 1:ncol(input_data)) {

  boxplot(mean_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(mean_sim[,j], .5)-.2,quantile(mean_sim[,j], .95)+.2))
  title(paste0("Mean-Site ", j), cex = 0.5)
  abline(h=mean(input_data[,j]), col = 'red')
  
}
#dev.off()

##Plotting the Standard Deviation###
#pdf(file = "plots/SETAR/Site Standard Deviation Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(sd_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(sd_sim[,j], .5)-.2,quantile(sd_sim[,j], .95)+.2))
  title(paste0("SD-Site ", j), cex = 0.5)
  abline(h=sd(input_data[,j]), col = 'red')
  
}
#dev.off()

##Plotting the Max###
#pdf(file = "plots/SETAR/Site Max Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(max_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(max_sim[,j], .5)-.2,quantile(max_sim[,j], .95)+.2))
  title(paste0("Max-Site ", j), cex = 0.5)
  abline(h=max(input_data[,j]), col = 'red')
  
}
#dev.off()

##Plotting the Min###
#pdf(file = "plots/SETAR/Site Min Simulations.pdf")
par(mfrow=c(1,5))
par(mar = c(4, 2, 1.5, 1))
for(j in 1:ncol(input_data)) {
  
  boxplot(min_sim[,j], 
          #ylim = c(quantile(mean_sim[,j], .5)-1,quantile(mean_sim[,j], .95))+1
          ylim = c(quantile(min_sim[,j], .5)-.2,quantile(min_sim[,j], .95)+.2))
  title(paste0("Min-Site ", j), cex = 0.5)
  abline(h=min(input_data[,j]), col = 'red')
  
}
#dev.off()

#####Plotting the Simulations by Sites######
#pdf(file = "plots/SETAR/Site Specific Simulation Moments.pdf")
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

#dev.off()
######################################


#pdf(file = "plots/SETAR/Spatial Distribution of Errors.pdf")

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
par(mfrow = c (1,1), mar = c(4,4,4,1))

N_Sims <- 1000 
wlt <- wavelet(x)
avg_pow_1 <- matrix(NA, nrow = N_Sims, ncol = length(wlt$p.avg))
avg_pow_2 <- matrix(NA, nrow = N_Sims, ncol = length(wlt$p.avg))
avg_pow_3 <- matrix(NA, nrow = N_Sims, ncol = length(wlt$p.avg))

for(j in 1:N_Sims) {
  #Setting up matrix to store the PC powers  
    
  #Simulating the three PCs
  sims.setar_1st <- setar.sim(B = B_1st,
                              lag = embd_1st,
                              nthresh = 1, 
                              Thresh = tail(mod.setar_1st$coefficients,1), 
                              type = "simul", 
                              thDelay = thDelay_1st,
                              n = 81)
  wlt <- wavelet(sims.setar_1st$serie)
  avg_pow_1[j,] <- wlt$p.avg
  
  sims.setar_2nd <- setar.sim(B = B_2nd,
                              lag = embd_2nd,
                              nthresh = 1, 
                              Thresh = tail(mod.setar_2nd$coefficients,1), 
                              type = "simul", 
                              thDelay = thDelay_2nd,
                              n = 81)
  wlt <- wavelet(sims.setar_2nd$serie)
  avg_pow_2[j,] <- wlt$p.avg
  
  sims.setar_3rd <- setar.sim(B = B_3rd,
                              lag = embd_3rd,
                              nthresh = 1, 
                              Thresh = tail(mod.setar_3rd$coefficients,1), 
                              type = "simul", 
                              thDelay = thDelay_3rd,
                              n = 81)
  wlt <- wavelet(sims.setar_3rd$serie)
  avg_pow_3[j,] <- wlt$p.avg
    
}
 
#pdf(file = "plots/SETAR/Global Wavelet Simulations.pdf")
for(i in 1:npcs) {
  x <- scale(pcs_sel[,i])
  Cw=CI(0.9,x,"w");
  C=CI(0.9,x,"r");
  wlt_og=wavelet(x)
  plot(wlt_og$period,wlt_og$p.avg,xlim=c(0,75),
       main=paste0("Global Wavelet Spectrum for PC",i),
       xlab="Period",ylab="Variance", col ='red',
       ylim = c((min(wlt_og$p.avg)-1),(max(wlt_og$p.avg)+7.5))) 
  lines(wlt_og$period,wlt_og$p.avg, col ='red');
  lines(wlt$period,Cw$sig, lty = 2);
  lines(wlt$period,C$sig,col="red", lty = 2);
  avg_pow_matrix <- matrix(NA, nrow = N_Sims, ncol = length(wlt$p.avg))
  if (i ==  1) {
    avg_pow_matrix <- avg_pow_1
  } else if ( i == 2) {
    avg_pow_matrix <- avg_pow_2
  } else  {
    avg_pow_matrix <- avg_pow_3
  }

  lower_percentile <- apply(avg_pow_matrix, 2, function(x) quantile(x, probs=.05))
  upper_percentile <- apply(avg_pow_matrix, 2, function(x) quantile(x, probs=.95))
  median_percentile <- apply(avg_pow_matrix, 2, median)
  lines(wlt_og$period, lower_percentile)
  lines(wlt_og$period, upper_percentile)
  lines(wlt_og$period, median_percentile, lwd = 2)
}
#dev.off()






##########Notes##############
'''For the First PC
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
'''



##################################################









