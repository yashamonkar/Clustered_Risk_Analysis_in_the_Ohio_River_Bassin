#####
#1. The objective of this file is to make predictions/simulations in PC space. 
#2. Use a dimension lowering algorithm - linear PCs
#3. Make predictions and simulations on those PCs. 
#4. Check if the relevant indices and the spectral structure has been replicated. 



###Loading Libraries##########
library(ggplot2)
library(e1071)


###Reading the Data and Initial Manipulation###
input_data <- read.table("data/Max_Annual_Streamflow.txt", sep="", header = TRUE)
######################################2


####Principal Component Analysis#############
data_pca <- prcomp(input_data, scale = TRUE)
var <- cumsum(data_pca$sdev^2)
plot(var/max(var),pch=19, main = "Variance explained by PCs",xlab = "PC's",ylab="Fraction Variance explained")
npcs <- 3 #This is the selected number of PC'S
abline(h = var[npcs]/max(var), lty = 2, col='red')
pcs_sel <- data_pca$x[,1:npcs]
#################################################2


pdf(file = "plots/KNN Method.pdf")
####################PART 2##################

###Provide the Info####
for(jks in 1:ncol(pcs_sel)) {
y <- scale(pcs_sel[,jks]) 
#y <- 1:100
ns <- 1 #One step ahead simulations. 
w = NULL #If NULL the model assumes the weights. 
k <- 5 #Number of Neighbors
lag_loc <- c(1,2,3,4) #These are the lag terms. 
N_Sims <- 500
t_sims <- rep(NA, length(y))
Sims_consolid <- matrix(NA, ncol = N_Sims, nrow = length(t_sims))


for(sim in 1:N_Sims) {
  
  #Converting to the data. 
  y=as.matrix(y)
  #Step 1: Dimension of the feature vector. 
  feature_vectors <- matrix(NA, ncol = length(lag_loc), nrow = length(y)-max(lag_loc))
  #Feature Vectors are the lag terms. 
  for(i in 1:(ncol(feature_vectors)-1)) {
    feature_vectors[,i] <-  tail(head(y,-lag_loc[i]),-max(lag_loc)+lag_loc[i])
  }
  feature_vectors[,ncol(feature_vectors)] <- head(y,-max(lag_loc))
  
  for(jk in 1:length(t_sims)) {
    
    #Step 2: Determine the k-nearest neighbours. 
    if(jk == 1) {
      st <- sample(1:nrow(feature_vectors),1, replace = T)
      xtest <- feature_vectors[st,]}
    
    #Finding and ranking the nearest neighbors. 
    a <- matrix(NA, ncol = ncol(feature_vectors), nrow = nrow(feature_vectors))
    for(i in 1:nrow(a)) a[i,] <- (xtest-feature_vectors[i,])^2
    c <- rowSums(a,na.rm=T)
    
    
    na=rep(NA,1)
    yknn=matrix(NA,nrow=ns,ncol(y))
    yk=na
    yr=seq(1,nrow(feature_vectors))  #This is the row corresponding to the sequence number
    yk[rank(c,ties.method='first')]=yr #yk has the nearest neighbors
    
    #Probabilities
    j=rank(c)		#equal prob for equidistant neighbours
    sj=sum(j[1:k]^(-1))
    pj=(j[1:k]^(-1))/sj #Probability of sampling
    
    
    ynp=sample(yk[1:k],ns,replace=T,prob=pj) #yk[1:k] are the closet neighbours. 
    #ynp <- feature_vectors[ynp,1]+1
    for(p in 1:ncol(y)) yknn[1,p] <- y[ynp+max(lag_loc),p]
      #yknn[,p] <- y[ynp,p]
    #t_sims[jk] <- yknn
    Sims_consolid[jk,sim] <- yknn
    
    #Resetting the xtest
    xtest <- c(yknn, head(xtest,-1))
  }
  
  #Sims_consolid[,sim] <- t_sims
}


###Density Histograms
og_density <- density(y)
consolid_sims <- as.list(1)
consolid_points <- as.list(1)
for(j in 1:ncol(Sims_consolid)) {
  sims_dens <- density(Sims_consolid[,j])
  consolid_sims[[j]] <- sims_dens$y
  consolid_points[[j]] <- sims_dens$x
}
consolid_points <- unlist(consolid_points)
consolid_sims <- unlist(consolid_sims)
consolid_df <- data.frame(x=consolid_points,y=consolid_sims)
consolid_df$x <- cut(consolid_df$x, breaks = seq(-4,4,.1))
mid_points <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", consolid_df$x) ),
                    upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", consolid_df$x) ))
consolid_df$x <- rowMeans(mid_points)
og_df <- data.frame(x1=og_density$x, y1= og_density$y)

plo <- ggplot(og_df, aes(x=x1,y=y1))+
  geom_line(size=1.25)+
  scale_x_continuous(limits=c(-4,4)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_boxplot(consolid_df, mapping = aes(y = y, x = x,group = x), outlier.shape = NA,outlier.colour = NA)+ 
  ggtitle(paste0("Simulated PDF \n for PC ", jks)) +
  xlab(" ") + 
  ylab("Probability Density")
print(plo)


##Computing Quantiles
mean_sim <- colMeans(Sims_consolid)
sd_sim <- apply(Sims_consolid, 2, sd)
skew_sim <- apply(Sims_consolid,2,skewness)




par(mfrow=c(1,3))
boxplot(mean_sim)
title(paste0("mean PC", jks), cex = 0.5)
abline(h=mean(y), col = 'red')

boxplot(sd_sim)
title(paste0("sd PC", jks), cex = 0.5)
abline(h=sd(y), col = 'red')

boxplot(skew_sim)
title(paste0("skew PC", jks), cex = 0.5)
abline(h=skewness(y), col = 'red')
par(mfrow=c(1,1)) 






####Checking the Spectral Signature######
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

wlt <- wavelet(y)
avg_pow <- matrix(NA, nrow = N_Sims, ncol = length(wlt$p.avg))

for(j in 1:N_Sims) {
  wlt_sim <- wavelet(Sims_consolid[,j])
  avg_pow[j,] <- wlt_sim$p.avg
  
}

#pdf(file = "plots/SETAR/Global Wavelet Simulations.pdf")
par(mar=c(4,4,3,1))
  Cw=CI(0.9,y,"w");
  C=CI(0.9,y,"r");
  wlt_og=wavelet(y)
  plot(wlt_og$period,wlt_og$p.avg,xlim=c(0,75),
       main=paste0("Global Wavelet Spectrum(Simulated) for PC ",jks),
       xlab="Period",ylab="Variance", col ='red',
       ylim = c((min(wlt_og$p.avg)-1),8)) 
  
  lower_percentile <- apply(avg_pow, 2, function(x) quantile(x, probs=.05))
  upper_percentile <- apply(avg_pow, 2, function(x) quantile(x, probs=.95))
  median_percentile <- apply(avg_pow, 2, median)
  low <- cbind(lower_percentile,wlt_og$period)
  high <- cbind(upper_percentile,wlt_og$period)
  sims <- as.data.frame(rbind(low,high))
  colnames(sims) <- c("Variance","Period")
  
  
  polygon(c(wlt_og$period,rev(wlt_og$period)),c(lower_percentile,rev(upper_percentile)),col="gray")
  lines(wlt_og$period, lower_percentile)
  lines(wlt_og$period, upper_percentile)
  lines(wlt_og$period,wlt_og$p.avg, col ='red', lwd = 2)
  points(wlt_og$period,wlt_og$p.avg, col ='red')
  lines(wlt_og$period, median_percentile, lwd = 2)

#dev.off()


###########Cumulative Distribution Functions#########

#pdf(file = "plots/SETAR/CDF Simulations.pdf")
x_eval <- seq(-3,3,0.01)
avg_cdf <- matrix(NA, nrow = N_Sims, ncol = length(x_eval))

for(j in 1:N_Sims) {
  #Setting up matrix to store the PC powers  
  p_t <- rep(NA,length(x_eval))
  cdf_eval <- ecdf(Sims_consolid[,j])
  for(k in 1:length(x_eval)) p_t[k] <- cdf_eval(x_eval[k])
  avg_cdf[j,] <- p_t
  
}


  cdf_og <- ecdf(y)
  for(k in 1:length(x_eval)) p_t[k] <- cdf_og(x_eval[k])
  
  
  
  #Getting the percentiles
  lower_percentile <- apply(avg_cdf, 2, function(x) quantile(x, probs=.05))
  upper_percentile <- apply(avg_cdf, 2, function(x) quantile(x, probs=.95))
  median_percentile <- apply(avg_cdf, 2, median)
  
  #Plotting the CDF
  plot(x_eval, p_t, type='l',col='red',
       lwd = 2, main = paste0("Simulated CDF PC-",jks), 
       xlab = "x", ylab = "CDF")
  polygon(c(x_eval,rev(x_eval)),c(lower_percentile,rev(upper_percentile)),col="gray")
  lines(x_eval, median_percentile, lwd = 2)
  lines(x_eval, p_t, col='red', lwd = 2)
  legend('bottomright', legend = c("Median Simulations", "True CDF","90% Confidence Region"),
         lty = 1, col = c('black','red','grey'), lwd = 2, cex = 0.75)
  
}
dev.off()
################################################3




















