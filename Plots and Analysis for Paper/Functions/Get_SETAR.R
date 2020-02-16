#Objective of the file is to write a function for the SETAR Model. 

#Input:- 1. Streamflow Field
#2. Site Information.
#3. The Embedding dimensions.  


#Output: Site Diagnostics. 


get_SETAR <- function(Max_Flow, Sites, np, embd, Num_Sims){

  #Install the libraries
  library(dplyr)
  detach("package:dplyr", unload=TRUE)
  #https://stackoverflow.com/questions/36983947/how-to-use-tsdynstar-while-dplyr-is-loaded  
  library(tsDyn)
  
  
  #Read the Data.
  Ann_Max <- Max_Flow
  start_year <- min(Ann_Max$Year);end_year <- max(Ann_Max$Year);
  Ann_Max$Year <- NULL
  Ann_Max <- log(Ann_Max) #Log Transformation to better fit the multi-normal conditions. 
  Site_info <- Sites
  
  #PCA
  data_pca=prcomp(Ann_Max, scale = TRUE, center = TRUE)
  npcs <- np
  pcs_sel <- data_pca$x[,1:npcs]
  
  #SETAR Model.
  embd_dim <- embd #Getting the embedding dimension. 
  N_Sims <- Num_Sims #Simulations for each Dimension
  Sims_Array <- array(NA, dim=c(npcs, N_Sims, dim(Ann_Max)[1])) #Array for storing the simulations
  
  for(i in 1:npcs) {
    x <- scale(pcs_sel[,i])
    selections <- selectSETAR(x, m=embd_dim[i], thDelay=0:(embd_dim[i]-1),
                              d=1,steps=1,plot = FALSE)
    sel_model <- selections$firstBests
    mod.setar <-  setar(x, m = embd_dim[i], thDelay = sel_model[1],
                        th = sel_model[2], mL = sel_model[4],
                        mH = sel_model[5])
    print(summary(mod.setar))
    plot(x,type='l', 
         main = paste0("SETAR Model Fit PC-",i),
         xlab = "Years")
    lines((embd_dim[i]+1):81,mod.setar$fitted.values,type='l',col='red')
    lines((embd_dim[i]+1):81, mod.setar$residuals, col='blue')
    legend(c('bottomright'), legend = c("Fitted","Residials"),
           col=c("red","blue"), lty=1, cex = 0.55)
  
  #Simulating the Time Series
  ml <- as.vector(head(mod.setar$coefficients,sel_model[4]+1)) #Getting the lower regime coefficients
  mh <- as.vector(head(mod.setar$coefficients,-1));mh <- tail(mh,-sel_model[4]-1)  
  ML <- MH <- rep(0,embd_dim[i]+1)
  ML[1:length(ml)] <- ml;MH[1:length(mh)] <- mh
  B <- c(ML,MH) 
  
  for(j in 1:N_Sims){
  sims.setar <- setar.sim(B = B,
                          lag = embd_dim[i],
                          nthresh = 1, 
                          Thresh = tail(mod.setar$coefficients,1), 
                          type = "simul", 
                          thDelay = sel_model[1],
                          n = length(x)+5)
  Sims_Array[i,j,] <-  tail(sims.setar$serie,81)}
  }
  
  #Plotting the Moments. 
  for(i in 1:npcs){
    x <- scale(pcs_sel[,i])
    sims <- as.data.frame(Sims_Array[i,,])
    mean_sims <- apply(sims, 1, mean)
    sd_sims <- apply(sims, 1, sd)
    max_sims <- apply(sims, 1, max)
    min_sims <- apply(sims, 1, min)

    par(mfrow=c(1,4),mar=c(2,2,3,1),oma=c(0.5,0.5,2,1))
    boxplot(mean_sims, main = 'Mean');abline(h = mean(x), lwd = 1.5)
    boxplot(sd_sims, main = 'SD');abline(h = sd(x), lwd = 1.5)
    boxplot(max_sims, main = 'Min');abline(h = max(x), lwd = 1.5)
    boxplot(min_sims, main = 'Max');abline(h = min(x), lwd = 1.5)
    mtext(paste0("PC-",i), outer=TRUE,  cex=1, line=-0.5)
    }
  par(mfrow=c(1,1),mar=c(2,2,2,2),oma=c(0.5,0.5,0.5,0.5))
    
  
  #Plotting the Wavelet Spectrums. 
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
  }  #Functions for the Wavelet Spectrum and Confidence Intervals
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
    
  }  ### Confidence level function
  
  for(i in 1:npcs) {
    x <- scale(pcs_sel[,i]);sims <- as.data.frame(Sims_Array[i,,])
    wlt_og <- wavelet(x)
    Cw=CI(0.9,x,"w");C=CI(0.9,x,"r")
    sim_power <- matrix(NA, nrow = N_Sims, ncol = length(wlt_og$p.avg))
    for(j in 1:N_Sims) {
      x_sim <- t(sims[j,])
      wlt <- wavelet(x_sim)
      sim_power[j,] <- wlt$p.avg } #Saving the power for each simulation
    lower_percentile <- apply(sim_power, 2, function(x) quantile(x, probs=.10))
    upper_percentile <- apply(sim_power, 2, function(x) quantile(x, probs=.9))
    median_percentile <- apply(sim_power, 2, median)
    low <- cbind(lower_percentile,wlt_og$period)
    high <- cbind(upper_percentile,wlt_og$period)
    polys <- as.data.frame(rbind(low,high));colnames(polys) <- c("Variance","Period")
    
    #Plotting
    par(mfrow=c(1,1))
    plot(wlt_og$period,wlt_og$p.avg,xlim=c(0,75),
         main=paste0("Global Wavelet Spectrum(Simulated) for PC",i),
         xlab="Period",ylab="Variance", col ='red',
         ylim = c((min(wlt_og$p.avg)-1),6)) 
    polygon(c(wlt_og$period,rev(wlt_og$period)),c(lower_percentile,rev(upper_percentile)),col="gray")
    lines(wlt_og$period,wlt_og$p.avg, col ='red', lwd = 2)
    lines(wlt_og$period,C$sig, col ='red', lwd = 2, lty=2)
    points(wlt_og$period,wlt_og$p.avg, col ='red')
    lines(wlt_og$period, median_percentile, lwd = 2)

  }

  #Plotting the CDF's
  for(i in 1:npcs) {
   x <- scale(pcs_sel[,i]);sims <- as.data.frame(Sims_Array[i,,])
   x_eval <- as.matrix(seq(-3,3,0.01))
   cdf_og <- ecdf(x);og_cdf <- apply(x_eval, 1, cdf_og)
   avg_cdf <- matrix(NA, nrow = N_Sims, ncol = length(x_eval))
   for(j in 1:N_Sims) {
     cdf_sim <- ecdf(sims[j,])
     avg_cdf[j,] <- apply(x_eval, 1, cdf_sim)}
   lower_percentile <- apply(avg_cdf, 2, function(x) quantile(x, probs=.05))
   upper_percentile <- apply(avg_cdf, 2, function(x) quantile(x, probs=.95))
   median_percentile <- apply(avg_cdf, 2, median)
   
   #Plotting the CDF
   plot(x_eval, og_cdf, type='l',col='red',
        lwd = 2, main = paste0("Simulated CDF Site PC-",i), 
        xlab = "x", ylab = "CDF")
   polygon(c(x_eval,rev(x_eval)),c(lower_percentile,rev(upper_percentile)),col="gray")
   lines(x_eval, median_percentile, lwd = 2)
   lines(x_eval, og_cdf, col='red', lwd = 2)
   legend('bottomright', legend = c("Median Simulations", "True CDF","90% Confidence Region"),
          lty = 1, col = c('black','red','grey'), lwd = 2, cex = 0.75)
   }
  
  
  #Generating the Lag Plots. 
  
  for(i in 1:npcs) {
    library(ggplot2)
    x <- scale(pcs_sel[,i]);sims <- as.data.frame(Sims_Array[i,,])
    for(lg in 1:4){
    lag_og <- cbind(x1=ts(x),x2=lag(ts(x),lg));lag_og <- as.data.frame(lag_og[complete.cases(lag_og), ])
    lag_t <- as.list(rep(NA,5));lead_t <- as.list(rep(NA,5)) 
    for(j in 1:Num_Sims){
      tt <- ts(t(sims[j,]))
      t_temp <- cbind(tt,lag(tt,lg))
      lag_t[[j]] <- t_temp[,1];lead_t[[j]] <- t_temp[,2]}
    lags <- cbind(unlist(lag_t),unlist(lead_t))
    lags <- as.data.frame(lags[complete.cases(lags), ])
    
    #Plotting the Contour Plots. 
    ggplot(lags, aes(x=V1, y=V2) ) +
      geom_density_2d()
    
    lag_plot <- ggplot(lags, aes(x=V1, y=V2)) +
      stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
      xlab("X") + 
      ylab("Lag X") +
      geom_point(lag_og, mapping = aes(x1,x2)) +
      ggtitle(paste0("PC-",i," Lag-",lg))
    print(lag_plot)
    }
  }
   


  
  
  
  
  
  
  
  
  
  
}