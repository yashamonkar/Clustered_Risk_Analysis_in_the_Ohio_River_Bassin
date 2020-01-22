#Function to compute the simulation skill.


#Input:- 1. Original DataSet of length T. 
#2. Matrix of NxT simulations.
#3. Time Series Name
#4. How many of the original simulations to be included. 

#Output:- 1. Moments -  Mean, SD, Max, Min. 
#Output:- 2. PDF.
#Output:- 3. CDF.
#Output:- 4. AutoCorrelation Structure.
#Output:- 5. Wavelet Analysis. 
#Output:- 6. Lag-1 and Lag-2 Structure.

get_SimSkill <- function(N_Sims, og_data, Sim_Mat, name, moments, prob_dens,
                         cumm_dens, auto, wavelt, lagged){
  
  #Libraries
  library(ggplot2)
  library(MASS)
  
  #Read the input data
  tx <- og_data
  sims <- Sim_Mat[,1:N_Sims]
  og_name <- name
  
  ####Functions_Needed###########
  #----------------------------------------------------------------------------
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
  
  
  ##Moments:- 
  if(moments == TRUE){
  mean_sim <- sd_sim <- max_sim <- min_sim <- rep(NA, ncol(sims))
  mean_sim <- apply(sims,2,mean)
  sd_sim <- apply(sims,2,sd)
  max_sim <- apply(sims,2,max)
  min_sim <- apply(sims, 2, min)
  
  par(mfrow = c(1,2), mar = c(2,3,4,2))
  boxplot(mean_sim, main = 'Mean', outline = FALSE)
  abline(h=mean(tx),col='red',lwd=2)
  boxplot(sd_sim, main = 'Standard Deviation', outline = FALSE)
  abline(h=sd(tx),col='red',lwd=2)
  mtext(paste0("Moments ", og_name),  side = 3, line = -1.5, outer = TRUE)
  
  boxplot(max_sim, main = 'Maximum', outline = FALSE)
  abline(h=max(tx),col='red',lwd=2)
  boxplot(min_sim, main = 'Minimum', outline = FALSE)
  abline(h=min(tx),col='red',lwd=2)
  mtext(paste0(og_name),  side = 3, line = -1.5, outer = TRUE)
  par(mfrow = c(1,1), mar = c(4,3,4,2))
  }
  
  ##Probability Density
  if(prob_dens == TRUE){
  og_density <- density(tx)
  consolid_sims <- list()
  consolid_points <- list()
  for(j in 1:ncol(sims)) {
    sims_dens <- density(sims[,j])
    consolid_sims[[j]] <- sims_dens$y
    consolid_points[[j]] <- sims_dens$x}
  consolid_points <- unlist(consolid_points)
  consolid_sims <- unlist(consolid_sims)
  consolid_df <- data.frame(x=consolid_points,y=consolid_sims)
  consolid_df$x <- cut(consolid_df$x, breaks = seq(min(tx)-5*sd(tx),max(tx)+5*sd(tx),(max(tx)-min(tx))/20))
  mid_points <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", consolid_df$x) ),
                        upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", consolid_df$x) ))
  consolid_df$x <- rowMeans(mid_points)
  og_df <- data.frame(x1=og_density$x, y1= og_density$y)
    
  plo <- ggplot(og_df, aes(x=x1,y=y1))+
      geom_line(size=1.25)+
      scale_x_continuous(limits=c(min(tx)-2*sd(tx),max(tx)+2*sd(tx))) +
      geom_boxplot(consolid_df, mapping = aes(y = y, x = x,group = x), outlier.shape = NA,outlier.colour = NA)+ 
      ggtitle(paste0("Probability Density of ", og_name)) +
      ylab("Probability Density") +
      theme(plot.title = element_text(size = 10, face = "bold"))
  print(plo)
  #Note:- The warning arises if the values are outside the cut values. 
  }
  
  ##Cumulative Density
  if(cumm_dens == TRUE){
  cdf_og <- ts_eval <- seq(min(tx)-2*sd(tx),max(tx)+2*sd(tx),(max(tx)-min(tx))/20)
  og_cdf <- ecdf(tx)
  for(j in 1:length(ts_eval)) cdf_og[j] <- og_cdf(ts_eval[j])
  sim_cdf <- matrix(NA, ncol = ncol(sims), nrow = length(ts_eval)) #Storing the Simulated CDF's
  for(j in 1:ncol(sims)){
    #Computing each CDF
    cdf_sim <- ecdf(sims[,j])
    for(i in 1:length(ts_eval)) {sim_cdf[i,j] <- cdf_sim(ts_eval[i])}
  }
  #Getting the percentiles
  lower_percentile <- apply(sim_cdf, 1, function(x) quantile(x, probs=.05))
  upper_percentile <- apply(sim_cdf, 1, function(x) quantile(x, probs=.95))
  median_percentile <- apply(sim_cdf, 1, median)
  par(mfrow=c(1,1), mar = c(4,4,3,1))
  plot(ts_eval, cdf_og, type='l',col='red',
       lwd = 2, main = paste0("Simulated CDF for", og_name), 
       xlab = "x", ylab = "CDF")
  polygon(c(ts_eval,rev(ts_eval)),c(lower_percentile,rev(upper_percentile)),col="gray")
  lines(ts_eval, median_percentile, lwd = 2)
  lines(ts_eval, cdf_og, col='red', lwd = 2)
  legend('bottomright', legend = c("Median Simulations", "True CDF","90% Confidence Region"),
         lty = 1, col = c('black','red','grey'), lwd = 2, cex = 0.75)
  }
  
  ##Auto-Correlation Structure
  #Lag-1 and Lag2 structure. 
  if(auto == TRUE){
  og_acf <- acf(tx, plot = FALSE)
  consolid_sims <- list()
  consolid_points <- list()
  acf_significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(tx)))
  
  for(j in 1:ncol(sims)) {
    t_acf <- acf(sims[,j], plot=FALSE)
    consolid_sims[[j]] <- t_acf$acf
    consolid_points[[j]] <- t_acf$lag}
  consolid_points <- unlist(consolid_points)
  consolid_sims <- unlist(consolid_sims)
  consolid_df <- data.frame(x=consolid_points,y=consolid_sims)
  consolid_df$x <- cut(consolid_df$x, breaks = seq(-1,20,1))
  mid_points <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", consolid_df$x) ),
                      upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", consolid_df$x) ))
  consolid_df$x <- rowMeans(mid_points)
  og_df <- data.frame(x1=og_acf$lag, y1= og_acf$acf)
  
  plo <- ggplot(og_df, aes(x=x1,y=y1))+
    geom_point(size=1.25)+
    scale_x_continuous(limits=c(-1,20)) +
    geom_boxplot(consolid_df, mapping = aes(y = y, x = x+0.5,group = x), outlier.shape = NA,outlier.colour = NA)+ 
    ggtitle(paste0("ACF Sims of ", og_name)) +
    ylab("ACF") +
    geom_point(size=1.25)+
    xlab(paste0(og_name)) +
    geom_hline(yintercept=acf_significance_level, linetype="dashed", color = "blue")+
    geom_hline(yintercept=-acf_significance_level, linetype="dashed", color = "blue")+
    theme(plot.title = element_text(size = 10, face = "bold"))
  print(plo)
  }
  
  
  ##Wavelet Strucuture - Low Frequency Structure
  if(wavelt == TRUE){
  wlt_og <- wavelet(tx)
  sim_pow <- matrix(NA, nrow = ncol(sims), ncol = length(wlt_og$p.avg))
  for(j in 1:ncol(sims)) {
    wlt <- wavelet(sims[,j])
    sim_pow[j,] <- wlt$p.avg
  }
  
  Cw=CI(0.9,tx,"w");
  C=CI(0.9,tx,"r");
  par(mar=c(4,4,4,2))
  plot(wlt_og$period,wlt_og$p.avg,xlim=c(0,length(tx)*0.33),
       main=paste0("Global Wavelet Spectrum(Simulated)"),
       xlab="Period",ylab="Variance", col ='red')
  
  lower_percentile <- apply(sim_pow, 2, function(x) quantile(x, probs=.05))
  upper_percentile <- apply(sim_pow, 2, function(x) quantile(x, probs=.95))
  median_percentile <- apply(sim_pow, 2, median)
  low <- cbind(lower_percentile,wlt_og$period)
  high <- cbind(upper_percentile,wlt_og$period)

  polygon(c(wlt_og$period,rev(wlt_og$period)),c(lower_percentile,rev(upper_percentile)),col="gray")
  lines(wlt_og$period, lower_percentile)
  lines(wlt_og$period, upper_percentile)
  lines(wlt_og$period,C$sig, lty=2, col ='red')
  lines(wlt_og$period,Cw$sig, lty=2, col ='black')
  lines(wlt_og$period,wlt_og$p.avg, col ='red', lwd = 2)
  points(wlt_og$period,wlt_og$p.avg, col ='red')
  lines(wlt_og$period, median_percentile, lwd = 2)
  }
  
  
  ##Lagged Plots
  #Lag-1
  #Color Scheme
  if(lagged == TRUE){
  library(RColorBrewer)
  k <- 10
  my.cols <- rev(brewer.pal(k, "RdYlBu"))
  og_lag  <- embed(tx,2)
  t_lat_mat <- matrix(NA, nrow = 1, ncol = 2)
  for(i in 1:ncol(sims)){
    t_lag <- embed(sims[,i],2)
    t_lat_mat <- rbind(t_lat_mat,t_lag)}
  t_lat_mat <- t_lat_mat[complete.cases(t_lat_mat), ]
  z <- kde2d(t_lat_mat[,1],t_lat_mat[,2])
  plot(og_lag,pch=19, cex=.6,main = paste0("Lagged Plot for ", og_name),
       ylab = "Lag-0", xlab = "Lag-1")
  contour(z,  nlevels=k, col=my.cols, add=TRUE)
  
  #Lag-2
  #Color Scheme
  og_lag  <- embed(tx,3); og_lag <- cbind(og_lag[,1],og_lag[,3])
  t_lat_mat <- matrix(NA, nrow = 1, ncol = 2)
  for(i in 1:ncol(sims)){
    t_lag <- embed(sims[,i],3)
    t_lag <- cbind(t_lag[,1],t_lag[,3])
    t_lat_mat <- rbind(t_lat_mat,t_lag)}
  t_lat_mat <- t_lat_mat[complete.cases(t_lat_mat), ]
  z <- kde2d(t_lat_mat[,1],t_lat_mat[,2])
  plot(og_lag,pch=19, cex=.6,main = paste0("Lagged Plot for ", og_name),
       ylab = "Lag-0", xlab = "Lag-2")
  contour(z,  nlevels=k, col=my.cols, add=TRUE)
  }
}
