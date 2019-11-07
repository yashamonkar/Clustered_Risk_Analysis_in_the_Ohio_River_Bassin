#Function to compute the WaveClust in ORB. 


#Input:- 1. Site Specific Annual Maximum
#2. Site Information

#Output:- 1. Visualizations of the WaveletClusters. 


get_WaveClust <- function(Max_Flow, Sites) {
  
  #Read the site data.
  Ann_Max <- Max_Flow
  start_year <- min(Ann_Max$Year);end_year <- max(Ann_Max$Year)
  Ann_Max$Year <- NULL
  Ann_Max <- log(Ann_Max)
  Site_info <- Sites
  
  ####Functions_Needed###########
  wavclust=function(xo,yr1,coord,nam){
    # nam  = name for series, e.g. ANnual Max flow
    #coord are lat long coords as a data frame and need to be in coords$Lon and coords$Lat
    #if maps are not desired then please modify this code to drop that part
    #if maps outside the USA are desired then modify the map function below
    #name is a character string that gives the name of the series being analyzed
    #yr1 is the start year of the data
    nr=nrow(xo)
    nc=ncol(xo)
    yr2=yr1+nr-1
    library("biwavelet")
    library("plotrix")
    library("maps")
    library("stats")
    np=length(wavelet(xo[,1])$p.avg)
    #above gets number of wavelet frequencies
    w.arr=array(NA,dim=c(nc,np))
    #creates place to store global wavelet for each series
    #standardize all the series so the global wavelet is comparable
    x=xo
    for (i in 1:nc)x[,i]=(xo[,i]-mean(xo[,i]))/sd(xo[,i])
    #compute global wavelet spectra for standardized series
    for ( i in 1:nc)w.arr[i,]=wavelet(x[,i])$p.avg
    #plot clustering
    par(mfrow=c(1,1))
    clus=hclust(dist(w.arr),method="ward.D2")
    plot(clus)
    (cls <- identify(clus))
    nclus=length(cls)
    ju=rep(NA,nc)
    for (i in 1:nclus)ju[cls[[i]]]=rep(i,length(cls[[i]]))
    par(mfrow=c(1,1))
    map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
    points(coord$lon,coord$lat,pch=19,cex=1,col=ju)
    title("Spatial Distribution of Clusters")
    
    #npcs is number of pcs to look at, x is the data matrix (time*spatial locs), 
    #coords are lat long coords as a data frame and need to be in coords$Lon and coords$Lat
    #if maps are not desired then please modify this code to drop that part
    #if maps outside the USA are desired then modify the map function below
    #name is a character string that gives the name of the series being analyzed
    #yr1 is the start year of the data
    #clp=matrix(NA,nrow=nr,ncol=nclus)
    
    for (i in 1:nclus){
      tpd=xo[,cls[[i]]];
      #change above to put the cluster members into a temp array
      # then do pca on the temp array and clp[,i] = pc$scores for the cluster i
      pcw=prcomp(tpd, scale = TRUE)
      p=pcw$x[,1]
      #get wavelet of first pc            
      wlt=wavelet(p);
      Cw=CI(0.9,p,"w");
      C=CI(0.9,p,"r");
      nam1=paste("Wavelet Cluster",nam);
      #nam1=paste("Wavelet Cluster",as.character(i),nam);
      par(mfrow=c(2,2))
      par(mar=c(4,2,3,0.5))
      plot(yr1:yr2,p,xlab="Year",ylab=nam,main="Cluster PC 1");
      lines(lowess(yr1:yr2,p,f=1/9),lwd=2,col="red");
      plot(wlt$period,wlt$p.avg,xlim=c(0,32),main="Global Wavelet Spectrum",xlab="Period",ylab="Variance"); lines(wlt$period,wlt$p.avg);
      lines(wlt$period,Cw$sig);
      lines(wlt$period,C$sig,col="red");
      wt1=wt(cbind(yr1:yr2,p));plot(wt1, type="power.corr.norm", xlab="Year",main=nam1);
      
      jj=rep(1,nc)
      jj[ju==i]=19
      drainage_areas <- (Site_info$drain_area_va/max(Site_info$drain_area_va))+1
      jc=rep(0,nc);jc[cls[[i]]]=abs(pcw$rotation[,1])
      par(mar=c(3,3,3,0))
      map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
      title("Stream Gauges with Loadings")
      points(coord$lon,coord$lat,pch=jj,cex=drainage_areas,col=color.scale(jc,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
      
      
    }
    par(mfrow=c(1,1))
  }
  
  
  
  ### Wavelet Function
  
  #WAVELET  1D Wavelet transform with optional singificance testing
  #
  #   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
  #
  #   Computes the wavelet transform of the vector Y (length N),
  #   with sampling rate DT.
  #
  #   By default, the Morlet wavelet (k0=6) is used.
  #   The wavelet basis is normalized to have total energy=1 at all scales.
  #
  #
  # INPUTS:
  #
  #    Y = the time series of length N.
  #    DT = amount of time between each Y value, i.e. the sampling time.
  #
  # OUTPUTS:
  #
  #    WAVE is the WAVELET transform of Y. This is a complex array
  #    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
  #    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
  #    The WAVELET power spectrum is ABS(WAVE)^2.
  #    Its units are sigma^2 (the time series variance).
  #
  #
  # OPTIONAL INPUTS:
  # 
  # *** Note *** setting any of the following to -1 will cause the default
  #               value to be used.
  #
  #    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
  #         N up to the next higher power of 2. This prevents wraparound
  #         from the end of the time series to the beginning, and also
  #         speeds up the FFT's used to do the wavelet transform.
  #         This will not eliminate all edge effects (see COI below).
  #
  #    DJ = the spacing between discrete scales. Default is 0.25.
  #         A smaller # will give better scale resolution, but be slower to plot.
  #
  #    S0 = the smallest scale of the wavelet.  Default is 2*DT.
  #
  #    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
  #        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
  #
  #    MOTHER = the mother wavelet function.
  #             The choices are 'MORLET', 'PAUL', or 'DOG'
  #
  #    PARAM = the mother wavelet parameter.
  #            For 'MORLET' this is k0 (wavenumber), default is 6.
  #            For 'PAUL' this is m (order), default is 4.
  #            For 'DOG' this is m (m-th derivative), default is 2.
  #
  #
  # OPTIONAL OUTPUTS:
  #
  #    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
  #           to the SCALEs.
  #
  #    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
  #            where J1+1 is the total # of scales.
  #
  #    COI = if specified, then return the Cone-of-Influence, which is a vector
  #        of N points that contains the maximum period of useful information
  #        at that particular time.
  #        Periods greater than this are subject to edge effects.
  #        This can be used to plot COI lines on a contour plot by doing:
  #
  #              contour(time,log(period),log(power))
  #              plot(time,log(coi),'k')
  #
  #----------------------------------------------------------------------------
  #   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
  #
  #   This software may be used, copied, or redistributed as long as it is not
  #   sold and this copyright notice is reproduced on each copy made. This
  #   routine is provided as is without any express or implied warranties
  #   whatsoever.
  #
  # Notice: Please acknowledge the use of the above software in any publications:
  #    ``Wavelet software was provided by C. Torrence and G. Compo,
  #      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
  #
  # Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
  #            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
  #
  # Please send a copy of such publications to either C. Torrence or G. Compo:
  #  Dr. Christopher Torrence               Dr. Gilbert P. Compo
  #  Research Systems, Inc.                 Climate Diagnostics Center
  #  4990 Pearl East Circle                 325 Broadway R/CDC1
  #  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
  #  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
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
  
  ####Calling the function####
  coords <- data.frame(lat = Site_info$dec_lat_va, lon = Site_info$dec_long_va)
  nam <- c("Ann_Max")
  wavclust(xo=Ann_Max, yr1=start_year, coord = coords, nam = nam)  
  
  
  
}
