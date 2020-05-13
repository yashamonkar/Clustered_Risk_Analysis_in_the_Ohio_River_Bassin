#Script to generate free Simulations. 
#Lenght - 1000 yr. 
#Check - Global Wavelet Spectrum. 

#______________________________________________________________________________#
#Setting Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#Load Pacakges


#Load Data
pois_par <- read.table("Emission_Parameters.txt", sep=" ", header = TRUE)
transition_matrix <- read.table('Transition_Parameters.txt', sep=" ", 
                                header = TRUE)
input_data <- read.table("~/Correlated Risk Analysis/Plots and Analysis for Paper/data/Annual_Maximum_Streamflows.txt", 
                         sep=" ", header = TRUE)



#----------------------------------------------------------------------------#
#Data Wrangling
input_data$Year <- NULL
thres_quantile <- apply( input_data, 2, quantile, 
                         probs = c(0.9), na.rm = TRUE)
Ann_Max_Proxy <- matrix(NA, nrow = nrow(input_data), ncol=ncol(input_data))
for(i in 1:ncol(Ann_Max_Proxy)){
  Ann_Max_Proxy[,i] <- ifelse(input_data[,i]>=thres_quantile[i], 1, 0)
}
ann_count <- rowSums(Ann_Max_Proxy)
Count_prox <- as.data.frame(ann_count)
colnames(Count_prox) <- c("Count")


#______________________________________________________________________________#
###Generate Simulations Using Transition Parameters. 

#Hyperparameters
N_Sim <- nrow(pois_par)/20
N_Length <- 1000 #Years. 
K <- ncol(pois_par)

#Local Storage
Sim_Mat <- matrix(NA,ncol=N_Sim,nrow=N_Length)

for(i in 1:N_Sim){
  emission_parameters <- pois_par[i,]
  trans_matrix <- matrix(transition_matrix[i,], ncol = K) 
  
  z <- rep(NA,N_Length)
  z[1] <- sample(1:K, size=1) 
  for(j in 2:N_Length) {
    z[j] <- sample(1:K,size = 1, prob = trans_matrix[z[j-1],])
  }
  
  psi_seq <- sapply(z, function(x){emission_parameters[,x]})
  Sim_Mat[,i] <- rpois(length(psi_seq), psi_seq)
}


#______________________________________________________________________________#
#Wavelet Function. 
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

##Setup Local Storage
wlt <- wavelet(Sim_Mat[,1])
N_scales <- length(wlt$p.avg)
Wt_Mat <- matrix(NA, nrow = N_scales, ncol=N_Sim)

pb = txtProgressBar(min = 1, max = N_Sim, initial = 1) 
for(i in 1:N_Sim){
  setTxtProgressBar(pb,i)
  wlt <- wavelet(Sim_Mat[,i])
  Wt_Mat[,i] <- wlt$p.avg
}

#______________________________________________________________________________#
#Computing the Confidence Interval over red and white noise over the true data. 
conf <- 0.9
#For the red noise spectrum the AR is fit to the real data. 

#Red Noise Spectrum.
zz=arima(Count_prox$Count, order = c(1, 0, 0))
alpha=zz$coef[1]
print(alpha)

dat <- Count_prox$Count
ps=wlt$period
LP= length(ps)
freq = 1/ps
na=length(Sim_Mat[,i])
CI_red=1:LP    ## confidence interval

for(i in 1:LP){
  
  P=(1-(alpha^2))/(1+(alpha^2)-(2*alpha*cos(2*pi*freq[i])))    # descrete fourier power spectrum page 69 [qn 16] ( torrence and compo)... null hypothesis test
  df=2*sqrt(1+((na/(2.32*ps[i]))^2))
  CI_red[i] =P*(qchisq(conf, df)/df)*var(dat)          #divide by 2 removes power of 2.....for mean no chi dist.[ eqn 17]
}

#White Noise Spectrum
alpha=0

dat <- Count_prox$Count
ps=wlt$period
LP= length(ps)
freq = 1/ps
na=length(Sim_Mat[,i])
CI=1:LP    ## confidence interval

for(i in 1:LP){
  
  P=(1-(alpha^2))/(1+(alpha^2)-(2*alpha*cos(2*pi*freq[i])))    # descrete fourier power spectrum page 69 [qn 16] ( torrence and compo)... null hypothesis test
  df=2*sqrt(1+((na/(2.32*ps[i]))^2))
  CI[i] =P*(qchisq(conf, df)/df)*var(dat)          #divide by 2 removes power of 2.....for mean no chi dist.[ eqn 17]
}


#______________________________________________________________________________#
#Plotting the results
lower_percentile <- apply(Wt_Mat, 1, function(x) quantile(x, probs=.25))
upper_percentile <- apply(Wt_Mat, 1, function(x) quantile(x, probs=.75))
median_percentile <- apply(Wt_Mat, 1, median)

pdf("Spectral Signature.pdf")
par(mar=c(4,4,4,2))
plot(wlt$period,wlt$p.avg,xlim=c(0,N_Length*0.33),
     main=paste0("Spectral Signature of Exteded Simulations - 1000 yr"),
     xlab="Period",ylab="Variance", col ='red',
     ylim = c(min(lower_percentile)*0.5, max(upper_percentile)*1.5),
     type="n")
polygon(c(wlt$period,rev(wlt$period)),c(lower_percentile,rev(upper_percentile)),col="lightgray")
lines(wlt$period, lower_percentile,lwd=2)
lines(wlt$period, upper_percentile,lwd=2)
lines(wlt$period, median_percentile, lwd = 2)
lines(wlt$period, CI_red, lty=2, lwd=3, col = "red")
lines(wlt$period, CI, lty=2, lwd=3)
legend('topright', legend = c("Median Simulations", "50% Confidence Region", "White Noise", "Red Noise"),
       lty = c(1,1,2,2), col = c('black','grey','black', 'red'), 
       lwd = 2, cex = 1.1)
dev.off()