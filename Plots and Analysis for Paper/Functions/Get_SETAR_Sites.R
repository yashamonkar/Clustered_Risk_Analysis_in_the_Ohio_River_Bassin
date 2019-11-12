#The purpose of this file is to fit SETAR models and compute the
#site specific SETAR performance

#Input:- 1. Annual Maximum Time Series
#2. Site Information
#3. Number of PCs to be included. 
#4. Embedding Dimension.
#5. Number of Simulations


get_SETAR_Sites <- function(Max_Flow, Sites, np, embd, Num_Sims) {
  
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
  data_pca=prcomp(Ann_Max, scale = TRUE)
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
    plot(x,type='l', 
         main = paste0("PC-",i))
    lines((embd_dim[i]+1):81,mod.setar$fitted.values,type='l',col='red')
    lines((embd_dim[i]+1):81, mod.setar$residuals, col='blue')
    
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
  
  #Converting the PC-Simulations to Sites-Simulations. 
  Num_Sites <- dim(Site_Info)[1]
  Site_Sims <- array(NA, dim=c(Num_Sites, N_Sims, dim(Ann_Max)[1]))
  pc_rotations <- data_pca$rotation[,1:npcs] #Sites x Np
  for(j in 1:N_Sims) {
    Single_Sim <- t(Sims_Array[,j,]) #Dimension Years x Np.
    Single_Sites <- t(Single_Sim %*% t(pc_rotations))* data_pca$scale + data_pca$center # Years x Sites
    Site_Sims[,j,] <- Single_Sites 
  }
  
  
  #Plotting the Site Specific Moments
  #Plotting the Moments. 
  for(i in 1:Num_Sites){
    x <- Ann_Max[,i]
    sims <- as.data.frame(Site_Sims[i,,])
    mean_sims <- apply(sims, 1, mean)
    sd_sims <- apply(sims, 1, sd)
    max_sims <- apply(sims, 1, max)
    min_sims <- apply(sims, 1, min)
    
    par(mfrow=c(1,4),mar=c(2,2,3,1),oma=c(0.5,0.5,2,1))
    boxplot(mean_sims, main = 'Mean');abline(h = mean(x), lwd = 1.5)
    boxplot(sd_sims, main = 'SD');abline(h = sd(x), lwd = 1.5)
    boxplot(max_sims, main = 'Min');abline(h = max(x), lwd = 1.5)
    boxplot(min_sims, main = 'Max');abline(h = min(x), lwd = 1.5)
    mtext(paste0("Site-",i), outer=TRUE,  cex=1, line=-0.5)
  }
  par(mfrow=c(1,1),mar=c(2,2,2,2),oma=c(0.5,0.5,0.5,0.5))
  
  
  #Plotting the CDF's
  for(i in 1:Num_Sites) {
    x <- Ann_Max[,i];sims <- as.data.frame(Site_Sims[i,,])
    x_eval <- as.matrix(seq(min(x),max(x),0.1))
    cdf_og <- ecdf(x);og_cdf <- apply(x_eval, 1, cdf_og)
    avg_cdf <- matrix(NA, nrow = N_Sims, ncol = length(x_eval))
    for(j in 1:N_Sims) {
      cdf_sim <- ecdf(sims[j,])
      avg_cdf[j,] <- apply(x_eval, 1, cdf_sim)}
    lower_percentile <- apply(avg_cdf, 2, function(x) quantile(x, probs=.05))
    upper_percentile <- apply(avg_cdf, 2, function(x) quantile(x, probs=.95))
    median_percentile <- apply(avg_cdf, 2, median)
    
    plot(x_eval, og_cdf, type='l',col='red',
         lwd = 2, main = paste0("Simulated CDF Site PC-",i), 
         xlab = "x", ylab = "CDF")
    polygon(c(x_eval,rev(x_eval)),c(lower_percentile,rev(upper_percentile)),col="gray")
    lines(x_eval, median_percentile, lwd = 2)
    lines(x_eval, og_cdf, col='red', lwd = 2)
    legend('bottomright', legend = c("Median Simulations", "True CDF","90% Confidence Region"),
           lty = 1, col = c('black','red','grey'), lwd = 2, cex = 0.75)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}