setwd("~/Correlated Risk Analysis/Plots and Analysis for Paper/Bayesian Principal Components HMM/Visualizations")

library(maps)


#Read the Simulations Data
sim_PC1 <- read.table("Simulations_PC_1.txt", sep=" ", header = TRUE)
sim_PC2 <- read.table("Simulations_PC_2.txt", sep=" ", header = TRUE)
sim_PC3 <- read.table("Simulations_PC_3.txt", sep=" ", header = TRUE)

#Read the true Data.
setwd("~/Correlated Risk Analysis/Plots and Analysis for Paper/Bayesian Principal Components HMM")
input_data <- read.table("data/Annual_Maximum_Streamflows.txt", sep=" ", header = TRUE)
Site_info <- read.table("data/Site_Information.txt", sep =" ", header = TRUE)
years <- input_data$Year;input_data$Year <- NULL 
input_data <- log(input_data) 

#Principal Component Analysis
pcs <- prcomp(input_data, scale = TRUE)
var <- cumsum(pcs$sdev^2)
p_true <- pcs$x[,1:3]

#Check Individual Skill
for(i in 1:3){
  p <- p_true[,i]
  lst <- get(paste0("sim_PC",i))
  cn <- paste0("PC-",i)
  get_SimSkill(og_data = p,
              Sim_Mat = lst,
              name = cn, 
              moments = TRUE,
              prob_dens = TRUE, 
              cumm_dens = TRUE,
              auto = TRUE,
              wavelt = TRUE,
              lagged = TRUE,
              N_Sims = 2000)
}


#Converting the PCs to N-Sites
for(site in 1:ncol(input_data)) {
 t_site <- list()
  for(i in 1:ncol(sim_PC1)) { 
    t_sim <- cbind(sim_PC1[,i],sim_PC2[,i],sim_PC3[,i]) 
    Site_Simulations <- t_sim %*% t(pcs$rotation[,1:3])
    for(j in 1:ncol(Site_Simulations)) {Site_Simulations[,j] <- scale(Site_Simulations[,j], center = FALSE , scale=1/pcs$scale[j]) }
    for(j in 1:ncol(Site_Simulations)) {   Site_Simulations[,j] <- scale(Site_Simulations[,j], center = -1 * pcs$center[j], scale=FALSE)}
    t_site[[i]] <- Site_Simulations[,site] 
  }
 t_site <- do.call(rbind, lapply(t_site, function(i) as.data.frame(t(unlist(i)))))
 t_site <- t(t_site)
 t_site <- exp(t_site)
 p <- input_data[,site]
 cn <- paste0("Site-",site)
 get_SimSkill(og_data = exp(p),
              Sim_Mat = t_site,
              name = cn, 
              moments = FALSE,
              prob_dens = TRUE, 
              cumm_dens = FALSE,
              auto = FALSE,
              wavelt = FALSE,
              lagged = FALSE,
              N_Sims = 2000)
 #plot(1:10, type ='n')
 }



#Visualizing the Good Sites and Bad Sites.
bad_sites <-       #This has to be manually done. 

map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
points(Site_info$dec_long_va,Site_info$dec_lat_va,
       pch=19,
       cex=0.5)
       #col=color.scale(final_sites$drain_area_va,c(1,0.5,0),c(0,0,1),color.spec="rgb"))
title("Error in Streamflow Prediction")
