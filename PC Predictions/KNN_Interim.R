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



####################PART 2##################
#Here the innovation includes a varying tau. 
#xt = f(xt-1,xt-tau) and so on. 

###Provide the Info####
y <- scale(pcs_sel[,3])
ns <- 1 #One step ahead simulations. 
w = NULL
k <- 40 #Number of Neighbors
lag_loc <- c(1,2,3,4) #These are the lag terms. 
N_Sims <- 100
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
      xtest <- feature_vectors[st,]*1.1}
    
    #Finding and ranking the nearest neighbors. 
    a <- matrix(NA, ncol = ncol(feature_vectors), nrow = nrow(feature_vectors))
    for(i in 1:nrow(a)) a[i,] <- (xtest-feature_vectors[i,])^2
    c <- rowSums(a,na.rm=T)
    
    
    na=rep(NA,1)
    yknn=matrix(NA,nrow=ns,ncol(y))
    yk=na
    yr=seq(1,nrow(feature_vectors))
    yk[rank(c,ties.method='first')]=yr
    
    #Probabilities
    j=rank(c)		#equal prob for equidistant neighbours
    sj=sum(j[1:k]^(-1))
    pj=(j[1:k]^(-1))/sj #Probability of sampling
    
    
    ynp=sample(yk[1:k],ns,replace=T,prob=pj) #yk[1:k] are the closet neighbours. 
    for(p in 1:ncol(y)) 
      yknn[,p] <- y[ynp,p]
    t_sims[jk] <- yknn
    
    #Resetting the xtest
    xtest <- c(yknn, head(xtest,-1))
  }
  
  Sims_consolid[,sim] <- t_sims
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
  ggtitle(paste0("Simulated PDF \n for PC ")) +
  xlab(" ") + 
  ylab("Probability Density")
print(plo)


##Computing Quantiles
mean_sim <- colMeans(Sims_consolid)
sd_sim <- apply(Sims_consolid, 2, sd)
max_sim <- apply(Sims_consolid,2,max)
min_sim <- apply(Sims_consolid,2,min)



par(mfrow=c(1,4))
boxplot(mean_sim)
title(paste0("mean"), cex = 0.5)
abline(h=mean(y), col = 'red')

boxplot(sd_sim)
title(paste0("SD"), cex = 0.5)
abline(h=sd(y), col = 'red')

boxplot(max_sim)
title(paste0("Max"), cex = 0.5)
abline(h=max(y), col = 'red')

boxplot(min_sim)
title(paste0("Min"), cex = 0.5)
abline(h=min(y), col = 'red')
par(mfrow=c(1,1))