###Script for PC-KNN Methods

#______________________________________________________________________________#
#Setting Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading Packages and Functions.
source("Functions/Get_Simulation_Skill.R")

#Reading Data
input_data <- read.table("data/Annual_Maximum_Streamflows.txt", 
                         sep="", header = TRUE)

#______________________________________________________________________________#
#Data Cleaning
str_yr <- min(input_data$Year);end_yr <-max(input_data$Year)
Ann_Max <- input_data
Ann_Max$Year <- NULL
Ann_Max <- log(Ann_Max)

#Principal Component Analysis
pca <- prcomp(input_data, scale = TRUE)
var <- cumsum(pca$sdev^2)
plot(var/max(var),pch=19, main = "Variance explained by PCs",
     xlab = "PC's",ylab="Fraction Variance explained")


#______________________________________________________________________________#
###Functions
#Embeddings
get_embedding <- function(ts,embd){
  Z <- embed(ts,embd)
  return(Z)
}

#KNN - non parametric nearest neighbor resampling based on training and test data.
knn <- function(y,x,xtest,k,ns,pj,w=NULL){
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
    
    c=rowSums(a,na.rm=T)
    
    yk[rank(c,ties.method='first')]=yr
    
    ynp=sample(yk[1:k],ns,replace=T,prob=pj)
    for(p in 1:ncol(y)) 
      yknn[,p]=y[ynp,p]
    
  }
  return(yknn)
}

#Overall KNN Simulator for a Single Time Series with embeddings
knn_sim <- function(ts,embd,N_Sim,ns,pj,nneib){
  #Local Variables
  N <- length(ts)
 
  #Getting the Data
  Z <- get_embedding(ts,embd+1) #One is added to get Y
  Y <- as.matrix(Z[,1],ncol=1)
  X <- as.matrix(Z[,-1],ncol=embd)
  
  #Weighting Kernel - Parabolic
  pr <- ((1:embd)-embd*0.5)^2
  sj=sum(pr)
  w=((pr))/sj
  w = matrix(w,nrow=1)
  
  #Setting up Global Storage
  tsSim <- matrix(NA, nrow=N, ncol=N_Sim)
  
  for(j in 1:N_Sim){
    #Set up Local Storage
    tsnew <- matrix(NA, nrow=N, ncol=1)
    str_pt <- sample(1:nrow(X),1)
    xtest <- matrix(jitter(X[1,]), nrow = 1)
    tsnew[1,1] <- xtest[1,1]
  
    for(i in 2:N){
      #Getting the Simulations
      yknn <- knn(y=Y,x=X,xtest,k=nneib,ns=ns,pj,w=w)
      tsnew[i,1] <- jitter(sample(yknn,1))
      xtest <- matrix(c(tsnew[i,1],xtest[,-ncol(xtest)]),nrow=1)
      #xtest <- jitter(xtest)
    }
    tsSim[,j] <- tsnew
  }
  return(tsSim)
}

#Hyper-Parameters and Kernels
nneib <- 15 #Number of Neighbors

#Resampling Kernel
sj=sum((1:nneib)^(-1))
pj=((1:nneib)^(-1))/sj



#Running the Simulation.
Sim_data <- knn_sim(ts=pca$x[,1],
        embd=7,
        N_Sim=1000,
        ns=40,
        pj=pj,
        nneib=nneib)

#Removing Simulations which did not enter required trajectory. 
max_col <- apply(Sim_data, 2, max)
sub_sim <- Sim_data[,which(max_col>9.45)]
dim(sub_sim)

pdf("Figures/KNN_Simulation.pdf")
#Checking Simulation Skill
get_SimSkill(og_data = pca$x[,1],
             Sim_Mat = sub_sim,
             name = "PC-1", 
             moments = TRUE,
             prob_dens = TRUE, 
             cumm_dens = TRUE,
             auto = TRUE,
             wavelt = TRUE,
             lagged = TRUE,
             N_Sims = ncol(sub_sim))
dev.off()
