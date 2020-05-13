#----------------------------------------------------------------------------#
#The Objective of this file is to generate simulations using the
# Emission Parameters and the Transition Matrix Simulations.
#Unlike the Viterbi Sequence which is the globally most likely sequence
#Here we generate sequences using only the Transition Matrix. 

#Setting Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#Load Modules/Packages/Functions
library(dplyr)
library(ggplot2)
library(gridExtra)
library(zoo)
source("Get_Simulation_Skill.R")

#----------------------------------------------------------------------------#
#Reading the Data. 
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


#----------------------------------------------------------------------------#
#Generating the Simulations from the Bayesian HMM
Sim_Length <- 30 #Number of Years to Simulate the Dataset. 
N_Sims <- nrow(pois_par) #Number of Simulations.
Hidden_States <- Sim_Matrix <- matrix(NA, ncol = N_Sims, nrow = Sim_Length)
K <- ncol(pois_par) #Number of Hidden States

for(i in 1:N_Sims){
  emission_parameters <- pois_par[i,]
  trans_matrix <- matrix(transition_matrix[i,], ncol = K) 
  #Generating the Hidden and Observed Sequence
  z <- rep(NA,Sim_Length)
  z[1] <- sample(1:K, size=1) #Randomly Sample the first state. 
  for(j in 2:Sim_Length){
    z[j] <- sample(1:K,size = 1, prob = trans_matrix[z[j-1],])
  }
  Hidden_States[,i] <- z  
  psi_seq <- sapply(z, function(x){emission_parameters[,x]})
  Sim_Matrix[,i] <- rpois(length(psi_seq), psi_seq)

}

#----------------------------------------------------------------------------#
#Simulations from an I.I.D
IID_Matrix <- matrix(NA, ncol = N_Sims, nrow = Sim_Length)
for(i in 1:N_Sims){
  IID_Matrix[,i] <- rbinom(Sim_Length,28,prob = 0.1)
}


#----------------------------------------------------------------------------#
Data_Matrix <- matrix(NA, ncol = N_Sims, nrow = Sim_Length)
N_sam <- nrow(Count_prox)-Sim_Length
for(i in 1:N_Sims){
  ct <- sample(1:N_sam, 1, replace = TRUE)
  Data_Matrix[,i] <- Count_prox$Count[ct:(ct+Sim_Length-1)]
}

#----------------------------------------------------------------------------#
#Developing Idea for Tests
#1. 20-30 year periods without Presence of State States
#2. Maximum Persistence in States. - Max and Min. 

#This is the unique number of windows for each time series. 
N_window <- 10
Windows <- Sim_Length-N_window+1 



#----------------------------------------------------------------------------#
###Average Wet Year Occurence
#HMM-Model
cur_state <- 3 #Specify the State
avg_occ <- max_occ <- min_occ <- list() 
for(i in 1:N_Sims){ 
  z <- Hidden_States[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
   temp <- z[j:(j+N_window-1)] 
   temp_occ[[j]] <- sum(temp == cur_state )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_matrix <- data.frame(unlist(avg_occ))
colnames(hist_matrix) <- c("Years")
hist_matrix$Model <- "HMM"

#IID-Model
avg_occ <- max_occ <- min_occ <- list() 
thresh <- 9
for(i in 1:N_Sims){ 
  z <- IID_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp > thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_iid_matrix <- data.frame(unlist(avg_occ))
colnames(hist_iid_matrix) <- c("Years")
hist_iid_matrix$Model <- "I.I.D"

#Data
avg_occ <- max_occ <- min_occ <- list() 
thresh <- 9
for(i in 1:N_Sims){ 
  z <- Data_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp > thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_data_matrix <- data.frame(unlist(avg_occ))
colnames(hist_data_matrix) <- c("Years")
hist_data_matrix$Model <- "Obs. Data"


hist_matrix <- rbind(hist_iid_matrix,hist_matrix, hist_data_matrix)

p_wet <- ggplot(hist_matrix, aes(x=Years, fill = Model)) +
  geom_histogram(position = "dodge2", binwidth=0.5,
                 aes(y = 3*(..count..)/sum(..count..)))+
  labs(title = "Average Wet Year Occurence \n for 10-years", 
       x = "Years", y = "Fractional Occurrence") +
  theme(plot.title = element_text(hjust = 0.25))


#----------------------------------------------------------------------------#
###Average Dry Year Occurence
#HMM-Model
cur_state <- 1 
avg_occ <- max_occ <- min_occ <- list()
for(i in 1:N_Sims){ 
  z <- Hidden_States[,i]
  temp_occ <- list()
  for(j in 1:Windows) { #For each window
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp == cur_state)
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_matrix <- data.frame(unlist(avg_occ))
colnames(hist_matrix) <- c("Years")
hist_matrix$Model <- "HMM"

#IID-Model
avg_occ <- max_occ <- min_occ <- list() 
thresh <- 4
for(i in 1:N_Sims){ 
  z <- IID_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp < thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_iid_matrix <- data.frame(unlist(avg_occ))
colnames(hist_iid_matrix) <- c("Years")
hist_iid_matrix$Model <- "I.I.D"
hist_matrix <- rbind(hist_matrix, hist_iid_matrix)


#Data
avg_occ <- max_occ <- min_occ <- list() 
thresh <- 4
for(i in 1:N_Sims){ 
  z <- Data_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp < thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_data_matrix <- data.frame(unlist(avg_occ))
colnames(hist_data_matrix) <- c("Years")
hist_data_matrix$Model <- "Obs. Data"

hist_matrix <- rbind(hist_matrix, hist_iid_matrix, hist_data_matrix)

p_dry <- ggplot(hist_matrix, aes(x=Years, fill = Model)) +
  geom_histogram(position = "dodge2", binwidth=0.5,
                 aes(y = 3*(..count..)/sum(..count..)))+
  labs(title = "Average Dry Year Occurence \n for 10-years", 
       x = "Years", y = "Fractional Occurrence") +
  theme(plot.title = element_text(hjust = 0.25))

#----------------------------------------------------------------------------#
###Computing Persistence for Wet States
#HMM
cur_state <- 3 
max_pers <- list() 
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- Hidden_States[,i]
  temp <- z;temp[temp!=cur_state] <- 0;temp[temp==cur_state] <- 1
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
  pers <- list()
  for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
  max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_matrix <- data.frame(unlist(max_pers))
colnames(hist_matrix) <- c("Years")
hist_matrix$Model <- "HMM"

#I.I.D
max_pers <- list() 
thresh <- 9
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- IID_Matrix[,i]
  temp <- z
  temp[temp <= thresh] <- 0
  temp[temp > thresh] <- 1
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_iid_matrix <- data.frame(unlist(max_pers))
colnames(hist_iid_matrix) <- c("Years")
hist_iid_matrix$Model <- "I.I.D"

#Data
max_pers <- list() 
thresh <- 9
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- Data_Matrix[,i]
  temp <- z
  temp[temp <= thresh] <- 0
  temp[temp > thresh] <- 1
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_data_matrix <- data.frame(unlist(max_pers))
colnames(hist_data_matrix) <- c("Years")
hist_data_matrix$Model <- "Obs. Data"


cons_matrix <- rbind(hist_matrix,hist_iid_matrix, hist_data_matrix)


p_wet_persist <- ggplot(cons_matrix, aes(x=Years, fill = Model)) +
  xlim(-1,7) +
  geom_histogram(position = "dodge2", binwidth=0.5,
                 aes(y = 3*(..count..)/sum(..count..)))+
  labs(title = "Maximum Wet Year Persistence \n for 30-years", 
       x = "Spell Length - Years", y = "Fractional Occurrence") +
  theme(plot.title = element_text(hjust = 0.25))


#----------------------------------------------------------------------------#
###Computing Persistence for Dry States
#HMM
cur_state <- 1 
max_pers <- list() 
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- Hidden_States[,i]
  temp <- z;temp[temp!=cur_state] <- 0;temp[temp==cur_state] <- 1
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_matrix <- data.frame(unlist(max_pers))
colnames(hist_matrix) <- c("Years")
hist_matrix$Model <- "HMM"

#I.I.D
max_pers <- list() 
thresh <- 4
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- IID_Matrix[,i]
  temp <- z
  temp[temp < thresh] <- 1
  temp[temp >= thresh] <- 0
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_iid_matrix <- data.frame(unlist(max_pers))
colnames(hist_iid_matrix) <- c("Years")
hist_iid_matrix$Model <- "I.I.D"


#Model
max_pers <- list() 
thresh <- 4
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- Data_Matrix[,i]
  temp <- z
  temp[temp < thresh] <- 1
  temp[temp >= thresh] <- 0
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_data_matrix <- data.frame(unlist(max_pers))
colnames(hist_data_matrix) <- c("Years")
hist_data_matrix$Model <- "Obs. Data"



cons_matrix <- rbind(hist_matrix,hist_iid_matrix, hist_data_matrix)


p_dry_persist <- ggplot(cons_matrix, aes(x=Years, fill = Model)) +
  xlim(-1,19) +
  geom_histogram(position = "dodge2", binwidth=0.5,
                 aes(y = 3*(..count..)/sum(..count..)))+
  labs(title = "Maximum Dry Year Persistence \n for 30-years", 
       x = "Spell Length - Years", y = "Fractional Occurrence") +
  theme(plot.title = element_text(hjust = 0.25))

pdf("Models_Test_Simulation.pdf")
grid.arrange(p_wet, p_dry,
             p_wet_persist,p_dry_persist,
             ncol=2)
dev.off()




#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#Code on Counts and not States. 


#----------------------------------------------------------------------------#
#Generating the Simulations from the Bayesian HMM
Sim_Length <- 30 
N_Sims <- nrow(pois_par) 
Hidden_States <- Sim_Matrix <- matrix(NA, ncol = N_Sims, nrow = Sim_Length)
K <- ncol(pois_par) 

for(i in 1:N_Sims){
  emission_parameters <- pois_par[i,]
  trans_matrix <- matrix(transition_matrix[i,], ncol = K) 
  z <- rep(NA,Sim_Length)
  z[1] <- sample(1:K, size=1) #Randomly Sample the first state. 
  for(j in 2:Sim_Length){
    z[j] <- sample(1:K,size = 1, prob = trans_matrix[z[j-1],])
  }
  Hidden_States[,i] <- z  
  psi_seq <- sapply(z, function(x){emission_parameters[,x]})
  Sim_Matrix[,i] <- rpois(length(psi_seq), psi_seq)
  
}

#Simulations from an I.I.D
IID_Matrix <- matrix(NA, ncol = N_Sims, nrow = Sim_Length)
for(i in 1:N_Sims){
  IID_Matrix[,i] <- rbinom(Sim_Length,28,prob = 0.1)
}

#Boot Strap from Data Matrix
Data_Matrix <- matrix(NA, ncol = N_Sims, nrow = Sim_Length)
N_sam <- nrow(Count_prox)-Sim_Length
for(i in 1:N_Sims){
  ct <- sample(1:N_sam, 1, replace = TRUE)
  Data_Matrix[,i] <- Count_prox$Count[ct:(ct+Sim_Length-1)]
}

#----------------------------------------------------------------------------#
###Average Wet Year Occurence
#HMM-Model
thresh <- 9
avg_occ <- max_occ <- min_occ <- list() 
for(i in 1:N_Sims){ 
  z <- Sim_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp > thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_matrix <- data.frame(unlist(avg_occ))
colnames(hist_matrix) <- c("Years")
hist_matrix$Model <- "HMM"

#IID-Model
avg_occ <- max_occ <- min_occ <- list() 
for(i in 1:N_Sims){ 
  z <- IID_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp > thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_iid_matrix <- data.frame(unlist(avg_occ))
colnames(hist_iid_matrix) <- c("Years")
hist_iid_matrix$Model <- "I.I.D"

#Data
avg_occ <- max_occ <- min_occ <- list() 
for(i in 1:N_Sims){ 
  z <- Data_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp > thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_data_matrix <- data.frame(unlist(avg_occ))
colnames(hist_data_matrix) <- c("Years")
hist_data_matrix$Model <- "Obs. Data"


hist_matrix <- rbind(hist_iid_matrix,hist_matrix, hist_data_matrix)

p_wet <- ggplot(hist_matrix, aes(x=Years, fill = Model)) +
  geom_histogram(position = "dodge2", binwidth=0.5,
                 aes(y = 3*(..count..)/sum(..count..)))+
  labs(title = "Occurence of 30% Sites crossing exceedances \n across a 10-yr horizon", 
       x = "Years", y = "Fractional Occurrence") +
  theme(plot.title = element_text(hjust = 0.15))



#----------------------------------------------------------------------------#
###Average Dry Year Occurence
#HMM-Model
thresh <- 4
avg_occ <- max_occ <- min_occ <- list() 
for(i in 1:N_Sims){ 
  z <- Sim_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp < thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_matrix <- data.frame(unlist(avg_occ))
colnames(hist_matrix) <- c("Years")
hist_matrix$Model <- "HMM"

#IID-Model
avg_occ <- max_occ <- min_occ <- list() 
for(i in 1:N_Sims){ 
  z <- IID_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp < thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_iid_matrix <- data.frame(unlist(avg_occ))
colnames(hist_iid_matrix) <- c("Years")
hist_iid_matrix$Model <- "I.I.D"
hist_matrix <- rbind(hist_matrix, hist_iid_matrix)


#Data
avg_occ <- max_occ <- min_occ <- list() 
for(i in 1:N_Sims){ 
  z <- Data_Matrix[,i]
  temp_occ <- list()
  for(j in 1:Windows) {
    temp <- z[j:(j+N_window-1)] 
    temp_occ[[j]] <- sum(temp < thresh )
  }
  temp_occ <- unlist(temp_occ)
  max_occ[[i]] <- max(temp_occ)
  min_occ[[i]] <- min(temp_occ)
  avg_occ[[i]] <- temp_occ
}
hist_data_matrix <- data.frame(unlist(avg_occ))
colnames(hist_data_matrix) <- c("Years")
hist_data_matrix$Model <- "Obs. Data"

hist_matrix <- rbind(hist_matrix, hist_iid_matrix, hist_data_matrix)

p_dry <- ggplot(hist_matrix, aes(x=Years, fill = Model)) +
  geom_histogram(position = "dodge2", binwidth=0.5,
                 aes(y = 3*(..count..)/sum(..count..)))+
  labs(title = "Occurence of 10% Sites crossing exceedances \n across a 10-yr horizon", 
       x = "Years", y = "Fractional Occurrence") +
  theme(plot.title = element_text(hjust = 0.15))



#----------------------------------------------------------------------------#
###Computing Persistence for Wet States
#HMM
thresh <- 9
max_pers <- list() 
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- Sim_Matrix[,i]
  temp <- z
  temp[temp <= thresh] <- 0
  temp[temp > thresh] <- 1
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_matrix <- data.frame(unlist(max_pers))
colnames(hist_matrix) <- c("Years")
hist_matrix$Model <- "HMM"

#I.I.D
max_pers <- list() 
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- IID_Matrix[,i]
  temp <- z
  temp[temp <= thresh] <- 0
  temp[temp > thresh] <- 1
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_iid_matrix <- data.frame(unlist(max_pers))
colnames(hist_iid_matrix) <- c("Years")
hist_iid_matrix$Model <- "I.I.D"

#Data
max_pers <- list() 
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- Data_Matrix[,i]
  temp <- z
  temp[temp <= thresh] <- 0
  temp[temp > thresh] <- 1
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_data_matrix <- data.frame(unlist(max_pers))
colnames(hist_data_matrix) <- c("Years")
hist_data_matrix$Model <- "Obs. Data"


cons_matrix1 <- rbind(hist_matrix,hist_iid_matrix, hist_data_matrix)


p_wet_persist <- ggplot(cons_matrix1, aes(x=Years, fill = Model)) +
  xlim(-1,7) +
  geom_histogram(position = "dodge2", binwidth=0.5,
                 aes(y = 3*(..count..)/sum(..count..)))+
  labs(title = "Persistence of more than 10% Sites  \n crossing exceedances for 30-yr horizon", 
       x = "Spell Length - Years", y = "Fractional Occurrence") +
  theme(plot.title = element_text(hjust = 0.15))



#----------------------------------------------------------------------------#
###Computing Persistence for Dry States
#HMM
thresh <- 4
max_pers <- list() 
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- Sim_Matrix[,i]
  temp <- z
  temp[temp < thresh] <- 1
  temp[temp >= thresh] <- 0
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_matrix <- data.frame(unlist(max_pers))
colnames(hist_matrix) <- c("Years")
hist_matrix$Model <- "HMM"

#I.I.D
max_pers <- list() 
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- IID_Matrix[,i]
  temp <- z
  temp[temp < thresh] <- 1
  temp[temp >= thresh] <- 0
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_iid_matrix <- data.frame(unlist(max_pers))
colnames(hist_iid_matrix) <- c("Years")
hist_iid_matrix$Model <- "I.I.D"


#Model
max_pers <- list() 
pb = txtProgressBar(min = 1, max = N_Sims, initial = 1) 
for(i in 1:N_Sims){ 
  setTxtProgressBar(pb,i)
  z <- Data_Matrix[,i]
  temp <- z
  temp[temp < thresh] <- 1
  temp[temp >= thresh] <- 0
  sum_temp <- sum(temp) #This is the maximum possible persistence
  if(sum_temp==0){
    max_pers[[i]] <- 0
  } else { 
    pers <- list()
    for(k in 1:sum_temp){pers[[k]] <- max(rollapply(temp, k, sum))}
    max_pers[[i]] <- max(which(1:sum_temp==unlist(pers)))
  }
}
hist_data_matrix <- data.frame(unlist(max_pers))
colnames(hist_data_matrix) <- c("Years")
hist_data_matrix$Model <- "Obs. Data"



cons_matrix <- rbind(hist_matrix,hist_iid_matrix, hist_data_matrix)


p_dry_persist <- ggplot(cons_matrix, aes(x=Years, fill = Model)) +
  xlim(-1,19) +
  geom_histogram(position = "dodge2", binwidth=0.5,
                 aes(y = 3*(..count..)/sum(..count..)))+
  labs(title = "Persistence of less than 10% Sites  \n crossing exceedances over a 30-yr horizon", 
       x = "Spell Length - Years", y = "Fractional Occurrence") +
  theme(plot.title = element_text(hjust = 0.15))

pdf("Models_Test_Simulation_Counts.pdf")
grid.arrange(p_wet, p_dry,
             p_wet_persist,p_dry_persist,
             ncol=2)
dev.off()