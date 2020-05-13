##Objective:- Relationship between HMM and the Hidden States. 

#Setting the Path
setwd("~/Correlated Risk Analysis/Plots and Analysis for Paper")

#Load Packages
library(ggplot2)
library(gridExtra)
library(dplyr)
library(zoo)



#Reading the Data.
hidden_states <- read.table("Bayesian_Poisson_HMM/Visualizations/Hidden_States.txt", header = TRUE, sep=" ")
amo <- read.table("data/AMO.txt", header = TRUE, sep=" ")
nao <- read.table("data/NAO.txt", header = TRUE, sep=" ")
pdo <-read.table("data/PDO.txt", header = TRUE, sep=" ")
enso <- read.table("data/Nino_34.txt", header = TRUE, sep=" ")


#Converting to Annual Index
amo <- amo %>% group_by(Water_Year) %>% summarise(AMO=mean(AMO))
amo <- amo %>% filter(Water_Year > 1934 & Water_Year < 2020)


#AMO and States
hs1 <- hs2 <- hs3 <- list()
for(i in 1:nrow(hidden_states)){
  temp <- data.frame(AMO=amo$AMO, hs=t(hidden_states[i,]))
  colnames(temp) <- c("AMO","hs")
  temp_seg <- temp %>% filter(hs==1);hs1[[i]] <- temp_seg$AMO
  temp_seg <- temp %>% filter(hs==2);hs2[[i]] <- temp_seg$AMO
  temp_seg <- temp %>% filter(hs==3);hs3[[i]] <- temp_seg$AMO
}


#Plotting the results. 
violin_dataset <- data.frame(S1 = unlist(hs1), ind = 1)
p_dry <- ggplot(violin_dataset, aes(x=ind, y=S1)) +
  geom_violin(trim = TRUE) + labs(y = "AMO", x = " ") +
  ggtitle(c("Dry")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 20, face = "bold"))

violin_dataset <- data.frame(S2 = unlist(hs2), ind = 1)
p_int <- ggplot(violin_dataset, aes(x=ind, y=S2)) +
  geom_violin() + labs(y = "AMO", x = " ") +
  ggtitle(c("Intermediate")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 20, face = "bold"))


violin_dataset <- data.frame(S3 = unlist(hs3), ind = 1)
p_wet <- ggplot(violin_dataset, aes(x=ind, y=S3)) +
  geom_violin() + labs(y = "AMO", x = " ",bw = 4) +
  ggtitle(c("Wet")) + geom_boxplot(width=0.1)+
  theme(plot.title = element_text(size = 20, face = "bold"))

grid.arrange(p_dry, p_int, p_wet,ncol=3)




pdf("Figures/AMO and Hidden States.pdf")
X <- data.frame(X1=unlist(hs3))

g <- ggplot()
g <- g + geom_density(data = X, aes(X1, colour= "Wet"), adjust = 5)

X <- data.frame(X1=unlist(hs2))
g <- g + geom_density(data = X, aes(X1, colour= "Intermediate"), adjust = 5)

X <- data.frame(X1=unlist(hs1))
g <- g + geom_density(data = X, aes(X1, colour= "Dry"), adjust = 5)
g <- g + labs(x = "AMO", x = " ") + ggtitle(c("AMO and Hidden States"))
g
dev.off()










#######################################################################################################################
#Computing the Rolling Means 
library(zoo)
from <- as.Date("1935-01-01")
to <- as.Date("2019-12-31")
years <- seq.Date(from=from,to=to,by="year")
amo.ts <- zoo(amo$AMO, years)

amo.ts <- rollmean(amo.ts, 5)


#AMO and States
hs1 <- hs2 <- hs3 <- list()
pb = txtProgressBar(min = 1, max = nrow(hidden_states), initial = 1) #Counter
for(i in 1:nrow(hidden_states)){
  setTxtProgressBar(pb,i)
  hs <- t(hidden_states[i,])
  hs <- hs[3:83]
  temp <- data.frame(AMO=amo.ts, hs=hs)
  colnames(temp) <- c("AMO","hs")
  temp_seg <- temp %>% filter(hs==1);hs1[[i]] <- temp_seg$AMO
  temp_seg <- temp %>% filter(hs==2);hs2[[i]] <- temp_seg$AMO
  temp_seg <- temp %>% filter(hs==3);hs3[[i]] <- temp_seg$AMO
}


pdf("Figures/AMO and Hidden States.pdf")
X <- data.frame(X1=unlist(hs3))
g <- ggplot()
g <- g + geom_density(data = X, aes(X1, colour= "Wet"), adjust = 7)

X <- data.frame(X1=unlist(hs2))
g <- g + geom_density(data = X, aes(X1, colour= "Intermediate"), adjust = 7)

X <- data.frame(X1=unlist(hs1))
g <- g + geom_density(data = X, aes(X1, colour= "Dry"), adjust = 7)
g <- g + labs(x = "AMO(5-yr Running Mean)", x = " ") + ggtitle(c("AMO and Hidden States")) + 
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size=15, face="bold"),
        axis.title.y = element_text(size=15, face="bold"))
g
dev.off()