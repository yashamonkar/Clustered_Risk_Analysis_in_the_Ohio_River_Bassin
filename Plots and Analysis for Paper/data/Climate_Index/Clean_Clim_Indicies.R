#The objective of this file is to clean the climate indices datasets. 


  
###Loading the requiste directories####
library(dplyr)
library(corrplot)
library(biwavelet)
  
#Reading in the climate indices. 
nino_34 <- read.table("ihadisst1_nino3.4a.dat.txt",  sep ="")
pdo <- read.table("ipdo_ersst.dat.txt", header = TRUE, sep ="", dec=".")
nao <- read.table("inao.dat.txt", header = TRUE, sep ="", dec=".")
Months <- c("Jan","Feb","Mar","April","May","June", "Jul","Aug","Sept","Oct","Nov","Dec")


#Nino
temp <- do.call(rbind, strsplit(as.character(nino_34[,1]),"\\."))
nino_34$V1 <- NULL
colnames(nino_34) <- c("ENSO")
nino_34$Year <- nino_34$Month <- temp[,1] 
nino_34$Month[1:1788] <- rep(Months,149)
nino_34$Month[1789:1798] <- Months[1:10]
write.table(nino_34, "Nino_34.txt", sep = " ")

#NAO
nao_temp <- list()
for(i in 1:nrow(nao)){
  nao_temp[[i]] <- nao[i,2:13]}
nao <- data.frame(Year = rep(1822:2019,each = 12), Month = rep(Months, 198), NAO=unlist(nao_temp))
nao <- head(nao,-3) #Removing the months where we do not have data. 
nao$NAO[nao$NAO == -999.9] <- NA
nao <- nao[complete.cases(nao), ]
write.table(nao, "NAO.txt", sep = " ")

#PDO
pdo_temp <- list()
for(i in 1:nrow(pdo)){
  pdo_temp[[i]] <- pdo[i,2:13]}
pdo <- data.frame(Year = rep(1881:2020,each=12), Month = rep(Months, 140), PDO = unlist(pdo_temp))
pdo$PDO[pdo$PDO == -999.9] <- NA
pdo <- pdo[complete.cases(pdo), ]
write.table(pdo, "PDO.txt", sep = " ")