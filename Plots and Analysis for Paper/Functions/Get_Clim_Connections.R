#The objective of this file is to get the climate teleconnections. 


get_Clim_Connections <- function(Max_Flow, Sites, npcs, target) {
  
  ###Loading the requiste directories####
  library(dplyr)
  library(corrplot)
  library(biwavelet)
  
  #Read the site data.
  Ann_Max <- Max_Flow
  Site_info <- Sites
  #Ann_Max <- Ann_Max_Streamflow #Delete this line
  #Site_info <- Site_Info
  start_year <- min(Ann_Max$Year);end_year <- max(Ann_Max$Year)
  Ann_Max$Year <- NULL
  Ann_Max <- log(Ann_Max)
  
  ####Principal Component Analysis####
  data_pca <- prcomp(Ann_Max, scale = TRUE)
  npcs <- npcs #This is the selected number of PC'S
  pcs_sel <- data_pca$x[,1:npcs] 
  
  #Reading in the climate indices. 
  nino_34 <- read.table("data/Climate_Index/iersst_nino3.4a.dat.txt", header = TRUE, sep ="", dec=".")
  pdo <- read.table("data/Climate_Index/ipdo_a.txt", header = TRUE, sep ="", dec=".")
  nao <- read.table("data/Climate_Index/inao_a.txt", header = TRUE, sep ="", dec=".")
  amo <- read.table("data/Climate_Index/iamo_hadsst_ts_a.txt", header = TRUE, sep ="", dec=".")
  
  #Renaming. 
  colnames(pdo) <- c("Time","PDO_Index")
  colnames(nao) <- c("Time","NAO_Index")
  colnames(nino_34) <- c("Time","Nino_Index")
  colnames(amo) <- c("Time","AMO_Index")
  Months <- c("Jan","Feb","Mar","April","May","June", "Jul","Aug","Sept","Oct","Nov","Dec")
  target <- target
  years <- c(1936, 2018)
  
  #Getting the combined PDO-NAO Index.
  nao_temp <- nao[23:dim(nao)[1],] #Starting from 1824. 
  nao_temp$Month <- rep(Months,195)
  nao_temp$Year <- rep(1824:2018,each= 12)
  nao_temp <- nao_temp %>% filter(Month %in% target)
  nao_temp <- nao_temp %>% group_by(Year) %>% summarise(NAO_Index = mean(NAO_Index))
  nao_temp$NAO_detrended <- nao_temp$NAO_Index - mean(nao_temp$NAO_Index)
  nao_temp <- nao_temp %>% filter(
    Year > years[1],
    Year < years[2]
  )
  
  
  pdo_temp <- pdo[12:dim(pdo)[1],]
  pdo_temp$Month <- rep(Months,117)
  pdo_temp$Year <- rep(1901:2017,each= 12)
  pdo_temp <- pdo_temp %>% filter(Month %in% target)
  pdo_temp <- pdo_temp %>% group_by(Year) %>% summarise(PDO_Index = mean(PDO_Index))
  pdo_temp$PDO_detrended <- pdo_temp$PDO_Index - mean(pdo_temp$PDO_Index)
  pdo_temp <- pdo_temp %>% filter(
    Year > years[1],
    Year < years[2]
  )
  
  nino_temp <- nino_34[12:1967,] #Start from 1855
  nino_temp$Month <- rep(Months,163)
  nino_temp$Year <- rep(1855:2017,each= 12)
  nino_temp <- nino_temp %>% filter(Month %in% target)
  nino_temp <- nino_temp %>% group_by(Year) %>% summarise(ENSO_Index = mean(Nino_Index))
  nino_temp$ENSO_detrended <- nino_temp$ENSO_Index - mean(nino_temp$ENSO_Index)
  nino_temp <- nino_temp %>% filter(
    Year > years[1],
    Year < years[2]
  )
  
  
  
  #Combining them together. 
  climate_indices <- cbind(nino_temp,pdo_temp[,2:3],nao_temp[,2:3])
  climate_indices$PDO_NAO <- climate_indices$NAO_detrended*climate_indices$PDO_Index
  climate_indices$ENSO_NAO <- climate_indices$ENSO_detrended*climate_indices$NAO_detrended
  climate_indices$ENSO_PDO <- climate_indices$ENSO_detrended*climate_indices$PDO_Index
  climate_indices$PDO_detrended <- NULL
  climate_indices$NAO_detrended <- NULL
  climate_indices$Year <- NULL
  climate_indices$ENSO_detrended <- NULL
  colnames(climate_indices) <- c("ENSO","PDO","NAO","PDO_NAO","ENSO_NAO","ENSO_PDO")
  
  climate_indices$ENSO_NAO <- scale(climate_indices$ENSO_NAO)
  climate_indices$PDO_NAO <- scale(climate_indices$PDO_NAO)
  
  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  
  
  for(i in 1:npcs){
    temp <-  cbind(pcs_sel[,i],climate_indices)
    colnames(temp)[1] <- paste0("PC_",i)
    M<-cor(temp)
    p.mat <- cor.mtest(temp)
    corrplot(M, type="upper", order="hclust", 
             p.mat = p.mat, sig.level = 0.05)
  }
  
  
  par(mfrow=c(3,2))
  for(i in 1:npcs) {
    pc <- cbind(1937:2017, pcs_sel[,i])
    for(j in 1:ncol(climate_indices)) {
      clim <- cbind(1937:2017, climate_indices[,j])
      wtc.plt <- wtc(pc, clim) 
      par(mar = c(4,4,2,6))
      plot(wtc.plt, main = paste0("Coherence: PC-",i," & ", colnames(climate_indices)[j]),
           plot.cb = TRUE,plot.phase = TRUE,
           xlab = c("Year"), ylab = c("Period(years)"))
    }
  }
par(mfrow=c(1,1))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
}