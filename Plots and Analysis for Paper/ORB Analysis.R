setwd("~/Correlated Risk Analysis/Plots and Analysis for Paper")

#Get Calendar Year Streamflow Data
get_max_streamflow(lat_min=35, 
                   lat_max=42, 
                   lon_min=-90, 
                   lon_max=-78, 
                   drain_area=5791*0.25,
                   start_date= "1937-01-01",
                   end_date = "2017-12-31",
                   missing_data = 0.001 #Interms of Percent
                   )



#Get Water Year Streamflow Data
source("Functions/Get_WaterYear_Streamflow.R")
get_max_wateryear_streamflow(lat_min=35, 
                   lat_max=42, 
                   lon_min=-90, 
                   lon_max=-78, 
                   drain_area=5791*0.25,
                   start_date= "1934-10-01",
                   end_date = "2019-09-30",
                   missing_data = 0.1 #Interms of Percent
)

#Reading the Input
Ann_Max_Streamflow <- read.table('Water_Year_Data/Annual_Maximum_Streamflows.txt', sep=" ", header = TRUE)
Site_Info <- read.table('Water_Year_Data/Site_Information.txt', sep=" ", header = TRUE)



#Getting the Space-Time and Time-Space Domains. 
source("Functions/Get_PCWavelet.R")
get_PCWavelet(Max_Flow=Ann_Max_Streamflow, 
              Sites=Site_Info, 
              npcs=3)

source("Functions/Get_WaveClust.R")
get_WaveClust(Max_Flow=Ann_Max_Streamflow, 
              Sites=Site_Info)

pdf("Climate_Connections.pdf")
source("Functions/Get_Clim_Connections.R")
get_Clim_Connections(Max_Flow=Ann_Max_Streamflow,
                     Sites=Site_Info,
                     npcs=3, 
                     target = c("Feb","Mar","April")) #Just check how the months are written. E.g Mar instead of March. 
dev.off()
























##Fitting the SETAR Models. 
get_SETAR(Max_Flow=Ann_Max_Streamflow,
          Sites=Site_Info,
          np=3,
          embd <- c(4,1,1),
          Num_Sims <- 1000)

##Fitting SETAR and getting Site-Replications##
get_SETAR_Sites(Max_Flow=Ann_Max_Streamflow,
                Sites=Site_Info,
                np=3,
                embd <- c(4,4,1),
                Num_Sims <- 1000) 