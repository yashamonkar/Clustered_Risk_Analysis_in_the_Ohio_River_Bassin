setwd("~/Correlated Risk Analysis/Plots and Analysis for Paper")


get_max_streamflow(lat_min=35, 
                   lat_max=42, 
                   lon_min=-90, 
                   lon_max=-78, 
                   drain_area=5791*0.25,
                   start_date= "1937-01-01",
                   end_date = "2017-12-31",
                   missing_data = 0.001 #Interms of Percent
                   )

#Reading the Input
Ann_Max_Streamflow <- read.table('data/Annual_Maximum_Streamflows.txt', sep=" ", header = TRUE)
Site_Info <- read.table('data/Site_Information.txt', sep=" ", header = TRUE)



#Getting the Space-Time and Time-Space Domains. 
get_PCWavelet(Max_Flow=Ann_Max_Streamflow, 
              Sites=Site_Info, 
              npcs=3)

get_WaveClust(Max_Flow=Ann_Max_Streamflow, 
              Sites=Site_Info)

get_Clim_Connections(Max_Flow=Ann_Max_Streamflow,
                     Sites=Site_Info,
                     npcs=3, 
                     target = c("Feb", "Mar","April"))


##Fitting the SETAR Models. 
get_SETAR(Max_Flow=Ann_Max_Streamflow,
          Sites=Site_Info,
          np=3,
          embd <- c(4,4,1),
          Num_Sims <- 1000)

##Fitting SETAR and getting Site-Replications##
get_SETAR_Sites(Max_Flow=Ann_Max_Streamflow,
                Sites=Site_Info,
                np=3,
                embd <- c(4,4,1),
                Num_Sims <- 1000) 