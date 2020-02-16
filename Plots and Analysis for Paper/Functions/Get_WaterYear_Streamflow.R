#Function to get the data from Data-Retrival and to transform it

#Input:- Lat-Lon Boxes. 
#2. Drainage Area

get_max_wateryear_streamflow <- function(lat_min, lat_max, lon_min, lon_max, drain_area,
                               start_date,end_date,missing_data){

#Output:- 1. Site-Specific Annual Maximums
#2. Streamflows selected. 
#3. Monthly Distribution of Annual Maximums. 


#Libraries
library(dplyr)
library(ggplot2)
library(dataRetrieval)
library(maps)


####Getting the Sites### 
  
#Divinding into boxes of smaller lat-lon dimension. (Done since query size can be only so big)
lat_min <- lat_min #lat_min
lat_max <- lat_max #lat_max
lon_min <- lon_min#lon_min
lon_max <- lon_max #lon_max
lats <- lat_max-lat_min 
lons <- lon_max-lon_min

#Dates
start_date <- start_date
end_date <- end_date    
#ParameterCD is discharge in cubic ft/sec and is daily mean. Source:- https://help.waterdata.usgs.gov/code/parameter_cd_nm_query?parm_nm_cd=%25discharge%25&fmt=html
#hasDataTypeCd tells how is the parameter of interest measured. dv means daily. Source:- https://waterservices.usgs.gov/rest/Site-Service.html#outputDataTypeCd
sites_info <- list()
for(i in 1:lats) {
 lat_mi <- lat_min+i-1;lat_ma <- lat_mi+1
  for(j in 1:lons) {
   lon_mi <- lon_min+j-1;lon_ma <- lon_mi+1
   ct <- lons*(i-1)+j
   sites_info[[ct]] <- whatNWISsites(bBox = c(lon_mi, lat_mi, lon_ma, lat_ma), parameterCd = c("00060"), hasDataTypeCd = "dv")
}}
sites_info <- do.call(rbind, sites_info)
print("All Sites in region of interest identified")

#Getting site-data and segregating by drainage area
site_INFO <- readNWISsite(sites_info$site_no)
drain_area <- drain_area #drain_area #This is in miles squared. 
req_sites <- site_INFO %>% filter(drain_area_va > (drain_area)) #Criteria in Dave's Paper greater than 
req_sites$Record_Length <- NA 
req_sites$Flow_Change <- NA
req_sites$SUM_NA <- NA

#Consistency Check One. 
map('state')
title("Consistency Check - Sites with needed Drainage Area")
points(req_sites$dec_long_va,req_sites$dec_lat_va,pch=19,cex=0.5)
lines(rep(lon_min,lats+1), lat_max:lat_min,  col ='red', lwd = 2)
lines(rep(lon_max,lats+1), lat_max:lat_min,  col ='red', lwd = 2)
lines(lon_max:lon_min, rep(lat_min,lons+1),  col ='red', lwd = 2)
lines(lon_max:lon_min, rep(lat_max,lons+1),  col ='red', lwd = 2)
box()


#####Subsetting based on data-availability####

#Getting the Record Length for period of interest.
#Notes on reading data from readNWISdv:- It gives the daily value. 
#Flow_Cd gives a measure of the quality of the data. A -> Approved. A_e - > Approved Estimated.
#Soure:- https://owi.usgs.gov/R/training-curriculum/intro-curriculum/Analyze/
#Source:- https://help.waterdata.usgs.gov/codes-and-parameters/instantaneous-value-qualification-code-uv_rmk_cd
startDate <- start_date
endDate <- end_date
pb = txtProgressBar(min = 1, max = dim(req_sites)[1], initial = 1)
for( i in 1:dim(req_sites)[1]) {
  setTxtProgressBar(pb,i)
  parameterCd <- "00060"  # Discharge. Refer above for source. Daily Mean discharge in cfs
  startDate <- start_date
  endDate <- end_date
  statCd = "00003" #Deafult Implies Daily Mean. Here we use the readNWISdv data. 
  discharge <- readNWISdv(req_sites$site_no[i], parameterCd, startDate, endDate, statCd)
  discharge <- renameNWISColumns(discharge)
  req_sites$Record_Length[i] <- dim(discharge)[1]
  req_sites$Flow_Change[i] <- ncol(discharge)
  if(ncol(discharge) < 4) {  
    req_sites$SUM_NA[i] <- 100
  } else {
    req_sites$SUM_NA[i] <- sum(is.na(discharge[,4]))
  }
}

##Subsetting if they have the data. 
max_days <- as.numeric(as.Date(endDate,"%Y-%m-%d")-as.Date(startDate,"%Y-%m-%d"))+1
days_needed <- max_days*(100-missing_data)/100
final_sites <- req_sites %>% filter(Record_Length > days_needed)
print("The Sites with the appropritate Data have been selected")

#Consistency Check Two. 
map('state')
title("Consistency Check - Sites with needed Data")
points(final_sites$dec_long_va,final_sites$dec_lat_va,pch=19,cex=0.5)
lines(rep(lon_min,lats+1), lat_max:lat_min,  col ='red', lwd = 2)
lines(rep(lon_max,lats+1), lat_max:lat_min,  col ='red', lwd = 2)
lines(lon_max:lon_min, rep(lat_min,lons+1),  col ='red', lwd = 2)
lines(lon_max:lon_min, rep(lat_max,lons+1),  col ='red', lwd = 2)
box() 
  

#Plotting it on the Watershed
library(ggplot2)
library(sf)
library(rgeos)
library(sp)
library(rgdal)
orb=st_read("Water_Year_Data/ORBShape/ORBShape/WBDHU6.shp")%>%st_transform(crs="+proj=longlat +datum=NAD83")

ORB_gauges <- ggplot() + 
  geom_sf(data = orb, size = 0.1, color = "black", fill = "cyan1") + 
  ggtitle("Ohio River Basin all sites") + 
  coord_sf() +
  xlab("Longitude")+
  ylab("Latitude")+
  geom_point(data=final_sites, aes(x=final_sites$dec_long_va,y=final_sites$dec_lat_va))
print(ORB_gauges) 
  
####Finding the intersection in the dataframe###
gauges=data.frame(
  lon=final_sites$dec_long_va,
  lat=final_sites$dec_lat_va,
  site=final_sites$site_no)
coordinates(gauges)=~lon+lat
  
#Reading the RiverBasin Map  
ORB.map <- readOGR("Water_Year_Data/ORBShape/ORBShape/WBDHU6.shp")  
proj4string(ORB.map) <- proj4string(gauges) #Getting them to the same projections

sites_present <- over(gauges, ORB.map)
final_sites$Presence <- sites_present$OBJECTID
final_sites <- final_sites[!is.na(final_sites$Presence), ]

ORB_gauges <- ggplot() + 
  geom_sf(data = orb, size = 0.1, color = "black") + 
  ggtitle("Ohio River Basin") + 
  coord_sf() +
  xlab("Longitude")+
  ylab("Latitude")+
  geom_point(data=final_sites, aes(x=final_sites$dec_long_va,y=final_sites$dec_lat_va))
print(ORB_gauges) 
print(paste0("The number of gauges in the site are ", dim(final_sites)[1]))

#######Extracting the Annual Maximum at each site########
num_sites <- dim(final_sites)[1]
Dataset_Years <- format(as.Date(startDate, format="%Y-%m-%d"),"%Y"):format(as.Date(endDate, format="%Y-%m-%d"),"%Y")
Dataset_Years <- tail(Dataset_Years,-1)
Ann_Max_Streamflow <- as.data.frame(matrix(NA,nrow=length(Dataset_Years),ncol = num_sites))

pb = txtProgressBar(min = 1, max = num_sites, initial = 1) 
for(i in 1:num_sites){
  setTxtProgressBar(pb,i)
  parameterCd <- "00060"  # Discharge ft3/sec
  startDate <- start_date
  endDate <- end_date
  statCd = "00003"
  discharge <- readNWISdv(final_sites$site_no[i], parameterCd, startDate, endDate, statCd)
  discharge <- renameNWISColumns(discharge)
  discharge$Year <- as.numeric(format(as.Date(discharge$Date, format="%Y-%m-%d"),"%Y"))
  discharge$Month <- as.numeric(format(as.Date(discharge$Date, format="%Y-%m-%d"),"%m"))
  
  #Function converting to Water Year based on Month. 
  water_year <- function(x){
    if(x > 9) {1
    } else {0}
  }
  discharge$Water_Year <- discharge$Year + apply(as.array(discharge$Month), MARGIN = 1, FUN = water_year)
  discharge <-  discharge %>% group_by(Water_Year) %>% summarise(Ann_Max = max(Flow))
  Ann_Max_Streamflow[,i] <- discharge$Ann_Max
}
Ann_Max_Streamflow$Year <- Dataset_Years

write.table(final_sites, 'Water_Year_Data/Site_Information.txt', sep=" ")
write.table(Ann_Max_Streamflow, 'Water_Year_Data/Annual_Maximum_Streamflows.txt', sep=" ")


return(final_sites)
return(Ann_Max_Streamflow)
}