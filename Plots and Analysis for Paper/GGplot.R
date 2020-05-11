library(usmap)
library(ggplot2)
library(dplyr)


us <- map_data("state")

site_info <- Site_Info
site_info$Eigenvector <- rnorm(nrow(site_info),0,10)



us_subset <- us %>% filter(region %in% c("ohio","indiana","illinois",'pennsylvania',
                                         'west virginia','virginia','kentucky'))
gg <- ggplot() + 
  geom_map(data=us_subset, map=us,aes(x=long, y=lat, map_id=region),
           fill="#ffffff", color="#000000", size=0.15) +
  geom_point(site_info, mapping = aes(x=dec_long_va, y=dec_lat_va, color = Eigenvector, size = 5)) +
  scale_color_gradient(low="red", high="blue") +
  ggtitle("Stream Gauges with Eigenvectors") +
  xlab("Latitude") + ylab("Longitude")








