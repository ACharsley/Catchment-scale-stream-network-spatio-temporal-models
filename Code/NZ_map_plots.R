
rm(list=ls())

################
## Directories
################

NZdata_dir <- file.path(getwd(), "NZ_data")

fig_dir <- file.path(NZdata_dir, "Figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

library(tidyverse)
library(proj4)
library(akima)
library(grid)
library(gridExtra)

################
## Load data
################

## Network data ##
netfull <- readRDS(file.path(NZdata_dir, "NZ_network.rds"))

#Waikato network
net_waikato <- readRDS(file.path(file.path(getwd(), "Waikato", "Data"), "Greater_waikato_network_AC.rds"))
net_waikato$Study_area <- "Greater Waikato region"


net_waikato_connected <- data.frame('parent_s'=net_waikato$parent_s, 'child_s'=net_waikato$child_s, "Lat"=net_waikato$lat, "Lon"=net_waikato$long)
net_waikato_connected <- lapply(1:nrow(net_waikato_connected), function(x){
  parent <- net_waikato_connected$parent_s[x]
  find <- net_waikato_connected %>% filter(child_s == parent)
  if(nrow(find)>0) out <- cbind.data.frame(net_waikato_connected[x,], 'Lon2'=find$Lon, 'Lat2'=find$Lat)
  if(nrow(find)==0) out <- cbind.data.frame(net_waikato_connected[x,], 'Lon2'=NA, 'Lat2'=NA)
  return(out)
})
net_waikato_connected <- do.call(rbind, net_waikato_connected)




#Waitaki network
load(file.path(getwd(), "Waitaki", "Data", "nz_waitaki_longfin_eel.rda"))
net_waitaki <- nz_waitaki_longfin_eel[["network"]]
net_waitaki$Study_area <- "Waitaki catchment"


net_waitaki_connected <- data.frame('parent_s'=net_waitaki$parent_s, 'child_s'=net_waitaki$child_s, "Lat"=net_waitaki$lat, "Lon"=net_waitaki$long)
net_waitaki_connected <- lapply(1:nrow(net_waitaki_connected), function(x){
  parent <- net_waitaki_connected$parent_s[x]
  find <- net_waitaki_connected %>% filter(child_s == parent)
  if(nrow(find)>0) out <- cbind.data.frame(net_waitaki_connected[x,], 'Lon2'=find$Lon, 'Lat2'=find$Lat)
  if(nrow(find)==0) out <- cbind.data.frame(net_waitaki_connected[x,], 'Lon2'=NA, 'Lat2'=NA)
  return(out)
})
net_waitaki_connected <- do.call(rbind, net_waitaki_connected)


#Join
net_waikato_waitaki <- rbind(net_waikato, net_waitaki)
####



## Observations ##
obsfull <- readRDS(file.path(NZdata_dir, "NZ_observations.rds"))

#Waikato data
obs_waikato <- obsfull %>% filter(data_type == "count")
obs_waikato$Study_area <- "Greater Waikato region"

#Waitaki data
obs_waitaki <- obsfull %>% filter(grepl("aitaki", CatName))
obs_waitaki$Study_area <- "Waitaki catchment"

#Join 
obs_waikato_waitaki <- rbind(obs_waikato, obs_waitaki)


#Full NZFFD data
obs_NZFFD <- obsfull %>% 
  filter(data_type == "encounter")

obs_NZFFD$e_ne <- ifelse(obs_NZFFD$data_value==1, "Encounter", "Non-encounter")
####



# ###########################
# ## Extract waikato network
# ###########################
# 
# net_waikato <- netfull %>% filter(CatName %in% obs_waikato$CatName)
# 
# 
# ###########################
# ## Extract waitaki network
# ###########################
# 
# net_waitaki <- netfull %>% filter(CatName %in% obs_waitaki$CatName)



############
#  NZ Map  #
############

#xlab <- "Temperature (°C)"

nzmap <- ggplot(netfull) +
  geom_point(aes(x = long, y = lat), cex=0.2) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  theme_bw(base_size = 14)

##########################################
# Maps of Waitaki and Waikato catchments #
##########################################

##Waikato
Waikato_catchment_map <- ggplot(obs_waikato)

Waikato_catchment_map <- Waikato_catchment_map + geom_segment(data=net_waikato_connected, aes(x = Lon,y = Lat, xend = Lon2, yend = Lat2), col="#F8766D")

Waikato_catchment_map <- Waikato_catchment_map +
  geom_point(data = obs_waikato, aes(x = long, y = lat, fill = "Sites"), pch = 21, alpha = 0.8) +
  scale_fill_manual(values="white") +
  #guides(fill=guide_legend(title="")) +
  #theme(legend.position = "none") +
  xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
ggsave(file.path(fig_dir, "Waikato_catchment_map.png"), Waikato_catchment_map)


##Waitaki
waitaki_catchment_map <- ggplot(obs_waitaki)

waitaki_catchment_map <- waitaki_catchment_map + geom_segment(data=net_waitaki_connected, aes(x = Lon,y = Lat, xend = Lon2, yend = Lat2), col="#00BFC4")

waitaki_catchment_map <- waitaki_catchment_map +
  geom_point(data = obs_waitaki, aes(x = long, y = lat, fill = "Sites"), pch = 21, alpha = 0.8) +
  scale_fill_manual(values="white") +
  #guides(fill=guide_legend(title="")) +
  #theme(legend.position = "none") +
  xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
ggsave(file.path(fig_dir, "Waitaki_catchment_map.png"), waitaki_catchment_map)



##############################################
# NZ Map with Waitaki and Waikato catchments #
##############################################

nzmap_waikato_waitaki <- nzmap +
  geom_point(data = net_waikato_waitaki, aes(x = long, y = lat, color = Study_area)) +
  geom_point(data = obs_waikato_waitaki, aes(x = long, y = lat, fill = "Sites"), pch = 21, alpha = 0.8) +
  scale_fill_manual(values="white") +
  guides(fill=guide_legend(title=""), color=guide_legend(title="Study area")) +
  theme(legend.position = c(0.3,0.8),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
ggsave(file.path(fig_dir, "Waikato_waitaki_on_NZmap.png"), nzmap_waikato_waitaki)



nzmap_waikato <- nzmap + 
  geom_point(data = net_waikato, aes(x = long, y = lat), col="grey")
ggsave(file.path(fig_dir, "Waikato_on_NZmap.png"), nzmap_waikato)


nzmap_waitaki <- nzmap + 
  geom_point(data = net_waitaki, aes(x = long, y = lat), col="grey")
ggsave(file.path(fig_dir, "Waitaki_on_NZmap.png"), nzmap_waitaki)



nzmap_waikato_waitaki2 <- nzmap_waikato + 
  geom_point(data = net_waitaki, aes(x = long, y = lat), col="grey")
ggsave(file.path(fig_dir, "Waikato_waitaki_on_NZmap2.png"), nzmap_waikato_waitaki2)



##########################
#  Arrange maps on grid  #
##########################

NZ_catchment_grid <- grid.arrange(nzmap_waikato_waitaki, arrangeGrob(Waikato_catchment_map, waitaki_catchment_map, nrow=2), ncol = 2,
                                  widths = c(0.6,0.4),
                                  bottom = textGrob("Longitude (\u00B0E)"),
                                  left = textGrob("Latitude (\u00B0N)", rot = 90, vjust = 1))
ggsave(file.path(fig_dir, "NZ_catchment_grid.png"), NZ_catchment_grid)



##################################
#  NZ Map of NZFFD observations  #
##################################

catchmap <- nzmap +
  geom_point(data=obs_NZFFD, aes(x = long, y = lat, col = e_ne), pch=19, alpha=0.6) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  #ggtitle("Encounter/non-encounter longfin eel observations") +
  guides(col = guide_legend(title = "")) +
  scale_colour_manual(values = c("#377EB8", "#E41A1C")) +
  coord_fixed(1.1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = rel(0.5)),
        axis.text.x = element_text(angle = 90))
ggsave(file.path(fig_dir, "NZ_map_with_obs.png"), catchmap)


table(obs_NZFFD$agency, obs_NZFFD$data_value)

table(obs_NZFFD$fishmethod, obs_NZFFD$data_value)

############################################
#  NZ Map of NZFFD observations by decade  #
############################################

obs_NZFFD2 <- obs_NZFFD %>%
  filter(data_type == "encounter" & year >= 1960 & fishmethod == "Electric fishing")
obs_NZFFD2$Decade <- ifelse(obs_NZFFD2$year <= 1969, "1960-1969", 
                          ifelse(obs_NZFFD2$year >= 1970 & obs_NZFFD2$year <= 1979, "1970-1979", 
                                 ifelse(obs_NZFFD2$year >= 1980 & obs_NZFFD2$year <= 1989, "1980-1989", 
                                        ifelse(obs_NZFFD2$year >= 1990 & obs_NZFFD2$year <= 1999, "1990-1999", 
                                               ifelse(obs_NZFFD2$year >= 2000 & obs_NZFFD2$year <= 2009, "2000-2009", 
                                                      ifelse(obs_NZFFD2$year >= 2010 & obs_NZFFD2$year <= 2019, "2010-2019", NA))))))


table(obs_NZFFD2$data_value ,obs_NZFFD2$Decade, useNA = "ifany")



# obs_NZFFD2 <- obs_NZFFD %>% 
#   filter(data_type == "encounter" & year >= 1960 & fishmethod == "Electric fishing")
# 
# 
# tab_year <- table(obs_NZFFD2$year)
# obs_NZFFD2 <- obs_NZFFD2[obs_NZFFD2$year %in% names(tab_year[tab_year>=30]),]



catchmap2 <- nzmap +
  geom_point(data=obs_NZFFD2, aes(x = long, y = lat, col = e_ne), pch=19, alpha=0.6) +
  facet_wrap(.~Decade) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  #ggtitle("Encounter/non-encounter longfin eel electric fishing observations") +
  guides(col = guide_legend(title = "")) +
  scale_colour_manual(values = c("#377EB8", "#E41A1C")) +
  coord_fixed(1.1) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = rel(0.5)),
        axis.text.x = element_text(angle = 90))
ggsave(file.path(fig_dir, "NZ_map_with_obs_eachyear.png"), catchmap2)
