
rm(list=ls())

################
## Directories
################

sub_dir <- file.path(getwd(), "Waikato")
dir.create(sub_dir, showWarnings = FALSE)

NZdata_dir <- file.path(getwd(), "NZ_data")


data_dir2 <- file.path(sub_dir, "Data")
dir.create(data_dir2, showWarnings = FALSE)

fig_dir <- file.path(sub_dir, "Figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

library(tidyverse)
library(proj4)
library(akima)

#################
## Load data
################

network_full <- readRDS(file.path(NZdata_dir, "NZ_network.rds"))
waikato_len <- readRDS(file.path(data_dir2, "Waikato_obs_lf_length_AC.rds")) # sourced from "Assessing_data_needs" project


################################################
# Subset NZ network to just Waikato catchments #
################################################


#Select important info from network
netfull2 <- network_full #%>% select(nzsegment, CatName, parent_s, child_s, dist_s, long, lat, easting, northing)

nettojoin <- network_full %>% select(nzsegment, CatName, parent_s, child_s, dist_s, long, lat, easting, northing)

## bring important info from network to observations, removing root nodes
obslen2 <- left_join(waikato_len, nettojoin, by="nzsegment") %>% 
  filter(parent_s != 0) 



#netsub <- netfull2 %>% filter(grepl("aikato", CatName))

#Find network with same catchment names as observations
netsub <- netfull2 %>% filter(CatName %in% obslen2$CatName)

all(obslen2$nzsegment %in% netfull2$nzsegment)

nzmap <- ggplot(network_full) +
  geom_point(aes(x = long, y = lat), cex=0.2) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 14)


#Map of the catchment
catmap <- nzmap + 
  geom_point(data = netsub, aes(x = long, y = lat, color = CatName)) +
  geom_point(data = obslen2, aes(x = long, y = lat), pch = 21, fill = "white", alpha = 0.8) +
  guides(color = FALSE)
ggsave(file.path(fig_dir, "Waikato_obs_network.png"), catmap, height = 6, width = 7)



#######################################
#Find nodes around the observations
#######################################

## NOTE:
### Going up and down along the network is necessary as there are weird locations outside
### of the Waikato region. See Waikato_obs_network.png in fig_dir


#Network with all observations (from observations down)
obs_child <- unique(obslen2$child_s) #Unique child observation nodes
net_obs <- netsub %>% filter(child_s %in% obs_child) #Filter these nodes in the network
nextdown <- netsub %>% filter(child_s %in% net_obs$parent_s) #1 down from this
save <- rbind.data.frame(net_obs,nextdown) #save

for(i in 1:200){ #Finds the next 200 nodes down
  nextdown <- netsub %>% filter(child_s %in% nextdown$parent_s) #1 down
  save <- unique(rbind.data.frame(save, nextdown)) #save
  print(nrow(save))
}
netsub2 <- save


#Finds from roots to observations
roots <- netsub2 %>% filter(parent_s == 0) #Filter root nodes (start from the ocean) from the network created above


nextup <- netsub %>% filter(parent_s %in% roots$child_s) #Filter the parents of the roots from the entire network
save <- rbind.data.frame(net_obs,nextup) #Save the roots + 1 up
for(i in 1:200){ #repeat this 200 times
  nextup <- netsub %>% filter(parent_s %in% nextup$child_s) #next node up
  save <- unique(rbind.data.frame(save, nextup)) #save
  print(nrow(save))
}
netsub3 <- save

netsub_toUse <- unique(rbind.data.frame(netsub2, netsub3, roots)) #Network of observations going up and down

netsub_toUse <- netsub #use entire network


#Check both contain the nzsegments we need
all(obslen2$nzsegment %in% netsub_toUse$nzsegment)
#all(obslen2$nzsegment %in% netsub_long$nzsegment)

#Both contain all roots necessary
length(which(netsub_toUse$parent_s == 0))
#length(which(netsub_long$parent_s == 0))

catmap2 <- nzmap + 
  geom_point(data = netsub_toUse, aes(x = long, y = lat, color = CatName)) +
  geom_point(data = obslen2, aes(x = long, y = lat), pch = 21, fill = "white", alpha = 0.8) +
  guides(color = FALSE) 
# facet_wrap(.~CatName) +
# geom_segment(data=l2, aes(x = long2,y = lat2, xend = long, yend = lat), col="gray") +
# theme_bw(base_size = 14)
ggsave(file.path(fig_dir, "Length_network.png"), catmap2, height = 6, width = 7)



## rename nodes
nodes <- unique(c(netsub_toUse$child_s, netsub_toUse$parent_s))
inodes <- seq_along(nodes)

net_parents <- sapply(1:nrow(netsub_toUse), function(x){
  if(netsub_toUse$parent_s[x] != 0) new_node <- inodes[which(nodes == netsub_toUse$parent_s[x])]
  if(netsub_toUse$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
net_children <- sapply(1:nrow(netsub_toUse), function(x) inodes[which(nodes == netsub_toUse$child_s[x])])

netsub_toUse$parent_s <- net_parents
netsub_toUse$child_s <- net_children

obs_parents <- sapply(1:nrow(obslen2), function(x){
  if(obslen2$parent_s[x] != 0) new_node <- inodes[which(nodes == obslen2$parent_s[x])]
  if(obslen2$parent_s[x] == 0) new_node <- 0
  return(new_node)  
})
obs_children <- sapply(1:nrow(obslen2), function(x) inodes[which(nodes == obslen2$child_s[x])])

obslen3 <- obslen2
obslen3$parent_i <- obs_parents
obslen3$child_i <- obs_children

#Remove parent_s and child_s from obs (not needed anymore)
obslen3 <- obslen3 %>%
  select(-c(parent_s, child_s)) %>%
  rename("dist_i"=dist_s, "Year"=year)


saveRDS(netsub_toUse, file.path(data_dir2, "Waikato_updownnetwork.rds"))
saveRDS(obslen3, file.path(data_dir2, "Waikato_observations_lf.rds"))
