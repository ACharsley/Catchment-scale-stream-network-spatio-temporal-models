########################################################
#            Waitaki Longfin eel VAST models           #
#                  Anthony Charsley                    #
#                 Updated April 2021                   #
########################################################


rm(list=ls())


#################
#  Packages     #
#################
##remember to pull upstream development branches
#devtools::install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
#devtools::install_github("merrillrudd/VASTPlotUtils")

## Rtools won't work without doing this (PATH isn't set right!)
library(devtools)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
##

library(VAST)
library(VASTPlotUtils)
library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(foreach)
library(doParallel)


################
## Directories
################

res_dir <- file.path(getwd(), "Waitaki")
data_dir <- file.path(res_dir, "Data")
fig_dir <- file.path(res_dir, "Figures")



########################################################
########################################################



#################
# Read in data  #       
#################
load(file.path(data_dir, "nz_waitaki_longfin_eel.rda"))
network <- nz_waitaki_longfin_eel[["network"]]

## format network data (data about the stream network)
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat)

## make sure to use only encounter data
obs <- nz_waitaki_longfin_eel[["observations"]] %>%
  dplyr::filter(data_type=="encounter") %>% 
  select(-data_type) %>%
  rename('present' = data_value) %>%
  mutate('method_agency' = paste0(fishmethod, "_", agency)) %>%
  rename("Year"=year) %>%
  filter(fishmethod == "Electric fishing")


##########
# Catchment plot

l2 <- lapply(1:nrow(network), function(x){
  parent <- network$parent_s[x]
  find <- network %>% filter(child_s == parent)
  if(nrow(find)>0) out <- cbind.data.frame(network[x,], 'Lon2'=find$long, 'Lat2'=find$lat)
  if(nrow(find)==0) out <- cbind.data.frame(network[x,], 'Lon2'=NA, 'Lat2'=NA)
  return(out)
})
l2 <- do.call(rbind, l2)

obs$e_ne <- ifelse(obs$present==1, "Encounter", "Non-encounter")
#obs$p_a <- ifelse(obs$present==1, "Present", "Absent")

catchmap <- ggplot() +
  geom_point(data=network, aes(x = long, y = lat), col="gray") +
  geom_segment(data=l2, aes(x = Lon2,y = Lat2, xend = long, yend = lat), col="gray") +
  geom_point(data=obs, aes(x = long, y = lat, fill = e_ne), pch=22, alpha=0.6) +
  xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  #ggtitle("Presence/absence longfin eel observations") +
  coord_fixed() +
  guides(fill = guide_legend(title = "")) +
  scale_colour_manual(values = c("#377EB8", "#E41A1C"), aesthetics = "fill") +
  theme_bw()
ggsave(file.path(fig_dir, "Waitaki_map_with_obs.png"), catchmap)


catchmap2 <- ggplot() +
  geom_point(data=network, aes(x = long, y = lat), col="gray") +
  geom_segment(data=l2, aes(x = Lon2,y = Lat2, xend = long, yend = lat), col="gray") +
  geom_point(data=obs, aes(x = long, y = lat, fill = e_ne), pch=22, alpha=0.6) +
  facet_wrap(.~Year) +
  xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  #ggtitle("Presence/absence longfin eel observations") +
  coord_fixed() +
  guides(fill = guide_legend(title = "")) +
  scale_colour_manual(values = c("#377EB8", "#E41A1C"), aesthetics = "fill") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = rel(0.5)),
        axis.text.x = element_text(angle = 90))
ggsave(file.path(fig_dir, "Waitaki_map_with_obs_eachyear.png"), catchmap2)


##########


##### add small value to encounter observations
set.seed(1234)
present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
                            "Year" = as.numeric(obs$Year),
                            "Method" = obs$fishmethod,
                            "Agency" = obs$agency,
                            "Method_Agency" = obs$method_agency, 
                            "AreaSwept_km2" = 150/1000, 
                            "Area_km2" = obs$length,
                            "Lat" = obs$lat, 
                            "Lon" = obs$long, 
                            "Pass" = 0,
                            "Knot" = obs$child_i,
                            "Category" = "Longfin_eels")

#When the areaswept is greater than the segment area, set areaswept as the segment area
Data_Geostat[Data_Geostat$AreaSwept_km2 > Data_Geostat$Area_km2,"AreaSwept_km2"] <- Data_Geostat[Data_Geostat$AreaSwept_km2 > Data_Geostat$Area_km2,"Area_km2"]
any(Data_Geostat$AreaSwept_km2 > Data_Geostat$Area_km2)

## habitat data
hab <- nz_waitaki_longfin_eel[['habitat']]

covar_toUse <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm",'DamAffected')
hab <- hab %>% filter(covariate %in% covar_toUse)



#######################################
#        Habitat information          #
# NOTE: Treated as density covariates #
#######################################

nodes <- network$child_s[order(network$child_s)]
years <- min(obs$Year):max(obs$Year)
covar <- unique(hab$covariate)
n_x <- length(nodes)
n_t <- length(years)
n_p <- length(covar)
n_i <- nrow(obs)

###############
#Look at distribution of raw covariates

covs_to_plot <- covar[covar!="DamAffected"]

lapply(1:length(covs_to_plot), function(x){
  sub <- hab %>% 
    filter(covariate == covs_to_plot[x])
  
  
  
  jpeg(paste0(fig_dir,"/", covs_to_plot[x], "_rawcov.jpeg"), height=8, width=8,units="in", res=600)
  boxplot(sub$value, main=paste0("Distribution of ", covs_to_plot[x]))
  dev.off()
  
})

p <- ggplot(hab %>% filter(covariate == "MeanFlowCumecs"), aes(y=log(value))) +
  geom_boxplot()
ggsave(file.path(fig_dir, paste0("Log_dist_MeanFlowCumecs.png")),p, width=8, height=8)


## LOG TRANSFORM MeanFlowCumecs ##

log_MeanFlowCumecs <- hab %>%
  filter(covariate == "MeanFlowCumecs") %>%
  mutate(value = log(value))


#Take out MeanFlowCumecs and add it logged

hab <- hab %>%
  filter(covariate != "MeanFlowCumecs") #Taken out

hab <- rbind(hab, log_MeanFlowCumecs) #Add in

###############


hab_std <- lapply(1:length(covar), function(x){
  sub <- hab %>% 
    filter(covariate == covar[x]) %>%
    mutate(value_std = (value - mean(value))/sd(value))
  return(sub)
})
hab_std <- do.call(rbind, hab_std)

p_habstd <- ggplot(hab_std) +
  geom_point(aes(x = easting, y = northing, color = value_std)) +
  facet_wrap(.~covariate, nrow = 2) +
  scale_color_distiller(palette = "Spectral") +
  scale_x_continuous(breaks = quantile(hab$easting, prob = c(0.1, 0.5, 0.95)), labels = round(quantile(hab$easting, prob = c(0.1,0.5,0.95)),0)) +
  guides(color = guide_legend(title = "Standardised\nvalue")) +
  theme_minimal()
ggsave(file.path(fig_dir, "Habitat_covariates_standardised.png"), p_habstd, width = 15, height = 8)

for(i in 1:length(covar)){
  p <- ggplot(hab_std %>% filter(covariate == covar[i])) +
  geom_point(aes(x = easting, y = northing, color = value)) +
  guides(color=guide_legend(title=covar[i])) +
  scale_color_distiller(palette = "Spectral") +
  theme_minimal()
  ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar[i],".png")),p, width=10, height=8)
}



#Distribution of std. covariates
for(i in 1:length(covar)){
  p <- ggplot(hab_std %>% filter(covariate == covar[i]), aes(y=value_std)) +
    geom_boxplot()
  ggsave(file.path(fig_dir, paste0("Std_covariate_dist_", covar[i],".png")),p, width=8, height=8)
}



#################

X_gtp_input1 <- array(0, dim=c(n_x, n_t, n_p))
for(p in 1:n_p){
  psub <- hab %>% filter(covariate == covar[p])
  mat <- matrix(0, nrow=n_x, ncol = 1)
  mat[psub$child_s,1] <- psub$value
  if(covar[p]=="DamAffected"){
    X_gtp_input1[,,p] <- mat
  } else {
    mat_sd <- (mat - mean(mat, na.rm=TRUE))/sd(mat, na.rm=TRUE)
    X_gtp_input1[,,p] <- mat_sd
  }
}

## years since dam impact
X_choose <- X_gtp_input1[,,which(covar == "DamAffected")]
X_gtp1 <- sapply(1:length(years), function(x){
  sub <- X_choose[,x]
  sub[which(sub == 1)] <- years[x] - 1935
  return(sub)
})
X_gtp1_sd <- (X_gtp1 - mean(X_gtp1))/sd(X_gtp1)

## years since impact squared
X_gtp2 <- sapply(1:length(years), function(x){
  sub <- X_choose[,x]
  sub[which(sub == 1)] <- (years[x] - 1935)^2
  return(sub)
})
X_gtp2_sd <- (X_gtp2 - mean(X_gtp2))/sd(X_gtp2)

covar2 <- c(covar, "YearsSinceDam","YearsSinceDam2")[-which(covar=="DamAffected")]
n_p <- length(covar2)
X_gtp_input <- array(0, dim=c(n_x,n_t,n_p))
for(p in 1:(n_p)){
  ## skip dam affected
  if(p < length(covar)) X_gtp_input[,,p] <- X_gtp_input1[,,p]
  
  ## in place of dam affected, years since dam
  if(p == length(covar)) X_gtp_input[,,p] <- X_gtp1_sd
  
  ## additional covariate, years since dam squared
  if(p ==length(covar)+1) X_gtp_input[,,p] <- X_gtp2_sd
}

## match habitat covariates to observations
## double check the indices will match up properly
X_itp_input <- array(0, dim=c(n_i,n_t,n_p))
for(i in 1:n_i){
  for(p in 1:n_p){
    child_i <- obs$child_i[i]
    index <- which(nodes == child_i)
    X_itp_input[i,,p] <- X_gtp_input[index,,p]
  }
}


##### Plot final distributions

for(p in 1:length(covar2)){
  
  jpeg(paste0(fig_dir,"/Final_cov_dist_", covar2[p], ".jpeg"), height=8, width=8,units="in", res=600)
  boxplot((as.vector(X_gtp_input[,,which(covar2==covar2[p])])), main=paste0("Distribution of ", covar2[p]))
  dev.off()
  
}








####################################
## sampling information
## treated as catchability covariates
###################################

# #Trap, visual and angling sampling methods are grouped
# method_info <- Data_Geostat %>%
#   group_by(Method) %>%
#   summarise(Samples = length(Method), 
#             Prop_samples = length(Method)/nrow(Data_Geostat),
#             Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Method))
# method_info <- method_info[rev(order(method_info$Samples)),]
# Data_Geostat$Method2 <- as.character(Data_Geostat$Method)
# Data_Geostat$Method2[which(Data_Geostat$Method %in% c("Trap","Visual","Angling"))] <- "Trap_Visual_Angling"
# Q_ik_method <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Method2"])[,-1,drop=FALSE]


#Fish&game and consultants are grouped
agency_info <- Data_Geostat %>%
  group_by(Agency) %>%
  summarise(Samples = length(Agency), 
            Prop_samples = length(Agency)/nrow(Data_Geostat),
            Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Agency))
agency_info <- agency_info[rev(order(agency_info$Samples)),]
Data_Geostat$Agency2 <- as.character(Data_Geostat$Agency)
Data_Geostat$Agency2[which(Data_Geostat$Agency %in% c("fish&game","consultants"))] <- "fish&game_consultants"
Q_ik_agency <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Agency2"])[,-2,drop=FALSE]


# method_agency_info <- Data_Geostat %>%
#   group_by(Method_Agency) %>%
#   summarise(Samples = length(Method_Agency), 
#             Prop_samples = length(Method_Agency)/nrow(Data_Geostat),
#             Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Method_Agency))
# method_agency_info <- method_agency_info[rev(order(method_agency_info$Samples)),]
# Data_Geostat$Method_Agency2 <- sapply(1:nrow(Data_Geostat), function(x){
#   if(grepl("Electric", Data_Geostat$Method_Agency[x]) == FALSE){ #figure out what this does
#     out <- as.character(Data_Geostat$Method[x])
#   } else{
#     out <- as.character(Data_Geostat$Method_Agency[x])
#   }
#   if(out %in% c("Trap","Visual","Angling")) out <- "Method_other"
#   if(grepl("fish&game",out) | grepl("consultants",out)) out <- "Electric fishing_agency_other"
#   return(out)
# })
# method_agency_info2 <- Data_Geostat %>%
#   group_by(Method_Agency2) %>%
#   summarise(Samples = length(Method_Agency2), 
#             Prop_samples = length(Method_Agency2)/nrow(Data_Geostat),
#             Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Method_Agency2))
# method_agency_info2 <- method_agency_info2[rev(order(method_agency_info2$Samples)),]
# Q_ik_methodagency <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Method_Agency2"])[,-3,drop=FALSE]


## Save inputs
save.image(file.path(data_dir, "general_inputs_Waitaki.Rdata"))


########################################################
########################################################

##################################
#            Models              #
##################################

##########################
## Model1
## Spatial
## Temporal (beta1 IID)
## Spatiotemporal (independent among years)
## Lognormal dist
## Habitat covariates
## Agency catchability covariate
###########################

rm(list=ls())

load("./Waitaki/Data/general_inputs_Waitaki.Rdata")

path <- file.path(res_dir, "model1")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.dll"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0

Q_ik_inp <- Q_ik_agency
#Q_ik_inp <- NULL

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v12_0_0",
                          n_x = nrow(Network_sz),
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE,
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp,
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- 5
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp,
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 c_i = rep(0, nrow(Data_inp)),
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 working_dir = path,
                 extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp,
                 X_gtp = X_gtp_inp,
                 X_itp = X_itp_inp,
                 Q_ik = Q_ik_inp,
                 model_args = list(Map = Map, Par = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1)
plot_maps(plot_set = c(6), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1)

plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")
plot_maps(plot_set = c(6), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

##########################
## model1a
## Spatial
## Temporal (beta1 IID)
## Spatiotemporal (Random walk)
## Lognormal dist
## Habitat covariates
## Agency catchability covariate
###########################

rm(list=ls())

load("./Waitaki/Data/general_inputs_Waitaki.Rdata")

path <- file.path(res_dir, "model1a")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.dll"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0
Q_ik_inp <- Q_ik_method

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 2, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0",
                          n_x = nrow(Network_sz),
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE,
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp,
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp,
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 c_i = rep(0, nrow(Data_inp)),
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 working_dir = path,
                 extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp,
                 X_gtp = X_gtp_inp,
                 X_itp = X_itp_inp,
                 Q_ik = Q_ik_inp,
                 model_args = list(Map = Map, Par = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1)
plot_maps(plot_set = c(6), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1)

plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")
plot_maps(plot_set = c(6), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


# ##########################
# ## model2
# ## spatial, temporal, 
# ## beta1 IID, lognormal dist, 
# ## habitat covariate, fishing method covariate
# ###########################
# path <- file.path(res_dir, "model2_AC")
# dir.create(path, showWarnings = FALSE)
# fig <- file.path(path, "figures")
# dir.create(fig)
# 
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)
# 
# Data_inp <- Data_Geostat
# X_gtp_inp <- X_gtp_input
# X_itp_inp <- X_itp_input
# Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
# Xconfig_zcp_inp[2,,] <- 0
# Q_ik_inp <- Q_ik_method
# 
# FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
# RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
# ObsModel <- cbind("PosDist" = 1, "Link" = 0)
# OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
# Options <- c("Calculate_range" = 1,
#              "Calculate_effective_area" = 1)
# 
# settings <- make_settings(Version = "VAST_v8_2_0", 
#                           n_x = nrow(Network_sz), 
#                           Region = "Stream_network",
#                           FieldConfig = FieldConfig,
#                           RhoConfig = RhoConfig,
#                           ObsModel = ObsModel,
#                           OverdispersionConfig = OverdispersionConfig,
#                           Options = Options,
#                           purpose = "index",
#                           fine_scale = FALSE, 
#                           bias.correct = FALSE)
# settings$Method <- "Stream_network"
# settings$grid_size_km <- 1
# 
# 
# fit0 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   run_model = FALSE)
# 
# Par <- fit0$tmb_list$Parameters
# Map <- fit0$tmb_list$Map
# Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
# Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))
# 
# fit1 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   model_args = list(Map = Map, Par = Par),
#                   optimize_args = list(getsd = FALSE, newtonsteps = 0),
#                   test_fit = FALSE)
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# 
# fit <- fit_model(settings = settings,
#                  Lat_i = Data_inp[,"Lat"],
#                  Lon_i = Data_inp[,"Lon"],
#                  t_i = Data_inp[,"Year"],
#                  c_i = rep(0, nrow(Data_inp)),
#                  b_i = Data_inp[,"Catch_KG"],
#                  a_i = Data_inp[,"AreaSwept_km2"],
#                  working_dir = path,
#                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                  spatial_args = list(Network_sz_LL = Network_sz_LL),
#                  Network_sz = Network_sz,
#                  Xconfig_zcp = Xconfig_zcp_inp,
#                  X_gtp = X_gtp_inp, 
#                  X_itp = X_itp_inp,
#                  Q_ik = Q_ik_inp,
#                  model_args = list(Map = Map, Par = Par),
#                  test_fit = FALSE)
# saveRDS(fit, file.path(path, "Fit.rds"))
# 
# fit <- readRDS(file.path(path, "Fit.rds"))
# 
# 
# plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
# plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")
# 
# map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
# Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
# plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
# plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

##########################
## model3
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## fishing method covariate
###########################

rm(list=ls())

load("./Waitaki/Data/general_inputs_Waitaki.Rdata")

#path <- file.path(res_dir, "model3_AC")
path <- file.path(res_dir, "model3")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- NULL
X_itp_inp <- NULL
Xconfig_zcp_inp <- NULL
Q_ik_inp <- Q_ik_method

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
TMBhelper::check_estimability(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 c_i = rep(0, nrow(Data_inp)),
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 working_dir = path,
                 extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp,
                 X_gtp = X_gtp_inp, 
                 X_itp = X_itp_inp,
                 Q_ik = Q_ik_inp,
                 model_args = list(Map = Map, Par = Par),
                 test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

# ##########################
# ## model3b
# ## spatiotemporal, spatial, temporal, 
# ## beta1 IID, lognormal dist, 
# ## fishing method covariate
# ## ST smoother -- random walk
# ###########################
# path <- file.path(res_dir, "model3b_AC")
# dir.create(path, showWarnings = FALSE)
# fig <- file.path(path, "figures")
# dir.create(fig)
# 
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)
# 
# Data_inp <- Data_Geostat
# X_gtp_inp <- NULL
# X_itp_inp <- NULL
# Xconfig_zcp_inp <- NULL
# Q_ik_inp <- Q_ik_method
# 
# FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
# RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 2, "Epsilon2" = 0)
# ObsModel <- cbind("PosDist" = 1, "Link" = 0)
# OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
# Options <- c("Calculate_range" = 1,
#              "Calculate_effective_area" = 1)
# 
# settings <- make_settings(Version = "VAST_v8_2_0", 
#                           n_x = nrow(Network_sz), 
#                           Region = "Stream_network",
#                           FieldConfig = FieldConfig,
#                           RhoConfig = RhoConfig,
#                           ObsModel = ObsModel,
#                           OverdispersionConfig = OverdispersionConfig,
#                           Options = Options,
#                           purpose = "index",
#                           fine_scale = FALSE, 
#                           bias.correct = FALSE)
# settings$Method <- "Stream_network"
# settings$grid_size_km <- 1
# 
# 
# fit0 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   run_model = FALSE)
# 
# Par <- fit0$tmb_list$Parameters
# Map <- fit0$tmb_list$Map
# Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
# Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))
# 
# fit1 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   model_args = list(Map = Map, Par = Par),
#                   optimize_args = list(getsd = FALSE, newtonsteps = 0),
#                   test_fit = FALSE)
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# 
# fit <- fit_model(settings = settings,
#                  Lat_i = Data_inp[,"Lat"],
#                  Lon_i = Data_inp[,"Lon"],
#                  t_i = Data_inp[,"Year"],
#                  c_i = rep(0, nrow(Data_inp)),
#                  b_i = Data_inp[,"Catch_KG"],
#                  a_i = Data_inp[,"AreaSwept_km2"],
#                  working_dir = path,
#                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                  spatial_args = list(Network_sz_LL = Network_sz_LL),
#                  Network_sz = Network_sz,
#                  Xconfig_zcp = Xconfig_zcp_inp,
#                  X_gtp = X_gtp_inp, 
#                  X_itp = X_itp_inp,
#                  Q_ik = Q_ik_inp,
#                  model_args = list(Map = Map, Par = Par),
#                  test_fit = FALSE)
# saveRDS(fit, file.path(path, "Fit.rds"))
# 
# fit <- readRDS(file.path(path, "Fit.rds"))
# 
# plot_maps(plot_set = c(1,6, 13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
# plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")
# 
# map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
# Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
# plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
# plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


# ##########################
# ## model4
# ## spatial, temporal, 
# ## beta1 IID, lognormal dist, 
# ## habitat covariate, agency covariate
# ###########################
# path <- file.path(res_dir, "model4_AC")
# dir.create(path, showWarnings = FALSE)
# fig <- file.path(path, "figures")
# dir.create(fig)
# 
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)
# 
# Data_inp <- Data_Geostat
# X_gtp_inp <- X_gtp_input
# X_itp_inp <- X_itp_input
# Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
# Xconfig_zcp_inp[2,,] <- 0
# Q_ik_inp <- Q_ik_agency
# 
# FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
# RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
# ObsModel <- cbind("PosDist" = 1, "Link" = 0)
# OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
# Options <- c("Calculate_range" = 1,
#              "Calculate_effective_area" = 1)
# 
# settings <- make_settings(Version = "VAST_v8_2_0", 
#                           n_x = nrow(Network_sz), 
#                           Region = "Stream_network",
#                           FieldConfig = FieldConfig,
#                           RhoConfig = RhoConfig,
#                           ObsModel = ObsModel,
#                           OverdispersionConfig = OverdispersionConfig,
#                           Options = Options,
#                           purpose = "index",
#                           fine_scale = FALSE, 
#                           bias.correct = FALSE)
# settings$Method <- "Stream_network"
# settings$grid_size_km <- 1
# 
# 
# fit0 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   run_model = FALSE)
# 
# Par <- fit0$tmb_list$Parameters
# Map <- fit0$tmb_list$Map
# Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
# Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))
# 
# fit1 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   model_args = list(Map = Map, Par = Par),
#                   optimize_args = list(getsd = FALSE, newtonsteps = 0),
#                   test_fit = FALSE)
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# 
# fit <- fit_model(settings = settings,
#                  Lat_i = Data_inp[,"Lat"],
#                  Lon_i = Data_inp[,"Lon"],
#                  t_i = Data_inp[,"Year"],
#                  c_i = rep(0, nrow(Data_inp)),
#                  b_i = Data_inp[,"Catch_KG"],
#                  a_i = Data_inp[,"AreaSwept_km2"],
#                  working_dir = path,
#                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                  spatial_args = list(Network_sz_LL = Network_sz_LL),
#                  Network_sz = Network_sz,
#                  Xconfig_zcp = Xconfig_zcp_inp,
#                  X_gtp = X_gtp_inp, 
#                  X_itp = X_itp_inp,
#                  Q_ik = Q_ik_inp,
#                  model_args = list(Map = Map, Par = Par),
#                  test_fit = FALSE)
# saveRDS(fit, file.path(path, "Fit.rds"))
# 
# fit <- readRDS(file.path(path, "Fit.rds"))
# 
# plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
# plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")
# 
# map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
# Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
# plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
# plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model5 (selected model)
## Spatial
## Temporal (beta1 random walk)
## Spatiotemporal OFF
## Lognormal dist
## Habitat covariates
## Agency catchability covariate
###########################


rm(list=ls())

load("./Waitaki/Data/general_inputs_Waitaki.Rdata")

path <- file.path(res_dir, "model5")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.dll"), to = path)


Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0

#Take out yearssincedam^2
Xconfig_zcp_inp[,,which(covar2=="YearsSinceDam2")] <- 0

formula_inp <- paste0("~",(paste0(covar2, collapse = "+")))

#Q_ik_inp <- Q_ik_methodagency
Q_ik_inp <- Q_ik_agency #Changing to this as method is only EF now

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 2, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v12_0_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index2",
                          fine_scale = FALSE, 
                          bias.correct = T)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"Area_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  formula = formula_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)


Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
#Par[["logkappa1"]] <- 5
#Par[["logkappa1"]] <- -5
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"Area_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
TMBhelper::check_estimability(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 c_i = rep(0, nrow(Data_inp)),
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 working_dir = path,
                 extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"Area_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp,
                 X_gtp = X_gtp_inp, 
                 X_itp = X_itp_inp,
                 Q_ik = Q_ik_inp,
                 model_args = list(Map = Map, Par = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

#fit <- readRDS(file.path(path, "Fit.rds"))

source("./Code/funcs.R")


network <- network %>% rename("Lat"=lat, "Lon"=long)


## Check model ##
#plot_anisotropy(FileName=file.path(fig,"Aniso.png"), Obj=fit$tmb_list$Obj)
dharmaRes = summary(fit, what="residuals", working_dir=paste0(fig,"/"), type=1)
#dharmaRes <- readRDS(file.path(path, "dharmaRes.rds"))
plot_residuals(residuals=dharmaRes$scaledResiduals, fit=fit, Data_inp = Data_inp, network=network, 
               coords="lat_long", save_dir = fig)

saveRDS(dharmaRes, file.path(path, "dharmaRes.rds"))

####
## Arnaud's residual plots ##

######## Check whether observed encounter frequencies for either low or high probability samples 
######## are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_Geostat, DirName = paste0(fig, '/'))


######## Run diagnostics for the positive-catch-rate component of the delta-Lognormal model
######## We can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
######## A good Q-Q plot will have residuals along the one-to-one line
# source("./Code/plotQuantileDiagnostic.r")
# Q = plotQuantileDiagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive",
#                            FileName_Phist = "Posterior_Predictive-Histogram", save_dir = paste0(fig, "/QQ_Fn/" ),
#                            FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")

####




## dharma plots ##
# Histogram of residuals #
val = dharmaRes$scaledResiduals
val[val == 0] = -0.01
val[val == 1] = 1.01

jpeg(file.path(fig, "Resid_hist.jpg"), width = 600, height = 600)
hist(val, 
     breaks = seq(-0.02, 1.02, len = 53),
     col = c("red",rep("lightgrey",50), "red"),
     #main = "Hist of DHARMa residuals",
     main = "",
     xlab = "Residuals (outliers are marked red)",
     cex.axis=1.5, cex.lab=1.5)
dev.off()
##
## QQ plot ##
jpeg(file.path(fig, "QQplot.jpg"), width = 600, height = 600)
gap::qqunif(dharmaRes$scaledResiduals,
            pch=2,
            bty="n", 
            logscale = F, 
            col = "black", 
            cex = 0.6, 
            cex.axis=1.5, cex.lab=1.5, 
            #main = "QQ plot residuals", 
            cex.main = 1)
dev.off()
##
####



## Model comparisons ##
fit$parameter_estimates$AIC


#####
## Plot key model results ##

# Map of lf probability of capture across time
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category", 
                  PlotName = "POC_lf_yearly",
                  #PlotTitle = "Longfin eel yearly probability of capture in Waitaki, NZ",
                  PlotTitle = "",
                  cex = 0.5, 
                  Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)

# Map of lf probability of capture for each year
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel = "Year",
                  PlotName = "POC_lf",
                  PlotTitle = "Longfin eel P.O.C in Waitaki, NZ",
                  cex = 0.75,  
                  #Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)

# Map of spatio-temporal variation of longfin eel probability of capture
plot_maps_network(plot_set = 5, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category",
                  PlotName = "Epsilon_lf",
                  PlotTitle = "Spatio-temporal variation of longfin eel probability of capture",
                  cex = 0.5, 
                  arrows=F)



# #Plot individual covariate effects
# #covar_names_to_use <- covars_all #Add names here
# covar_names_to_use <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
#                         "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam")
# 
# n_p = length(covar_names_to_use)
# 
# 
# ##Covariate effects across time
# cov_effects_array <- plot_maps_network(plot_set = 10, 
#                                        fit = fit, 
#                                        Sdreport = fit$parameter_estimates$SD, 
#                                        TmbData = fit$data_list, 
#                                        spatial_list = fit$spatial_list, 
#                                        DirName = fig, 
#                                        Panel = "category",
#                                        category_names="Longfin eel",
#                                        covar_names=covar_names_to_use,
#                                        PlotName="Covariate_effects",
#                                        PlotTitle = "Individual covariate effect",
#                                        cex = 0.5, 
#                                        #Zlim = c(0,1), 
#                                        arrows=F, 
#                                        pch=15,
#                                        n_p = n_p,
#                                        which_np_touse = which(covar_names_to_use == covar_names_to_use),
#                                        which_cat_cov_toplot=1)
# 
# saveRDS(cov_effects_array, file.path(fig, "cov_effects_array.rds"))
# #cov_effects_array <- readRDS(file.path(fig, "cov_effects_array.rds"))
# 
# 
# 
# 
# ## Plot raw (standardised) covariate values across all time ##
# 
# Network_sz_EN <- data.frame('parent_s'=fit$data_list$parent_s, 'child_s'=fit$data_list$child_s, fit$spatial_list$latlon_g)
# l2 <- lapply(1:nrow(Network_sz_EN), function(x){
#   parent <- Network_sz_EN$parent_s[x]
#   find <- Network_sz_EN %>% filter(child_s == parent)
#   # if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$E_km, 'N2'=find$N_km)
#   # if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
#   if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$Lon, 'N2'=find$Lat)
#   if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
#   return(out)
# })
# l2 <- do.call(rbind, l2)
# 
# 
# for(covar in which(covar_names_to_use %in% covar_names_to_use[c(1:6)])){ #change 'which' depending on what I want to plot - don't include years since dam
#   
#   
#   #Which covariate to plot
#   print(paste0("Plotting: ", covar_names_to_use[covar], " covariate effect"))
#   
#   #Create data frame for plotting
#   xct <- data.frame('value'=cov_effects_array[,covar,1], fit$spatial_list$latlon_g, "category"="Longfin") #'1' is first year - years are the same!
#   
#   Xlim = c(min(xct$Lon),max(xct$Lon))
#   Ylim = c(min(xct$Lat),max(xct$Lat))
#   
#   #Z limits should be set for consistency between years
#   inp_Zlim = quantile(cov_effects_array[,covar,1], prob = c(0,1), na.rm=TRUE)
#   
#   p <- ggplot(xct)
#   
#   p <- p + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = E2, yend = N2), col="gray92")  
#   
#   p <- p +
#     geom_point(aes(x = Lon, y = Lat, color = value), cex = 0.5, pch=19) +
#     scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
#     coord_cartesian(xlim = Xlim, ylim = Ylim) +
#     scale_x_continuous(breaks=quantile(xct$Lon, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$Lon, prob=c(0.1,0.5,0.9)),1)) +
#     # guides(color=guide_legend(title=plot_codes[plot_num])) +
#     mytheme() +
#     xlab("Longitude") + ylab("Latitude")
#   
#   p <- p + ggtitle(wrapper(paste0("Individual covariate effect - ", covar_names_to_use[covar]), width = 50))
#   ggsave(file.path(fig, paste0("Covariate_effects_", covar, ".png")), p, width=8,height=8)
#   
#   
# }
# ####


# Covariate names in plotting
covar_names_to_use <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
                        "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam", "Years since dam squared")

Raw_cov_array = plot_maps_network(plot_set = 9,
                                  fit = fit,
                                  Sdreport = fit$parameter_estimates$SD,
                                  TmbData = fit$data_list,
                                  spatial_list = fit$spatial_list,
                                  DirName = fig,
                                  Panel = "category",
                                  PlotName = "Raw_covariate_values",
                                  PlotTitle = "Standardised covariate values",
                                  covar_names = covar_names_to_use,
                                  cex = 0.5,
                                  arrows=F,
                                  pch=15)

saveRDS(Raw_cov_array, file.path(fig, "Raw_cov_array.rds"))


Network_sz_EN <- data.frame('parent_s'=fit$data_list$parent_s, 'child_s'=fit$data_list$child_s, fit$spatial_list$latlon_g)
l2 <- lapply(1:nrow(Network_sz_EN), function(x){
  parent <- Network_sz_EN$parent_s[x]
  find <- Network_sz_EN %>% filter(child_s == parent)
  # if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$E_km, 'N2'=find$N_km)
  # if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
  if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$Lon, 'N2'=find$Lat)
  if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
  return(out)
})
l2 <- do.call(rbind, l2)


covs_to_plot <- covar2[!(covar2 %in% c("YearsSinceDam", "YearsSinceDam2"))]

for(covar in 1:length(covs_to_plot)){ #change 'which' depending on what I want to plot - don't include years since dam


  #Which covariate to plot
  print(paste0("Plotting: ", covar_names_to_use[covar], " raw covariate"))

  #Create data frame for plotting
  xct <- data.frame('value'=Raw_cov_array[,covar,1], fit$spatial_list$latlon_g) #'1' is first year - years are the same!

  Xlim = c(min(xct$Lon),max(xct$Lon))
  Ylim = c(min(xct$Lat),max(xct$Lat))

  #Z limits
  inp_Zlim = quantile(Raw_cov_array[,covar,1], prob = c(0,1), na.rm=TRUE)

  p <- ggplot(xct)

  p <- p + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = E2, yend = N2), col="gray92")

  p <- p +
    geom_point(aes(x = Lon, y = Lat, color = value), cex = 0.5, pch=19) +
    scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
    coord_cartesian(xlim = Xlim, ylim = Ylim) +
    #scale_x_continuous(breaks=quantile(xct$Lon, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$Lon, prob=c(0.1,0.5,0.9)),1)) +
    # guides(color=guide_legend(title=plot_codes[plot_num])) +
    mytheme() +
    xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)")

  p <- p + ggtitle(wrapper(covar_names_to_use[covar], width = 50))
  ggsave(file.path(fig, paste0(covar, "_raw_values", ".png")), p, width=8,height=8)


}



waitaki_river_length = (sum(network$length))

plot_range_index_riverlength(Report = fit$Report, 
                             TmbData = fit$data_list, 
                             Sdreport = fit$parameter_estimates$SD, 
                             Znames = colnames(fit$data_list$Z_gm), 
                             PlotDir = fig, 
                             year_labels = as.character(fit$year_labels), 
                             #years_to_plot = fit$years_to_plot,
                             FileName_EffArea = paste0(fig,"/River_length_occupied.png"),
                             use_biascorr = TRUE, 
                             category_names = "Longfin_eels",
                             total_river_length = waitaki_river_length)



plot_index(Index_ctl = fit$Report$Index_ctl,
           year_labels = as.character(c(fit$year_labels)),
           DirName = paste0(fig, "/"),
           Yrange = c(0, 5))


#####













##########################
## model5b
## Spatial
## Temporal (beta1 random walk)
## Spatiotemporal OFF
## Lognormal dist
## Habitat covariates
## Agency catchability covariate OFF
###########################


rm(list=ls())

load("./Waitaki/Data/general_inputs_Waitaki.Rdata")

path <- file.path(res_dir, "model5_woagency")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.dll"), to = path)


Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0

#Take out yearssincedam^2
Xconfig_zcp_inp[,,which(covar2=="YearsSinceDam2")] <- 0

formula_inp <- paste0("~",(paste0(covar2, collapse = "+")))

#Q_ik_inp <- Q_ik_methodagency
#Q_ik_inp <- Q_ik_agency #Changing to this as method is only EF now
Q_ik_inp <- NULL

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 2, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v12_0_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index2",
                          fine_scale = FALSE, 
                          bias.correct = T)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"Area_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  formula = formula_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)


Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
#Par[["logkappa1"]] <- 5
#Par[["logkappa1"]] <- -5
#Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"Area_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
TMBhelper::check_estimability(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 c_i = rep(0, nrow(Data_inp)),
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 working_dir = path,
                 extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"Area_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp,
                 X_gtp = X_gtp_inp, 
                 X_itp = X_itp_inp,
                 Q_ik = Q_ik_inp,
                 model_args = list(Map = Map, Par = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

#fit <- readRDS(file.path(path, "Fit.rds"))

source("./Code/funcs.R")


network <- network %>% rename("Lat"=lat, "Lon"=long)


## Check model ##
#plot_anisotropy(FileName=file.path(fig,"Aniso.png"), Obj=fit$tmb_list$Obj)
dharmaRes = summary(fit, what="residuals", working_dir=paste0(fig,"/"), type=1)
#dharmaRes <- readRDS(file.path(path, "dharmaRes.rds"))
plot_residuals(residuals=dharmaRes$scaledResiduals, fit=fit, Data_inp = Data_inp, network=network, 
               coords="lat_long", save_dir = fig)

saveRDS(dharmaRes, file.path(path, "dharmaRes.rds"))

####
## Arnaud's residual plots ##

######## Check whether observed encounter frequencies for either low or high probability samples 
######## are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_Geostat, DirName = paste0(fig, '/'))


######## Run diagnostics for the positive-catch-rate component of the delta-Lognormal model
######## We can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
######## A good Q-Q plot will have residuals along the one-to-one line
# source("./Code/plotQuantileDiagnostic.r")
# Q = plotQuantileDiagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive",
#                            FileName_Phist = "Posterior_Predictive-Histogram", save_dir = paste0(fig, "/QQ_Fn/" ),
#                            FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")

####




## dharma plots ##
# Histogram of residuals #
val = dharmaRes$scaledResiduals
val[val == 0] = -0.01
val[val == 1] = 1.01

jpeg(file.path(fig, "Resid_hist.jpg"), width = 600, height = 600)
hist(val, 
     breaks = seq(-0.02, 1.02, len = 53),
     col = c("red",rep("lightgrey",50), "red"),
     #main = "Hist of DHARMa residuals",
     main = "",
     xlab = "Residuals (outliers are marked red)",
     cex.axis=1.5, cex.lab=1.5)
dev.off()
##
## QQ plot ##
jpeg(file.path(fig, "QQplot.jpg"), width = 600, height = 600)
gap::qqunif(dharmaRes$scaledResiduals,
            pch=2,
            bty="n", 
            logscale = F, 
            col = "black", 
            cex = 0.6, 
            cex.axis=1.5, cex.lab=1.5, 
            #main = "QQ plot residuals", 
            cex.main = 1)
dev.off()
##
####



## Model comparisons ##
fit$parameter_estimates$AIC


#####
## Plot key model results ##

# Map of lf probability of capture across time
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category", 
                  PlotName = "POC_lf_yearly",
                  #PlotTitle = "Longfin eel yearly probability of capture in Waitaki, NZ",
                  PlotTitle = "",
                  cex = 0.5, 
                  Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)

# Map of lf probability of capture for each year
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel = "Year",
                  PlotName = "POC_lf",
                  PlotTitle = "Longfin eel P.O.C in Waitaki, NZ",
                  cex = 0.75,  
                  Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)

# Map of spatio-temporal variation of longfin eel probability of capture
plot_maps_network(plot_set = 5, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category",
                  PlotName = "Epsilon_lf",
                  PlotTitle = "Spatio-temporal variation of longfin eel probability of capture",
                  cex = 0.5, 
                  arrows=F)



# #Plot individual covariate effects
# #covar_names_to_use <- covars_all #Add names here
# covar_names_to_use <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
#                         "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam")
# 
# n_p = length(covar_names_to_use)
# 
# 
# ##Covariate effects across time
# cov_effects_array <- plot_maps_network(plot_set = 10, 
#                                        fit = fit, 
#                                        Sdreport = fit$parameter_estimates$SD, 
#                                        TmbData = fit$data_list, 
#                                        spatial_list = fit$spatial_list, 
#                                        DirName = fig, 
#                                        Panel = "category",
#                                        category_names="Longfin eel",
#                                        covar_names=covar_names_to_use,
#                                        PlotName="Covariate_effects",
#                                        PlotTitle = "Individual covariate effect",
#                                        cex = 0.5, 
#                                        #Zlim = c(0,1), 
#                                        arrows=F, 
#                                        pch=15,
#                                        n_p = n_p,
#                                        which_np_touse = which(covar_names_to_use == covar_names_to_use),
#                                        which_cat_cov_toplot=1)
# 
# saveRDS(cov_effects_array, file.path(fig, "cov_effects_array.rds"))
# #cov_effects_array <- readRDS(file.path(fig, "cov_effects_array.rds"))
# 
# 
# 
# 
# ## Plot raw (standardised) covariate values across all time ##
# 
# Network_sz_EN <- data.frame('parent_s'=fit$data_list$parent_s, 'child_s'=fit$data_list$child_s, fit$spatial_list$latlon_g)
# l2 <- lapply(1:nrow(Network_sz_EN), function(x){
#   parent <- Network_sz_EN$parent_s[x]
#   find <- Network_sz_EN %>% filter(child_s == parent)
#   # if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$E_km, 'N2'=find$N_km)
#   # if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
#   if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$Lon, 'N2'=find$Lat)
#   if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
#   return(out)
# })
# l2 <- do.call(rbind, l2)
# 
# 
# for(covar in which(covar_names_to_use %in% covar_names_to_use[c(1:6)])){ #change 'which' depending on what I want to plot - don't include years since dam
#   
#   
#   #Which covariate to plot
#   print(paste0("Plotting: ", covar_names_to_use[covar], " covariate effect"))
#   
#   #Create data frame for plotting
#   xct <- data.frame('value'=cov_effects_array[,covar,1], fit$spatial_list$latlon_g, "category"="Longfin") #'1' is first year - years are the same!
#   
#   Xlim = c(min(xct$Lon),max(xct$Lon))
#   Ylim = c(min(xct$Lat),max(xct$Lat))
#   
#   #Z limits should be set for consistency between years
#   inp_Zlim = quantile(cov_effects_array[,covar,1], prob = c(0,1), na.rm=TRUE)
#   
#   p <- ggplot(xct)
#   
#   p <- p + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = E2, yend = N2), col="gray92")  
#   
#   p <- p +
#     geom_point(aes(x = Lon, y = Lat, color = value), cex = 0.5, pch=19) +
#     scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
#     coord_cartesian(xlim = Xlim, ylim = Ylim) +
#     scale_x_continuous(breaks=quantile(xct$Lon, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$Lon, prob=c(0.1,0.5,0.9)),1)) +
#     # guides(color=guide_legend(title=plot_codes[plot_num])) +
#     mytheme() +
#     xlab("Longitude") + ylab("Latitude")
#   
#   p <- p + ggtitle(wrapper(paste0("Individual covariate effect - ", covar_names_to_use[covar]), width = 50))
#   ggsave(file.path(fig, paste0("Covariate_effects_", covar, ".png")), p, width=8,height=8)
#   
#   
# }
# ####


# Covariate names in plotting
covar_names_to_use <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
                        "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam", "Years since dam squared")

Raw_cov_array = plot_maps_network(plot_set = 9,
                                  fit = fit,
                                  Sdreport = fit$parameter_estimates$SD,
                                  TmbData = fit$data_list,
                                  spatial_list = fit$spatial_list,
                                  DirName = fig,
                                  Panel = "category",
                                  PlotName = "Raw_covariate_values",
                                  PlotTitle = "Standardised covariate values",
                                  covar_names = covar_names_to_use,
                                  cex = 0.5,
                                  arrows=F,
                                  pch=15)

saveRDS(Raw_cov_array, file.path(fig, "Raw_cov_array.rds"))


Network_sz_EN <- data.frame('parent_s'=fit$data_list$parent_s, 'child_s'=fit$data_list$child_s, fit$spatial_list$latlon_g)
l2 <- lapply(1:nrow(Network_sz_EN), function(x){
  parent <- Network_sz_EN$parent_s[x]
  find <- Network_sz_EN %>% filter(child_s == parent)
  # if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$E_km, 'N2'=find$N_km)
  # if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
  if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$Lon, 'N2'=find$Lat)
  if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
  return(out)
})
l2 <- do.call(rbind, l2)


covs_to_plot <- covar2[!(covar2 %in% c("YearsSinceDam", "YearsSinceDam2"))]

for(covar in 1:length(covs_to_plot)){ #change 'which' depending on what I want to plot - don't include years since dam
  
  
  #Which covariate to plot
  print(paste0("Plotting: ", covar_names_to_use[covar], " raw covariate"))
  
  #Create data frame for plotting
  xct <- data.frame('value'=Raw_cov_array[,covar,1], fit$spatial_list$latlon_g) #'1' is first year - years are the same!
  
  Xlim = c(min(xct$Lon),max(xct$Lon))
  Ylim = c(min(xct$Lat),max(xct$Lat))
  
  #Z limits
  inp_Zlim = quantile(Raw_cov_array[,covar,1], prob = c(0,1), na.rm=TRUE)
  
  p <- ggplot(xct)
  
  p <- p + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = E2, yend = N2), col="gray92")
  
  p <- p +
    geom_point(aes(x = Lon, y = Lat, color = value), cex = 0.5, pch=19) +
    scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
    coord_cartesian(xlim = Xlim, ylim = Ylim) +
    #scale_x_continuous(breaks=quantile(xct$Lon, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$Lon, prob=c(0.1,0.5,0.9)),1)) +
    # guides(color=guide_legend(title=plot_codes[plot_num])) +
    mytheme() +
    xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)")
  
  p <- p + ggtitle(wrapper(covar_names_to_use[covar], width = 50))
  ggsave(file.path(fig, paste0(covar, "_raw_values", ".png")), p, width=8,height=8)
  
  
}



waitaki_river_length = (sum(network$length))

plot_range_index_riverlength(Report = fit$Report, 
                             TmbData = fit$data_list, 
                             Sdreport = fit$parameter_estimates$SD, 
                             Znames = colnames(fit$data_list$Z_gm), 
                             PlotDir = fig, 
                             year_labels = as.character(c(fit$year_labels)), 
                             #years_to_plot = fit$years_to_plot,
                             FileName_EffArea = paste0(fig,"/River_length_occupied.png"),
                             use_biascorr = TRUE, 
                             category_names = "Longfin_eels",
                             total_river_length = waitaki_river_length,
                             Yrange = c(0,85),
                             Xrange = c(1934, 2018))

plot_index(Index_ctl = fit$Report$Index_ctl,
           year_labels = as.character(c(fit$year_labels)),
           DirName = paste0(fig, "/"),
           Yrange = c(0, 8))

#####














# ##########################
# ## model6
# ## spatial, temporal, 
# ## beta1 IID, lognormal dist, 
# ## method/agency covariate
# ###########################
# path <- file.path(res_dir, "model6_AC")
# dir.create(path, showWarnings = FALSE)
# fig <- file.path(path, "figures")
# dir.create(fig)
# 
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)
# 
# Data_inp <- Data_Geostat
# X_gtp_inp <- NULL
# X_itp_inp <- NULL
# Xconfig_zcp_inp <- NULL
# Q_ik_inp <- Q_ik_methodagency
# 
# FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
# RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
# ObsModel <- cbind("PosDist" = 1, "Link" = 0)
# OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
# Options <- c("Calculate_range" = 1,
#              "Calculate_effective_area" = 1)
# 
# settings <- make_settings(Version = "VAST_v8_2_0", 
#                           n_x = nrow(Network_sz), 
#                           Region = "Stream_network",
#                           FieldConfig = FieldConfig,
#                           RhoConfig = RhoConfig,
#                           ObsModel = ObsModel,
#                           OverdispersionConfig = OverdispersionConfig,
#                           Options = Options,
#                           purpose = "index",
#                           fine_scale = FALSE, 
#                           bias.correct = FALSE)
# settings$Method <- "Stream_network"
# settings$grid_size_km <- 1
# 
# 
# fit0 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   run_model = FALSE)
# 
# Par <- fit0$tmb_list$Parameters
# Map <- fit0$tmb_list$Map
# Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
# Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))
# 
# fit1 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   model_args = list(Map = Map, Par = Par),
#                   optimize_args = list(getsd = FALSE, newtonsteps = 0),
#                   test_fit = FALSE)
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# 
# fit <- fit_model(settings = settings,
#                  Lat_i = Data_inp[,"Lat"],
#                  Lon_i = Data_inp[,"Lon"],
#                  t_i = Data_inp[,"Year"],
#                  c_i = rep(0, nrow(Data_inp)),
#                  b_i = Data_inp[,"Catch_KG"],
#                  a_i = Data_inp[,"AreaSwept_km2"],
#                  working_dir = path,
#                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                  spatial_args = list(Network_sz_LL = Network_sz_LL),
#                  Network_sz = Network_sz,
#                  Xconfig_zcp = Xconfig_zcp_inp,
#                  X_gtp = X_gtp_inp, 
#                  X_itp = X_itp_inp,
#                  Q_ik = Q_ik_inp,
#                  model_args = list(Map = Map, Par = Par),
#                  test_fit = FALSE)
# saveRDS(fit, file.path(path, "Fit.rds"))
# 
# fit <- readRDS(file.path(path, "Fit.rds"))
# 
# plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
# plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")
# 
# map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
# Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
# plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
# plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

##########################
## model7
## Spatial
## Temporal (beta1 IID)
## Spatiotemporal (independent among years)
## Lognormal dist
## Habitat covariates
## Agency catchability covariate
###########################

rm(list=ls())

load("./Waitaki/Data/general_inputs_Waitaki.Rdata")

#path <- file.path(res_dir, "model7_AC")
path <- file.path(res_dir, "model7")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)

ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.dll"), to = path)


Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0

#Take out yearssincedam^2
Xconfig_zcp_inp[,,which(covar2=="YearsSinceDam2")] <- 0

formula_inp <- paste0("~",(paste0(covar2, collapse = "+")))

#Q_ik_inp <- Q_ik_methodagency
Q_ik_inp <- Q_ik_agency #Changing to this as method is only EF now

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1)

settings <- make_settings(#Version = "VAST_v8_2_0", 
                          Version = "VAST_v12_0_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index2",
                          fine_scale = FALSE, 
                          bias.correct = T)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  formula = formula_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- 5
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
TMBhelper::check_estimability(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 c_i = rep(0, nrow(Data_inp)),
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 working_dir = path,
                 extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp,
                 X_gtp = X_gtp_inp, 
                 X_itp = X_itp_inp,
                 Q_ik = Q_ik_inp,
                 model_args = list(Map = Map, Par = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

#fit <- readRDS(file.path(path, "Fit.rds"))



## OLD ##
# plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
# plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")
# 
# map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
# Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
# plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
# plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)
####




source("./Code/funcs.R")


network <- network %>% rename("Lat"=lat, "Lon"=long)


## Check model ##
#plot_anisotropy(FileName=file.path(fig,"Aniso.png"), Obj=fit$tmb_list$Obj)
dharmaRes = summary(fit, what="residuals", working_dir=paste0(fig,"/"), type=1)
#dharmaRes <- readRDS(file.path(path, "dharmaRes.rds"))
plot_residuals(residuals=dharmaRes$scaledResiduals, fit=fit, Data_inp = Data_inp, network=network, 
               coords="lat_long", save_dir = fig)

saveRDS(dharmaRes, file.path(path, "dharmaRes.rds"))

####
## Arnaud's residual plots ##

######## Check whether observed encounter frequencies for either low or high probability samples 
######## are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_Geostat, DirName = paste0(fig, '/'))


######## Run diagnostics for the positive-catch-rate component of the delta-Lognormal model
######## We can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
######## A good Q-Q plot will have residuals along the one-to-one line
source("./Code/plotQuantileDiagnostic.r")
Q = plotQuantileDiagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive",
                           FileName_Phist = "Posterior_Predictive-Histogram", save_dir = paste0(fig, "/QQ_Fn/" ),
                           FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")

####


## Model comparisons ##
fit$parameter_estimates$AIC


#####
## Plot key model results ##

# Map of lf probability of capture across time
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category", 
                  PlotName = "POC_lf_yearly",
                  PlotTitle = "Longfin eel yearly probability of capture in Waitaki, NZ",
                  cex = 0.5, 
                  Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)

# Map of lf probability of capture for each year
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel = "Year",
                  PlotName = "POC_lf",
                  PlotTitle = "Longfin eel P.O.C in Waitaki, NZ",
                  cex = 0.75,  
                  #Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)

# Map of spatio-temporal variation of longfin eel probability of capture
plot_maps_network(plot_set = 5, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category",
                  PlotName = "Epsilon_lf",
                  PlotTitle = "Spatio-temporal variation of longfin eel probability of capture",
                  cex = 0.5, 
                  arrows=F)



# Covariate names in plotting
covar_full_names_all <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
                          "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam", "Years since dam squared")



Raw_cov_array = plot_maps_network(plot_set = 9,
                                  fit = fit,
                                  Sdreport = fit$parameter_estimates$SD,
                                  TmbData = fit$data_list,
                                  spatial_list = fit$spatial_list,
                                  DirName = fig,
                                  Panel = "category",
                                  PlotName = "Raw_covariate_values",
                                  PlotTitle = "Standardised covariate values",
                                  covar_names = covar_full_names_all,
                                  cex = 0.5,
                                  arrows=F,
                                  pch=15)



#Plot individual covariate effects
#covar_names_to_use <- covars_all #Add names here
covar_names_to_use <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
                        "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam")

n_p = length(covar_names_to_use)


##Covariate effects across time
cov_effects_array <- plot_maps_network(plot_set = 10, 
                                       fit = fit, 
                                       Sdreport = fit$parameter_estimates$SD, 
                                       TmbData = fit$data_list, 
                                       spatial_list = fit$spatial_list, 
                                       DirName = fig, 
                                       Panel = "category",
                                       category_names="Longfin eel",
                                       covar_names=covar_names_to_use,
                                       PlotName="Covariate_effects",
                                       PlotTitle = "Individual covariate effect",
                                       cex = 0.5, 
                                       #Zlim = c(0,1), 
                                       arrows=F, 
                                       pch=15,
                                       n_p = length(covar_names_to_use),
                                       which_np_touse = which(covar_names_to_use == covar_names_to_use),
                                       which_cat_cov_toplot=1)



plot_range_index(Report = fit$Report, 
                 TmbData = fit$data_list, 
                 Sdreport = fit$parameter_estimates$SD, 
                 Znames = colnames(fit$data_list$Z_gm), 
                 PlotDir = fig, 
                 Year_Set = fit$year_labels, 
                 use_biascorr = TRUE, 
                 category_names = "Longfin_eels")
#####






##########################
## model7b
## Spatial
## Temporal (beta1 random walk)
## Spatiotemporal  (independent among years)
## Lognormal dist
## Habitat covariates
## Agency catchability covariate
###########################

rm(list=ls())

load("./Waitaki/Data/general_inputs_Waitaki.Rdata")

#path <- file.path(res_dir, "model7_AC")
path <- file.path(res_dir, "model7b")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)

ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v12_0_0.dll"), to = path)


Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0

#Take out yearssincedam^2
Xconfig_zcp_inp[,,which(covar2=="YearsSinceDam2")] <- 0

formula_inp <- paste0("~",(paste0(covar2, collapse = "+")))

#Q_ik_inp <- Q_ik_methodagency
Q_ik_inp <- Q_ik_agency #Changing to this as method is only EF now

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 2, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1)

settings <- make_settings(#Version = "VAST_v8_2_0", 
  Version = "VAST_v12_0_0", 
  n_x = nrow(Network_sz), 
  Region = "Stream_network",
  FieldConfig = FieldConfig,
  RhoConfig = RhoConfig,
  ObsModel = ObsModel,
  OverdispersionConfig = OverdispersionConfig,
  Options = Options,
  purpose = "index2",
  fine_scale = FALSE, 
  bias.correct = T)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  formula = formula_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- 5
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
TMBhelper::check_estimability(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 c_i = rep(0, nrow(Data_inp)),
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 working_dir = path,
                 extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp,
                 X_gtp = X_gtp_inp, 
                 X_itp = X_itp_inp,
                 Q_ik = Q_ik_inp,
                 model_args = list(Map = Map, Par = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

#fit <- readRDS(file.path(path, "Fit.rds"))



## OLD ##
# plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
# plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")
# 
# map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
# Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
# plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
# plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)
####




source("./Code/funcs.R")


network <- network %>% rename("Lat"=lat, "Lon"=long)


## Check model ##
#plot_anisotropy(FileName=file.path(fig,"Aniso.png"), Obj=fit$tmb_list$Obj)
dharmaRes = summary(fit, what="residuals", working_dir=paste0(fig,"/"), type=1)
#dharmaRes <- readRDS(file.path(path, "dharmaRes.rds"))
plot_residuals(residuals=dharmaRes$scaledResiduals, fit=fit, Data_inp = Data_inp, network=network, 
               coords="lat_long", save_dir = fig)

saveRDS(dharmaRes, file.path(path, "dharmaRes.rds"))

####
## Arnaud's residual plots ##

######## Check whether observed encounter frequencies for either low or high probability samples 
######## are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_Geostat, DirName = paste0(fig, '/'))


######## Run diagnostics for the positive-catch-rate component of the delta-Lognormal model
######## We can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
######## A good Q-Q plot will have residuals along the one-to-one line
source("./Code/plotQuantileDiagnostic.r")
Q = plotQuantileDiagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive",
                           FileName_Phist = "Posterior_Predictive-Histogram", save_dir = paste0(fig, "/QQ_Fn/" ),
                           FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")

####


## Model comparisons ##
fit$parameter_estimates$AIC


#####
## Plot key model results ##

# Map of lf probability of capture across time
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category", 
                  PlotName = "POC_lf_yearly",
                  PlotTitle = "Longfin eel yearly probability of capture in Waitaki, NZ",
                  cex = 0.5, 
                  Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)

# Map of lf probability of capture for each year
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel = "Year",
                  PlotName = "POC_lf",
                  PlotTitle = "Longfin eel P.O.C in Waitaki, NZ",
                  cex = 0.75,  
                  #Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)

# Map of spatio-temporal variation of longfin eel probability of capture
plot_maps_network(plot_set = 5, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category",
                  PlotName = "Epsilon_lf",
                  PlotTitle = "Spatio-temporal variation of longfin eel probability of capture",
                  cex = 0.5, 
                  arrows=F)



# Covariate names in plotting
covar_full_names_all <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
                          "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam", "Years since dam squared")



Raw_cov_array = plot_maps_network(plot_set = 9,
                                  fit = fit,
                                  Sdreport = fit$parameter_estimates$SD,
                                  TmbData = fit$data_list,
                                  spatial_list = fit$spatial_list,
                                  DirName = fig,
                                  Panel = "category",
                                  PlotName = "Raw_covariate_values",
                                  PlotTitle = "Standardised covariate values",
                                  covar_names = covar_full_names_all,
                                  cex = 0.5,
                                  arrows=F,
                                  pch=15)



#Plot individual covariate effects
#covar_names_to_use <- covars_all #Add names here
covar_names_to_use <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
                        "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam")

n_p = length(covar_names_to_use)


##Covariate effects across time
cov_effects_array <- plot_maps_network(plot_set = 10, 
                                       fit = fit, 
                                       Sdreport = fit$parameter_estimates$SD, 
                                       TmbData = fit$data_list, 
                                       spatial_list = fit$spatial_list, 
                                       DirName = fig, 
                                       Panel = "category",
                                       category_names="Longfin eel",
                                       covar_names=covar_names_to_use,
                                       PlotName="Covariate_effects",
                                       PlotTitle = "Individual covariate effect",
                                       cex = 0.5, 
                                       #Zlim = c(0,1), 
                                       arrows=F, 
                                       pch=15,
                                       n_p = length(covar_names_to_use),
                                       which_np_touse = which(covar_names_to_use == covar_names_to_use),
                                       which_cat_cov_toplot=1)



plot_range_index(Report = fit$Report, 
                 TmbData = fit$data_list, 
                 Sdreport = fit$parameter_estimates$SD, 
                 Znames = colnames(fit$data_list$Z_gm), 
                 PlotDir = fig, 
                 Year_Set = fit$year_labels, 
                 use_biascorr = TRUE, 
                 category_names = "Longfin_eels")
#####









##########################
## model8
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## method/agency covariate
###########################
path <- file.path(res_dir, "model8_AC")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- NULL
X_itp_inp <- NULL
Xconfig_zcp_inp <- NULL
Q_ik_inp <- Q_ik_methodagency

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- 8
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 c_i = rep(0, nrow(Data_inp)),
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 working_dir = path,
                 extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp,
                 X_gtp = X_gtp_inp, 
                 X_itp = X_itp_inp,
                 Q_ik = Q_ik_inp,
                 model_args = list(Map = Map, Par = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.01, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


# ##########################
# ## model9
# ## spatial, temporal, 
# ## beta1 IID, lognormal dist, 
# ## habitat covariate
# ###########################
# path <- file.path(res_dir, "model9_AC")
# dir.create(path, showWarnings = FALSE)
# fig <- file.path(path, "figures")
# dir.create(fig)
# 
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
# ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.dll"), to = path)
# 
# Data_inp <- Data_Geostat
# X_gtp_inp <- X_gtp_input
# X_itp_inp <- X_itp_input
# Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
# Xconfig_zcp_inp[2,,] <- 0
# Q_ik_inp <- NULL
# 
# FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
# RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
# ObsModel <- cbind("PosDist" = 1, "Link" = 0)
# OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
# Options <- c("Calculate_range" = 1,
#              "Calculate_effective_area" = 1)
# 
# settings <- make_settings(Version = "VAST_v8_2_0", 
#                           n_x = nrow(Network_sz), 
#                           Region = "Stream_network",
#                           FieldConfig = FieldConfig,
#                           RhoConfig = RhoConfig,
#                           ObsModel = ObsModel,
#                           OverdispersionConfig = OverdispersionConfig,
#                           Options = Options,
#                           purpose = "index",
#                           fine_scale = FALSE, 
#                           bias.correct = FALSE)
# settings$Method <- "Stream_network"
# settings$grid_size_km <- 1
# 
# 
# fit0 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   run_model = FALSE)
# 
# Par <- fit0$tmb_list$Parameters
# Map <- fit0$tmb_list$Map
# Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
# Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))
# 
# fit1 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   c_i = rep(0, nrow(Data_inp)),
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   working_dir = path,
#                   extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz,
#                   Xconfig_zcp = Xconfig_zcp_inp,
#                   X_gtp = X_gtp_inp, 
#                   X_itp = X_itp_inp,
#                   Q_ik = Q_ik_inp,
#                   model_args = list(Map = Map, Par = Par),
#                   optimize_args = list(getsd = FALSE, newtonsteps = 0),
#                   test_fit = FALSE)
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# 
# fit <- fit_model(settings = settings,
#                  Lat_i = Data_inp[,"Lat"],
#                  Lon_i = Data_inp[,"Lon"],
#                  t_i = Data_inp[,"Year"],
#                  c_i = rep(0, nrow(Data_inp)),
#                  b_i = Data_inp[,"Catch_KG"],
#                  a_i = Data_inp[,"AreaSwept_km2"],
#                  working_dir = path,
#                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                  spatial_args = list(Network_sz_LL = Network_sz_LL),
#                  Network_sz = Network_sz,
#                  Xconfig_zcp = Xconfig_zcp_inp,
#                  X_gtp = X_gtp_inp, 
#                  X_itp = X_itp_inp,
#                  Q_ik = Q_ik_inp,
#                  model_args = list(Map = Map, Par = Par),
#                  optimize_args = list(startpar = fit1$parameter_estimates$par),
#                  test_fit = FALSE)
# saveRDS(fit, file.path(path, "Fit.rds"))
# 
# fit <- readRDS(file.path(path, "Fit.rds"))
# 
# plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
# plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")
# 
# map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
# Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
# plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
# plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)



df <- data.frame("Model" = c(#"model1",
  "model2",
  "model3",
  "model3b",
  "model4",
  "model5",
  "model6",
  "model7",
  "model8",
  "model9"))
df$AIC <- NA
for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df$deltaAIC <- df$AIC - min(df$AIC)
df[order(df$AIC),]
#     Model         AIC   deltaAIC
# 2  model3 -159.979788   0.000000
# 8  model8 -158.511960   1.467829
# 7  model7 -148.587055  11.392733 # Lowest AIC whilst still retaining covariates
# 5  model5  -17.984667 141.995121
# 1  model2  -17.571630 142.408159
# 4  model4    4.240237 164.220025
# 9  model9    9.113004 169.092792
# 6  model6   38.112503 198.092291
# 3 model3b   54.638451 214.618239


### MODEL A: SST, no habitat, gear
modelA <- readRDS(file.path(res_dir, "model3_AC", "Fit.rds"))

### MODEL B: SST, no habitat, gear
modelB <- readRDS(file.path(res_dir, "model3b_AC", "Fit.rds"))

### MODEL C: ST, habitat, gear
modelC <- readRDS(file.path(res_dir, "model2_AC", "Fit.rds"))


## compare maps
ep_byModel <- lapply(1:3, function(x){
  if(x == 1){
    Report <- modelA$Report
    year_labels = modelA$year_labels
    years_to_plot = modelA$years_to_plot
    spatial_list <- modelA$spatial_list
    name <- "Spatiotemporal variation,\n Method covariates"
  }
  if(x == 2){
    Report <- modelB$Report
    year_labels = modelB$year_labels
    years_to_plot = modelB$years_to_plot
    spatial_list <- modelB$spatial_list
    name <- "Spatiotemporal variation w/smoother,\nMethod covariates"
  }
  if(x == 3){
    Report <- modelC$Report
    year_labels = modelC$year_labels
    years_to_plot = modelC$years_to_plot
    spatial_list <- modelC$spatial_list
    name <- "No spatiotemporal variation,\n Habitat & Method covariates"
  }
  Array_xct = Report$R1_gcy
  dimnames(Array_xct) <- list(Node = 1:dim(Array_xct)[1], Category = "Longfin eels", Year = year_labels)
  xct <- reshape2::melt(Array_xct) %>% mutate(Model = name)
  xctll <- full_join(xct, cbind.data.frame("Node" = 1:spatial_list$n_g,spatial_list$latlon_g))
  return(xctll)
})
ep <- do.call(rbind, ep_byModel)

plot_ep <- ep %>% filter(Year %in% c(1965, 1995, 2018))

p <- ggplot(plot_ep) +
  geom_point(aes(x = Lon, y = Lat, color = value), cex = 0.5, alpha = 0.75) +
  scale_color_distiller(palette = "Spectral") +
  facet_grid(Year ~ Model) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_dir, "Compare_encounter_prob_maps.png"), p, height = 10, width = 12)

## effective area occupied
eao_byModel <- lapply(1:3, function(x){
  if(x == 1){
    SD <- TMB::summary.sdreport(modelA$parameter_estimates$SD)
    TmbData <- modelA$data_list
    year_labels = modelA$year_labels
    years_to_plot = modelA$years_to_plot
    spatial_list <- modelA$spatial_list
    name <- "Spatiotemporal variation,\n Method covariates"
  }
  if(x == 2){
    SD <- TMB::summary.sdreport(modelB$parameter_estimates$SD)
    TmbData <- modelB$data_list
    year_labels = modelB$year_labels
    years_to_plot = modelB$years_to_plot
    spatial_list <- modelB$spatial_list
    name <- "Spatiotemporal variation w/smoother,\nMethod covariates"
  }
  if(x == 3){
    SD <- TMB::summary.sdreport(modelC$parameter_estimates$SD)
    TmbData <- modelC$data_list
    year_labels = modelC$year_labels
    years_to_plot = modelC$years_to_plot
    spatial_list <- modelC$spatial_list
    name <- "No spatiotemporal variation,\n Habitat & Method covariates"
  }
  EffectiveName = "effective_area_cyl"
  SD_effective_area_ctl = SD_log_effective_area_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
  use_biascorr = TRUE
  # Extract estimates
  SD_effective_area_ctl = SD_log_effective_area_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
  # Effective area
  if( use_biascorr==TRUE && "unbiased"%in%names(SD) ){
    SD_effective_area_ctl[] = SD[which(rownames(SD)==EffectiveName),c('Est. (bias.correct)','Std. Error')]
  }
  if( !any(is.na(SD_effective_area_ctl)) ){
    message("Using bias-corrected estimates for effective area occupied (natural scale)...")
  }else{
    message("Not using bias-corrected estimates for effective area occupied (natural scale)...")
    SD_effective_area_ctl[] = SD[which(rownames(SD)==EffectiveName),c('Estimate','Std. Error')]
  }
  # Log-Effective area
  if( use_biascorr==TRUE && "unbiased"%in%names(SD) ){
    SD_log_effective_area_ctl[] = SD[which(rownames(SD)==paste0("log_",EffectiveName)),c('Est. (bias.correct)','Std. Error')]
  }
  if( !any(is.na(SD_log_effective_area_ctl)) ){
    message("Using bias-corrected estimates for effective area occupied (log scale)...")
  }else{
    message("Not using bias-corrected estimates for effective area occupied (log scale)...")
    SD_log_effective_area_ctl[] = SD[which(rownames(SD)==paste0("log_",EffectiveName)),c('Estimate','Std. Error')]
  }
  
  Index_ctl=array(SD_log_effective_area_ctl[,,,'Estimate'],dim(SD_log_effective_area_ctl)[1:3])
  dimnames(Index_ctl) <- list(Category = "Longfin eel", Year = year_labels, Stratum = NA)
  
  sd_Index_ctl=array(SD_log_effective_area_ctl[,,,'Std. Error'],dim(SD_log_effective_area_ctl)[1:3])
  dimnames(sd_Index_ctl) <- list(Category = "Longfin eel", Year = year_labels, Stratum = NA)
  
  df1 <- reshape2::melt(Index_ctl) %>% rename("Estimate" = value)
  df2 <- reshape2::melt(sd_Index_ctl) %>% rename("SD" = value)
  df <- full_join(df1, df2) %>% mutate(Model = name)
  return(df)
})
eao <- do.call(rbind, eao_byModel)

p <- ggplot(eao) +
  geom_segment(aes(x = Year, xend = Year, y = Estimate - 1.96 * SD, yend = Estimate + 1.96 * SD), color = "red", lwd = 1.2) +
  geom_point(aes(x = Year, y = Estimate), color = "red", cex = 3) +
  geom_line(aes(x = Year, y = Estimate), color = "red") +
  coord_cartesian(ylim = c(0,max(eao$Estimate + 1.96 * eao$SD)*1.01)) +
  facet_grid(~Model) +
  ylab("Effective area occupied (km^2)") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_dir, "Compare_effective_area_occupied.png"), p, height = 6, width = 14)

### encounter diagnostic
ep_diag_byModel <- lapply(1:3, function(x){
  if(x == 1){
    Report <- modelA$Report
    name <- "Spatiotemporal variation,\n Method covariates"
  }
  if(x == 2){
    Report <- modelB$Report
    name <- "Spatiotemporal variation w/smoother,\nMethod covariates"
  }
  if(x == 3){
    Report <- modelC$Report
    name <- "No spatiotemporal variation,\n Habitat & Method covariates"
  }
  # Get bin for each datum
  cutpoints_z=seq(0,1,length=21)
  z_i = cut( Report$R1_i, breaks=cutpoints_z, include.lowest=TRUE )
  midpoints_z = rowMeans( cbind(cutpoints_z[-1],cutpoints_z[-length(cutpoints_z)]) )
  
  # Get encounter frequency for each bin
  freq_z = tapply( ifelse(Data_Geostat[,'Catch_KG']>0,1,0), INDEX=z_i, FUN=mean )
  
  # Get expectation given model
  num_z = tapply( Report$R1_i, INDEX=z_i, FUN=length )
  mean_z = tapply( Report$R1_i, INDEX=z_i, FUN=mean )
  var_z = tapply( Report$R1_i, INDEX=z_i, FUN=function(vec){sum(vec*(1-vec))} )
  sd_mean_z = sqrt(var_z / num_z^2)
  
  df_z <- data.frame('midpoint' = midpoints_z, 
                     'frequency' = freq_z, 
                     "num" = num_z,
                     'mean' = mean_z,
                     'var' = var_z,
                     'sd_mean' = sd_mean_z,
                     'Model' = name)
  return(df_z)
})
ep_diag <- do.call(rbind, ep_diag_byModel)

p <- ggplot(ep_diag %>% filter(is.na(frequency) == FALSE)) +
  geom_ribbon(aes(x = midpoint, ymin = mean - 1.96 * sd_mean, ymax = mean + 1.96 * sd_mean), fill = "red", color = NA, alpha = 0.25) +
  geom_line(aes(x = midpoint, y = mean), lwd = 1.5, color = "red") +
  geom_point(aes(x = midpoint, y = frequency), cex = 3) +
  geom_abline(aes(slope = 1, intercept = 0), lty = 2) +
  facet_grid(~Model) + 
  xlab("Observed encounter probability") + ylab("Predicted encounter probability") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_dir, "Compare_encounter_diagnostic.png"), p, height = 6, width = 14)
