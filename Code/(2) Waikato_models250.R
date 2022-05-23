########################################## NOTE ##########################################
#                                                                                        #
# Data created in "Assessing_data_needs" project.                                        #
# This is coded in "(a) Create_Waikato_data.R" and "(b) Create_Waikato_hab_data"         #
# I source 
##########################################################################################


###############################################################
#             VAST model using count and size data            #
#              for the Longfin eel in Waikato, NZ             #
#                                                             #
#                         April, 2021                         #
#                      Anthony Charsley                       #
###############################################################


##Setup
rm(list=ls())

#################
#  Packages     #
#################

## IMPORTANT NOTE: This is version 2. Updated 30/11/20 to change habitat covariates and 
##                 hopefully improve model.

# Install packages
#devtools::install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
#devtools::install_github("merrillrudd/VASTPlotUtils")

Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

library(devtools)
library(VAST)
library(VASTPlotUtils)
library(tidyverse)
library(RColorBrewer)
#library(proj4)
#library(foreach)
#library(doParallel)
##


#################
#  Directories  #
#################

#Directories to work from and save to
res_dir <- file.path(getwd(), "Waikato")
data_dir <- file.path(res_dir, "Data")

fig_dir <- file.path(res_dir, "Figures")
dir.create(fig_dir, showWarnings=FALSE)

#########################
# Form the network data #       
#########################
#This is a stream network and defines the spatial domain

####
##Network
network <- readRDS(file.path(data_dir, "Greater_waikato_network_AC.rds")) #%>% mutate('AreaSwept_km2'=dist_s*width)

## format network data
Network_sz = network %>% select(c('parent_s','child_s','dist_s')) # Data frame of nodes and distance between nodes
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat) # Data frame of nodes, distance between nodes and latitude/longitude information
####



###############################
# Form the observational data #       
###############################
#This is the data frame of longfin observations

obs <- readRDS(file.path(data_dir, "Greater_waikato_observations_lf_AC.rds")) 
####



#########################
#Form the habitat data #
#########################

#covar_all <- c("Dist2Coast_FromMid", "upElev", "DamAffected")
covar_all <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm",'DamAffected')

hab_dat <- network %>% select(lat, long, all_of(covar_all)) %>%
  rename("Lat"=lat, "Lon"=long) 



data_list <- list()

for(i in min(obs$Year):max(obs$Year)){
  
  sub1 <- hab_dat %>% mutate(Year=i)
  
  damvar <- sub1[,"DamAffected"]
  damvar[which(damvar==1)] <- i - 1947
  
  sub1 <- sub1 %>% 
    mutate("YearsSinceDam"=damvar) %>%
    select(-c("DamAffected"))
  
  
  data_list[[(i-2008)]] <- sub1
  
}

hab_dat <- do.call(rbind, data_list) 

covar_all <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm","YearsSinceDam")


### loc_elev SHOULD BE TRANFORMED BUT HAS ZEROS AND NEGATIVE VALUES SO I'LL LEAVE IT
sum(hab_dat$loc_elev <= 0)

# ## loc_elev needs to be log transformed ##
# 
# boxplot(log((as.vector(hab_dat[,which(covar_all=="loc_elev")]))))
# 
# hab_dat$loc_elev <- log(hab_dat$loc_elev+0.01)
# 
# ##



#Standardise covariates (models fit more easily when covariates are standardised)
hab_std <- lapply(1:length(covar_all), function(x){
  sub <- hab_dat %>%
    select(covar_all[x])
  sub <- scale(sub)
  return(sub)
})
hab_std <- do.call(cbind, hab_std)

hab_dat[,covar_all] <- hab_std

nodes <- network$child_s[order(network$child_s)]
years <- unique(obs$Year)[order(unique(obs$Year))]

sapply(1:ncol(hab_dat), function(x) length(which(is.na(hab_dat[,x])))/nrow(hab_dat))

hab_dat <- hab_dat[complete.cases(hab_dat),]

###

##### Plot final distributions

for(p in 1:length(covar_all)){
  
  jpeg(paste0(fig_dir,"/Final_cov_dist_", covar_all[p], ".jpeg"), height=8, width=8,units="in", res=600)
  boxplot((as.vector(hab_dat[,which(covar_all==covar_all[p])])), main=paste0("Distribution of ", covar_all[p]))
  dev.off()
  
}


########################
# Form the length data #       
########################
#This part categorises the length data and forms the count data.
#Before this, the data has counts in bins of 10mm

###
# TWO categories ("<250", ">=250")
absent <- obs %>% filter(Category == "Absent") #All the absent data
present <- obs %>% filter(Category != "Absent") #All the present data
present$Category <- as.numeric(present$Category)
quantile(present$Category)

present$Cat250 <- sapply(1:nrow(present), function(x) ifelse(present$Category[x] < 250, "<250mm", ">=250mm")) #Categorise into two categories

#classifying absent data as both under250 and over250 so that <250 and >=250 have same impact of absences
absent1 <- absent %>% mutate(Cat250 = "<250mm") 
absent2 <- absent %>% mutate(Cat250 = ">=250mm")

#bind and group counts
obs_len250 <- rbind.data.frame(present, absent1, absent2) #bind all presences and absences
obs_len250 <- obs_len250 %>% select(Year, nzsegment, count, dist_i, long, lat, child_i, Cat250) %>% #select these variables
  group_by(Year, dist_i, long, lat, child_i, Cat250) %>% #group these together and then summarise by sum(count) in these groups
  summarise(count = sum(count))
obs_len250$CatNum <- sapply(1:nrow(obs_len250), function(x) ifelse(obs_len250$Cat250[x] == "<250mm", 1, 2)) #Number categories
table(obs_len250$Cat250, obs_len250$Year)
###


########################
# Format data for VAST #
########################

Data_len250 <- data.frame( "Catch_KG" = obs_len250$count, 
                           "Year" = as.numeric(obs_len250$Year),
                           "AreaSwept_km2" = (150/1000), 
                           "Lat" = obs_len250$lat, 
                           "Lon" = obs_len250$long, 
                           "Pass" = 0, 
                           "Knot" = obs_len250$child_i,
                           "Category" = obs_len250$Cat250,
                           "CategoryNum" = obs_len250$CatNum) %>% #include data in likelihood
  mutate(Density = Catch_KG /AreaSwept_km2) #density is the number of eels / area surveyed


###########################################
# Plot of the count data with the cut-off #
###########################################

#Plot
pnet <- ggplot(Network_sz_LL) + 
  geom_point(aes(x = Lon, y = Lat), col = "gray", alpha = 0.6) + 
  xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  theme_bw(base_size = 14)

pobs3 <- pnet +
  geom_point(data = Data_len250, aes(x = Lon, y = Lat, size = Density, color = Category), alpha = 0.6) +
  scale_color_brewer(palette = "Set1") +
  guides(size = guide_legend(title = "Catch rate")) +
  facet_wrap(.~Year)
ggsave(file.path(fig_dir, "Obs_length_cut250.png"), pobs3, height = 7, width = 8)




########
# Save #
########
save.image(file.path(data_dir, "general_inputs250_Waikato.Rdata"))







##################################
# Model 1
##################################


rm(list=ls())


#Load data
load("./Waikato/Data/general_inputs250_Waikato.Rdata")

FieldConfig <-  c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
RhoConfig <- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 2, "Epsilon2" = 2) #spatio-temporal terms a random effect autoregressive among years

#ObsModel <- c(2,0) #conventional delta model
ObsModel <- c(2,1) #gamma distribution with 'Poisson-link'
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0) #can ignore as vessel effects weren't used
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1) #additional calculations

#VAST model settings
settings <- make_settings(n_x = nrow(Network_sz),
                          Region = "Stream_network",
                          purpose = "index2",
                          fine_scale = FALSE,
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          OverdispersionConfig = OverdispersionConfig,
                          ObsModel = ObsModel,
                          bias.correct = T,
                          Options = Options,
                          Version = "VAST_v12_0_0")
settings$Method <- "Stream_network"
settings$grid_size_km <- 1
##

#Covariates to use. Using all covariates
covar_touse <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm","YearsSinceDam")

n_p <- length(covar_touse) #number of covariates



#### Data input ######
Data_inp <- Data_len250 #the input data
Data_inp$Catch_KG <- Data_inp$Density #Use density as the outcome
######################



cov_input <- hab_dat
X1config_inp <- array(0, dim = c(2,n_p)) #all turned off
X2config_inp <- array(0, dim = c(2,n_p)) #all turned off

#Now turn the ones I want on
X1config_inp[,which(covar_touse=="Dist2Coast_FromMid")] <- 1 
X2config_inp[,which(covar_touse=="Dist2Coast_FromMid")] <- 1 

X1config_inp[,which(covar_touse=="loc_elev")] <- 1 
X2config_inp[,which(covar_touse=="loc_elev")] <- 1 

X1config_inp[,which(covar_touse=="YearsSinceDam")] <- 1 
X2config_inp[,which(covar_touse=="YearsSinceDam")] <- 1 


X1_formula_inp = paste0("~",(paste0(covar_touse, collapse = "+"))) # a symbolic description of the model to be fitted
X2_formula_inp = paste0("~",(paste0(covar_touse, collapse = "+"))) # a symbolic description of the model to be fitted




#Saving directories
path <- file.path(getwd(), paste0("Waikato/Model1_250"))
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "Figures")
dir.create(fig, showWarnings = FALSE)


#VAST files needed
ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.o"), to = path)
ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.dll"), to = path)


start=Sys.time() #measure how long it takes

## fit0 - Compile the model and set up the model structure specified in settings.
fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
                  working_dir = path,
                  
                  covariate_data = cov_input,
                  X1config_cp = X1config_inp,
                  X2config_cp = X2config_inp,
                  X1_formula = X1_formula_inp,
                  X2_formula = X2_formula_inp,
                  
                  run_model = FALSE, #model isn't run here
                  
                  extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], 
                                                               "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map

## fit1 - check if the model parameters are identifiable.
fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
                  working_dir = path,
                  
                  covariate_data = cov_input,
                  X1config_cp = X1config_inp,
                  X2config_cp = X2config_inp,
                  X1_formula = X1_formula_inp,
                  X2_formula = X2_formula_inp,
                  
                  
                  extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz, 
                  model_args = list(Map = Map, Parameters = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0), #SD isn't estimated
                  test_fit = FALSE)
TMBhelper::check_estimability(fit1$tmb_list$Obj)
#check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) #Check model fit
saveRDS(fit1, file.path(path, "fit1.rds"))

## fit - Run the model, estimating standard errors (should converge if model is checked properly)
fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
                 working_dir = path,
                 
                 
                 covariate_data = cov_input,
                 X1config_cp = X1config_inp,
                 X2config_cp = X2config_inp,
                 X1_formula = X1_formula_inp,
                 X2_formula = X2_formula_inp,
                 
                 extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 model_args = list(Map = Map, Parameters = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)

saveRDS(fit, file.path(path, "Fit.rds"))
end=Sys.time()
time = end - start ; time


#fit <- readRDS(file.path(path, "Fit.rds"))
fit$parameter_estimates$diagnostics



## New plots based on updated functions (18/11/20) ##
source("./Code/funcs.r")


Network_sz_EN <- data.frame('parent_s'=fit$data_list$parent_s, 'child_s'=fit$data_list$child_s, fit$spatial_list$latlon_g)
l2 <- lapply(1:nrow(Network_sz_EN), function(x){
  parent <- Network_sz_EN$parent_s[x]
  find <- Network_sz_EN %>% filter(child_s == parent)
  if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'Lon2'=find$Lon, 'Lat2'=find$Lat)
  if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'Lon2'=NA, 'Lat2'=NA)
  return(out)
})
l2 <- do.call(rbind, l2)



#####
## Check model ##
dharmaRes = summary(fit, what="residuals", working_dir=paste0(fig,"/"), type=1)
#dharmaRes <- readRDS(file.path(path, "dharmaRes.rds"))
plot_residuals(residuals = dharmaRes$scaledResiduals, fit = fit, save_dir=fig,
               Data_inp = Data_inp, network=Network_sz_LL, coords="lat_long")

saveRDS(dharmaRes, file.path(path, "dharmaRes.rds"))

####
## Arnaud's residual plots ##

######## Check whether observed encounter frequencies for either low or high probability samples 
######## are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = paste0(fig, '/'))


######## Run diagnostics for the positive-catch-rate component of the delta-Lognormal model
######## We can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
######## A good Q-Q plot will have residuals along the one-to-one line
source("./Code/plotQuantileDiagnostic.r")
Q = plotQuantileDiagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive",
                           FileName_Phist = "Posterior_Predictive-Histogram", save_dir = paste0(fig, "/QQ_Fn/" ),
                           FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")

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
            cex.axis=1.5, cex.lab=1.5 
            #main = "QQ plot residuals", 
            #cex.main = 1
            )
dev.off()
##
####

######

dens <- quantile(log(fit$Report$D_gct))

# #Index for count
# plot_index(Index_ctl = fit$Report$Index_ctl,
#            #fit = fit, 
#            #Sdreport = fit$parameter_estimates$SD, 
#            DirName = fig, 
#            interval_width = 1.96,
#            category_names = c("<250",">=250")#, 
#            #use_biascorr=TRUE, 
#            #Plot_suffix = "Count"
#            )
# 
# #Index for biomass
# plot_index(fit = fit, 
#            Sdreport = fit$parameter_estimates$SD, 
#            DirName = fig, 
#            interval_width = 1.96,
#            category_names = c("<250",">=250"), 
#            use_biascorr=TRUE, 
#            Plot_suffix = "Biomass")

# Map of lf probability of capture across time
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category", 
                  PlotName = "POC_lf_yearly",
                  PlotTitle = "Longfin eel yearly probability of capture in Waikato, NZ",
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
                  PlotTitle = "Longfin eel P.O.C in Waikato, NZ",
                  cex = 0.75,  
                  Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)



# Log-predicted density by category (<250, >=250)
plot_maps_network(plot_set = 3, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel="Category",
                  category_names = c("<250mm",">=250mm"),
                  PlotName = "lflogdens",
                  PlotTitle = "Log-predicted density",
                  Zlim = c(min(dens),max(dens)),
                  legend=TRUE, 
                  cex=1, 
                  pch=19)


####
# Log-predicted density by year
log_dens_yr <- plot_maps_network(plot_set = 3, 
                                 fit = fit, 
                                 Sdreport = fit$parameter_estimates$SD,
                                 TmbData = fit$data_list, 
                                 spatial_list = fit$spatial_list, 
                                 DirName = fig,
                                 Panel="Year",
                                 category_names = c("<250mm",">=250mm"),
                                 PlotName = "lflogdens",
                                 PlotTitle = "Log-predicted density",
                                 Zlim = c(min(dens),max(dens)),
                                 legend=TRUE, 
                                 cex=0.5, 
                                 pch=19)


#Plot manuscript
# subsetting to extract the years 2009, 2012, 2015 and 2017 for both categories

#set zlims so all on same scale:
#Z limits
inp_Zlim = quantile(log_dens_yr, prob = c(0,1), na.rm=TRUE)


# <250mm #
log_dens_less250_df <- as.data.frame(log_dens_yr[,1,])
colnames(log_dens_less250_df) <- c(2009:2017)
log_dens_less250_df <- cbind(log_dens_less250_df, Network_sz_LL)

log_dens_less250_df <- reshape(log_dens_less250_df, direction = "long", varying = list(1:9), v.names = "Log_density", times = c(2009:2017))
log_dens_less250_df <- log_dens_less250_df %>% rename("Year"="time")

less250_sub1 <- log_dens_less250_df %>% filter(Year == 2009 | Year == 2012)
less250_sub2 <- log_dens_less250_df %>% filter(Year == 2015 | Year == 2017)

#Now plot_1
less250_map1 <- ggplot(less250_sub1)
less250_map1 <- less250_map1 + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = Lon2, yend = Lat2), col="gray90")

less250_map1 <- less250_map1 +
  geom_point(aes(x = Lon, y = Lat, color = Log_density), cex = 0.5,pch=19) +
  scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
  xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "grey90", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Year, ncol=1)
ggsave(file.path(fig, "Waikato_less250_lndensity1.png"), less250_map1, width = 8, height=12, units = "cm")


#Now plot_2
less250_map2 <- ggplot(less250_sub2)
less250_map2 <- less250_map2 + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = Lon2, yend = Lat2), col="gray90")

less250_map2 <- less250_map2 +
  geom_point(aes(x = Lon, y = Lat, color = Log_density), cex = 0.5,pch=19) +
  scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
  xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "grey90", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Year, ncol=1)
ggsave(file.path(fig, "Waikato_less250_lndensity2.png"), less250_map2, width = 8, height=12, units = "cm")


# NEED TO ADD FACET WRAP AND DO THE >=250 PLOTS

# <250mm #
log_dens_greater250_df <- as.data.frame(log_dens_yr[,2,])
colnames(log_dens_greater250_df) <- c(2009:2017)
log_dens_greater250_df <- cbind(log_dens_greater250_df, Network_sz_LL)

log_dens_greater250_df <- reshape(log_dens_greater250_df, direction = "long", varying = list(1:9), v.names = "Log_density", times = c(2009:2017))
log_dens_greater250_df <- log_dens_greater250_df %>% rename("Year"="time")

greater250_sub1 <- log_dens_greater250_df %>% filter(Year == 2009 | Year == 2012)
greater250_sub2 <- log_dens_greater250_df %>% filter(Year == 2015 | Year == 2017)

#Now plot
greater250_map1 <- ggplot(greater250_sub1)
greater250_map1 <- greater250_map1 + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = Lon2, yend = Lat2), col="gray90")

greater250_map1 <- greater250_map1 +
  geom_point(aes(x = Lon, y = Lat, color = Log_density), cex = 0.5,pch=19) +
  scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
  #xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "grey90", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Year, ncol=1)
ggsave(file.path(fig, "Waikato_greater250_lndensity1.png"), greater250_map1, width = 8, height=12, units = "cm")


greater250_map2 <- ggplot(greater250_sub2)
greater250_map2 <- greater250_map2 + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = Lon2, yend = Lat2), col="gray90")

greater250_map2 <- greater250_map2 +
  geom_point(aes(x = Lon, y = Lat, color = Log_density), cex = 0.5,pch=19) +
  scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
  #xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = 6),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "grey90", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Year, ncol=1)
ggsave(file.path(fig, "Waikato_greater250_lndensity2.png"), greater250_map2, width = 8, height=12, units = "cm")


ggplot() + geom_point() +
  scale_color_distiller(palette = "Spectral", limits = inp_Zlim)


#Put together in grid
waikato_logdens_grid <- grid.arrange(less250_map1, greater250_map1, less250_map2, greater250_map2, ncol = 2, nrow=2,
                                     widths = c(2,2,2,4),
                                  bottom = textGrob("Longitude (\u00B0E)"),
                                  left = textGrob("Latitude (\u00B0N)", rot = 90, vjust = 1))
ggsave(file.path(fig, "waikato_logdens_map.png"), waikato_logdens_grid, height = 10)
####




#Spread of longfin eel density over time
jpeg(file.path(fig, "Log_density_spread.jpg"), width = 600, height = 600)
plot(as.vector(log(fit$Report$D_gct[,2,]))~ rep(c(2009:2017), each=nrow(network)),
     xlab="Year", ylab="Log(density)", main="Spread of log(density) of longfin eels >=250mm in the Waikato")
dev.off()

jpeg(file.path(fig, "Density_spread.jpg"), width = 600, height = 600)
plot(as.vector((fit$Report$D_gct[,2,]))~ rep(c(2009:2017), each=nrow(network)),
     xlab="Year", ylab="Density", main="Spread of density of longfin eels >=250mm in the Waikato")
dev.off()


# Total biomass across all categories
plot_maps_network(plot_set = 4,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel="Category",
                  PlotName = "Totalbiomass",
                  PlotTitle = "Total biomass",
                  legend=TRUE,
                  cex=0.5,
                  pch=19)

# Total biomass across all categories
plot_maps_network(plot_set = 4,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel="Year",
                  PlotName = "Totalbiomass",
                  PlotTitle = "Total biomass",
                  legend=TRUE,
                  cex=0.5,
                  pch=19)

# Spatio-temporal variation in probability of capture
plot_maps_network(plot_set = 5,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  PlotName = "Epsilon_POC_lf",
                  PlotTitle = "Spatio-temporal variation of longfin eel P.O.C",
                  legend=TRUE,
                  cex=0.5,
                  pch=19)

# Spatio-temporal variation in positive catch rates
plot_maps_network(plot_set = 6, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  PlotName = "Epsilon_PCRs_lf",
                  PlotTitle = "Spatio-temporal variation of longfin eel P.C.Rs",
                  legend=TRUE, 
                  cex=0.5, 
                  pch=19)



covar_names_to_use <- c("Mean flow (cumecs)", "Distance to coast (km)", "Mean elevation (m)", "Mean slope (degrees)",
                        "CV of annual catchment rainfall (mm)", "Mean Jan. air temp (°C x 10)", "Years since dam")

# Covariates that are included in the model (measured values)
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
                                  cex = 1,
                                  pch=15)

saveRDS(Raw_cov_array, file.path(fig, "Raw_cov_array.rds"))


covs_to_plot <- covar_touse[!(covar_touse %in% c("YearsSinceDam"))]

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
    geom_point(aes(x = Lon, y = Lat, color = value), cex = 2, pch=19) +
    scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
    coord_cartesian(xlim = Xlim, ylim = Ylim) +
    #scale_x_continuous(breaks=quantile(xct$Lon, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$Lon, prob=c(0.1,0.5,0.9)),1)) +
    # guides(color=guide_legend(title=plot_codes[plot_num])) +
    mytheme() +
    xlab("Longitude (\u00B0E)") + ylab("Latitude (\u00B0N)")
  
  p <- p + ggtitle(wrapper(covar_names_to_use[covar], width = 50))
  ggsave(file.path(fig, paste0(covar, "_raw_values", ".png")), p, width=8,height=8)
  
  
}






#Plot individual covariate effects for P.O.C
covar_names_to_use <- c("Distance to coast (km)", "Upstream elevation (m)", "Years since dam") #Add names here
n_p = length(covar_names_to_use)
which_np_touse = which(covar_all %in% c("Dist2Coast_FromMid","loc_elev","YearsSinceDam")) 


plot_maps_network(plot_set = 10,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  covar_names=covar_names_to_use,
                  PlotName="Covariate_effects_poc_lessthan250",
                  PlotTitle = "Individual covariate effect P.O.C",
                  cex = 0.5,
                  pch=15,
                  n_p = n_p,
                  which_np_touse = which_np_touse,
                  which_cat_cov_toplot = 1) #1st category

plot_maps_network(plot_set = 10,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  covar_names=covar_names_to_use,
                  PlotName="Covariate_effects_poc_greaterthan250",
                  PlotTitle = "Individual covariate effect P.O.C",
                  cex = 0.5,
                  pch=15,
                  n_p = n_p,
                  which_np_touse = which_np_touse,
                  which_cat_cov_toplot = 2) #1st category


plot_maps_network(plot_set = 11,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  covar_names=covar_names_to_use,
                  PlotName="Covariate_effects_pcr_lessthan250",
                  PlotTitle = "Individual covariate effect P.C.Rs",
                  cex = 0.5,
                  pch=15,
                  n_p = n_p,
                  which_np_touse = which_np_touse,
                  which_cat_cov_toplot = 1)

plot_maps_network(plot_set = 11,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  covar_names=covar_names_to_use,
                  PlotName="Covariate_effects_pcr_greaterthan250",
                  PlotTitle = "Individual covariate effect P.C.Rs",
                  cex = 0.5,
                  pch=15,
                  n_p = n_p,
                  which_np_touse = which_np_touse,
                  which_cat_cov_toplot = 2)


# # Combined covariate effects for probability of capture
# plot_maps_network(plot_set = 12, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   PlotName = "Combined_covareffect_poc",
#                   PlotTitle = "Combined covariate effect of longfin eel P.O.C",
#                   legend=TRUE, 
#                   cex=0.5, 
#                   pch=19)
# 
# # Combined covariate effects for positive catch rates
# plot_maps_network(plot_set = 13, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   PlotName = "Combined_covareffect_pcr",
#                   PlotTitle = "Combined covariate effect of longfin eel P.C.Rs",
#                   legend=TRUE, 
#                   cex=0.5, 
#                   pch=19)

plot_range_index(Report = fit$Report, 
                 TmbData = fit$data_list, 
                 Sdreport = fit$parameter_estimates$SD, 
                 Znames = colnames(fit$data_list$Z_gm), 
                 PlotDir = fig, 
                 Year_Set = fit$year_labels, 
                 use_biascorr = TRUE, 
                 category_names = c("<250mm",">=250mm"))


waikato_river_length = (sum(network$length))

plot_range_index_riverlength(Report = fit$Report, 
                             TmbData = fit$data_list, 
                             Sdreport = fit$parameter_estimates$SD, 
                             Znames = colnames(fit$data_list$Z_gm), 
                             PlotDir = fig, 
                             year_labels = as.character(fit$year_labels), 
                             #years_to_plot = fit$years_to_plot,
                             FileName_EffArea = paste0(fig,"/River_length_occupied.png"),
                             use_biascorr = TRUE, 
                             category_names = c("<250mm",">=250mm"),
                             total_river_length = waikato_river_length)


####

###################################################
###################################################






# ##################################
# # Model 2 - the same but 2 factors
# ##################################
# 
# 
# rm(list=ls())
# 
# 
# #Load data
# load("./Waikato/Data/general_inputs250_Waikato.Rdata")
# 
# FieldConfig <-  c("Omega1" = 2, "Epsilon1" = 2, "Omega2" = 2, "Epsilon2" = 2)
# RhoConfig <- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 2, "Epsilon2" = 2) #spatio-temporal terms a random effect autoregressive among years
# 
# #ObsModel <- c(2,0) #conventional delta model
# ObsModel <- c(2,1) #gamma distribution with 'Poisson-link'
# OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0) #can ignore as vessel effects weren't used
# Options <- c("Calculate_range" = 1,
#              "Calculate_effective_area" = 1) #additional calculations
# 
# #VAST model settings
# settings <- make_settings(n_x = nrow(Network_sz),
#                           Region = "Stream_network",
#                           purpose = "index2",
#                           fine_scale = FALSE,
#                           FieldConfig = FieldConfig,
#                           RhoConfig = RhoConfig,
#                           OverdispersionConfig = OverdispersionConfig,
#                           ObsModel = ObsModel,
#                           bias.correct = T,
#                           Options = Options,
#                           Version = "VAST_v12_0_0")
# settings$Method <- "Stream_network"
# settings$grid_size_km <- 1
# ##
# 
# #Covariates to use. Using all covariates
# covar_touse <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm","YearsSinceDam")
# 
# n_p <- length(covar_touse) #number of covariates
# 
# 
# 
# #### Data input ######
# Data_inp <- Data_len250 #the input data
# Data_inp$Catch_KG <- Data_inp$Density #Use density as the outcome
# ######################
# 
# 
# 
# cov_input <- hab_dat
# X1config_inp <- array(0, dim = c(2,n_p)) #all turned off
# X2config_inp <- array(0, dim = c(2,n_p)) #all turned off
# 
# #Now turn the ones I want on
# X1config_inp[,which(covar_touse=="Dist2Coast_FromMid")] <- 1 
# X2config_inp[,which(covar_touse=="Dist2Coast_FromMid")] <- 1 
# 
# X1config_inp[,which(covar_touse=="loc_elev")] <- 1 
# X2config_inp[,which(covar_touse=="loc_elev")] <- 1 
# 
# X1config_inp[,which(covar_touse=="YearsSinceDam")] <- 1 
# X2config_inp[,which(covar_touse=="YearsSinceDam")] <- 1 
# 
# 
# X1_formula_inp = paste0("~",(paste0(covar_touse, collapse = "+"))) # a symbolic description of the model to be fitted
# X2_formula_inp = paste0("~",(paste0(covar_touse, collapse = "+"))) # a symbolic description of the model to be fitted
# 
# 
# 
# 
# #Saving directories
# path <- file.path(getwd(), paste0("Waikato/Model2_250"))
# dir.create(path, showWarnings = FALSE)
# fig <- file.path(path, "Figures")
# dir.create(fig, showWarnings = FALSE)
# 
# 
# #VAST files needed
# ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.cpp"), to = path)
# ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.o"), to = path)
# ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.dll"), to = path)
# 
# 
# start=Sys.time() #measure how long it takes
# 
# ## fit0 - Compile the model and set up the model structure specified in settings.
# fit0 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
#                   working_dir = path,
#                   
#                   covariate_data = cov_input,
#                   X1config_cp = X1config_inp,
#                   X2config_cp = X2config_inp,
#                   X1_formula = X1_formula_inp,
#                   X2_formula = X2_formula_inp,
#                   
#                   run_model = FALSE, #model isn't run here
#                   
#                   extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], 
#                                                                "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz)
# 
# Par <- fit0$tmb_list$Parameters
# Map <- fit0$tmb_list$Map
# 
# ## fit1 - check if the model parameters are identifiable.
# fit1 <- fit_model(settings = settings,
#                   Lat_i = Data_inp[,"Lat"],
#                   Lon_i = Data_inp[,"Lon"],
#                   t_i = Data_inp[,"Year"],
#                   b_i = Data_inp[,"Catch_KG"],
#                   a_i = Data_inp[,"AreaSwept_km2"],
#                   c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
#                   working_dir = path,
#                   
#                   covariate_data = cov_input,
#                   X1config_cp = X1config_inp,
#                   X2config_cp = X2config_inp,
#                   X1_formula = X1_formula_inp,
#                   X2_formula = X2_formula_inp,
#                   
#                   
#                   extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                   spatial_args = list(Network_sz_LL = Network_sz_LL),
#                   Network_sz = Network_sz, 
#                   model_args = list(Map = Map, Parameters = Par),
#                   optimize_args = list(getsd = FALSE, newtonsteps = 0), #SD isn't estimated
#                   test_fit = FALSE)
# TMBhelper::check_estimability(fit1$tmb_list$Obj)
# #check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) #Check model fit
# saveRDS(fit1, file.path(path, "fit1.rds"))
# 
# ## fit - Run the model, estimating standard errors (should converge if model is checked properly)
# fit <- fit_model(settings = settings,
#                  Lat_i = Data_inp[,"Lat"],
#                  Lon_i = Data_inp[,"Lon"],
#                  t_i = Data_inp[,"Year"],
#                  b_i = Data_inp[,"Catch_KG"],
#                  a_i = Data_inp[,"AreaSwept_km2"],
#                  c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
#                  working_dir = path,
#                  
#                  
#                  covariate_data = cov_input,
#                  X1config_cp = X1config_inp,
#                  X2config_cp = X2config_inp,
#                  X1_formula = X1_formula_inp,
#                  X2_formula = X2_formula_inp,
#                  
#                  extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
#                  spatial_args = list(Network_sz_LL = Network_sz_LL),
#                  Network_sz = Network_sz,
#                  model_args = list(Map = Map, Parameters = Par),
#                  optimize_args = list(startpar = fit1$parameter_estimates$par),
#                  test_fit = FALSE)
# 
# saveRDS(fit, file.path(path, "Fit.rds"))
# end=Sys.time()
# time = end - start ; time
# 
# 
# #fit <- readRDS(file.path(path, "Fit.rds"))
# fit$parameter_estimates$diagnostics
# 
# 
# 
# ## New plots based on updated functions (18/11/20) ##
# source("./Code/funcs.r")
# 
# #####
# ## Check model ##
# dharmaRes = summary(fit, what="residuals", working_dir=paste0(fig,"/"), type=1)
# #dharmaRes <- readRDS(file.path(path, "dharmaRes.rds"))
# plot_residuals(residuals = dharmaRes$scaledResiduals, fit = fit, save_dir=fig,
#                Data_inp = Data_inp, network=Network_sz_LL, coords="lat_long")
# 
# saveRDS(dharmaRes, file.path(path, "dharmaRes.rds"))
# 
# ####
# ## Arnaud's residual plots ##
# 
# ######## Check whether observed encounter frequencies for either low or high probability samples 
# ######## are within the 95% predictive interval for predicted encounter probability
# Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = paste0(fig, '/'))
# 
# 
# ######## Run diagnostics for the positive-catch-rate component of the delta-Lognormal model
# ######## We can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
# ######## A good Q-Q plot will have residuals along the one-to-one line
# source("./Code/plotQuantileDiagnostic.r")
# Q = plotQuantileDiagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive",
#                            FileName_Phist = "Posterior_Predictive-Histogram", save_dir = paste0(fig, "/QQ_Fn/" ),
#                            FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")
# 
# ####
# 
# 
# ## dharma plots ##
# # Histogram of residuals #
# val = dharmaRes$scaledResiduals
# val[val == 0] = -0.01
# val[val == 1] = 1.01
# 
# jpeg(file.path(fig, "Resid_hist.jpg"), width = 600, height = 600)
# hist(val, 
#      breaks = seq(-0.02, 1.02, len = 53),
#      col = c("red",rep("lightgrey",50), "red"),
#      main = "Hist of DHARMa residuals",
#      xlab = "Residuals (outliers are marked red)",
#      cex.main = 1)
# dev.off()
# ##
# ## QQ plot ##
# jpeg(file.path(fig, "QQplot.jpg"), width = 600, height = 600)
# gap::qqunif(dharmaRes$scaledResiduals,
#             pch=2,
#             bty="n", 
#             logscale = F, 
#             col = "black", 
#             cex = 0.6, 
#             main = "QQ plot residuals", 
#             cex.main = 1)
# dev.off()
# ##
# ####
# 
# ######
# 
# dens <- quantile(log(fit$Report$D_gct))
# 
# # #Index for count
# # plot_index(Index_ctl = fit$Report$Index_ctl,
# #            #fit = fit, 
# #            #Sdreport = fit$parameter_estimates$SD, 
# #            DirName = fig, 
# #            interval_width = 1.96,
# #            category_names = c("<250",">=250")#, 
# #            #use_biascorr=TRUE, 
# #            #Plot_suffix = "Count"
# #            )
# # 
# # #Index for biomass
# # plot_index(fit = fit, 
# #            Sdreport = fit$parameter_estimates$SD, 
# #            DirName = fig, 
# #            interval_width = 1.96,
# #            category_names = c("<250",">=250"), 
# #            use_biascorr=TRUE, 
# #            Plot_suffix = "Biomass")
# 
# # Map of lf probability of capture across time
# plot_maps_network(plot_set = 1, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD, 
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig, 
#                   Panel = "category", 
#                   PlotName = "POC_lf_yearly",
#                   PlotTitle = "Longfin eel yearly probability of capture in Waikato, NZ",
#                   cex = 0.5, 
#                   Zlim = c(0,1), 
#                   arrows=F, 
#                   pch=15)
# 
# # Map of lf probability of capture for each year
# plot_maps_network(plot_set = 1, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD, 
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig,
#                   Panel = "Year",
#                   PlotName = "POC_lf",
#                   PlotTitle = "Longfin eel P.O.C in Waikato, NZ",
#                   cex = 0.75,  
#                   Zlim = c(0,1), 
#                   arrows=F, 
#                   pch=15)
# 
# 
# 
# # Log-predicted density by category (<250, >=250)
# plot_maps_network(plot_set = 3, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig,
#                   Panel="Category",
#                   category_names = c("<250mm",">=250mm"),
#                   PlotName = "lflogdens",
#                   PlotTitle = "Log-predicted density",
#                   Zlim = c(min(dens),max(dens)),
#                   legend=TRUE, 
#                   cex=0.5, 
#                   pch=19)
# 
# # Log-predicted density by year
# plot_maps_network(plot_set = 3, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig,
#                   Panel="Year",
#                   category_names = c("<250mm",">=250mm"),
#                   PlotName = "lflogdens",
#                   PlotTitle = "Log-predicted density",
#                   Zlim = c(min(dens),max(dens)),
#                   legend=TRUE, 
#                   cex=0.5, 
#                   pch=19)
# 
# #Spread of longfin eel density over time
# jpeg(file.path(fig, "Log_density_spread.jpg"), width = 600, height = 600)
# plot(as.vector(log(fit$Report$D_gct[,2,]))~ rep(c(2009:2017), each=nrow(network)),
#      xlab="Year", ylab="Log(density)", main="Spread of log(density) of longfin eels >=250mm in the Waikato")
# dev.off()
# 
# jpeg(file.path(fig, "Density_spread.jpg"), width = 600, height = 600)
# plot(as.vector((fit$Report$D_gct[,2,]))~ rep(c(2009:2017), each=nrow(network)),
#      xlab="Year", ylab="Density", main="Spread of density of longfin eels >=250mm in the Waikato")
# dev.off()
# 
# 
# # Total biomass across all categories
# plot_maps_network(plot_set = 4,
#                   fit = fit,
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list,
#                   spatial_list = fit$spatial_list,
#                   DirName = fig,
#                   Panel="Category",
#                   PlotName = "Totalbiomass",
#                   PlotTitle = "Total biomass",
#                   legend=TRUE,
#                   cex=0.5,
#                   pch=19)
# 
# # Total biomass across all categories
# plot_maps_network(plot_set = 4,
#                   fit = fit,
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list,
#                   spatial_list = fit$spatial_list,
#                   DirName = fig,
#                   Panel="Year",
#                   PlotName = "Totalbiomass",
#                   PlotTitle = "Total biomass",
#                   legend=TRUE,
#                   cex=0.5,
#                   pch=19)
# 
# # Spatio-temporal variation in probability of capture
# plot_maps_network(plot_set = 5,
#                   fit = fit,
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list,
#                   spatial_list = fit$spatial_list,
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   PlotName = "Epsilon_POC_lf",
#                   PlotTitle = "Spatio-temporal variation of longfin eel P.O.C",
#                   legend=TRUE,
#                   cex=0.5,
#                   pch=19)
# 
# # Spatio-temporal variation in positive catch rates
# plot_maps_network(plot_set = 6, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   PlotName = "Epsilon_PCRs_lf",
#                   PlotTitle = "Spatio-temporal variation of longfin eel P.C.Rs",
#                   legend=TRUE, 
#                   cex=0.5, 
#                   pch=19)
# 
# 
# 
# 
# # Covariates that are included in the model (measured values)
# plot_maps_network(plot_set = 9,
#                   fit = fit,
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list,
#                   spatial_list = fit$spatial_list,
#                   DirName = fig,
#                   Panel = "category",
#                   PlotName = "Raw_covariate_values",
#                   PlotTitle = "Standardised covariate values",
#                   covar_names = covar_touse,
#                   cex = 0.5,
#                   pch=15)
# 
# 
# 
# 
# #Plot individual covariate effects for P.O.C
# covar_names_to_use <- c("Distance to coast (km)", "Upstream elevation (m)", "Years since dam") #Add names here
# n_p = length(covar_names_to_use)
# which_np_touse = which(covar_all %in% c("Dist2Coast_FromMid","loc_elev","YearsSinceDam")) 
# 
# 
# plot_maps_network(plot_set = 10,
#                   fit = fit,
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list,
#                   spatial_list = fit$spatial_list,
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   covar_names=covar_names_to_use,
#                   PlotName="Covariate_effects_poc_lessthan250",
#                   PlotTitle = "Individual covariate effect P.O.C",
#                   cex = 0.5,
#                   pch=15,
#                   n_p = n_p,
#                   which_np_touse = which_np_touse,
#                   which_cat_cov_toplot = 1) #1st category
# 
# plot_maps_network(plot_set = 10,
#                   fit = fit,
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list,
#                   spatial_list = fit$spatial_list,
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   covar_names=covar_names_to_use,
#                   PlotName="Covariate_effects_poc_greaterthan250",
#                   PlotTitle = "Individual covariate effect P.O.C",
#                   cex = 0.5,
#                   pch=15,
#                   n_p = n_p,
#                   which_np_touse = which_np_touse,
#                   which_cat_cov_toplot = 2) #1st category
# 
# 
# plot_maps_network(plot_set = 11,
#                   fit = fit,
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list,
#                   spatial_list = fit$spatial_list,
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   covar_names=covar_names_to_use,
#                   PlotName="Covariate_effects_pcr_lessthan250",
#                   PlotTitle = "Individual covariate effect P.C.Rs",
#                   cex = 0.5,
#                   pch=15,
#                   n_p = n_p,
#                   which_np_touse = which_np_touse,
#                   which_cat_cov_toplot = 1)
# 
# plot_maps_network(plot_set = 11,
#                   fit = fit,
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list,
#                   spatial_list = fit$spatial_list,
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   covar_names=covar_names_to_use,
#                   PlotName="Covariate_effects_pcr_greaterthan250",
#                   PlotTitle = "Individual covariate effect P.C.Rs",
#                   cex = 0.5,
#                   pch=15,
#                   n_p = n_p,
#                   which_np_touse = which_np_touse,
#                   which_cat_cov_toplot = 2)
# 
# 
# # # Combined covariate effects for probability of capture
# # plot_maps_network(plot_set = 12, 
# #                   fit = fit, 
# #                   Sdreport = fit$parameter_estimates$SD,
# #                   TmbData = fit$data_list, 
# #                   spatial_list = fit$spatial_list, 
# #                   DirName = fig,
# #                   Panel = "category",
# #                   category_names = c("<250mm",">=250mm"),
# #                   PlotName = "Combined_covareffect_poc",
# #                   PlotTitle = "Combined covariate effect of longfin eel P.O.C",
# #                   legend=TRUE, 
# #                   cex=0.5, 
# #                   pch=19)
# # 
# # # Combined covariate effects for positive catch rates
# # plot_maps_network(plot_set = 13, 
# #                   fit = fit, 
# #                   Sdreport = fit$parameter_estimates$SD,
# #                   TmbData = fit$data_list, 
# #                   spatial_list = fit$spatial_list, 
# #                   DirName = fig,
# #                   Panel = "category",
# #                   category_names = c("<250mm",">=250mm"),
# #                   PlotName = "Combined_covareffect_pcr",
# #                   PlotTitle = "Combined covariate effect of longfin eel P.C.Rs",
# #                   legend=TRUE, 
# #                   cex=0.5, 
# #                   pch=19)
# 
# plot_range_index(Report = fit$Report, 
#                  TmbData = fit$data_list, 
#                  Sdreport = fit$parameter_estimates$SD, 
#                  Znames = colnames(fit$data_list$Z_gm), 
#                  PlotDir = fig, 
#                  Year_Set = fit$year_labels, 
#                  use_biascorr = TRUE, 
#                  category_names = c("<250mm",">=250mm"))
# 
# ####
# 
# ###################################################
# ###################################################
# 



##################################
# Model 2
##################################


rm(list=ls())


#Load data
load("./Waikato/Data/general_inputs250_Waikato.Rdata")

FieldConfig <-  c("Omega1" = 2, "Epsilon1" = 2, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 2, "Epsilon2" = 2) #spatio-temporal terms a random effect autoregressive among years

#ObsModel <- c(2,0) #conventional delta model
ObsModel <- c(2,1) #gamma distribution with 'Poisson-link'
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0) #can ignore as vessel effects weren't used
Options <- c("Calculate_range" = 1,
             "Calculate_effective_area" = 1) #additional calculations

#VAST model settings
settings <- make_settings(n_x = nrow(Network_sz),
                          Region = "Stream_network",
                          purpose = "index2",
                          fine_scale = FALSE,
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          OverdispersionConfig = OverdispersionConfig,
                          ObsModel = ObsModel,
                          bias.correct = T,
                          Options = Options,
                          Version = "VAST_v12_0_0")
settings$Method <- "Stream_network"
settings$grid_size_km <- 1
##

#Covariates to use. Using all covariates
covar_touse <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm","YearsSinceDam")

n_p <- length(covar_touse) #number of covariates



#### Data input ######
Data_inp <- Data_len250 #the input data
Data_inp$Catch_KG <- Data_inp$Density #Use density as the outcome
######################



cov_input <- hab_dat
X1config_inp <- array(0, dim = c(2,n_p)) #all turned off
X2config_inp <- array(0, dim = c(2,n_p)) #all turned off

#Now turn the ones I want on
X1config_inp[,which(covar_touse=="Dist2Coast_FromMid")] <- 1 
X2config_inp[,which(covar_touse=="Dist2Coast_FromMid")] <- 1 

X1config_inp[,which(covar_touse=="loc_elev")] <- 1 
X2config_inp[,which(covar_touse=="loc_elev")] <- 1 

X1config_inp[,which(covar_touse=="YearsSinceDam")] <- 1 
X2config_inp[,which(covar_touse=="YearsSinceDam")] <- 1 


X1_formula_inp = paste0("~",(paste0(covar_touse, collapse = "+"))) # a symbolic description of the model to be fitted
X2_formula_inp = paste0("~",(paste0(covar_touse, collapse = "+"))) # a symbolic description of the model to be fitted




#Saving directories
path <- file.path(getwd(), paste0("Waikato/Model2_250"))
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "Figures")
dir.create(fig, showWarnings = FALSE)


#VAST files needed
ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.o"), to = path)
ignore <- file.copy(from = file.path(paste0(getwd(), "/Waikato"), "VAST_v12_0_0.dll"), to = path)


start=Sys.time() #measure how long it takes

## fit0 - Compile the model and set up the model structure specified in settings.
fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
                  working_dir = path,
                  
                  covariate_data = cov_input,
                  X1config_cp = X1config_inp,
                  X2config_cp = X2config_inp,
                  X1_formula = X1_formula_inp,
                  X2_formula = X2_formula_inp,
                  
                  run_model = FALSE, #model isn't run here
                  
                  extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], 
                                                               "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map

## fit1 - check if the model parameters are identifiable.
fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_i = Data_inp[,"Year"],
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
                  working_dir = path,
                  
                  covariate_data = cov_input,
                  X1config_cp = X1config_inp,
                  X2config_cp = X2config_inp,
                  X1_formula = X1_formula_inp,
                  X2_formula = X2_formula_inp,
                  
                  
                  extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz, 
                  model_args = list(Map = Map, Parameters = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0), #SD isn't estimated
                  test_fit = FALSE)
TMBhelper::check_estimability(fit1$tmb_list$Obj)
#check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) #Check model fit
saveRDS(fit1, file.path(path, "fit1.rds"))

## fit - Run the model, estimating standard errors (should converge if model is checked properly)
fit <- fit_model(settings = settings,
                 Lat_i = Data_inp[,"Lat"],
                 Lon_i = Data_inp[,"Lon"],
                 t_i = Data_inp[,"Year"],
                 b_i = Data_inp[,"Catch_KG"],
                 a_i = Data_inp[,"AreaSwept_km2"],
                 c_iz = as.numeric(Data_inp[,"CategoryNum"]) - 1,
                 working_dir = path,
                 
                 
                 covariate_data = cov_input,
                 X1config_cp = X1config_inp,
                 X2config_cp = X2config_inp,
                 X1_formula = X1_formula_inp,
                 X2_formula = X2_formula_inp,
                 
                 extrapolation_args = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                 spatial_args = list(Network_sz_LL = Network_sz_LL),
                 Network_sz = Network_sz,
                 model_args = list(Map = Map, Parameters = Par),
                 optimize_args = list(startpar = fit1$parameter_estimates$par),
                 test_fit = FALSE)

saveRDS(fit, file.path(path, "Fit.rds"))
end=Sys.time()
time = end - start ; time


#fit <- readRDS(file.path(path, "Fit.rds"))
fit$parameter_estimates$diagnostics



## New plots based on updated functions (18/11/20) ##
source("./Code/funcs.r")

#####
## Check model ##
dharmaRes = summary(fit, what="residuals", working_dir=paste0(fig,"/"), type=1)
#dharmaRes <- readRDS(file.path(path, "dharmaRes.rds"))
plot_residuals(residuals = dharmaRes$scaledResiduals, fit = fit, save_dir=fig,
               Data_inp = Data_inp, network=Network_sz_LL, coords="lat_long")

saveRDS(dharmaRes, file.path(path, "dharmaRes.rds"))

####
## Arnaud's residual plots ##

######## Check whether observed encounter frequencies for either low or high probability samples 
######## are within the 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = paste0(fig, '/'))


######## Run diagnostics for the positive-catch-rate component of the delta-Lognormal model
######## We can visualize fit to residuals of catch-rates given encounters using a Q-Q plot.  
######## A good Q-Q plot will have residuals along the one-to-one line
source("./Code/plotQuantileDiagnostic.r")
Q = plotQuantileDiagnostic(TmbData = fit$data_list, Report = fit$Report, FileName_PP = "Posterior_Predictive",
                           FileName_Phist = "Posterior_Predictive-Histogram", save_dir = paste0(fig, "/QQ_Fn/" ),
                           FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist")

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
     main = "Hist of DHARMa residuals",
     xlab = "Residuals (outliers are marked red)",
     cex.main = 1)
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
            main = "QQ plot residuals", 
            cex.main = 1)
dev.off()
##
####

######

dens <- quantile(log(fit$Report$D_gct))

# #Index for count
# plot_index(Index_ctl = fit$Report$Index_ctl,
#            #fit = fit, 
#            #Sdreport = fit$parameter_estimates$SD, 
#            DirName = fig, 
#            interval_width = 1.96,
#            category_names = c("<250",">=250")#, 
#            #use_biascorr=TRUE, 
#            #Plot_suffix = "Count"
#            )
# 
# #Index for biomass
# plot_index(fit = fit, 
#            Sdreport = fit$parameter_estimates$SD, 
#            DirName = fig, 
#            interval_width = 1.96,
#            category_names = c("<250",">=250"), 
#            use_biascorr=TRUE, 
#            Plot_suffix = "Biomass")

# Map of lf probability of capture across time
plot_maps_network(plot_set = 1, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD, 
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig, 
                  Panel = "category", 
                  PlotName = "POC_lf_yearly",
                  PlotTitle = "Longfin eel yearly probability of capture in Waikato, NZ",
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
                  PlotTitle = "Longfin eel P.O.C in Waikato, NZ",
                  cex = 0.75,  
                  Zlim = c(0,1), 
                  arrows=F, 
                  pch=15)



# Log-predicted density by category (<250, >=250)
plot_maps_network(plot_set = 3, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel="Category",
                  category_names = c("<250mm",">=250mm"),
                  PlotName = "lflogdens",
                  PlotTitle = "Log-predicted density",
                  Zlim = c(min(dens),max(dens)),
                  legend=TRUE, 
                  cex=0.5, 
                  pch=19)

# Log-predicted density by year
plot_maps_network(plot_set = 3, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel="Year",
                  category_names = c("<250mm",">=250mm"),
                  PlotName = "lflogdens",
                  PlotTitle = "Log-predicted density",
                  Zlim = c(min(dens),max(dens)),
                  legend=TRUE, 
                  cex=0.5, 
                  pch=19)

#Spread of longfin eel density over time
jpeg(file.path(fig, "Log_density_spread.jpg"), width = 600, height = 600)
plot(as.vector(log(fit$Report$D_gct[,2,]))~ rep(c(2009:2017), each=nrow(network)),
     xlab="Year", ylab="Log(density)", main="Spread of log(density) of longfin eels >=250mm in the Waikato")
dev.off()

jpeg(file.path(fig, "Density_spread.jpg"), width = 600, height = 600)
plot(as.vector((fit$Report$D_gct[,2,]))~ rep(c(2009:2017), each=nrow(network)),
     xlab="Year", ylab="Density", main="Spread of density of longfin eels >=250mm in the Waikato")
dev.off()


# Total biomass across all categories
plot_maps_network(plot_set = 4,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel="Category",
                  PlotName = "Totalbiomass",
                  PlotTitle = "Total biomass",
                  legend=TRUE,
                  cex=0.5,
                  pch=19)

# Total biomass across all categories
plot_maps_network(plot_set = 4,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel="Year",
                  PlotName = "Totalbiomass",
                  PlotTitle = "Total biomass",
                  legend=TRUE,
                  cex=0.5,
                  pch=19)

# Spatio-temporal variation in probability of capture
plot_maps_network(plot_set = 5,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  PlotName = "Epsilon_POC_lf",
                  PlotTitle = "Spatio-temporal variation of longfin eel P.O.C",
                  legend=TRUE,
                  cex=0.5,
                  pch=19)

# Spatio-temporal variation in positive catch rates
plot_maps_network(plot_set = 6, 
                  fit = fit, 
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list, 
                  spatial_list = fit$spatial_list, 
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  PlotName = "Epsilon_PCRs_lf",
                  PlotTitle = "Spatio-temporal variation of longfin eel P.C.Rs",
                  legend=TRUE, 
                  cex=0.5, 
                  pch=19)




# Covariates that are included in the model (measured values)
plot_maps_network(plot_set = 9,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  PlotName = "Raw_covariate_values",
                  PlotTitle = "Standardised covariate values",
                  covar_names = covar_touse,
                  cex = 0.5,
                  pch=15)




#Plot individual covariate effects for P.O.C
covar_names_to_use <- c("Distance to coast (km)", "Upstream elevation (m)", "Years since dam") #Add names here
n_p = length(covar_names_to_use)
which_np_touse = which(covar_all %in% c("Dist2Coast_FromMid","loc_elev","YearsSinceDam")) 


plot_maps_network(plot_set = 10,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  covar_names=covar_names_to_use,
                  PlotName="Covariate_effects_poc_lessthan250",
                  PlotTitle = "Individual covariate effect P.O.C",
                  cex = 0.5,
                  pch=15,
                  n_p = n_p,
                  which_np_touse = which_np_touse,
                  which_cat_cov_toplot = 1) #1st category

plot_maps_network(plot_set = 10,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  covar_names=covar_names_to_use,
                  PlotName="Covariate_effects_poc_greaterthan250",
                  PlotTitle = "Individual covariate effect P.O.C",
                  cex = 0.5,
                  pch=15,
                  n_p = n_p,
                  which_np_touse = which_np_touse,
                  which_cat_cov_toplot = 2) #1st category


plot_maps_network(plot_set = 11,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  covar_names=covar_names_to_use,
                  PlotName="Covariate_effects_pcr_lessthan250",
                  PlotTitle = "Individual covariate effect P.C.Rs",
                  cex = 0.5,
                  pch=15,
                  n_p = n_p,
                  which_np_touse = which_np_touse,
                  which_cat_cov_toplot = 1)

plot_maps_network(plot_set = 11,
                  fit = fit,
                  Sdreport = fit$parameter_estimates$SD,
                  TmbData = fit$data_list,
                  spatial_list = fit$spatial_list,
                  DirName = fig,
                  Panel = "category",
                  category_names = c("<250mm",">=250mm"),
                  covar_names=covar_names_to_use,
                  PlotName="Covariate_effects_pcr_greaterthan250",
                  PlotTitle = "Individual covariate effect P.C.Rs",
                  cex = 0.5,
                  pch=15,
                  n_p = n_p,
                  which_np_touse = which_np_touse,
                  which_cat_cov_toplot = 2)


# # Combined covariate effects for probability of capture
# plot_maps_network(plot_set = 12, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   PlotName = "Combined_covareffect_poc",
#                   PlotTitle = "Combined covariate effect of longfin eel P.O.C",
#                   legend=TRUE, 
#                   cex=0.5, 
#                   pch=19)
# 
# # Combined covariate effects for positive catch rates
# plot_maps_network(plot_set = 13, 
#                   fit = fit, 
#                   Sdreport = fit$parameter_estimates$SD,
#                   TmbData = fit$data_list, 
#                   spatial_list = fit$spatial_list, 
#                   DirName = fig,
#                   Panel = "category",
#                   category_names = c("<250mm",">=250mm"),
#                   PlotName = "Combined_covareffect_pcr",
#                   PlotTitle = "Combined covariate effect of longfin eel P.C.Rs",
#                   legend=TRUE, 
#                   cex=0.5, 
#                   pch=19)

plot_range_index(Report = fit$Report, 
                 TmbData = fit$data_list, 
                 Sdreport = fit$parameter_estimates$SD, 
                 Znames = colnames(fit$data_list$Z_gm), 
                 PlotDir = fig, 
                 Year_Set = fit$year_labels, 
                 use_biascorr = TRUE, 
                 category_names = c("<250mm",">=250mm"))


waikato_river_length = (sum(network$length))

plot_range_index_riverlength(Report = fit$Report, 
                             TmbData = fit$data_list, 
                             Sdreport = fit$parameter_estimates$SD, 
                             Znames = colnames(fit$data_list$Z_gm), 
                             PlotDir = fig, 
                             year_labels = as.character(fit$year_labels), 
                             #years_to_plot = fit$years_to_plot,
                             FileName_EffArea = paste0(fig,"/River_length_occupied.png"),
                             use_biascorr = TRUE, 
                             category_names = c("<250mm",">=250mm"),
                             total_river_length = waikato_river_length)

####

###################################################
###################################################

