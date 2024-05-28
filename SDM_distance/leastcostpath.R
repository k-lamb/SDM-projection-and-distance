library(raster) #for converting df to raster and used in cost distance
library(igraph) #for cost distance
library(gdistance) #for transition and cost distance 
library(SpatialEpi) #for converting latlong to km grid
library(geosphere)
library(tidyr)
library(plyr)
library(dplyr)
library(devtools)

set.seed(1066)
setwd("~/Desktop/Documents/Research/Q3/SDM/")

# bio_curr_pred_all_old <- read.csv("./SDM_dist/data/bioclim_spatial_predictions_all_autocorr10km_curr_newtest.csv")[,-1]
# preds <- bio_curr_pred_all[,c("x", "y", "prediction")] #restrict to just long, lat, and predictions
bio_curr_pred_all <- read.csv("SDM_dist/data/100rf_seed_avg.csv")[,-1] # averages from 100 rf models

#binning predictions by 0.1 and ensuring disconnected areas can still be routed through (thus adding 0.01)
preds.temp <- bio_curr_pred_all %>% 
  dplyr::select(x, y, prediction) %>% 
  mutate_at(vars(prediction), ~round(., 1)) %>%
  # mutate_at(vars(x), ~round(., 3)) %>%
  # mutate_at(vars(y), ~round(., 3)) %>%
  mutate(prediction = if_else(prediction==0, 0.01, prediction))


preds.temp <- cbind(bio_curr_pred_all_old[,c("x", "y")], preds.temp[,"prediction"])

# preds.temp$prediction[preds.temp$prediction == 0] <- 0.01 #makes sure dead zones can still be traversed

#make raster and assign coordinate system of data
r <- rasterFromXYZ(preds.temp) #set up raster from predictions
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
raster::plot(r) #just checking it out

#create transition layer
tr1x <- transition(r, mean, directions=16, symm=F)
tr1Cx <- geoCorrection(tr1x, type="c") #geo-correct the matrix

#make points set
# populations <- read.csv("data/populations.csv")
populations <- read.csv("~/Downloads/West_seq.txt", sep="\t")
populations <- populations %>% drop_na(longitude)
# populations <- populations[-c(263),]
IDs <- populations$ID

#create empty data frame and check the script is working appropriately

loc1 <- cbind(populations[1,][1], populations[1,][2], populations[1,][4]) #long, lat, pop
loc2 <- cbind(populations[2,][1], populations[2,][2], populations[2,][4]) #long, lat, pop

#calculate 1 run of distances
dist <- shortestPath(tr1Cx, as.matrix(loc1[c(1,2)]), as.matrix(loc2[c(1,2)]), output="SpatialLines")
dist <- SpatialLinesLengths(dist)

temp <- cbind(dist, distm(loc1[c(1,2)], loc2[c(1,2)], fun = distVincentyEllipsoid)/1000, loc1[3], loc2[3]) #calculates distance (km) between 2 pops using geo-corrected cost map tr1C
colnames(temp) <- c("env_dist", "dir_dist", "pop1", "pop2")
cost.dist <- temp[-c(1), ] #make new empty df without first row of data

#for loops that calculates distance (km) between pops using geo-corrected cost map tr1C

for (i in 1:length(IDs)) {
  loc1 <- cbind(populations[i,][1], populations[i,][2], populations[i,][4]) #lat, long, pop for loc1
  
  for (x in 1:length(IDs)) {
    loc2 <- cbind(populations[x,][1], populations[x,][2], populations[x,][4]) #lat, long, pop for loc2
    dist <- shortestPath(tr1Cx,as.matrix(loc1[c(1,2)]),as.matrix(loc2[c(1,2)]), 
                         output="SpatialLines")
    dist <- SpatialLinesLengths(dist)
    temp <- cbind(dist, distm(loc1[c(1,2)], loc2[c(1,2)], fun = distVincentyEllipsoid)/1000,
                  loc1[3], loc2[3]) #calculate distance between loc1 and loc2 using geo-corrected cost map tr1C
    colnames(temp) <- c("env_dist", "dir_dist", "pop1", "pop2")
    cost.dist <- rbind(cost.dist, temp) #append to temp df 
  }
}

write.csv(cost.dist, "./SDM_dist/data/SEQ_eco_dist_path.csv")
# cost.dist <- read.csv("./SDM_dist/data/eco_dist_path.csv")

ggplot(cost.dist, aes(x=dir_dist, y=env_dist))+
  geom_point()+
  geom_smooth(method="lm")+
  ylab("migration distance (km)")+
  xlab("direct distance (km)")+
  geom_abline(color="red")+
  xlim(c(0,1500))+
  ylim(c(0,1500))+
  theme_bw()

#plot shortest path

#examples 
A <- c(-80.61300, 37.27800) #VA5
B <- c(-93.192, 44.901) #MN117

C <- c(-77.8865, 34.1777) #Wilmington
D <- c(-82.5496, 28.1907) #Tampa-- completely separate

AtoB <- shortestPath(tr1Cx, AL19, AL23, output="SpatialLines")
# AtoB <- shortestPath(tr1Cx, A, B, output="SpatialLines")
# CtoD <- shortestPath(tr1Cx, C, D, output="SpatialLines")

jpeg("./SDM_dist/plots/eco_path_plots/example_cost_plot_multi.jpg", width=10, height=10, units="in", res=300)
plot(r)
plot(AtoB, add=T)
# plot(CtoD, add=T)
dev.off()

#plot all distances of interest
res_pops <- read.csv("./SDM_dist/data/res_pops.csv") #shortened list of pops of interest
res_IDs <- res_pops$alt_ID #used for for loop length dim

#for loops that plot shortest distance and saves to file
for (i in 1:length(res_IDs)) {
  
  loc1 <- c(as.numeric(res_pops[i,][1]), as.numeric(res_pops[i,][2]))
  loc1 <- c(as.numeric(loc1[1]), as.numeric(loc1[2])) #sets up coords as vector
  loc1_name <- as.character(res_pops[i,][5]) #gets loc1 name
  
  for (x in 1:length(res_IDs)) {
    
    loc2 <- c(as.numeric(res_pops[x,][1]), as.numeric(res_pops[x,][2]))
    loc2 <- c(as.numeric(loc2[1]), as.numeric(loc2[2])) #sets up coords as vector
    loc2_name <- as.character(res_pops[x,][5]) #gets loc2 name
    
    temp <- shortestPath(tr1Cx, loc1, loc2, output="SpatialLines") #calculates shortest path as graph-able object
    
    jpeg(sprintf("./SDM_dist/plots/eco_path_plots/populations of interest/%s_to_%s_plot.jpg", loc1_name, loc2_name), width=10, height=10, units="in", res=300)
    plot(r) #plots raster
    plot(temp, add=T) #plots line onto raster
    dev.off()
    
  }
}
