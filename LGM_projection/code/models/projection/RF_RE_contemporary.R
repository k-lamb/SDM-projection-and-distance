
# libraries

# modeling
library(randomForest) # random forest modeling 
library(splitTools) # splitting data into testing and training
library(mecofun) # predictions from random forest
library(dplyr) # pipes

# map projection
library(geodata) # worldclim spatrastering
library(rbioclim) #from github: "MoisesExpositoAlonso/rbioclim" ... needed for LGM download
library(sp) # cropping/clipping
library(sf) # cropping/clipping
library(raster) # extent clipping
library(terra) # dealing with SpatRasters
library(rnaturalearth) # setting up coordinate grid for prediction

setwd("~/Desktop/Documents/Research/Paper Code/LGM_AP-KL/")

# define variables
pca.perc = 0.95 # percent cut off for number of principal components to retain for modeling
n.model = 100 # number of random forest models to average
size.t = 16 # text size for plots
### RUNNING THE RANDOM FOREST MODELS

## Actions performed in this code:
# 1 - read in data and transform with PCA
# 2 - compute random forest model using contemporary climate data (1970-2000)
# 3 - predict N American presence/absence and save prediction raster
# 4 - predict N American LGM using contemporary model and save prediction raster


## 1. read in data and transform with PCA (data generated in dat_prep.R)
c.americana <- read.csv("data/RearEdge_pa2m30s.csv") %>% na.omit() # 2 rows of NA's (coordinates are in ocean by NYC which escaped mask filter (resolution issues?)
cam.pca <- prcomp(c.americana[,c(paste0("bio", 1:19))], 
                  center=T, scale.=T)

# extract the number of PC which explain at least N% of the data
cumsum <- summary(cam.pca)$importance["Cumulative Proportion",]
pc.num <- which(cumsum >= pca.perc)[1] %>% as.numeric()

# extract and combine data
cam.pcN <- cam.pca[["x"]]
c.americana <- cbind(c.americana[,c("presence", "latitude", "longitude")], cam.pcN)

## 2. compute random forest model using data obtained in dat_prep.R
# i - split data in testing and training to evaluate model accuracy/type I & II errors
# ii - run model on total aggregate data set (100 model)
# iii - set up worldclim data for broader prediction
# iv - predict full raster from each model and average across cells
# v - plot predictions and save contemporary prediction as raster

# i. split data into training (70%) and testing (30%) data sets
inds <- partition(c.americana$presence, p=c(train = 0.7, test = 0.3))
c.am.training <- c.americana[inds$train, ]
c.am.testing <- c.americana[inds$test, ]

# baseline model with 70% training set
set.seed(216) # set random set so model is replicable 

# run model on the number of PCA identified earlier as explaining 95% of climate variance (pc.num)
m_rf <- randomForest(x=c.am.training[,paste0("PC", 1:pc.num)],
                     y=as.factor(c.am.training[,"presence"]), 
                     ntree=5000, importance =T)

# test error rates
prediction.table <- predict(m_rf, c.am.testing[,paste0("PC", 1:pc.num)]) #create a prediction table
con_matrix <- table(observed=c.am.testing[,"presence"],predicted=prediction.table) %>% as.data.frame()
false.neg <- con_matrix %>% filter(observed == 1 & predicted == 0) %>% dplyr::select(Freq) %>% as.numeric() / con_matrix %>% filter(predicted == 0) %>% dplyr::select(Freq) %>% sum()
false.pos <- con_matrix %>% filter(observed == 0 & predicted == 1) %>% dplyr::select(Freq) %>% as.numeric / con_matrix %>% filter(predicted == 1) %>% dplyr::select(Freq) %>% sum()

false.neg
false.pos

# ii. run full model with 100 seeds to get average prediction per cell. takes ~6 minutes to run
m_rf <- randomForest(x=c.americana[,paste0("PC", 1:pc.num)], y=as.factor(c.americana[,"presence"]), 
                     ntree=500, importance=T)
temp <- list(m_rf)

# run 10 seeds and save models... 100 took forever to predict because of higher resolution
for (i in 2:n.model){
  set.seed(i)
  t <- randomForest(x=c.americana[,paste0("PC", 1:pc.num)], y=as.factor(c.americana[,"presence"]), 
                    ntree=500, importance=T)
  x <- list(t)
  temp <- c(temp, x)
  print(i)
}

# iii. set up worldclim data for broader prediction
saveme <- temp # allow easy restarting without having to re-run models

# download world data to crop extent of full bioclim data
world.sp <- ne_countries(scale = "large", returnclass = "sp")
world.crop <- crop(world.sp, extent(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-8, 
                                    max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+8, 
                                    min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-8, 
                                    max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+8))
crop.sf <- st_as_sf(world.crop)
coord.sf <- st_coordinates(crop.sf) %>% as.data.frame() %>% dplyr::select(X,Y)
names(coord.sf) <- c("lon", "lat")

# retrieve bioclim data for world.crop... using resolution of 2.5

# use an exist check to cut down on computation time if re-running
if (exists("bio.crop") == FALSE) {
  bio.world <- worldclim_global(var="bio", res=2.5, path="~/Desktop/Documents/Research/Q3/SDM/geodata/", version="2.1")
  # bio.raster <- stack(bio.world) # convert to workable format
  # bio.extract <- extract(bio.raster, world.crop)
  
  bio.crop <- terra::extract(bio.world, vect(world.crop), cells=T, xy=T)
  bio.crop <- bio.crop[,-1] # gets rid of useless ID column
  names(bio.crop) <- c(paste0("bio", 1:19), "cell", "longitude", "latitude") # make names match cam.pca
}

bio.bio <- bio.crop[,paste0("bio", 1:19)]
bio.xy <- bio.crop[,c("longitude", "latitude")]
bio.pca <- predict(cam.pca, bio.bio) # predict values of world data in existing PCA so as to not change PCA values. takes a while (~1-2 min.)
bio.pcN <- bio.pca[,1:pc.num]  # limit to same number of PC as before

# determine which cells have missing data and remove
bio.NAomit <- cbind(bio.xy, bio.pcN)
bio.NAomit <- bio.NAomit %>% na.omit()
coord.sf <- bio.NAomit %>% dplyr::select(longitude, latitude) %>% as.data.frame()
bio.pcN <- bio.NAomit %>% dplyr::select(paste0("PC", 1:pc.num)) %>% as.data.frame()

# iv. generate predictions from all models for range extent +/- 3 lat/long
pred.list <- lapply(saveme, mecofun::predictSDM, newdata=bio.pcN) # ~2-3 minutes run time
pred.df <- as.data.frame(pred.list)
pred.cell <- pred.df %>% rowMeans()

# combine predictions with spatial data
pred.grid <- cbind(pred.cell, coord.sf)

# v. plot check
pred.grid <- pred.grid %>% mutate(bins = cut(pred.cell, breaks = c(-Inf,0.2,0.4,0.6,0.8,Inf))) %>% as.data.frame()

world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

jpeg("plots/rf/m_rf_100avg_RearEdge_contemporary.jpeg", width=10, height=10, units="in", res=500)
ggplot()+
  geom_sf(data=world.sf, fill=NA)+
  geom_sf(data=states, fill=NA)+ #adds state boundaries
  geom_tile(data=pred.grid, aes(x=longitude, y=latitude, color=bins), size=1)+ #for crop
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  scale_color_manual(values = c("steelblue4",  "yellow3", "yellow", "chartreuse3", "green4"))+
  geom_sf(data=lakes, fill="gray10")+ #adds the great lakes
  coord_sf(xlim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-8, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+8), 
           ylim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-8, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+8), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  geom_point(data=subset(c.americana, presence==1), aes(x=longitude, y=latitude), colour="black", size=0.25)+ #adds lab collections
  theme_bw()+
  theme(text = element_text(size = size.t), legend.position = "none", # right
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = (size.t - 6)))
dev.off()

# save files for export
write.csv(pred.grid, "data/models/rf/rf_100avg_contemporary.csv", row.names = F)

# save raster
pred.raster <- rasterFromXYZ(pred.grid[,c("longitude", "latitude", "pred.cell")],
                             crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(pred.raster) # check that it worked
writeRaster(pred.raster, 'data/models/rf/rf_100avg_contemporary.tif', overwrite=T)


# 4 - predict LGM using contemporary model framework
# i - import and extract values from lgm
# ii - predict distribution from all RF models 
# iii - plot models

# i. import and extract climate data from LGM 
lgm <- recursive.getData(times="lgm", var="bio", res="2.5") # by default resolution 2.5m and the 19 bioclim variables. (but can be changed)
lgm <- lgm[["lgm"]]

# crop data to area of interest
bio.lgm.crop <- terra::extract(rast(lgm), vect(world.crop), cells=T, xy=T)
bio.lgm.crop <- bio.lgm.crop[,-1] # gets rid of useless ID column
names(bio.lgm.crop) <- c(paste0("bio", 1:19), "cell", "longitude", "latitude") # make names match cam.pca

bio.lgm.bio <- bio.lgm.crop[,paste0("bio", 1:19)]
bio.lgm.xy <- bio.lgm.crop[,c("longitude", "latitude")]
bio.lgm.pca <- predict(cam.pca, bio.lgm.bio) # predict values of world data in existing PCA so as to not change PCA values. takes a while (~1-2 min.)
bio.lgm.pcN <- bio.lgm.pca[,1:pc.num]  # limit to same number of PC as before

# determine which cells have missing data and remove
bio.lgm.NAomit <- cbind(bio.lgm.xy, bio.lgm.pcN)
bio.lgm.NAomit <- bio.lgm.NAomit %>% na.omit()
coord.lgm.sf <- bio.lgm.NAomit %>% dplyr::select(longitude, latitude) %>% as.data.frame()
bio.lgm.pcN <- bio.lgm.NAomit %>% dplyr::select(paste0("PC", 1:pc.num)) %>% as.data.frame()

# iv. generate predictions from all models for range extent +/- 3 lat/long
pred.lgm.list <- lapply(saveme, mecofun::predictSDM, newdata=bio.lgm.pcN) # ~2-3 minutes run time
pred.lgm.df <- as.data.frame(pred.lgm.list)
pred.lgm.cell <- pred.lgm.df %>% rowMeans()

# combine predictions with spatial data
pred.lgm.grid <- cbind(pred.lgm.cell, coord.lgm.sf)

# v. plot check
pred.lgm.grid <- pred.lgm.grid %>% mutate(bins = cut(pred.lgm.cell, breaks = c(-Inf,0.2,0.4,0.6,0.8,Inf))) %>% as.data.frame()

world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

jpeg("plots/rf/m_rf_100avg_RearEdge_LGM.jpeg", width=10, height=10, units="in", res=500)
ggplot()+
  geom_sf(data=world.sf, fill=NA)+
  geom_sf(data=states, fill=NA)+ #adds state boundaries
  geom_tile(data=pred.lgm.grid, aes(x=longitude, y=latitude, color=bins), size=1)+ #for crop
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  scale_color_manual(values = c("steelblue4",  "yellow3", "yellow", "chartreuse3", "green4"))+
  geom_sf(data=lakes, fill="gray10")+ #adds the great lakes
  coord_sf(xlim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-8, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+8), 
           ylim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-8, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+8), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  geom_point(data=subset(c.americana, presence==1), aes(x=longitude, y=latitude), colour="black", size=0.25)+ #adds lab collections
  theme_bw()+
  theme(text = element_text(size = size.t), legend.position = "none", # right
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = (size.t - 6)))
dev.off()

# save files for export
write.csv(pred.lgm.grid, "data/models/rf/rf_100avg_LGM.csv", row.names = F)

# save raster
pred.lgm.raster <- rasterFromXYZ(pred.lgm.grid[,c("longitude", "latitude", "pred.lgm.cell")],
                             crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(pred.lgm.raster) # check that it worked
writeRaster(pred.lgm.raster, 'data/models/rf/rf_100avg_LGM.tif', overwrite=T)

