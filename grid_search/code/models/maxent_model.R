
# libraries

# modeling
library(dismo)
library(dplyr) # pipes
library(SpatialEpi) # for latlong to km grid conversion

# map projection
library(maps)
library(ggplot2)
library(geodata) # worldclim spatrastering
library(sp) # cropping/clipping
library(sf) # cropping/clipping
library(raster) # extent clipping
library(terra) # dealing with SpatRasters
library(rnaturalearth) # setting up coordinate grid for prediction

setwd("~/Desktop/Documents/Research/Q3/SDM/")

# define variables
pca.perc = 0.95 # percent cut off for number of principal components to retain for modeling
n.model = 10 # number of random forest models to average
window.fac = 10 # number of windows on X and Y axis (multiply by itself to get total number of discrete non-overlapping windows)
pts.rep = 5 # number of times a point should be represented in a window in an X or Y plane (how far window should slide each model)

### RUNNING THE RANDOM FOREST MODELS

## Actions performed in this code:
# 1 - read in data and transform with PCA
# 2 - compute random forest model using contemporary climate data (1970-2000)
# 3 - predict N American presence/absence and save prediction raster


## 1. read in data and transform with PCA (data generated in dat_prep.R)
c.americana <- read.csv("code/grid_search/data/cam_pa2m30s_5x.csv") %>% na.omit() # 2 rows of NA's (coordinates are in ocean by NYC which escaped mask filter (resolution issues?)
cam.pca <- prcomp(c.americana[,c(paste0("bio", 1:19))], 
                  center=T, scale.=T)

# extract the number of PC which explain at least N% of the data
cumsum <- summary(cam.pca)$importance["Cumulative Proportion",]
pc.num <- which(cumsum >= pca.perc)[1] %>% as.numeric()

# extract and combine data
cam.pcN <- cam.pca[["x"]]
c.americana <- cbind(c.americana[,c("presence", "latitude", "longitude")], cam.pcN)

## 2. compute random forest model using data obtained in dat_prep.R
# i - convert lat/long of points to grid, then subset points by sliding window grid dimension. save models in list
# ii - for each model in list, generate a prediction raster. convert lat/long raster to km grid and find highest predictive value for cell
# v - plot predictions and save contemporary prediction as raster

# convert to km grid
grid <- round(latlong2grid(c.americana[,c("longitude","latitude")])) # gets km grid and rounds to nearest whole km
c.americana <- cbind(c.americana, grid)

# establish range for dimensionality of windows
lat.range <- (max(c.americana$y) - min(c.americana$y)) / window.fac
long.range <- (max(c.americana$x) - min(c.americana$x)) / window.fac

# how far window should slide each time (km)
lat.slide.interval <- lat.range / pts.rep
long.slide.interval <- long.range / pts.rep

# number of sliding windows on each axis
lat.win <- (lat.range*11)/lat.slide.interval
long.win <- (long.range*11)/long.slide.interval

# set up data for prediction
bio.world <- worldclim_global(var="bio", res=2.5, path="~/Desktop/Documents/Research/Q3/SDM/geodata/", version="2.1")
world.sp <- ne_countries(scale = "large", returnclass = "sp")
world.crop <- crop(world.sp, extent(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-3, 
                                    max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+3, 
                                    min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-3, 
                                    max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+3))

bio.crop <- terra::extract(bio.world, vect(world.crop), cells=T, xy=T)
bio.crop <- bio.crop[,-1] # gets rid of useless ID column
names(bio.crop) <- c(paste0("bio", 1:19), "cell", "longitude", "latitude") # make names match cam.pca

bio.bio <- bio.crop[,paste0("bio", 1:19)]
bio.xy <- bio.crop[,c("longitude", "latitude")]
bio.pca <- predict(cam.pca, bio.bio) # predict values of world data in existing PCA so as to not change PCA values. takes a while (~1-2 min.)
bio.pcN <- bio.pca[,1:pc.num]  # limit to same number of PC as before

# determine which cells have missing data and remove
bio.NAomit <- cbind(bio.xy, bio.pcN)
bio.NAomit <- bio.NAomit %>% na.omit()
bio.grid <- round(latlong2grid(bio.NAomit[,c("longitude","latitude")]))
bio.NAomit <- cbind(bio.NAomit, bio.grid)

# run models in sliding framework and save models to list
# predict data from models for specific windows
pred.df <- data.frame(pred = numeric(), longitude=numeric(), latitude=numeric(), x=numeric(), y=numeric())

for (i in 1:long.win) {
  # isolate longitude sliding window
  long.min = min(c.americana$x) + (long.slide.interval*(i-1))
  long.max = long.min + long.range
  sub.long <- c.americana %>% filter(x > long.min & x <= long.max) # presence/absence points with climate/location data
  bio.sub1 <- bio.NAomit %>% filter(x > long.min & x <= long.max) # broader climate data for prediction
  
  print("long window = ")
  print(i)
  
  for (n in 1:lat.win) {
    # isolate latitude sliding window
    lat.min = min(sub.long$y) + (lat.slide.interval*(n-1))
    lat.max = lat.min + lat.range
    sub <- sub.long %>% filter(y > lat.min & y <= lat.max)
    
    # predict window for model
    bio.sub2 <- bio.sub1 %>% filter(y > lat.min & y <= lat.max) # broader climate data for prediction
    
    # random forest model for window
    tryCatch(
     {
        if (length(unique(sub[,"presence"] )) == 2) { # want to avoid any windows where there are only absences or presences
          m_me <- maxent(sub[,paste0("PC", 1:pc.num)], as.factor(sub[,"presence"]), silent=T)
          pred <- predict(m_me, bio.sub2[,paste0("PC", 1:pc.num)])
          temp <- cbind(pred, bio.sub2[,c("longitude", "latitude", "x", "y")]) %>% as.data.frame()
          names(temp)[1] <- "pred"
          pred.df <- rbind(pred.df, temp) 
        }
     }, 
    error = function(e) {
      message("An error occurred: ", e)
    }
    )
  }
}

pred.grid <- pred.df %>%
  # mutate(x=round_any(x, 10)) %>%
  # mutate(y=round_any(y, 10)) %>%
  group_by(x,y) %>% # group by lat/long grid
  arrange(desc(pred)) %>%  # keep highest prediction value for each lat/long grid cell
  slice(1)

world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

pred.grid <- pred.grid %>% mutate(bins = cut(pred, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf)))

jpeg("code/grid_search/plots/maxent/maxent_cellmax_baseline.jpeg", width=10, height=10, units="in", res=500)
ggplot()+
  geom_sf(data=world.sf, fill=NA)+
  geom_sf(data=states, fill=NA)+ #adds state boundaries
  geom_tile(data=pred.grid, aes(x=longitude, y=latitude, color=bins), size=0.35)+ # for crop. 1.35=10k, 0.35=1km
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  scale_color_manual(values = c("steelblue4",  "yellow3", "yellow", "chartreuse3", "green4"))+
  geom_sf(data=lakes, fill="gray10")+ #adds the great lakes
  coord_sf(xlim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-3, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+3), 
           ylim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-3, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+3), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  geom_point(data=subset(c.americana, presence==1), aes(x=longitude, y=latitude), colour="black", size=0.25)+ #adds lab collections
  theme_bw()
dev.off()

# save files for export
# write.csv(pred.df, "code/grid_search/data/models/maxent/me_all_windows.csv", row.names = F)
write.csv(pred.grid, "code/grid_search/data/models/maxent/me_cellmax_windows.csv", row.names = F)

# save raster
pred.raster <- rasterFromXYZ(pred.grid[,c("longitude", "latitude", "pred")],
                             crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(pred.raster) # check that it worked
writeRaster(pred.raster, 'code/grid_search/data/models/maxent/me_cellmax.tif', overwrite=T)



# smoothing?
library(smoothr)
library(stars)
pred.sf <- st_as_stars(pred.raster, coords=c("longitude", "latitude"), crs=4326)
pred.sf <- st_as_sf(pred.sf, crs=4326)
pred.smooth <- smooth(pred.sf, method = "spline")

smooth.raster <- st_rasterize(pred.smooth %>% dplyr::select(pred, geometry)) # layer if starting from saved raster, pred if starting from saved df
write_stars(smooth.raster, 'code/grid_search/data/models/maxent/me_cellmax_smooth.tif')

smooth.raster <- raster("code/grid_search/data/models/maxent/me_cellmax_smooth.tif")
smooth.df <- rasterToPoints(smooth.raster) %>% as.data.frame()
names(smooth.df) <- c("longitude", "latitude", "pred")

smooth.df <- smooth.df %>% mutate(bins = cut(pred, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf)))

jpeg("code/grid_search/plots/maxent/me_cellmax_smooth.jpeg", width=10, height=10, units="in", res=500)
ggplot()+
  geom_sf(data=world.sf, fill=NA)+
  geom_sf(data=states, fill=NA)+ #adds state boundaries
  geom_tile(data=smooth.df, aes(x=longitude, y=latitude, color=bins), size=1.25)+ # for crop. 1.25 for smoothed
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  scale_color_manual(values = c("steelblue4",  "yellow3", "yellow", "chartreuse3", "green4"))+
  geom_sf(data=lakes, fill="gray10")+ #adds the great lakes
  coord_sf(xlim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))-3, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(longitude))+3), 
           ylim=c(min(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))-3, 
                  max(c.americana %>% filter(presence == 1) %>% dplyr::select(latitude))+3), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  geom_point(data=subset(c.americana, presence==1), aes(x=longitude, y=latitude), colour="black", size=0.25)+ #adds lab collections
  theme_bw()
dev.off()
