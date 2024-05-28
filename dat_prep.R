
# GIS-related packages
library(SpatialEpi) # for km lat/long grid
library(sp) # format conversion, general GIS
library(sf) # format conversion, general GIS
library(geodata) # Worldclim 2.1 bioclim download
library(geosphere) # outlier distance calculation
library(concaveman) # concave polygon for pseudo-absence generation now that rangemap is defunct
library(raster) # masking, general GIS
library(rgeos) # for clipping/cropping
library(maptools) # world map for simplest clipping/cropping
library(spatialEco) # pseudo-absence generation
library(terra) # bioclim rastering

# Data manipulation packages
library(dplyr) # filter
library(plyr) # round_any

# plotting
library(maps)
library(rnaturalearth)
library(ggplot2)

setwd("~/Desktop/Documents/Research/Q3/SDM/")
set.seed(1066) #set seed for ML random forest 

# Declare variables
grid_size = 10 # km resolution for grid cells
# threshold = 200 * 1e3 # first number is kilometer distance to nearest neighbor for a point to be declared an outlier and removed. 1e3 converts distance to meters
cluster.size = 5 # number of neighbors within given distance
buffer.dist = 200 # buffer distance for concave hull
resolution = 2.5 # worldclim resolution. finer resolution takes longer (0.5, 2.5, etc.)
nearest.neighbor.TF = T # determines whether to use cluster size to subsample presence
pa.rep = 1 # multiplier for creating pseudo-absence data
sanity.plot = T # whether to plot sanity checks

# function for cluster size
cluster.func <- function(x) {
  sorted_x <- sort(x)
  return(sorted_x[cluster.size])
}


### DATA PREPARATION

## Actions performed in this code:
# 1 - Import iNaturalist data
# 2 - Spatial down-sampling of presence data
# 3 - Generate pseudo-absence data bound by some reasonable proximity to known points
# 4 - Download/import WorldClim 2.1 data
# 5 - Save data


## 1. Import set of coordinates, with columns "id", "latitude" and "longitude"
coordinates.full <- read.csv("./data/iNaturalist_coords.csv", sep=',') # data from iNaturalist to May 2022


## 2. Spatial down-sample of presence data
# i - create spatial grids
# ii - remove outlier points (bad observations?)
# iii - remove duplicate points (down-sample to 1 point per grid cell)
# iv - sanity check plot

# i. create spatial grids (10km*10km cell size)
grid <- round(latlong2grid(coordinates.full[,c("longitude","latitude")])) # gets km grid and rounds to nearest whole km
grid$x <- round_any(grid$x, grid_size)
grid$y <- round_any(grid$y, grid_size)
coordinates.full <- cbind(coordinates.full, grid)

# ii. remove outlier points ≥200km away from the nearest neighboring point (<1min to run)
dist.mat <- distm(coordinates.full %>% dplyr::select("longitude", "latitude"), fun=distHaversine)
diag(dist.mat) <- NA # remove distances to self (dist=0)
nearest.neighbor <- apply(dist.mat, 1, cluster.func) # grabs nth nearest neighbor from matrix
coordinates.full$distance <- nearest.neighbor

# remove points further than threshold
# threshold is decided as the 5th percentile distance of the nth neighbor
if (nearest.neighbor.TF == T) {
  coordinates.full <- coordinates.full %>% 
    filter(distance <= (quantile(nearest.neighbor, 0.05) %>% as.numeric)*1e3) # distance is multiplied by 1e3 to convert to meters
}

# iii. remove duplicate observations in a single 10km*10km cell
coordinates <- coordinates.full %>%
  group_by(x, y) %>% 
  slice(1)

# iv. sanity check plotting
# set up map attributes
world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

if (sanity.plot == T) {
  ggplot(data=world.sf)+
    geom_sf()+
    geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
    geom_sf(data=lakes, fill="white")+ #adds the great lakes
    coord_sf(xlim=c(min(coordinates$longitude)-3, max(coordinates$longitude)+3), 
             ylim=c(min(coordinates$latitude)-3, max(coordinates$latitude)+3), expand=FALSE)+
    geom_point(data=coordinates.full, aes(x=longitude, y=latitude), size=0.5)+
    # geom_text(data=coordinates, aes(x=longitude, y=latitude, label=id))+
    xlab("Latitude")+
    ylab("Longitude")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}


## 3. Generate pseudo-absence data
# i - generate a minimum concave polygon and buffer to pre-determined extent
# ii - mask the species range to avoid generating points outside buffer, in ocean, or in lakes
# iii - generate pseudo-absence data within mask and away from ANY known presence observation
# iv - merge pseudo-absence and presence data

# i. generate a minimum concave polygon around known observed points with a spatial buffer
sf.df <- st_as_sf(coordinates, coords=c("longitude", "latitude"), crs=4326)
concave.hull <- concaveman(sf.df) # takes ~3 minutes to run. from package 'concaveman'
buffer.hull <- st_buffer(concave.hull, dist = buffer.dist*1e3) # takes ~2 minute to run

# ii. mask hulls. takes ~1 minute (plotting time not included)
# in case sanity check in 2.iv is skipped, load relevant data
world.sf <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
lakes <- st_as_sf(map("lakes", plot = FALSE, fill = TRUE))

# get data for north american countries of interest
world.sp <- ne_countries(scale = "large", returnclass = "sp")
world.crop <- crop(world.sp, extent(min(coordinates.full$longitude)-10, max(coordinates.full$longitude)+10,
                             min(coordinates.full$latitude)-10, max(coordinates.full$latitude)+10))

r <- raster(ncol=2e3, nrow=2e3) # create raster object to fill into. bigger gives better resolution but takes longer

# mask out oceans
r.mask <- rasterize(world.crop, r, mask=F) 

# mask out lakes
lakes.sp <- as(lakes, "Spatial")
r.mask <- raster::mask(r.mask, lakes.sp, inverse=T) 

# mask outside of buffer.hull
buffer.sp <- as(buffer.hull, "Spatial")
r.mask <- raster::mask(r.mask, buffer.sp)

## iii. generate pseudo-absence data
coord.sp <- SpatialPoints(coordinates.full[,c("longitude", "latitude")])

pa.df <- data.frame(x=numeric(), y=numeric())
for (i in 1:pa.rep) {
  pa <- pseudo.absence(coord.sp, n=length(coordinates$latitude), window='hull', Mask=r.mask, KDE=T) # prevents point overlap with ANY known C americana point (regardless of spatial down-sampling)
  temp <- as.data.frame(pa[["sample"]]@coords) #convert to df to plot
  pa.df <- rbind(pa.df, temp)
}

# restrict pseudo-absence data to 1 non-observation per km grid cell
pa.grid <- round(latlong2grid(pa.df[,c(1,2)])) #gets km grid and rounds to nearest whole km
pa.grid <- cbind(pa.df, pa.grid)
colnames(pa.grid) <- c("longitude", "latitude", "x", "y")

# make sure points aren't in the same grid cell
pa.grid$x <- round_any(pa.grid$x, 10)
pa.grid$y <- round_any(pa.grid$y, 10)
pa.grid <- pa.grid %>% 
  group_by(x, y) %>% 
  slice(1)

# sanity check plot
raster.pts <- rasterToPoints(r.mask)
r.sf <- st_as_sf(data.frame(raster.pts), coords = c("x", "y"), crs = st_crs(r.mask))

if (sanity.plot == T) {
  ggplot(data=world.sf)+
    geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
    geom_sf(data=lakes, fill="white")+ #adds the great lakes
    geom_sf(data=r.sf, fill=NA, color="gray", alpha=0.3)+ # full mask
    geom_sf(data=concave.hull, fill=NA, color="red")+ # concave hull
    geom_sf(data=buffer.hull, fill=NA, color="blue3")+ # buffer hull (used in masking)
    geom_point(data=pa.grid, aes(x=longitude,y=latitude), color="coral1")+
    geom_point(data=coordinates, aes(x=longitude, y=latitude), color="forestgreen")+
    coord_sf(xlim=c(min(coordinates$longitude)-3, max(coordinates$longitude)+3),
             ylim=c(min(coordinates$latitude)-3, max(coordinates$latitude)+3), expand=FALSE)+
    xlab("Latitude")+
    ylab("Longitude")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}

# iv. merge presence and pseudo-absence data
presence <- coordinates %>% 
  ungroup %>% 
  dplyr::select(latitude, longitude) %>% 
  mutate(presence = 1)

pseudo.absence <- pa.grid %>% 
  ungroup %>% 
  dplyr::select(latitude, longitude) %>% 
  mutate(presence = 0)

c.americana <- rbind(presence, pseudo.absence) %>% as.data.frame()


## 4. download and prep bioclim data from WorldClim 2.1
# i - download and extract data
# ii - combine data 

# convert points for extracting bioclim data
cam.sf <- st_as_sf(c.americana, coords=c("longitude", "latitude"), crs=4326)

# download spatial data (1km^2 resolution)
bio.world <- worldclim_global(var="bio", res=resolution, path="~/Desktop/Documents/Research/Q3/SDM/geodata/", version="2.1")
cam.bio <- terra::extract(bio.world, cam.sf) 

# rename variable
cam.bio <- cam.bio[,-1]
names(cam.bio) <- paste0("bio", 1:19)

# combine data
c.americana <- cbind(c.americana, cam.bio)

if (resolution == 0.5) {
  write.csv(c.americana, "code/final/data/cam_pa30s.csv", row.names = F)
}

if (resolution == 2.5) {
  write.csv(c.americana, "code/final/data/cam_pa2m30s.csv", row.names = F)
}

# save plot
jpeg("code/final/plots/data_prep2m30s.jpeg", res=5e2, units="in", height=5, width=6)
ggplot(data=world.sf)+
  geom_sf(data=states, fill=NA, alpha=0.5, color="black")+ #adds state boundaries
  geom_sf(data=lakes, fill="white")+ #adds the great lakes
  geom_sf(data=r.sf, fill=NA, color="gray", alpha=0.3)+ # full mask
  geom_sf(data=concave.hull, fill=NA, color="red4")+ # concave hull
  geom_sf(data=buffer.hull, fill=NA, color="blue3")+ # buffer hull (used in masking)
  geom_point(data=pa.grid, aes(x=longitude,y=latitude), color="coral1", size=0.25)+
  geom_point(data=coordinates, aes(x=longitude, y=latitude), color="forestgreen", size=0.25)+
  coord_sf(xlim=c(min(coordinates$longitude)-3, max(coordinates$longitude)+3), 
           ylim=c(min(coordinates$latitude)-3, max(coordinates$latitude)+3), expand=FALSE)+
  xlab("Latitude")+
  ylab("Longitude")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
