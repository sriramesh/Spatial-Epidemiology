# load libraries ----
library(raster)
library(ModelMetrics)
library(spaMM)
library(tidyverse)
library(rgdal)
library(spatstat)
library(plotrix)
library(fields)
library(leaflet)
library(maptools)
library(RColorBrewer)
library(lattice)
library(geoR)
library(geosphere)
library(ape)
library(spdep)
library(pgirmess)
library(sp)
library(gtools)
library(car)

# global options ----
options(stringsAsFactors = F) # import quality
par(mar=c(1,1,1,1)) # ensure plot margin sizes are always large enough

# load soil-transmitted helminths data ----
STH_kenya <- "Kenya_STH.csv" %>% read.csv(header = T) %>%
  select(full_name_paper, lat, long, tri_prev, tri_np)

STH_kenya <- STH_kenya %>%
  mutate(tri_prev = tri_prev/100,
         tri_examined = tri_np/tri_prev,
         tri_nneg = tri_examined - tri_np,
         tri_prev_logit = car::logit(tri_prev) + 0.0001) %>%
  drop_na(tri_prev_logit)

STH_kenya <- SpatialPointsDataFrame(coords = STH_kenya %>% select(long,lat),
                                    data = STH_kenya)

# load covariates ----
alt <- raster::getData('alt', country = 'KEN', mask = F)
bioclim <- raster::getData('worldclim', var='bio', res=2.5) # Bio12, higher resolution

roads <- readOGR("KEN_rds", "KEN_roads") #roads in Uganda
railroads <- readOGR("KEN_rrd", "KEN_rails") #railroads in Uganda
lakes <- readOGR("KEN_wat", "KEN_water_areas_dcw") #inland waterbodies in Uganda
canals <- readOGR("KEN_wat", "KEN_water_lines_dcw") #inland water lines in Uganda

# set crs---
crs <- "+proj=utm +zone=37 +south +datum=WGS84"
proj4string(STH_kenya) <- crs %>% CRS()
proj4string(roads) <- crs %>% CRS()
proj4string(railroads) <- crs %>% CRS()
proj4string(lakes) <- crs %>% CRS()
proj4string(canals) <- crs %>% CRS()

# resample rasters ----
KEN_Adm_0 <- raster::getData("GADM", country="KEN", level=0)
KEN_bioclim <- bioclim[[12]] %>% crop(KEN_Adm_0 %>% extent())
KEN_alt <- alt %>% crop(KEN_bioclim %>% extent())
# resample higher res raster using lower res raster as reference
KEN_bioclim_resamp <- KEN_bioclim %>% resample(KEN_alt, method = "ngb")
KEN_alt
KEN_bioclim_resamp

# calc distance to prevalence points ----
dist_to_lakes <- STH_kenya %>% dist2Line(lakes)
dist_to_canals <- STH_kenya %>% dist2Line(canals)
dist_to_railroads <- STH_kenya %>% dist2Line(railroads)
dist_to_roads <- STH_kenya %>% dist2Line(roads)

STH_kenya$dist_to_lakes <- dist_to_lakes[,1]
STH_kenya$dist_to_canals <- dist_to_canals[,1]
STH_kenya$dist_to_railroads <- dist_to_railroads[,1]
STH_kenya$dist_to_roads <- dist_to_roads[,1]

STH_kenya_df <- STH_kenya %>% as.data.frame() %>% select(-long.1, -lat.1)

STH_kenya_df <- STH_kenya_df %>%
  mutate(alt = raster::extract(KEN_alt, STH_kenya_df %>% select(long, lat)),
         bio12 = raster::extract(KEN_bioclim_resamp, STH_kenya_df %>% select(long, lat)))

STH_kenya_df %>% write.csv("STH_kenya_covariates.csv")
