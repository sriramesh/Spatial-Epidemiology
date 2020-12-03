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
KEN_Adm_0 <- raster::getData("GADM", country="KEN", level=0)
alt <- raster::getData('alt', country = 'KEN', mask = F)
bioclim <- raster::getData('worldclim', var='bio', res=2.5) # Bio12, higher resolution
roads <- readOGR("KEN_rds", "KEN_roads") #roads in Uganda
railroads <- readOGR("KEN_rrd", "KEN_rails") #railroads in Uganda
lakes <- readOGR("KEN_wat", "KEN_water_areas_dcw") #inland waterbodies in Uganda
canals <- readOGR("KEN_wat", "KEN_water_lines_dcw") #inland water lines in Uganda

# set CRS ----
crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
KEN_Adm_0 <- KEN_Adm_0 %>% spTransform(CRS(crs))
proj4string(STH_kenya) <- crs %>% CRS()

roads <- roads %>% spTransform(CRS(crs))
canals <- canals %>% spTransform(CRS(crs))
railroads <- railroads %>% spTransform(CRS(crs))
lakes <- lakes %>% spTransform(CRS(crs))

# rasterize roads ----
rs <- raster(extent(roads), crs=projection(crs))
rs[] <- 1:ncell(rs)

rsp <- rs %>% rasterToPolygons()
rp <- roads %>% raster::intersect(rsp)
rp$length <- rgeos::gLength(rp, byid=TRUE) / 1000
x <- tapply(rp$length, rp$layer, sum)

roads_r <- rs %>% raster()
roads_r[as.integer(names(x))] <- x

roads_r %>% plot()
roads %>% lines()
KEN_Adm_0 %>% lines()

# rasterize canals ----
rs <- raster(extent(canals), crs=projection(crs))
rs[] <- 1:ncell(rs)

rsp <- rs %>% rasterToPolygons()
rp <- roads %>% raster::intersect(rsp)
rp$length <- rgeos::gLength(rp, byid=TRUE) / 1000
x <- tapply(rp$length, rp$layer, sum)

canals_r <- rs %>% raster()
canals_r[as.integer(names(x))] <- x

canals_r %>% plot()
canals %>% lines()
KEN_Adm_0 %>% lines()



# rasterize railroads ----
rs <- raster(extent(railroads), crs=projection(crs))
rs[] <- 1:ncell(rs)

rsp <- rs %>% rasterToPolygons()
rp <- roads %>% raster::intersect(rsp)
rp$length <- rgeos::gLength(rp, byid=TRUE) / 1000
x <- tapply(rp$length, rp$layer, sum)

railroads_r <- rs %>% raster()
railroads_r[as.integer(names(x))] <- x

railroads_r %>% plot()
railroads %>% lines()
KEN_Adm_0 %>% lines()


# rasterize lakes ----
ext <-  extent(33.90959, 41.92622, -4.720417, 5.061166)
xy <- ext %>% bbox() %>% as.matrix() %>% apply(1, diff) %>% abs()
n <- 5
r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)
lakes_r <- lakes %>% rasterize(r, fun='first')

lakes_r %>% plot()
lakes %>% lines()
KEN_Adm_0 %>% lines()


# resample rasters----
KEN_alt <- alt %>% crop(KEN_Adm_0 %>% extent())
KEN_bioclim <- bioclim[[12]] %>% crop(KEN_Adm_0 %>% extent())
KEN_canals_r <- canals_r %>% crop(KEN_Adm_0 %>% extent())
KEN_roads_r <- roads_r %>% crop(KEN_Adm_0 %>% extent())
KEN_railroads_r <- railroads_r %>% crop(KEN_Adm_0 %>% extent())
KEN_lakes_r <- lakes_r %>% crop(KEN_Adm_0 %>% extent())

KEN_bioclim_resamp <- KEN_bioclim %>% resample(KEN_alt, method = "ngb")
KEN_canals_r_resamp <- KEN_canals_r %>% resample(KEN_alt, method = "ngb")
KEN_roads_r_resamp <- KEN_roads_r %>% resample(KEN_alt, method = "ngb")
KEN_railroads_r_resamp <- KEN_railroads_r %>% resample(KEN_alt, method = "ngb")
KEN_lakes_r_resamp <- KEN_lakes_r %>% resample(KEN_alt, method = "ngb")

distance_to_lakes_r <- distance(KEN_lakes_r_resamp)
distance_to_canals_r <- distance(KEN_canals_r_resamp)
distance_to_roads_r <- distance(KEN_roads_r_resamp)
distance_to_railroads_r <- distance(KEN_railroads_r_resamp)

STH_kenya_df <- STH_kenya %>% as.data.frame() %>% select(-long.1, -lat.1)

STH_kenya_df <- STH_kenya_df %>%
  mutate(alt = raster::extract(KEN_alt, STH_kenya_df %>% select(long, lat)),
         bio12 = raster::extract(KEN_bioclim_resamp, STH_kenya_df %>% select(long, lat)),
         lakes = raster::extract(distance_to_lakes_r, STH_kenya_df %>% select(long, lat)),
         canals = raster::extract(distance_to_canals_r, STH_kenya_df %>% select(long, lat)),
         roads = raster::extract(distance_to_roads_r, STH_kenya_df %>% select(long, lat)),
         railroads = raster::extract(distance_to_railroads_r, STH_kenya_df %>% select(long, lat)))

STH_kenya_df %>% write.csv("STH_kenya_covariates.csv", row.names = F)
