# load libraries ----
library(ggplot2)
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

# view vaccination prevalence data ----
# load vaccination data
BCG_vaccination_UGA <- read.table("BCG_vaccination_UGA.txt", header = T, sep = ",", dec = ".") # vaccination data

# calculate vaccination prevalence
BCG_vaccination_UGA <- BCG_vaccination_UGA %>%
  mutate(prevalence = number_positive / numer_examined)

# global options
options(stringsAsFactors = F)

# load covariates and set crs ----
ntl_2013_raster <- 'gp2013africa.tif' %>% raster() # Annual mighttime lights (NTL) images for Africa, 2013
ntl_2016_raster <- raster("SVDNB_npp_20160101-20161231_00N060E_v10_c201807311200/SVDNB_npp_20160101-20161231_00N060E_vcm-orm_v10_c201807311200.avg_rade9.tif")
roads <- readOGR("UGA_rds", "UGA_roads") #roads in Uganda
railroads <- readOGR("UGA_rrd", "UGA_rails") #railroads in Uganda
waterbodies <- readOGR("UGA_wat", "UGA_water_areas_dcw") #inland waterbodies in Uganda
canals <- readOGR("UGA_wat", "UGA_water_lines_dcw") #inland water lines in Uganda
health_facilities <- read.csv("UGA_00_SSA_MFL_(130219).csv") %>% drop_na(Lat)

# convert to spdf
BCG_vaccination_UGA_spdf <- SpatialPointsDataFrame(coords = BCG_vaccination_UGA %>% select(lng,lat),
                                                   data = BCG_vaccination_UGA %>%
                                                     select(ID, numer_examined, number_positive, prevalence))

health_facilities_spdf <- SpatialPointsDataFrame(coords = health_facilities %>% select(Long,Lat),
                                                 data = health_facilities %>%
                                                   select(Facility.name, Lat, Long))

crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj4string(BCG_vaccination_UGA_spdf) <- crs %>% CRS()
proj4string(roads) <- crs %>% CRS()
proj4string(railroads) <- crs %>% CRS()
proj4string(waterbodies) <- crs %>% CRS()
proj4string(canals) <- crs %>% CRS()
proj4string(health_facilities_spdf) <- crs %>% CRS()

# plot on leaflet ----
pal = colorNumeric("Oranges", BCG_vaccination_UGA$prevalence)

leaflet() %>% addProviderTiles("Stamen.TonerLite") %>%
  addCircleMarkers(data = BCG_vaccination_UGA,
                   fillOpacity=1,
                   fillColor= ~pal(prevalence),
                   radius=~prevalence*2, stroke=TRUE, weight=.1) %>%
  addPolylines(data = roads, weight=1.5, color = "dimgrey") %>%
  addPolylines(data = railroads, weight=2, color = "grey") %>%
  addPolygons(data = waterbodies, weight=1.5, color="dodgerblue", opacity = .5) %>%
  addLegend(pal = pal, values = BCG_vaccination_UGA$prevalence, title = "Vacc Prev (%)")

# calc distance between points and waterbodies and points and roads ----
# dist_to_water <- BCG_vaccination_UGA_spdf %>% dist2Line(waterbodies) #takes a while
# dist_to_roads <- BCG_vaccination_UGA_spdf %>% dist2Line(roads) #takes forever, see pre-written csv file
# dist_to_railroads <- BCG_vaccination_UGA_spdf %>% dist2Line(railroads)

# BCG_vaccination_UGA_spdf$dist_to_water <- dist_to_water[,1]
# BCG_vaccination_UGA_spdf$dist_to_roads <- dist_to_roads[,1]
# BCG_vaccination_UGA_spdf$dist_to_railroads <- dist_to_railroads[,1]

# convert spdf to df and append nighttime lights raster layer
BCG_vaccination_UGA <- BCG_vaccination_UGA_spdf %>% as.data.frame()

BCG_vaccination_UGA <- BCG_vaccination_UGA %>%
  mutate(ntl_2013 = raster::extract(ntl_2013_raster, BCG_vaccination_UGA %>% select(lng, lat)))

# ntl data
# BCG_vaccination_kampala <- BCG_vaccination_UGA_spdf[Kampala, ]
# railroads_kampala <- railroads[Kampala, ]
# waterbodies_kampala <- waterbodies[Kampala, ]
# canals_kampala <- canals[Kampala, ]
