#simplest data is point data (i.e. a table with coordinates) (lat-long or x-y?)
# we'll work with malaria point prevalence data in Ethiopia

# Working with Spatial Data in R ----
# global options
options(stringsAsFactors = F)
# first load libraries
library(sp)
library(raster)
library(rgdal)
library(leaflet)
library(tidyverse)

# import data
ETH_malaria_data <- read.csv("https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week1/Lab_files/Data/mal_data_eth_2009_no_dups.csv",
                             header=T)
# examine columns
head(ETH_malaria_data)
# pf_pr = infection prevalence, proportion infected (pf_pos/examined)
# examined: numbers tested
# pf_pos = of those tested, how many were positive for Pf malaria
# longitude = longitude of school in decimal degrees
# latitude = latitude of school in decimal degrees

hist(ETH_malaria_data$pf_pr, breaks=10)

# Plotting and mapping spatial data ----

plot(ETH_malaria_data$longitude, ETH_malaria_data$latitude,
     cex = ETH_malaria_data$pf_pr*10, #size of the circle is a function of infection prevalence rate
     ylab = "Latitude",
     xlab = "Longitude")

# Using the sp package, we can package spatial data up into a "Spatial" class of objects
# This makes it easier to work with and is often a requirement for other functions
# SpatialPoints or SpatialPolygons
# For data associated with each spatial feature, you can use:
# SpatialPointsDataFrames
# SpatialPolygonsDataFrames

# Spatial Points DataFrame ----
ETH_malaria_data_SPDF <- SpatialPointsDataFrame(
  coords = ETH_malaria_data[,c("longitude", "latitude")], # specify lat-longs via "coords"
  data = ETH_malaria_data[,c("examined", "pf_pos", "pf_pr")], #specify columns to include?
  proj4string = CRS("+init=epsg:4326") #specify CRS to WGS1984, optional but good to specify
)

# SPDFs partition data elements
# this means the coordinates are stored separately from the data
head(ETH_malaria_data_SPDF@coords)

head(ETH_malaria_data_SPDF@data) # the 3 columns specified earlier!

# accessing the data as a standard data frame:
head(ETH_malaria_data_SPDF$examined)
head(ETH_malaria_data_SPDF$pf_pr)

# you can use the plot or spplot function to get quick plots
plot(ETH_malaria_data_SPDF)
spplot(ETH_malaria_data_SPDF, zcol = "pf_pr")

# Spatial Polygons Data Frame ----
library(rgdal)
ETH_Adm_1 <- readOGR("ETH_Adm1_shapefile", "ETH_Adm_1")

# getData function from the raster package gives you boundary data (for Ethiopia? we'll see)

# First need IS03 codes for the country of interest, get these using "ccodes()"
# For Ethiopia, the IS03 is ETH
# The getData function allows you to retrieve the relevant admin level boundaries from GADM

ETH_Adm_1 <- raster:: getData("GADM", country = "ETH", level = 1)

plot(ETH_Adm_1) ## plot the country boundaries

## plot the malaria prevalence points
points(ETH_malaria_data$longitude, ETH_malaria_data$latitude,
       cex = ETH_malaria_data$pf_pr * 10,
       ylab = "Latitude", xlab = "Longitude",
       col="red")

## plotting data using WEB MAPS! :) ----
library(leaflet)
basemap <- leaflet() %>% addProviderTiles("Esri.OceanBasemap")
basemap

## adding the Ethiopia bounadaries to the Esri Ocean BaseMap
basemap %>% addPolygons(data=ETH_Adm_1) 

## stylin it!
basemap %>% addPolygons(data=ETH_Adm_1, color="purple",
                        weight = 1, fillOpacity = 0.2)
## adding popups!
basemap %>% addPolygons(data=ETH_Adm_1,
                        color = "purple",
                        weight = 1,
                        fillOpacity = 0.2,
                        popup = ETH_Adm_1$NAME_1)
## you can overlay the points on this as well
basemap %>% addPolygons(data=ETH_Adm_1,
                        color = "purple",
                        weight = 1,
                        fillOpacity = 0.2,
                        popup = ETH_Adm_1$NAME_1) %>%
  addCircleMarkers(
    ## SpatialPointsDataFrame containing malaria point-prevalence data
    data=ETH_malaria_data_SPDF,
    ## circles will appear red with a radius of 2 (units?)
                   color="darkblue", radius = 0.3)

## using color palettes from leaflet to turn some points into difference colors (maybe for cases v. controls?)
library(wesanderson)

colorPal <- colorNumeric(wes_palette("Zissou1")[1:5],
                         ETH_malaria_data_SPDF$pf_pr,
                         n = 5)
# colorPal is now a function you can apply to get the corresponding color for  value
colorPal(0.6)

basemap %>% addPolygons(data=ETH_Adm_1,
                        weight = 2,
                        fillOpacity = 0,
                        popup = ETH_Adm_1$NAME_1) %>%
  addCircleMarkers(data = ETH_malaria_data_SPDF,
                   # applying ColorPal to values of the malaria prevalence points
                   color = colorPal(ETH_malaria_data_SPDF$pf_pr),
                   radius = 2,
                   # applying a popup which is a label for this data
                   popup = as.character(ETH_malaria_data_SPDF$pf_pr)) %>%
  
  ## now add a legend
  
  addLegend(pal = colorPal,
            title = "Prevalence",
            values = ETH_malaria_data_SPDF$pf_pr)

## more complicated pop-ups using HTML

basemap %>% addPolygons(data=ETH_Adm_1,
                        weight = 2,
                        fillOpacity = 0,
                        popup = ETH_Adm_1$NAME_1) %>%
  addCircleMarkers(data = ETH_malaria_data_SPDF,
                   # applying ColorPal to values of the malaria prevalence points
                   color = colorPal(ETH_malaria_data_SPDF$pf_pr),
                   radius = 2,
                   # applying a popup which is a label for this data
                   popup = paste("<p>", "Prevalence:",
                                 round(ETH_malaria_data_SPDF$pf_pr,2),
                                 "<p>")) %>%
  
  ## now add a legend
  
  addLegend(pal = colorPal,
            title = "Prevalence",
            values = ETH_malaria_data_SPDF$pf_pr)

# Plotting Raster Data ----

## load raster of ELEVATION in Ethiopia into memory from Hugh's github:
elev <- raster(
  "https://github.com/HughSt/HughSt.github.io/raw/master/course_materials/week1/Lab_files/Data/elev_ETH.tif"
)

## getData function from raster package

elev <- raster::getData("alt", country="ETH")
elev

## plotting the elevation using standard plot function
plot(elev)

## ...or use leaflet which is dope
basemap %>% addRasterImage(elev)

## define color palette for the legend, then add legend
raster_colorPal <- colorNumeric(
  topo.colors(64),
  values(elev),
  na.color = NA)

# PLot
basemap %>% addRasterImage(elev,
                           color = raster_colorPal) %>%
  addLegend(values = values(elev),
            # pal or palette is where you input the legend palette you defined earlier
            pal = raster_colorPal)

# exporting the data
library(plotKML)
# see ?plotKML for more options
plotKML(ETH_malaria_data_SPDF)

## Notes ----

#sp package: important
# raster package: important
# sf: provides alternative ways of handing spatial data (not as important?)
# sf is good for larger datasets

# spatial epidemiology! WHERE diseases happen, and why they happen there












  



