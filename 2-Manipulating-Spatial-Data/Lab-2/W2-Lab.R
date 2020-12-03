## OBJECTIVES ----
# 1) Learning how to CROP and SUBSET Spatial Data!!
# 2) Go through the process of RESAMPLING RASTERS in R as well

## Load libraries ----
library(sp)
library(raster)
library(leaflet)
library(rgdal)
library(geosphere)
library(rgeos)
library(wesanderson)
library(stats)
library(RANN)
library(geosphere)

## Set global options ----
options(stringsAsFactors = F)




## Load data ----

ETH_Adm_1 <- raster::getData("GADM", country="ETH", level = 1)
ETH_Adm_1_cropped <- ETH_Adm_1[1,] ## just use the subset function

# what does this subset? the **FIRST ROW**

# Get a summary of the cropped data
ETH_Adm_1_cropped # remember this is the first row

# Plot over the top of the full dataset
plot(ETH_Adm_1)
lines(ETH_Adm_1_cropped, col="red", lwd=2)

# Subsetting by name
ETH_Adm_1_Amhara <- subset(ETH_Adm_1, ETH_Adm_1$NAME_1=="Amhara")
plot(ETH_Adm_1) #plot the main dataset first
lines(ETH_Adm_1_Amhara, col="blue", lwd=2) #then plot the subsetted part over it (overlay)


## Pop Quiz 1) How would you plot all provinces EXCEPT Amhara? ----
ETH_Adm_1_allbutAmhara <- subset(ETH_Adm_1, ETH_Adm_1$NAME_1 != "Amhara")
plot(ETH_Adm_1) #plot the main dataset first
lines(ETH_Adm_1_allbutAmhara, col="purple", lwd=2) 

## Pop Quiz 2) Try plotting all provinces using leaflet, with Amhara colored red and all others colored orange ----
basemap <- leaflet() %>% addProviderTiles("Esri.OceanBasemap")
basemap

## adding the Ethiopia boundaries to the Esri Ocean BaseMap
basemap %>% addPolygons(data=ETH_Adm_1_Amhara,
                        color="red",
                        weight = 1,
                        fillOpacity = 0.2) %>%
  addPolygons(data=ETH_Adm_1_allbutAmhara,
                        color="orange",
                        weight = 1,
                        fillOpacity = 0.2)

# Spatial Overlays ----
# Often we have point and polygon data and wish to relate them
# we might want to summarize point data over regions
# aggregating point data -> province level on malaria across Ethiopia

# Load data
ETH_malaria_data <- read.csv("https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week1/Lab_files/Data/mal_data_eth_2009_no_dups.csv",header=T)
class(ETH_malaria_data)
# Convert data to spatial data type (a SPDF or spatial points data frame) using sp package
ETH_malaria_data_SPDF <- SpatialPointsDataFrame(
  coords = ETH_malaria_data[,c("longitude", "latitude")], #let's use lat longs
  data = ETH_malaria_data[,c("examined", "pf_pos", "pf_pr")], # i only want these 3 cols
  proj4string = CRS("+init=epsg:4326")) # set project CRS

# this line of code will error out...
ETH_Adm_1_per_point <- over(ETH_malaria_data_SPDF, ETH_Adm_1)
# this is because the two datasets do not have the exact same CRS!! CRS must match.
crs(ETH_malaria_data_SPDF)
crs(ETH_Adm_1)

# reproject the same data using the spTransform function, also from the sp package BY....
# changing the CRS so that the spatial dataframes one matches the regular dataframe one
ETH_malaria_data_SPDF <- spTransform(ETH_malaria_data_SPDF, crs(ETH_Adm_1))

# check the CRS to make sure it worked
crs(ETH_malaria_data_SPDF)

# re-reun the over command that initially errored out...
ETH_Adm_1_per_point <- over(ETH_malaria_data_SPDF, ETH_Adm_1)
# we're combining POINT data with BOUNDARIES

# it worked now! wooooOooooOoo

# this gives us a table where:
# - each row represents a point from ETH_malaria_data_SPDF (point data)
# - each column represents data from ETH_Adm_1 (areal boundaries data)
head(ETH_Adm_1_per_point)

# Using this combined data to calculate STATS at the unit of aggregation (admin unit specific)
# Num of sites per admin unit:
table(ETH_Adm_1_per_point$NAME_1) # there are 201 in Oromia, 1 in the other 2

# apply a function across groups:
Nexamined_per_Adm1 <- tapply(ETH_malaria_data_SPDF$examined,
                             ETH_Adm_1_per_point$NAME_1,
                             sum) # a sum of "examined" per NAME_1
Nexamined_per_Adm1

Npositives_per_Adm1 <- tapply(ETH_malaria_data_SPDF$pf_pos,
                              ETH_Adm_1_per_point$NAME_1,
                              sum) # a sum of "pf_pos" per NAME_1
Npositives_per_Adm1

# Now divide the positives by the number examined to get the prevalence rate per Admin unit:
prev_per_Adm1 <- Npositives_per_Adm1 / Nexamined_per_Adm1
prev_per_Adm1 # it's made a vector

# turn this vector into a data frame
prev_per_Adm1_df <- data.frame(NAME_1 = names(prev_per_Adm1),
                               prevalence = prev_per_Adm1,
                               row.names = NULL)

# merge dataframe with the existing data so final product is an added COLUMN for prevalence
ETH_Adm_1 <- merge(ETH_Adm_1, prev_per_Adm1_df,
                   by = "NAME_1")
head(ETH_Adm_1)

## now we can plot the prevalence data on our pretty leaflet map!!

# remember, you must define your color palette/legend:

colorPal <- colorNumeric(wes_palette("Zissou1")[1:5], ETH_Adm_1$prevalence)

# create chloropleth:
leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(data = ETH_Adm_1,
              weight = 1,
              col=colorPal(ETH_Adm_1$prevalence), ##apply color palette to values of this column
              fillOpacity = 0.4) %>%
  
  # and then add a legend
  
  addLegend(pal = colorPal,
            values = ETH_Adm_1$prevalence,
            title = "Prevalence")

## Pop Quiz 1) Try generating the same plot using a different color palette ----

palette <- "Darjeeling1"

# copied code
colorPal <- colorNumeric(wes_palette(palette)[1:5], ETH_Adm_1$prevalence)

# create chloropleth:
leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(data = ETH_Adm_1,
              weight = 1,
              col=colorPal(ETH_Adm_1$prevalence), ##apply color palette to values of this column
              fillOpacity = 0.4) %>%
  
  # and then add a legend
  
  addLegend(pal = colorPal,
            values = ETH_Adm_1$prevalence,
            title = "Prevalence")


## Pop Quiz 2) How would you plot only the provinces for which you have prevalence estimates? ----

ETH_prev_Adm1_only <- subset(ETH_Adm_1, !is.na(ETH_Adm_1$prevalence))

leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(data = ETH_prev_Adm1_only,
              weight = 1,
              col=colorPal(ETH_Adm_1$prevalence), ##apply color palette to values of this column
              fillOpacity = 0.4) %>%
  
  # and then add a legend
  
  addLegend(pal = colorPal,
            values = ETH_Adm_1$prevalence,
            title = "Prevalence")













## Manipulating Raster Data ----

# Load Ethiopia elevation data:
ETH_elev <- raster::getData("alt", country="ETH")
plot(ETH_elev)

# Load Ethiopia land use data:
# Land use (# For information on land use classifications see:
# http://due.esrin.esa.int/files/GLOBCOVER2009_Validation_Report_2.2.pdf)
ETH_land_use <- raster(
  "https://github.com/HughSt/HughSt.github.io/blob/master/course_materials/week2/Lab_files/ETH_land_use.tif?raw=true"
  )
ETH_land_use

# Plot the land use raster:
plot(ETH_land_use)

# For a break down of the land use CLASSES in Ethiopia
# (aka how often each land use type occurs in Ethiopia)
# Note: this is just the number of pixels per land use type - NOT acres
table(ETH_land_use[])

# Resampling rasters -- do it!!
# remember -- this line of code will take TIME.... 
# We are just moving from ONE GRID (e.g. 9 by 9 pixels) to another (e.g. 3 by 3 pixels)
ETH_land_use_resampled <- resample(ETH_land_use, ## larger pixels (zoomed in)
                                   ETH_elev, ## smaller pixels (zommed out)
                                   method="ngb")

# Get summaries of both raster objects to check resolution and EXTENT!!
# and to see whether RESAMPLED values look right
ETH_land_use_resampled
ETH_elev

## Manipulating Rasters ----
res(ETH_elev) # in decimal degrees. 1 decimal degree = 111km at the equator

# Let's aggregate by a factor of 10 (i.e. making it lower resolution, less detailed)
ETH_elev_low_res <- aggregate(ETH_elev,
                              fact = 10) # by default, calculates mean
res(ETH_elev_low_res)
plot(ETH_elev_low_res)

# categorizing raster values:
ETH_elev_categorized <- cut(ETH_elev, 4) # range of colors is from 1 to 4
plot(ETH_elev_categorized)

# performing joint operations on rasters of the same resolution and extent
# like subtracting values from one another
new_raster <- ETH_elev - ETH_land_use_resampled # just for illustrative purposes
plot(new_raster)

## Extracting data from rasters ----
# appended an "elev" column to this data now:
ETH_malaria_data_SPDF$elev <- extract(ETH_elev,
                                      ETH_malaria_data_SPDF)
ETH_malaria_data_SPDF 

# you can also extract values using the polygons to get admin 1 level elevations
# you just have to define a function to apply, otherwise you get all the pixel values per polygon
# for very large rasters, check out the velox package

## extract values using POLYGONS to get admin 1 level elevations
# This will take a while...
ETH_Adm_1$elev <- extract(ETH_elev, ETH_Adm_1, fun=mean, na.rm=TRUE) # takes a little longer..

## Exploratory spatial analysis ----

# Analyze the relationship between prevalence and elevation:
ETH_malaria_data_SPDF$prevalence <-
  ETH_malaria_data_SPDF$pf_pos / ETH_malaria_data_SPDF$examined

plot(ETH_malaria_data_SPDF$elev,
     ETH_malaria_data_SPDF$prevalence)

# Overlay water bodies data
waterbodies <- readOGR(
  "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week2/Lab_files/ETH_waterbodies.geojson"
  )
waterbodies
plot(waterbodies)

basemap %>% addPolygons(data=waterbodies, #can add geojson items to leaflet
                        color="red",
                        weight = 1,
                        fillOpacity = 0.2)

# the "dist2Line" function from the geosphere package calculates distance in meters from
# spatial data recorded using decimal degrees
dist_to_water <- dist2Line(ETH_malaria_data_SPDF, waterbodies)
head(dist_to_water)

# Can add to your data frame by extracting the first column:
ETH_malaria_data_SPDF$dist_to_water <- dist_to_water[,1]

# must calculate the distance to EVERY POINT, then identify the minimum
# calculate the centroid (the center) of each polygon:
waterbodies_points <- gCentroid(waterbodies, byid=T) # polygons = waterbodies

# show distances between each malaria point and waterbody centroid using a distance matrix
# distm works OK if you don't have too many points
# 1 degree at the equator is a larger distance than nearer the poles
dist_matrix <- distm(ETH_malaria_data_SPDF, waterbodies_points) # must use meters as DEGREES

# use the apply function to apply the minimum function to each row
# as each row represents the distance of every waterbody point from our first obs
ETH_malaria_data_SPDF$dist_to_water_point <- apply(dist_matrix, 1, min)

# get the index of the waterbody points/centroids NEAREST to each malaria case point
nn <- nn2(waterbodies_points@coords,ETH_malaria_data_SPDF@coords, k=1)

# Calculate the distance in meters between each observations and its nearest waterbody point
ETH_malaria_data_SPDF$dist_to_water_point <-
  distGeo(ETH_malaria_data_SPDF@coords,
          waterbodies_points@coords[nn$nn.idx,])

# office hours Qs ----

# x 1. head of the table sufficient?
# for table:
# - similar to in the lab, how we generated table. PREVALENCE PER PROVINCE.
# - code chunk that has admin 2 name + prevalence below them

# x 1. it looks like this data only has 3 levels (admin level 1)? does this look right?

# 2 I only generated 1 plot and it looks like this. what other kinds of plots can i draw?

# 3 

# ETH_alt <- raster::getData("alt", country="ETH", level = 2)
# 
# colorPal <- colorNumeric(heat.colors(5, 0.6, 0.4), ETH_Adm_2$prevalence)
# 
# pal <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(ETH_alt),
#                     na.color = "transparent")
# 
# alt_map <- basemap %>% addRasterImage(ETH_alt, colors = pal, opacity = 0.8) 
# 
# alt_map %>%
#   addPolygons(data=ETH_prev_Adm2_only,
#               color="black",
#               weight = 1,
#               fillOpacity = 0) %>%
#   addPolygons(data = ETH_prev_Adm2_only,
#               weight = 1,
#               col=colorPal(ETH_prev_Adm2_only$prevalence),
#               fillOpacity = 0.5,
#               label = paste("Prevalence of Malaria:",
#                             round(ETH_prev_Adm2_only$prevalence,4),"%")) %>%
#   addLegend(pal = pal, values = values(ETH_alt),
#             title = "Elevation </br> in Ethiopia")
# ETH_malaria_data_SPDF$alt <- raster::extract(ETH_alt, ETH_malaria_data_SPDF)

# plot(ETH_malaria_data_SPDF$alt, ETH_malaria_data_SPDF$prevalence)


# ETH_Adm_2$landclass <- raster::extract(ETH_land_use, ETH_Adm_2)
## 

# visualize the land class raster, aggregate, and fact of 20

# only a few of the classes are actually observed
# A11 and A12 are observed 



# calculate mean per land use class
# calculate boxplot (5 number summary) of prevalence by land use class














