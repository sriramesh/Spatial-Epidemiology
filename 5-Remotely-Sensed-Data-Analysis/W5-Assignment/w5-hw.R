# load libraries----
library(tidyverse)
library(raster)
library(leaflet)

# India ----
IND_Adm0 <- raster::getData(name = "GADM", country = "IND", level = 0)

IND_prec_2.5 <- raster::getData(name = "worldclim", var = "prec",
                                res = 2.5)

# Restrict to Jul
IND_prec_2.5_Jul <- IND_prec_2.5[[7]]

# Crop and Mask to India extent 
IND_prec_2.5_Jul_Crop_Unmasked <- raster::crop(x = IND_prec_2.5_Jul,
                                             y = IND_Adm0)
IND_prec_2.5_Jul_Crop <- raster::mask(x = IND_prec_2.5_Jul_Crop_Unmasked,
                                    mask = IND_Adm0)

# Layer map
raster_colorPal_prec_JUL <- colorNumeric(palette = topo.colors(64),
                                         domain = values(IND_prec_2.5_Jul_Crop),
                                         na.color = NA)

leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addRasterImage(x = IND_prec_2.5_Jul_Crop,
                 color = raster_colorPal_prec_JUL) %>%
  addPolygons(data = IND_Adm0,
              popup = IND_Adm0$NAME_0,
              label = IND_Adm0$NAME_0,
              fillOpacity = 0,
              color = "purple",
              weight = 3) %>%
  addLegend(title = "Jul precipitation (mm)<br>(2.5' res)",
            values = values(IND_prec_2.5_Jul_Crop),
            pal = raster_colorPal_prec_JUL)


# Kenya ----
# Precipitation from worldclim at the 0.5 resolution (minutes of a degree)... 
KEN_Adm0 <- raster::getData(name = "GADM", country = "KEN", level = 0)

# bottom half of Kenya (below equator)
KEN_prec_0.5_bottom <- raster::getData(name = "worldclim", var = "prec",
                                       res = 0.5, lat = -3, lon = 39)

# top half of Kenya (above equator)
KEN_prec_0.5_top <- raster::getData(name = "worldclim", var = "prec",
                                    res = 0.5, lat = 3, lon = 35)
# merge two rasters
KEN_prec_0.5 <- raster::merge(x = KEN_prec_0.5_bottom,
                              y = KEN_prec_0.5_top)

# Restrict to Jan
KEN_prec_0.5_Jan <- KEN_prec_0.5[[1]]

# Aggregating data to 1 minute of a degree resolution
KEN_prec_0.5_Jan_inches <- KEN_prec_0.5_Jan * 0.0393701

KEN_prec_1_Jan_inches <- KEN_prec_0.5_Jan_inches %>%
  raster::aggregate(fact = 2, fun = sum)

# Crop to Kenya extent 
KEN_prec_1_Jan_Crop_Unmasked <- raster::crop(x = KEN_prec_1_Jan_inches,
                                               y = KEN_Adm0)

# Mask to Kenya
KEN_prec_1_Jan_Crop <- raster::mask(x = KEN_prec_1_Jan_Crop_Unmasked,
                                      mask = KEN_Adm0)

# Layer map
raster_colorPal_prec_JAN <- colorNumeric(palette = topo.colors(64),
                                         domain = values(KEN_prec_1_Jan_Crop),
                                         na.color = NA)

leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addRasterImage(x = KEN_prec_1_Jan_Crop,
                 color = raster_colorPal_prec_JAN) %>%
  addPolygons(data = KEN_Adm0,
              popup = KEN_Adm0$NAME_0,
              label = KEN_Adm0$NAME_0,
              fillOpacity = 0,
              color = "purple",
              weight = 3) %>%
  addLegend(title = "Jan precipitation (in)<br>(1' res)",
            values = values(KEN_prec_1_Jan_Crop),
            pal = raster_colorPal_prec_JAN)
