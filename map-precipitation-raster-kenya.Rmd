---
title: "Precipitation in Kenya"
author: "Sri Ramesh"
date: "4/7/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(raster)
library(leaflet)
```

In the following, I pulled WordClim data on the long-term average precipitation in Kenya each January from 1970 to 2000 to illustrate on a leaflet map. The time-series data can be found here: *worldclim.org*.

## Q1 Restricting precipitation data to January

```{r, warning=F, message=F}
KEN_Adm0 <- raster::getData(name = "GADM", country = "KEN", level = 0)
KEN_prec_0.5_bottom <- raster::getData(name = "worldclim", var = "prec", res = 0.5, lat = -3, lon = 39)
KEN_prec_0.5_top <- raster::getData(name = "worldclim", var = "prec", res = 0.5, lat = 3, lon = 35)

KEN_prec_0.5 <- raster::merge(x = KEN_prec_0.5_bottom, y = KEN_prec_0.5_top)

KEN_prec_0.5_Jan <- KEN_prec_0.5[[1]]
```

## Q2 Aggregating data to 1 minute of a degree resolution

```{r, warning=F, message=F}
KEN_prec_0.5_Jan_inches <- KEN_prec_0.5_Jan * 0.0393701

KEN_prec_1_Jan_inches <- KEN_prec_0.5_Jan_inches %>%
  raster::aggregate(fact = 2, fun = sum)
```

## Q3 Cropping and masking to Kenya

```{r, warning=F, message=F}
KEN_prec_1_Jan_Crop_Unmasked <- KEN_prec_1_Jan_inches %>% raster::crop(y = KEN_Adm0)
KEN_prec_1_Jan_Crop <- KEN_prec_1_Jan_Crop_Unmasked %>% raster::mask(mask = KEN_Adm0)
```

## Q4 Adding a leaflet base map

```{r, warning=F, message=F}
basemap <- leaflet() %>% addProviderTiles("CartoDB.Positron")
```

## Q5 Adding a legend with a title

```{r, warning=F, message=F}
raster_colorPal_prec_JAN <- colorNumeric(palette = topo.colors(64), domain = values(KEN_prec_1_Jan_Crop), 
                                         na.color = NA)

map_kenya <- basemap %>%
  addRasterImage(x = KEN_prec_1_Jan_Crop,
                 color = raster_colorPal_prec_JAN) %>%
  addPolygons(data = KEN_Adm0,
              fillOpacity = 0,
              color = "purple",
              weight = 3)

map_kenya <- map_kenya %>%
  addLegend(title = "Jan precipitation (in)<br>(1' res)",
            values = values(KEN_prec_1_Jan_Crop),
            pal = raster_colorPal_prec_JAN)

map_kenya
```
