---
title: "Children Underweight in Burundi"
author: "Sri Ramesh"
date: "3/16/2020"
output: html_document
---

```{r setup, include=F, eval=T, echo=F}
# set global options  ----
options(stringsAsFactors = F)

# load packages  ----
library(tidyverse)
library(leaflet)
library(rgdal)
library(htmltools)
library(htmlwidgets)

# read in shapefile of Burundi in 2016 DHS ----
dhs_2016_burundi <- readOGR("input/sdr_subnational_data_2020-03-16/shps", "sdr_subnational_data_dhs_2016")

```

This document outlines the 4 steps I took to produce a chloropleth map of the percentage of children underweight from the [2016 Demographic Health Survey (DHS) done in Burundi](https://spatialdata.dhsprogram.com/data/#/):

1. Load basemap
2. Add polygons of Burundi sub-national boundaries, 2016 DHS
3. Add polygons of the % of children underweight per unit of aggregation, 2016 DHS
4. Add legend with a title

*"Percentage of children underweight"* is the percentage of children underweight (below -2 SD of weight for age according to the WHO standard). This data is found in the data linked above.

### 1. Load basemap

I really enjoy these tiles of an oceanic basemap from Esri so I used those as my basemap:

```{r basemap, warning=F, message=F, echo=T}

basemap <- leaflet() %>% addProviderTiles("Esri.OceanBasemap")
basemap

```


### 2. Add polygons of Burundi sub-national boundaries, 2016 DHS

The DHS' Spatial Data Repository provides shapefiles of Burundi's sub-national boundaries. I loaded the boundaries with minimal opacity so that the names of these sub-national regions can be seen through from the underlying basemap itself:

```{r boundaries, warning=F, message=F, echo=T} 
# Overlay sub-national boundaries of Burundi per the 2016 DHS
dhs_burundi_map <- basemap %>% addPolygons(data=dhs_2016_burundi,
                        color = "purple",
                        weight = 1,
                        fillOpacity = 0.2)
dhs_burundi_map
```


### 3. Add polygons of the % of children underweight per unit of aggregation, 2016 DHS

I chose to make a chloropleth map of the % of children underweight across the sub-national regions in Burundi. This data was contained in the same shapefile that provided the sub-national boundaries in Burundi per the 2016 DHS:

```{r chloropleth, warning=F, message=F, echo=T} 
# Define color palette for chloropleth map
bins <- c(10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0)
pal <- colorBin("Spectral", domain = dhs_2016_burundi$CNNUTSCWA2, bins = bins)

# Shade in sub-national polygons per the color palette, and add various effects
dhs_burundi_map <- dhs_burundi_map %>%
  addPolygons(data=dhs_2016_burundi,
  fillColor = ~pal(dhs_2016_burundi$CNNUTSCWA2),
  weight = 0.7,
  opacity = 0,
  color = "white",
  dashArray = "3",
  fillOpacity = 0.35,
  highlight = highlightOptions(
    weight = 3,
    color = "#666",
    dashArray = "",
    fillOpacity = 0,
    bringToFront = TRUE),
  label = paste("% children underweight:", dhs_2016_burundi$CNNUTSCWA2),
  labelOptions = labelOptions(
    style = list("font-weight" = "normal", padding = "3px 8px"),
    textsize = "15px",
    direction = "auto"))
dhs_burundi_map

```




### 4. Add legend with a title

Added in a legend with a title for the final product:

```{r legend, warning=F, message=F, echo=T} 
## Add legend
dhs_burundi_map <- dhs_burundi_map %>% addLegend(pal = pal,
                     values = dhs_2016_burundi$CNNUTSCWA2,
                     opacity = 0.35,
                     title = "% children </br> underweight, </br> 2016 DHS in </br> Burundi",
                position = "topright")

dhs_burundi_map

```

# Citations
