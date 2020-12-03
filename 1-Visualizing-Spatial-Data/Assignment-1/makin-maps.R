# set global options  ----
options(stringsAsFactors = F)


# load packages  ----
library(tidyverse)
library(leaflet)
library(rgdal)
library(htmltools)
library(htmlwidgets)

# read in shapefile of Burundi in 2016 DHS ----
dhs_2016_burundi <- readOGR("sdr_subnational_data_2020-03-16/shps", "sdr_subnational_data_dhs_2016")


# define basemap ----
basemap <- leaflet() %>% addProviderTiles("Esri.OceanBasemap")
basemap


# overlay burundi DHS sub-national boundaries ----
m <- basemap %>% addPolygons(data=dhs_2016_burundi,
                        color = "purple",
                        weight = 1,
                        fillOpacity = 0.2)


# shade in chloropleth -----

# define color palette for legend
bins <- c(10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, Inf)
pal <- colorBin("Spectral", domain = dhs_2016_burundi$CNNUTSCWA2, bins = bins)

# make chloropleth
m <- m %>%
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
  label = dhs_2016_burundi$CNNUTSCWA2,
  labelOptions = labelOptions(
    style = list("font-weight" = "normal", padding = "3px 8px"),
    textsize = "15px",
    direction = "auto"))
m


# add legend and title ----

tag.map.title <- tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(-50%,20%);
    position: fixed !important;
    left: 50%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 28px;
  }
"))

title <- tags$div(
  tag.map.title, HTML("% of Children Underweight in Burundi, 2016 DHS")
)  

m <- m %>% addLegend(pal = pal,
                values = dhs_2016_burundi$CNNUTSCWA2,
                opacity = 0.35,
                title = "% children </br> underweight, </br> 2016 DHS in </br> Burundi",
                position = "topright") %>%
  addControl(title, position = "topleft", className="map-title")

m

