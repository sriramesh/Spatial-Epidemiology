library(leaflet)
library(tidyverse)

#Esri.WorldGrayCanvas
#Esri.OceanBasemap
#Esri.WorldImagery
#Esri.DeLorme
#Stamen.Watercolor
#Stamen.TonerLite
#Stamen.Terrain
#CartoDB.Positron (this is kind of standard?)

basemap <- leaflet() %>% addProviderTiles("Esri.WorldGrayCanvas")
basemap