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
library(ape)
library(spdep) # Spatial Dependence: Weighting Schemes, Statistics and Models
library(pgirmess) # Data Analysis in Ecology
library(sp)
library(gtools)
library(car)

# global options ----
options(stringsAsFactors = F) # import quality
par(mar=c(1,1,1,1)) # ensure plot margin sizes are always large enough

# load soil-transmitted helminths data ----
STH_kenya <- "Kenya_STH.csv" %>% read.csv(header = T)
STH_kenya <- STH_kenya %>% select(tri_prev)

STH_kenya <- STH_kenya %>%
  mutate(tri_prev = tri_prev/100,
         tri_prev_logit = car::logit(tri_prev) + 0.0001)

# plot histograms: Trichuris trichiura (Whipworm) ----
STH_kenya %>% ggplot(aes(x=tri_prev)) + 
  geom_histogram(aes(y=..density..), colour="darkgrey", fill="lightblue")+
  geom_density(alpha=.2, fill="#FF6666") +
  xlab("T. trichiura prevalence") + ylab("Density") +
  ggtitle("Distribution of T. trichiura prevalence (tri_prev)")

STH_kenya %>% ggplot(aes(x=tri_prev_logit)) + 
  geom_histogram(aes(y=..density..), colour="darkgrey", fill="lightblue")+
  geom_density(alpha=.2, fill="#FF6666") +
  xlab("T. trichiura prevalence") + ylab("Density") +
  ggtitle("Distribution of T. trichiura prevalence - logit transform")

# plot on leaflet - whipworm ----
pal = colorNumeric("Oranges", STH_kenya$tri_prev_logit) 

STH_kenya %>%
  leaflet() %>% addProviderTiles("CartoDB.DarkMatter") %>%
  addCircleMarkers(~long, ~lat, fillOpacity=1,
                   fillColor= ~pal(tri_prev_logit),
                   radius=3,
                   weight=0.1,
                   stroke=TRUE) %>%
  addLegend(pal = pal, values = ~tri_prev_logit, title = "Whipworm Prevalence (%)")
