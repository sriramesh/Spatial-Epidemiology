# Set global options ----
options(stringsAsFactors = F)

# Load libraries ----
library(tidyverse)
library(sp)
library(raster)
library(leaflet)
library(wesanderson)
library(knitr)

# Load data ----

# utility data
url_malaria <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week1/Lab_files/Data/mal_data_eth_2009_no_dups.csv"
url_landuse <- "https://github.com/HughSt/HughSt.github.io/blob/master/course_materials/week2/Lab_files/ETH_land_use.tif?raw=true"

# Ethiopia admin boundaries
ETH_Adm_2 <- raster::getData("GADM", country="ETH", level = 2)

# Ethiopia malaria data
ETH_malaria_data <- read.csv(url_malaria, header=T)

# Convert data to SPDF
ETH_malaria_data_SPDF <- SpatialPointsDataFrame(
  coords = ETH_malaria_data[,c("longitude", "latitude")],
  data = ETH_malaria_data[,c("examined", "pf_pos", "pf_pr", "rural_urban")],
  proj4string = CRS("+init=epsg:4326"))


# 1 overlay malaria point data with boundaries data ----

# overlay data after reprojecting CRS if needed
crs(ETH_malaria_data_SPDF)
crs(ETH_Adm_2)

ETH_malaria_data_SPDF <- spTransform(ETH_malaria_data_SPDF, crs(ETH_Adm_2))

crs(ETH_malaria_data_SPDF)
crs(ETH_Adm_2)

ETH_Adm_2_per_point <- over(ETH_malaria_data_SPDF, ETH_Adm_2)

# calculate prevalence
Nexamined_per_Adm2 <- tapply(ETH_malaria_data_SPDF$examined,
                             ETH_Adm_2_per_point$NAME_2, sum)

Npositives_per_Adm2 <- tapply(ETH_malaria_data_SPDF$pf_pos,
                              ETH_Adm_2_per_point$NAME_2, sum)

prev_per_Adm2 <- Npositives_per_Adm2 / Nexamined_per_Adm2

prev_per_Adm2_df <- data.frame(NAME_2 = names(prev_per_Adm2),
                               prevalence = prev_per_Adm2,
                               row.names = NULL)

ETH_Adm_2 <- merge(ETH_Adm_2, prev_per_Adm2_df, by = "NAME_2")

prev_per_prov <- ETH_Adm_2 %>%
  as.data.frame() %>%
  select(NAME_2, prevalence) %>%
  na.omit()

# basemap
basemap <- leaflet() %>% addProviderTiles("Esri.WorldGrayCanvas")

# layer the map
colorPal <- colorNumeric(wes_palette("Zissou1")[1:5], ETH_Adm_2$prevalence)
ETH_prev_Adm2_only <- subset(ETH_Adm_2, !is.na(ETH_Adm_2$prevalence))

map_prev <- basemap %>% addPolygons(data=ETH_prev_Adm2_only,
                        color="navy",
                        weight = 1,
                        fillOpacity = 0.2) %>%
  addPolygons(data = ETH_prev_Adm2_only,
                        weight = 1,
                        col=colorPal(ETH_prev_Adm2_only$prevalence),
                        fillOpacity = 0.4,
                        label = round(ETH_prev_Adm2_only$prevalence,5)) %>%
  addLegend(pal = colorPal, values = ETH_prev_Adm2_only$prevalence,
            title = "Prevalence </br> of Malaria")

map_prev


# 2 explore relationship between land class and infection prevalence ----
boxplot <- function(df, x, y, group, xlab, ylab) {
  
  df %>%
    ggplot(aes(x=x, y=y, group=group)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    scale_x_continuous(xlab,
                       labels = as.character(x),
                       breaks = x) +
    ylab(ylab) +
    theme_minimal()
  
}

# load land use data
ETH_land_use <- raster(url_landuse)

# map land use data to malaria data
ETH_malaria_data_SPDF$prevalence <- ETH_malaria_data_SPDF$pf_pos / ETH_malaria_data_SPDF$examined
ETH_malaria_data_SPDF$landuse <- raster::extract(ETH_land_use, ETH_malaria_data_SPDF)

# organize data to get prevalence rate by land use classifications
df_lu_prev <- ETH_malaria_data_SPDF %>%
  as.data.frame() %>%
  select(landuse, prevalence) %>%
  filter(prevalence != 0) %>%
  mutate(landuse = ifelse(landuse <= 30, "A11", "A12"))

avg_lu_prev <- df_lu_prev %>%
  group_by(landuse) %>%
  summarise(avg_prev = mean(prevalence))

# plot 1: Bar plot of average malaria prevalence by land class
bp1 <- avg_lu_prev %>%
  ggplot(aes(x=landuse, y=avg_prev)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Land Use Reclassified") +
  ylab("Infection Prevalence") +
  theme_minimal()
bp1

# plot 2: Stacked bar plot of malaria prevalence by land classification
bp2 <- df_lu_prev %>%
  ggplot(aes(x=landuse, y=prevalence, fill = prevalence)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Land Use Reclassified") +
  ylab("Infection Prevalence") +
  theme_minimal()
bp2

# plot 3: Box and whisker plot of malaria prevalence by land classification
bp3 <- df_lu_prev %>%
  ggplot(aes(x=landuse, y=prevalence, group=landuse)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Land Use Reclassified") +
  ylab("Infection Prevalence") +
  theme_minimal()
bp3 

# 3 get more data and analyze in tandem ----
make_map <- function(popup_label, column, rastercolor) {
  
  # define color ramp for legend
  colorPal <- colorNumeric("RdYlGn", ETH_prev_Adm2_only$prevalence)
  
  ETH_prev_Adm2_only$var <-
    raster::extract(column, ETH_prev_Adm2_only, fun=mean, na.rm=T)
  
  clim_map <- basemap %>% addRasterImage(column,
                                         colors = rastercolor,
                                         opacity = 2) %>%
    addPolygons(data=ETH_prev_Adm2_only,
                color="black",
                weight = 1,
                fillOpacity = 0) %>%
    addLegend(pal = colorPal, values = ETH_prev_Adm2_only$prevalence,
              title = "Malaria Prevalence")
  
  clim_map <- clim_map %>%
    addPolygons(data = ETH_prev_Adm2_only,
                weight = 2,
                col = colorPal(ETH_prev_Adm2_only$prevalence),
                fillOpacity = 0.3,
                label = paste(popup_label, round(ETH_prev_Adm2_only$var, 2)))
  clim_map
}

# load and crop world climate data to ethopia admin 2 boundaries
climate_eth <- getData('worldclim', var='bio', res=10) %>% crop(ETH_prev_Adm2_only)

# Map and plots of Mean of Annual Mean Temp
map_meanannualmeantemp <- make_map("Mean of Annual Mean Temp:", climate_eth$bio1, "Blues")
map_meanannualmeantemp

# Map and plots of Mean Annual Precipitation
map_meanannualprecip <- make_map("Mean Annual Precipitation:", climate_eth$bio12, "Purples")
map_meanannualprecip

ETH_malaria_data_SPDF$bio1 <- raster::extract(climate_eth$bio1,
                                          ETH_malaria_data_SPDF, fun=mean, na.rm=T)

ETH_malaria_data_SPDF$bio12 <- raster::extract(climate_eth$bio12,
                                               ETH_malaria_data_SPDF, fun=mean, na.rm=T)

df_plots <- ETH_malaria_data_SPDF %>%
  as.data.frame() %>%
  select(latitude, longitude, rural_urban, landuse, bio1, bio12, prevalence) %>%
  filter(prevalence != 0) %>%
  mutate(landuse = ifelse(landuse <= 30, "A11", "A12"))

barplot <- function(df, xvar, yvar, fill, xlab, ylab) {
  df %>%
    ggplot(aes(x=xvar, y=yvar, fill = fill)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    xlab(xlab) +
    ylab(ylab) +
    theme_minimal()
}

barplot(df_plots, df_plots$rural_urban, df_plots$bio1, df_plots$prevalence,
        "Rural v. Unknown", "Bio1")

barplot(df_plots, df_plots$rural_urban, df_plots$bio12, df_plots$prevalence,
        "Rural v. Unknown", "Bio12")

scatterplots <- function(df, xvar, yvar, xlab, ylab) {
  
  df_plots %>%
    ggplot(aes(x=xvar, y=yvar)) + 
    geom_point(shape=18, color="blue") +
    geom_smooth(method=lm, linetype="dashed",
                color="darkred", fill="blue") +
    xlab(xlab) +
    ylab(ylab)
  
}

scatterplots(df_plots, df_plots$bio1, df_plots$prevalence,
             ylab = "Infection Prevalence", xlab = "Bio1")

scatterplots(df_plots, df_plots$bio12, df_plots$prevalence,
             ylab = "Infection Prevalence", xlab = "Bio12")


