---
title: 'Analyzing Malaria Prevalence in Ethiopia'
author: "Sri Ramesh"
date: "3/22/2020"
output: html_document
always_allow_html: true
---

The following provides an analysis of Professor Hugh Sturrock's (UCSF) malaria point prevalence data taken from surveys done in Ethiopia in 2009. This data is found on Professor Sturrock's GitHub page:

* https://github.com/HughSt/HughSt.github.io/tree/master/course_materials/week1/Lab_files/Data

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(kableExtra)
library(leaflet)
library(raster)
library(sp)
library(tidyverse)
library(wesanderson)

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

# project data ----
crs(ETH_malaria_data_SPDF)
crs(ETH_Adm_2)

ETH_malaria_data_SPDF <- spTransform(ETH_malaria_data_SPDF, crs(ETH_Adm_2))

crs(ETH_malaria_data_SPDF)
crs(ETH_Adm_2)

ETH_Adm_2_per_point <- over(ETH_malaria_data_SPDF, ETH_Adm_2)

# calculate prevalence ----
Nexamined_per_Adm2 <- tapply(ETH_malaria_data_SPDF$examined,
                             ETH_Adm_2_per_point$NAME_2, sum)

Npositives_per_Adm2 <- tapply(ETH_malaria_data_SPDF$pf_pos,
                              ETH_Adm_2_per_point$NAME_2, sum)

prev_per_Adm2 <- Npositives_per_Adm2 / Nexamined_per_Adm2

prev_per_Adm2_df <- data.frame(NAME_2 = names(prev_per_Adm2),
                               prevalence = prev_per_Adm2,
                               row.names = NULL)

ETH_Adm_2 <- merge(ETH_Adm_2, prev_per_Adm2_df, by = "NAME_2")

ETH_malaria_data_SPDF$prevalence <- ETH_malaria_data_SPDF$pf_pos / ETH_malaria_data_SPDF$examined

```

## Q1: Table showing prevalence of infection at admin 2 level

This table provides infection prevalence for each admin 2 level in Ethiopia.

```{r q1table, echo=T, warning=F, message=F}
# generate table using admin 2 level Ethiopia data and calculated prevalence value ----
prev_per_prov <- ETH_Adm_2 %>% as.data.frame() %>% select(NAME_2, prevalence) %>% na.omit()

prev_per_prov %>% 
  kable(col.names = c("Admin 2 Level", "Prevalence of Malaria Infection"), row.names = F) %>% kable_styling()

```

## Q1: Map showing prevalence of infection at admin 2 level

This leaflet map shows infection prevalence at the admin 2 level using the columns *NAME_2*. **NA values were not included on the map.**

```{r q1map, echo=T, warning=F, message=F}
# generate basemap ----
basemap <- leaflet() %>% addProviderTiles("Esri.WorldGrayCanvas")

# define color ramp to use infection prevalence values as the legend ----
colorPal <- colorNumeric(wes_palette("Zissou1")[1:5], ETH_Adm_2$prevalence)

# drop admin 2 level boundaries that have no prevalence data (will use this throughout) ----
ETH_prev_Adm2_only <- subset(ETH_Adm_2, !is.na(ETH_Adm_2$prevalence))

# add admin 2 level boundaries layer ----
map_prev <- basemap %>% addPolygons(data=ETH_prev_Adm2_only, color="navy", weight = 1,
                        fillOpacity = 0.2) %>%
  
  # add chloropleth of infection prevalence
  addPolygons(data = ETH_prev_Adm2_only, weight = 1, col=colorPal(ETH_prev_Adm2_only$prevalence),
                        fillOpacity = 0.4, label = round(ETH_prev_Adm2_only$prevalence,5)) %>%
  
  # add legend for infection prevalence
  addLegend(pal = colorPal, values = ETH_prev_Adm2_only$prevalence,
            title = "Infection </br> Prevalence")

map_prev
```

## Q2: Land class is reclassified using the LCCS Entry column information

This is how I reclassified land use using the LCCS Entry column information:

```{r q2reclassify, echo=T, warning=F, message=F}

# load land use data ----
ETH_land_use <- raster(url_landuse)

# extract land use data to malaria data ----
ETH_malaria_data_SPDF$landuse <- raster::extract(ETH_land_use, ETH_malaria_data_SPDF)

# organize data into 2 tables to get prevalence rate by land use classifications ----
df_lu_prev <- ETH_malaria_data_SPDF %>%  as.data.frame() %>% select(landuse, prevalence) %>%
  filter(prevalence != 0) %>% mutate(landuse = ifelse(landuse <= 30, "A11", "A12")) %>%
  arrange(landuse, prevalence)

avg_lu_prev <- df_lu_prev %>% group_by(landuse) %>% summarise(avg_prev = mean(prevalence))

```

This table illustrates the *average* infection prevalence by the land classes found in the LCCS Entry column:

```{r q2table2, echo=T, warning=F, message=F}
avg_lu_prev %>% kable(col.names = c("LCCS Entry Data", "Avg. Infection Prevalence")) %>% kable_styling()

```

## Q2: Generate exploratory plots to look at the relationship between land class and infection prevalence

The following contains 3 plots to explore the relationship between the reclassified land classes and infection prevalence:

* a bar plot,
* a stacked barplot, and
* a boxplot.

**Plot 1: Bar plot of average infection prevalence by land class.** This bar plot shows each land class' *average* malaria infection prevalence.

```{r q2plot1, echo=T, warning=F, message=F}
bp1 <- avg_lu_prev %>% ggplot(aes(x=landuse, y=avg_prev)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Land Use Reclassified") + ylab("Infection Prevalence") + 
  ggtitle("Bar plot of average infection prevalence by land class") + theme_minimal()

bp1

```

**Plot 2: Stacked bar plot of infection prevalence by land class.** This stacked bar plot shows the distribution of infection prevalence by each land class. We see that the class A11 has a higher average infection prevalence than the class A12.

```{r q2plot2, echo=T, warning=F, message=F}
bp2 <- df_lu_prev %>% ggplot(aes(x=landuse, y=prevalence, fill = prevalence)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Land Use Reclassified") + ylab("Infection Prevalence") + 
  ggtitle("Stacked bar plot of infection prevalence by land class") + theme_minimal()

bp2
```

**Plot 3: Box and whisker plot of infection prevalence by land class.** This box and whisker plot shows five-number summaries of infection prevalence by land class. We see that the class A11 also has more variety in the infection prevalence values than the class A12. We see that the class A11 has a wider range than the infection prevalence than the class A12.

```{r q2plot3, echo=T, warning=F, message=F}
boxplot1 <- df_lu_prev %>% ggplot(aes(x=landuse, y=prevalence, group=landuse)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Land Use Reclassified") + ylab("Infection Prevalence") + 
  ggtitle("Box plot of infection prevalence by land class") + theme_minimal()

boxplot1 
```

## Q3: 1st variable to relate to the survey data: **Bio1, World Climate Data**

The following relates Bio1 which is the **Mean of Annual Mean Temperature** at the Admin 2 level in Ethiopia to the survey data in three ways:

* a leaflet map,
* a stacked bar plot, and
* a scatterplot.

Here are the functions I used to make each one:

```{r q3functions, echo=T, warning=F, message=F}
make_map <- function(popup_label, column, rastercolor) {

  # define color palette  
  colorPal <- colorNumeric("RdYlGn", ETH_prev_Adm2_only$prevalence)
  
  # extract values of dataset to append column to ETH_Adm_2 boundaries data that does not include NAs
  ETH_prev_Adm2_only$var <- raster::extract(column, ETH_prev_Adm2_only, fun=mean, na.rm=T)
  
  # add bioclimatic data raster image as first layer
  clim_map <- basemap %>% addRasterImage(column, colors = rastercolor, opacity = 2) %>%
    
    # add boundaries layer
    addPolygons(data=ETH_prev_Adm2_only, color="black", weight = 1, fillOpacity = 0)
  
    # add legend for infection prevalence
  clim_map <- clim_map %>% addLegend(pal = colorPal, values = ETH_prev_Adm2_only$prevalence, 
                                     title = "Infection Prevalence") %>%
    
    # add infection prevalence chloropleth layer
    addPolygons(data = ETH_prev_Adm2_only, weight = 2, col = colorPal(ETH_prev_Adm2_only$prevalence),
                fillOpacity = 0.3, label = paste(popup_label, round(ETH_prev_Adm2_only$var, 2)))
  clim_map
}

scatterplot <- function(df, xvar, yvar, xlab, ylab, ggtitle) {
  
  df_plots %>%
    ggplot(aes(x=xvar, y=yvar)) +  geom_point(shape=18, color="blue") +
    geom_smooth(method=lm, linetype="dashed", color="darkred", fill="blue") + xlab(xlab) + ylab(ylab) + ggtitle(ggtitle) + theme_minimal()
  
}

barplot <- function(df, xvar, yvar, fill, xlab, ylab, ggtitle) {
  
  df %>% ggplot(aes(x=xvar, y=yvar, fill = fill)) + geom_bar(stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + xlab(xlab) + ylab(ylab) +
    ggtitle(ggtitle) + theme_minimal()
  
}

```

**Map of Mean of Annual Mean Temperature vs. Infection Prevalence at Admin 2 levels.** The map seems to show no direct relationship between these two variables (Bio1 and infection prevalence). **Note: If you hover over the map, you will see the Bio1 value. You can compare this to the underlying chloropleth of the infection prevalence data.**

```{r q3layer1map, echo=T, warning=F, message=F}
# load and crop world climate data to ethopia admin 2 boundaries ----
climate_eth <- getData('worldclim', var='bio', res=10) %>% raster::crop(ETH_prev_Adm2_only)

# extract bioclimatic data to cropped ethiopia raster data ----
ETH_malaria_data_SPDF$bio1 <- raster::extract(climate_eth$bio1, ETH_malaria_data_SPDF, fun=mean, na.rm=T)

ETH_malaria_data_SPDF$bio12 <- raster::extract(climate_eth$bio12, ETH_malaria_data_SPDF, fun=mean, na.rm=T)

# convert data to data.frame for easier use ----
df_plots <- ETH_malaria_data_SPDF %>% as.data.frame() %>%
  select(latitude, longitude, rural_urban, landuse, bio1, bio12, prevalence) %>%
  filter(prevalence != 0) %>% mutate(landuse = ifelse(landuse <= 30, "A11", "A12"))

# layer map of bio1 variable ----
map_meanannualmeantemp <- make_map("Mean of Annual Mean Temp:", climate_eth$bio1, "Blues")
map_meanannualmeantemp

```

**Scatterplot of Mean of Annual Mean Temperature and Infection Prevalence.** With a wide confidence interval, we can see that infection prevalence does not seem to increase with the average of annual mean temperature.

```{r q3layer1scatter, echo=T, warning=F, message=F}
scatterplot(df = df_plots, xvar = df_plots$bio1, yvar = df_plots$prevalence,
            xlab = "Bio1", ylab = "Infection Prevalence",
            ggtitle = "Scatterplot of Bio1 and Infection Prevalence")

```

**Barplot of Mean of Annual Mean Temperature, Infection Prevalence, and Rural v. Unknown Sites.** Infection prevalence seems to be higher at rural sites vs. sites that were unknown (rural or urban). These rural sites also seemed to have more average annual mean temperature values. 

```{r q3layer1barplot, echo=T, warning=F, message=F}
barplot(df = df_plots, xvar = df_plots$rural_urban, yvar = df_plots$bio1, 
        fill = df_plots$prevalence, xlab = "Rural v. Unknown", ylab = "Bio1",
        ggtitle = "Barplot of Bio1, Infection Prevalence, and Rural v. Unknown Sites")
```

## Q3: 2nd variable to relate to the survey data: **Bio12, World Climate Data**

The following relates Bio12 which is the **Mean of Annual Precipitation** at the Admin 2 level in Ethiopia. 

**Map of Mean of Annual Precipitation and Infection Prevalence.** This map shows higher average annual precipitation in the northwestern portion of this region of Ethiopia, but not necessarily any kind of correlation between those values and infection prevalence.

```{r q3layer2map, echo=T, warning=F, message=F}
map_meanannualprecip <- make_map("Mean Annual Precipitation:", climate_eth$bio12, "Purples")
map_meanannualprecip
```

**Scatterplot of Mean of Annual Precipitation and Infection Prevalence.** With a wide range of confidence, we might say that there is a positive correlation between the mean of annual precipitation (Bio12) and infection prevalence.

````{r q3layer2scatter, echo=T, warning=F, message=F}
scatterplot(df = df_plots, xvar = df_plots$bio12, yvar = df_plots$prevalence,
            xlab = "Bio12", ylab = "Infection Prevalence",
            ggtitle="Scatterplot of Bio12 and Infection Prevalence")

```

**Barplot of Mean of Annual Precipitation, Infection Prevalence, and Rural v. Unknown Sites.** Similar to the Bio1 barplot, this barplot shows higher infection prevalence in rural sites as opposed to sites that were unknown to be either rual or urban. It also shows higher average annual precipitation at the rural sites.

````{r q3layer2barplot, echo=T, warning=F, message=F}
barplot(df = df_plots, xvar = df_plots$rural_urban, yvar = df_plots$bio12,
        fill = df_plots$prevalence, xlab = "Rural v. Unknown", ylab = "Bio12",
        ggtitle="Barplot of Bio12, Infection Prevalence, and Rural v. Unknown Sites")
```
