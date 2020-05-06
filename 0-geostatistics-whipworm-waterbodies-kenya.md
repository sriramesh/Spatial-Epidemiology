---
title: "Geostatistical analysis of the association between proximity to water bodies and Trichuris Trichiura prevalence in western Kenya"
author: "Sri Ramesh"
date: "5/1/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# load libraries ----
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
library(spdep)
library(pgirmess)
library(sp)
library(gtools)
library(car)

# global options ----
options(stringsAsFactors = F) # import quality
par(mar=c(1,1,1,1)) # ensure plot margin sizes are always large enough

# load soil-transmitted helminths data ----
STH_kenya <- "Kenya_STH.csv" %>% read.csv(header = T) %>%
  select(full_name_paper, lat, long, tri_prev, tri_np)

STH_kenya <- STH_kenya %>%
  mutate(tri_prev = tri_prev/100,
         tri_examined = tri_np/tri_prev,
         tri_nneg = tri_examined - tri_np,
         tri_prev_logit = car::logit(tri_prev) + 0.0001) %>%
  drop_na(tri_prev_logit)

STH_kenya <- SpatialPointsDataFrame(coords = STH_kenya %>% select(long,lat),
                                    data = STH_kenya %>%
                                      select(lat, long, tri_prev, tri_prev_logit))

STH_kenya_df <- read.csv("STH_kenya_covariates.csv")
```

## Introduction

In the following, I will use geostatistics to explore Professor Hugh Sturrock's data on soil-transmitted helminths prevalence in western Kenya with a focus on *T. Trichiura*, or the human whipworm. *T. Trichiura* is one of three well-documented soil-transmitted helminths common to Asia and Sub-Saharan Africa. In human beings, it causes a neglected tropical disease known as *trichuriasis*. Understanding the geographic distribution of soil-transmitted helminths like whipworm is important for developing anthelminthic treatment strategies in tropical countries.

A sizable body of evidence out of several tropical countries points to a potential correlation between altitude and *T. Trichiura*, as well as annual precipitation and *T. Trichiura*. Recent studies including Chaiyos et. al.'s 2018 study out of Thailand, as well as older studies like Appleton et. al.'s 1995 study in South Africa, find that prevalence of *T. trichiura* was significantly correlated with annual precipitation and annual rainfall. Chammartin et. al.'s 2013 study in Bolivia finds that altitude has a negative effect on *T. Trichiura*.

A smaller set of studies in East Africa point to a correlation between *T. Trichiura* prevalence and proximity to Lake Victoria. Handzel et. al.'s 2003 study out of western Kenya finds that closer proximity to Lake Victoria was associated with increased *T. trichiura* infection rates among school-aged children. A similar 2015 study in Tanzania by Siza et. al. finds high whipworm infection rates among residents of Tanzania's Lake Victoria Basin.

Based on this short literature review, I investigated the following two hypotheses using geostatistical methods:

1. Closer proximity to waterbodies (i.e. lakes, rivers, and canals) is associated with higher  prevalence of *T. Trichiura* in western Kenya. In other words, distance to lakes, distance to rivers/canals, and whipworm prevalence are negatively correlated.
2. This association is tempered by altitude and annual precipitation as well as transit accessibility. In other words, lower altitudes and lower annual precipitation are associated with lower whipworm prevalence, and transit accessibility, which will be measured by a closer proximity to roads and railroads, is negatively correlated with whipworm prevalence.

## Methods

In the following, I will discuss my methods which included pre-processing the whipworm prevalence data, cluster analyses, processing model covariates, spatial and non-spatial regression analyses, model cross-validation, and risk mapping.

### Pre-processing whipworm prevalence data

My first step to examine my hypothesis was to pre-process the whipworm prevalence data in preparation for a cluster analysis. The "raw" whipworm prevalence data had a sample size of 68 prevalence points. My data pre-processing steps included generating a histogram, which revealed that the data was not normally distributed. I conducted a logit transformation on this data to calculate log whipworm prevalence and dropped any NA values in the log whipworm prevalence variable. This reduced my sample size to 61 observations.

The logit transformation did not correct the skew in the data, as shown in the plot below.

```{r, echo=F, warning=F, message=F}
STH_kenya_df %>% ggplot(aes(x=tri_prev_logit)) + 
  geom_histogram(aes(y=..density..), colour="darkgrey", fill="lightblue")+
  geom_density(alpha=.2, fill="#FF6666") +
  xlab("T. trichiura prevalence") + ylab("Density") +
  ggtitle("Distribution of T. trichiura prevalence - logit transform (still not quite normal)")

```

The limitations of using prevalence points that are not normally distributed is detailed in the conclusion.

### Global and local cluster analysis

After pre-processing the whipworm prevalence data, I examined the data for instances of global and local clustering.

#### Global clustering

To examine instances of global clustering, I generated correlograms to determine the highest distance classes at which statistically significant positive spatial autocorrelation could be discerned. I then used a Monte Carlo simulation of the global Moran's I test to determine the Moran's I test statistic.

The neighbor structure used within my Moran's I tests followed a k-nearest neighbors approach with the parameters **k = 1** and row-standardized weights. With these parameters, only immediately adjacent points were defined as “neighbors,” and each neighbor linkage was assigned a spatial weight value equivalent to the 1 divided by the number of neighbors.

Correlograms calculated with 15 and 20 distance classes initially revealed that spatial autocorrelation exists up to about **0.17 decimal degrees** and I conducted a Monte Carlo simulation of the Moran's I test to confirm the value of the test statistic, as shown below.

```{r, echo=F, warning=F}
xy <- STH_kenya_df %>%
  select(long, lat) %>% bind_cols() %>% as.matrix()
coords <- xy %>% coordinates()
IDs <- coords %>% as.data.frame() %>% row.names()
Neigh_nb <- coords %>% knearneigh(k=1, longlat = TRUE) %>% knn2nb(row.names=IDs)
dsts <- Neigh_nb %>% nbdists(coords) %>% unlist()

# find maximum distance used to provide at least one neighbor to each point
max_1nn <- dsts %>% max()

# define distance-based neighbors within maximum distance:
Neigh_kd1 <- coords %>% dnearneigh(d1=0, d2=max_1nn, row.names=IDs) # when k=1 (immediate neighbors)
# define distance-based neighbors within 2 times the maximum distance:
Neigh_kd2 <- coords %>% dnearneigh(d1=0, d2=2*max_1nn, row.names=IDs) # when k=1
# assign weights to neighbor linkages
weights <- Neigh_kd1 %>% nb2listw(style="W")

# Monte Carlo simulation for Moran's I test
set.seed(1234)
bperm <- STH_kenya$tri_prev_logit %>% moran.mc(listw=weights,nsim=999)
bperm

# Plot of Moran's I
par(mfrow=c(1,1), mar= c(5, 4, 4, 2))
hist(bperm$res, freq=T, breaks=20, xlab="Moran's I statistics",
     main = "Distribution of Moran's I statistics over 999 simulations")
abline(v=0.5979, col="red")

```

The Monte Carlo simulation reveals that to a high degree of statistical significance, there is positive global clustering in our whipworm prevalence data with a Moran's I statistic of **0.597.**

However, running the checks for neighbor symmetry below, we find that neighbors are not symmetrical:

```{r, echo=F, warning=F}
# list of neighbor structures
nb_1 <- list(d1=Neigh_kd1, d2=Neigh_kd2) #when k=1

# Check for symmetry
nb_1 %>% sapply(function(x) is.symmetric.nb(x, verbose=F, force=T))
nb_1 %>% sapply(function(x) n.comp.nb(x)$nc)
```

This limitation in assigning neighbors and spatial weights to identify clusters is explored in the conclusion.

#### Local clustering

To examine instances of local clustering, I decomposed the global Moran's I statistic into its local parts by generating windows of Local Instances of Spatial Autocorrelation (LISAs) around each whipworm prevalence point. The LISAs for each point show several instances of local clustering, as shown below.

```{r, echo=F, warning=F}
# First calculate the local Moran's I around each point based on row standardized weights where k = 1
I <- STH_kenya$tri_prev_logit %>% localmoran(weights) # "spdep" package

# Print 'LISA' for each point
Coef <- I %>% as.data.frame(row.names=row.names(coords), check.names=FALSE) %>% printCoefmat()


```

However, the plot of observed vs. expected whipworm prevalence (log) shown below indicates that there are not many outliers in the data, only 7 outliers. This means that most of the data does not deviate from expectation, where expectation is defined as the weighted mean of each point’s neighbors. 

```{r, echo=F}
nci <- STH_kenya$tri_prev_logit %>% moran.plot(listw=weights,
                                               main = "Observed vs. Expected Log prevalence",
                                               xlab="Log prevalence",
                                               ylab="Spatially lagged log prev",
                                               labels=T, pch=16, col="grey")

# find which points are statistically significant outliers
infl <- apply(nci$is.inf,1,any) 
num_outliers <- sum(infl %in% T)

```

### Processing model covariates

Following a cluster analysis, I used geostatistics to discern the correlation between whipworm prevalence and distance to waterbodies, distance to transit, altitude, and annual precipitation. I will briefly describe 1) the repository from which each covariate was pulled, and 2) transformations applied on covariates.

#### Covariate data sources

The following describes each covariate, the data source it was pulled from, and relevant sample sizes.

* **Annual precipitation data (RasterLayer):** The annual precipitation data was available on the WorldClimate database. I pulled the variable *Bio12* which corresponds to annual precipitation at a resolution of 2.5 meters.
* **Altitude data (RasterLayer):** The altitude data was aggregated from NASA's Shuttle Radar Topography Mission (SRTM) in February 2000, which obtained 90 meter (3 arc-second) resolution between -60 and 60 latitude. This data provides altitude information at 90m resolution for many parts of the world. I pulled this data corresponding to the Kenya GADM administrative region. 
* **Inland lakes (Polygons):** The inland lakes data provides the locations of lakes across Kenya. This data comes from the Digital Chart of the World (DCW) inland waterbodies database. The DCW was commissioned by the United States Defense Mapping Agency around 1992 for global scientific research purposes and to support military analysis. I only downloaded a subset of the DCW inland waterbodies database corresponding to Kenya's national borders. This subset had a sample size of 255 inland lakes.
* **Inland rivers and canals (Lines):** The inland rivers and canals data provides the locations of inland rivers and canals in Kenya. This data comes from the DCW inland waterbodies database. I only downloaded a subset of this database which corresponds to Kenya's national borders. This subset had a sample size of 2569 inland lakes.
* **Roads (Lines):** The roads data provides the locations of roads in Kenya. This data comes from the DCW roads database. I only downloaded a subset of this database which corresponds to Kenya's national borders. This subset had a sample size of 2097 roads.
* **Railroads (Lines):** The railroads data provides the locations of railroads in Kenya. This data comes from the DCW railroads database. I only downloaded a subset of this database which corresponds to Kenya's national borders. This subset had a sample size of 42 railroads.

Here is a visualization of the vector-based covariates listed above: 

```{r, include=F}
# load covariates ----
KEN_Adm_0 <- raster::getData("GADM", country="KEN", level=0)
alt <- raster::getData('alt', country = 'KEN', mask = F)
bioclim <- raster::getData('worldclim', var='bio', res=2.5) # Bio12n
roads <- readOGR("KEN_rds", "KEN_roads") # roads
railroads <- readOGR("KEN_rrd", "KEN_rails") # railroads
lakes <- readOGR("KEN_wat", "KEN_water_areas_dcw") # inland lakes
canals <- readOGR("KEN_wat", "KEN_water_lines_dcw") #inland canals and rivers

# set crs ----
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj4string(STH_kenya) <- wgs84 %>% CRS()
proj4string(roads) <- wgs84 %>% CRS()
proj4string(railroads) <- wgs84 %>% CRS()
proj4string(lakes) <- wgs84 %>% CRS()
proj4string(canals) <- wgs84 %>% CRS()

roads_map <- roads %>% crop(STH_kenya %>% extent())
railroads_map <- railroads %>% crop(STH_kenya %>% extent())
lakes_map <- lakes %>% crop(STH_kenya %>% extent())
canals_map <- canals %>% crop(STH_kenya %>% extent())

pal = colorNumeric("Reds", STH_kenya_df$tri_prev)

map <- leaflet() %>% addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(data = STH_kenya_df, lng = ~long, lat = ~lat,
                   fillOpacity=1,
                   fillColor= ~pal(tri_prev),
                   radius=~5, stroke=TRUE, weight=.1) %>%
  addPolylines(data = roads_map, weight=1.5, color = "dimgrey") %>%
  addPolylines(data = railroads_map, weight=2, color = "grey") %>%
  addPolylines(data = canals_map, weight=2, color = "dodgerblue", opacity = 0.5) %>%
  addPolygons(data = lakes_map, weight=1.5, color="dodgerblue", opacity = .5) %>%
  addLegend(pal = pal, values = STH_kenya_df$tri_prev, title = "Whipworm Prevalence")

```

```{r, echo=F, message=F}
map
```

In the map above, **points in shades of red** are whipworm prevalence points, **blue lines and polygons** are inland waterbodies and **grey lines** are roads and railroads. For illustration purposes, I have clipped these features to the bounding box of the whipworm prevalence data, which is in and around western Kenya.

After loading the covariates, I rasterized all vector data, resampled the rasters, and calculated distances between each prevalence point and the nearest lake, rivers/canal, road, and railroad.

#### Rasterizing all vector data

I rasterized my four covariates containing vector data from the DCW database:

* lakes shapefile
* canals/lakes shapefile
* roads shapefile
* railroads shapefile

First, I re-projected all of my vector data into the WGS84 coordinate reference system:

```{r, echo=F}
KEN_Adm_0 <- KEN_Adm_0 %>% spTransform(CRS(wgs84))
proj4string(STH_kenya) <- wgs84 %>% CRS()

roads <- roads %>% spTransform(CRS(wgs84))
canals <- canals %>% spTransform(CRS(wgs84))
railroads <- railroads %>% spTransform(CRS(wgs84))
lakes <- lakes %>% spTransform(CRS(wgs84))
```

Next, I rasterized each covariate with a step-wise process. Here is the process used for rasterizing the roads shapefile, which contains spatial lines:

```{r, echo=F, warning=F}
# rasterize roads ----
rs <- raster(extent(roads), crs=projection(wgs84))
rs[] <- 1:ncell(rs)

rsp <- rs %>% rasterToPolygons()
rp <- roads %>% raster::intersect(rsp)
rp$length <- rgeos::gLength(rp, byid=TRUE) / 1000
x <- tapply(rp$length, rp$layer, sum)

roads_r <- rs %>% raster()
roads_r[as.integer(names(x))] <- x
```

A plot of the newly produced raster (in this case, roads) revealed that rasterization worked:

```{r, include=F}
raster_colorPal <- colorNumeric(palette = terrain.colors(64), reverse = T, domain = values(roads_r), na.color = NA)

basemap <- leaflet() %>% addProviderTiles("CartoDB.Positron")

map_kenya <- basemap %>%
  addRasterImage(x = roads_r, color = raster_colorPal, opacity = 0.5) %>%
    addPolylines(data = roads, weight=1.5, color = "dimgrey") %>%
    addLegend(title = "Density of roads, Kenya, 1992 DCW",
            values = values(roads_r), pal = raster_colorPal)
```

```{r, echo=F}
map_kenya
```

```{r, echo=F, warning=F}
# rasterize canals ----
rs <- raster(extent(canals), crs=projection(wgs84))
rs[] <- 1:ncell(rs)

rsp <- rs %>% rasterToPolygons()
rp <- roads %>% raster::intersect(rsp)
rp$length <- rgeos::gLength(rp, byid=TRUE) / 1000
x <- tapply(rp$length, rp$layer, sum)

canals_r <- rs %>% raster()
canals_r[as.integer(names(x))] <- x

# canals_r %>% plot()
# canals %>% lines()
# KEN_Adm_0 %>% lines()

# rasterize railroads ----
rs <- raster(extent(railroads), crs=projection(wgs84))
rs[] <- 1:ncell(rs)

rsp <- rs %>% rasterToPolygons()
rp <- roads %>% raster::intersect(rsp)
rp$length <- rgeos::gLength(rp, byid=TRUE) / 1000
x <- tapply(rp$length, rp$layer, sum)

railroads_r <- rs %>% raster()
railroads_r[as.integer(names(x))] <- x

# railroads_r %>% plot()
# railroads %>% lines()
# KEN_Adm_0 %>% lines()
```

To rasterize my lakes vector data, I applied a method for rasterizing spatial polygons as shown below:

```{r, echo=T}
# rasterize lakes ----
ext <-  extent(33.90959, 41.92622, -4.720417, 5.061166)
xy <- ext %>% bbox() %>% as.matrix() %>% apply(1, diff) %>% abs()
n <- 5
r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)
lakes_r <- lakes %>% rasterize(r, fun='first')
```

This also successfully rasterized my lakes data:

```{r, include=F}
raster_colorPal <- colorNumeric(palette = hcl.colors(64), reverse = T, domain = values(lakes_r), na.color = NA)

basemap <- leaflet() %>% addProviderTiles("CartoDB.Positron")

map_kenya <- basemap %>%
  addRasterImage(x = lakes_r, color = raster_colorPal, opacity = 0.5) %>%
    addPolygons(data = lakes, weight=1.5, color = "lightblue") %>%
    addLegend(title = "Density of lakes, Kenya, 1992 DCW",
            values = values(lakes_r), pal = raster_colorPal)
```

```{r, echo=F}
map_kenya
```

#### Resampling rasters

After rasterizing my vector data, I cropped each of the 6 raster layers to match the extent of the Kenya Admin 1 Level per the GADM, then resampled all 6 of my raster layers to the lowest resolution raster of the set: the altitude raster layer.

```{r, echo=F}
KEN_alt <- alt %>% crop(KEN_Adm_0 %>% extent())
KEN_bioclim <- bioclim[[12]] %>% crop(KEN_Adm_0 %>% extent())
KEN_canals_r <- canals_r %>% crop(KEN_Adm_0 %>% extent())
KEN_roads_r <- roads_r %>% crop(KEN_Adm_0 %>% extent())
KEN_railroads_r <- railroads_r %>% crop(KEN_Adm_0 %>% extent())
KEN_lakes_r <- lakes_r %>% crop(KEN_Adm_0 %>% extent())

KEN_bioclim_resamp <- KEN_bioclim %>% resample(KEN_alt, method = "ngb")
KEN_canals_r_resamp <- KEN_canals_r %>% resample(KEN_alt, method = "ngb")
KEN_roads_r_resamp <- KEN_roads_r %>% resample(KEN_alt, method = "ngb")
KEN_railroads_r_resamp <- KEN_railroads_r %>% resample(KEN_alt, method = "ngb")
KEN_lakes_r_resamp <- KEN_lakes_r %>% resample(KEN_alt, method = "ngb")

# takes quite a while
# distance_to_lakes_r <- KEN_lakes_r_resamp %>% distance()
# distance_to_canals_r <- KEN_canals_r_resamp %>% distance()
# distance_to_roads_r <- KEN_roads_r_resamp %>% distance()
# distance_to_railroads_r <- KEN_railroads_r_resamp %>% distance()
```

Because the annual precipitation data is continuous, I resampled using a nearest-neighbor approach, which assigns the value of the nearest cell to cells in the prediction grid. The resampling procedure worked.

We can see that the extent and resolution of all of raster layers now match up, and are ready for calculating distances:

```{r, echo=T}
KEN_alt
KEN_bioclim_resamp
KEN_canals_r_resamp
KEN_roads_r_resamp
KEN_railroads_r_resamp
KEN_lakes_r_resamp
```

#### Calculating distances

I calculated distance between prevalence points to each feature using the `distance` command. This command served to generate raster layers of nothing but distances. Finally, I extracted values from the distance rasters as columns into my data.frame.

```{r, echo=F}
# STH_kenya_df <- STH_kenya %>% as.data.frame() %>% select(-long.1, -lat.1)
# 
# STH_kenya_df <- STH_kenya_df %>%
#   mutate(alt = raster::extract(KEN_alt, STH_kenya_df %>% select(long, lat)),
#          bio12 = raster::extract(KEN_bioclim_resamp, STH_kenya_df %>% select(long, lat)),
#          lakes = raster::extract(distance_to_lakes_r, STH_kenya_df %>% select(long, lat)),
#          canals = raster::extract(distance_to_canals_r, STH_kenya_df %>% select(long, lat)),
#          roads = raster::extract(distance_to_roads_r, STH_kenya_df %>% select(long, lat)),
#          railroads = raster::extract(distance_to_railroads_r, STH_kenya_df %>% select(long, lat)))
```

The extraction revealed two immediate problems with data availability: 1) I could not run regressions using distance to canals/lakes because all extractions had a value of 0:

```{r, echo=F}
STH_kenya_df %>% select(canals) %>% summary()
```

And 2) I could not run regressions using distance to roads because only the following two values were non-zero extractions:

```{r, echo=F}
STH_kenya_df %>% select(roads) %>% head(2)
```

### Spatial and non-spatial regression analyses

I specified non-spatial generalized linear models (GLMs) with various combinations of the 4 covariates available (*bio12, altitude, distance to lakes,* and *distance to railroads*). All models specified included binomial outputs. Importantly, I determined that models that included *distance to roads* resulted in a negative coefficient. I dropped this variable from the models. 

I used two criteria to determine model fit:
* lowest residual deviances
* lowest cross-validated mean squared error.

The model that fit this criteria was the 3-covariate non-spatial GLM with *altitude, annual precipitation, and distance to lakes:*

```{r, echo = F}
glm_mod_1 <- glm(cbind(tri_np, round(tri_nneg)) ~
                   alt + bio12 + lakes,
                 data=STH_kenya_df, family=binomial())
glm_mod_1 %>% summary()
```

I also used the deviance residuals evaluate model robustness. The third quartile (3Q) of the deviance residuals has a value of 0.91. This means that the first 75 percent of our model observations are within an acceptable distance of our predicted values, with "acceptable" being defined as a value less than 2.

After choosing the 4-covariate GLM, I ran a few diagnostics to determine the extent of residual spatial autocorrelation. I first examined a residual plot for any linear patterns:

```{r, echo = F}
residual_plot <- ggplot() + geom_point(aes(glm_mod_1$fitted, STH_kenya_df$tri_prev_logit)) +
  xlab("Expected log prevalence values") +
  ylab("Observed log prevalence values") +
  ggtitle("Observed vs. expected log prevalence values, GLM")
residual_plot
```

The plot illustrates a lack of linearity. This indicates that these covariates, while statistically significant, still do not explain much of the spatial pattern underlying the whipworm prevalence data.

I then used a correlogram to examine clustering effects. I took note of the highest distance class at which statistically significant residual clustering can be detected.

```{r, echo=F}
nbc <- 10
cor_r <- pgirmess::correlog(coords=STH_kenya_df %>% select(long, lat),
                            z=glm_mod_1$residuals,
                            method="Moran", nbclass=nbc)
cor_r 

cor <- cor_r %>% as.data.frame
cor$variable <- "residuals_glm"

cor %>% subset(variable %in% "residuals_glm") %>% ggplot(aes(dist.class, coef)) + 
  geom_hline(yintercept = 0, col="grey") +
  geom_line(col="steelblue") + 
  geom_point(col="steelblue") +
  xlab("distance") + 
  ylab("Moran's coefficient")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Correlogram (GLM)")

```

It seems that the Moran's I statistic is not statistically significant at any of the distance classes. This indicates that we have **no** residual spatial autocorrelation. Our non-spatial GLM has independent and identically distributed residuals (iid), satisfying a key assumption of linear regression modeling.

However, I ran a spatial regression regardless to see how the inclusion of a spatially-correlated fixed effect might improve model fit.

I ran a Matern model with binomial outputs:

```{r, echo=F}
matern <- spaMM::fitme(cbind(tri_np, round(tri_nneg)) ~
                         alt + bio12 + lakes +
                         Matern(1|lat+long),
                       data=STH_kenya_df, family=binomial())
matern %>% summary()
```

I calculated 95% confidence intervals on each covariate and expressed them as odds ratios, as follows:

```{r, include=F}
coefs <- as.data.frame(summary(matern)$beta_table)
row <- row.names(coefs) %in% c('alt')
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, 'Cond. SE']
or_ci_alt <- c(lower, upper) %>% exp()

row <- row.names(coefs) %in% c('bio12')
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, 'Cond. SE']
or_ci_bio12 <- c(lower, upper) %>% exp()

row <- row.names(coefs) %in% c('lakes')
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, 'Cond. SE']
or_ci_dist2lakes <- c(lower, upper) %>% exp()

```

* Altitude: **`r or_ci_alt`** 
* Annual Precipitation: **`r or_ci_bio12`** 
* Distance to lakes: **`r or_ci_dist2lakes`** 

Next, I checked the residual plot and correlogram for the spatial regression to determine the effect of including the spatially-correlated fixed effect.

Residual plot:

```{r, echo=F}
plot <- ggplot() + geom_point(aes(matern$fv, STH_kenya_df$tri_prev_logit)) +
  xlab("Expected log prevalence values") +
  ylab("Observed log prevalence values") +
  ggtitle("Observed vs. expected log prevalence values, Matern")
plot
```

Correlogram:

```{r, echo=F}
nbc <- 10
cor_r <- pgirmess::correlog(coords = STH_kenya_df %>% select(long, lat),
                            z = residuals(matern),
                            method="Moran", nbclass=nbc)
cor_r

cor <- cor_r %>% as.data.frame
cor$variable <- "residuals_glm" 

cor %>% subset(variable %in% "residuals_glm") %>% ggplot(aes(dist.class, coef)) + 
  geom_hline(yintercept = 0, col="grey") +
  geom_line(col="steelblue") + 
  geom_point(col="steelblue") +
  xlab("distance") + 
  ylab("Moran's coefficient")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Correlogram (Matern)")
```

The spatial regression residual plot indicates linearity between observed and expected values. This means that the inclusion of the spatially-correlated fixed effect improved our model fit. 

### Model cross-validation

```{r, include=F}
folds_list <- caret::createFolds(STH_kenya_df$tri_prev)
folds_list[[1]]

cross_validated_prediction <- NULL
observed <- NULL

for(fold in 1:length(folds_list)){
  
  training_data <- STH_kenya_df[-folds_list[[fold]], ]
  validation_data <- STH_kenya_df[folds_list[[fold]], ]
  
  fold_mod_spatial <- spaMM::fitme(cbind(tri_np, round(tri_nneg)) ~
                           alt + bio12 + lakes +
                           Matern(1|lat+long),
                         data=STH_kenya_df, family=binomial())
  
  x_valid_pred <- fold_mod_spatial %>% predict(validation_data)
  cross_validated_prediction <- c(cross_validated_prediction, x_valid_pred)
  observed <- c(observed, validation_data$tri_prev)
}

```

Next, I validated the spatial regression using a 10-fold cross validation technique. This is a plot of cross-validated predictions vs. observed values:

```{r, echo=F}
ggplot() + geom_point(aes(cross_validated_prediction, observed)) +
  xlab("Cross-validated predicted values") +
  ylab("Observed values") +
  ggtitle("Cross-validated predicted vs. observed values of whipworm prevalence")
```

This plot also follows a linear pattern, indicating no residual spatial autocorrelation. The mean squared error of the cross-validated spatial model with 4 covariates was the smallest of the other model specifications explored:

```{r, echo=F}
cross_validated_prediction %>% mse(observed)
```

### Risk mapping

Finally, I used the spatial regression to generate a risk map of whipworm prevalence across Kenya, as shown below. To generate the risk map, I rasterized model covariates, created a rasterstack with column names that matched model covariate names, and used the Matern model specification to predict across the raster stack. The risk map is shown below, masked to Kenya Admin 1 bounds.

```{r, include=F}
latitude_raster <- longitude_raster <- raster(nrows = KEN_alt %>% nrow(),
                                              ncols = KEN_alt %>% ncol(),
                                              ext = KEN_alt %>% extent())

longitude_raster[] <- coordinates(longitude_raster)[,1]
latitude_raster[] <- coordinates(latitude_raster)[,2]

pred_stack <- stack(KEN_alt, KEN_bioclim_resamp, KEN_lakes_r_resamp,
                    longitude_raster, latitude_raster)
names(pred_stack) <- c('alt', 'bio12', 'lakes', 'long', 'lat')

matern2 <- spaMM::fitme(cbind(tri_np, round(tri_nneg)) ~
               alt + bio12 + lakes +
               Matern(1|lat+long),
             data=STH_kenya_df, family=binomial())

predicted_risk <- raster::predict(pred_stack, matern2, type='response')

basemap <- leaflet() %>% addProviderTiles("CartoDB.Positron")

raster_colorPal <- colorNumeric(palette = topo.colors(64), domain = values(predicted_risk),
                                na.color = NA)

map_kenya <- basemap %>%
  addRasterImage(x = predicted_risk, color = raster_colorPal) %>%
  addPolygons(data = KEN_Adm_0, fillOpacity = 0, color = "grey", weight = 3)

map_kenya <- map_kenya %>%
  addLegend(title = "Predicted risk, whipworm prevalence",
            values = values(predicted_risk),
            pal = raster_colorPal)
```

```{r, echo=F}
map_kenya
```

## Results

Two key findings can be delineated from this analysis:

1. Altitude, annual precipitation, and distance to lakes are correlated with whipworm prevalence to a statistically significant degree, with p-values of `<6.75e-11`, `< 2e-16`, and `< 2e-16` respectively. Our non-spatial GLM also does *not* exhibit global clustering.

```{r, include=F}
or_alt = matern$fixef[[2]] %>% exp()
or_bio12 = matern$fixef[[3]] %>% exp()
or_dist2lakes = matern$fixef[[4]] %>% exp()
```

2. The inclusion of a spatially-correlated fixed effect improves model fit. However, the marginal impact of each of these factors on whipworm prevalence in western Kenya is extremely small. Here are the odds ratios:

* Altitude: **`r or_alt`**
* Bio12: **`r or_bio12`**
* Distance to lakes: **`r or_dist2lakes`**

3. The risk map generated with our spatial regression actually indicates *lower risk* of whipworm prevalence close to waterbodies, particularly Lake Victoria. This contradicts findings that closer proximity to Lake Victoria increases whipworm prevalence. However, it is important to note that this conclusion is limited by the fact that this region is primarily where our sample was taken. Thus, our risk estimation suffers somewhat from the lack of availability of whipworm prevalence data from regions outside of western Kenya.

## Conclusion

The findings from this study indicated statistically significant marginal impacts of altitude, annual precipitation, distance to waterbodies, or transit accessibility on whipworm prevalence in Western Kenya, however, those impacts were found to be negligible. Thus, this study neither evidenced nor contradicted the primary and secondary hypothesis.

Potential reasons for these results include the two factors that affected our analysis:

1. The data is not normally distributed even with a logit transformation, which is problematic for conducting a hypothesis test using the Moran's I.
2. There does not seem to be symmetry in the neighbor structure, which is problematic for defining spatial weights.

The first point limits our conclusions around the existence of clustering because the null hypothesis of the Moran's I is a normal distribution. Our outcomes data should be normally distributed for the hypothesis test to be conclusive. The second point is problematic because it points to a logical fallacy; it means that if k is a neighbor of j, then j is *not* a neighbor of k. Constructing neighbor structures in a cluster analysis of the whipworm prevalence data yields asymmetrical results and this should be investigated.

Other potential reasons include the outdatedness and incompleteness of the DCW data. The DCW database was compiled only in 1992, and Langaas et. al. (1995) has pointed to the challenges of using this data for statistical analysis due its outdatedness and incompleteness.

## Citations

Appleton, C. C., and E. Gouws. "The distribution of common intestinal nematodes along an altitudinal transect in KwaZulu-Natal, South Africa." Annals of Tropical Medicine & Parasitology 90.2 (1996): 181-188.

Appleton, C. C., M. Maurihungirire, and E. Gouws. "The distribution of helminth infections along the coastal plain of Kwazulu-Natal province, South Africa." Annals of Tropical Medicine & Parasitology 93.8 (1999): 859-868.

Chaiyos, J., et al. "MaxEnt modeling of soil-transmitted helminth infection distributions in Thailand." Parasitology research 117.11 (2018): 3507-3517.

Chammartin, Frédérique, et al. "Modelling the geographical distribution of soil-transmitted helminth infections in Bolivia." Parasites & vectors 6.1 (2013): 152.

Defense Mapping Agency (DMA), 1992. Digital Chart of the World. Defense Mapping Agency, Fairfax, Virginia. (Four CD-ROMs.)

Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37 (12): 4302-4315.

Global Administrative Areas ( 2012). GADM database of Global Administrative Areas, version 2.0. [online] URL: www.gadm.org.

Handzel, Thomas, et al. "Geographic distribution of schistosomiasis and soil-transmitted helminths in Western Kenya: implications for anthelminthic mass treatment." The American journal of tropical medicine and hygiene 69.3 (2003): 318-323.

Langaas, Sindre. Completeness of the Digital Chart of the World (DCW) database. UNEP/GRID-Arendal, 1995.

NASA Shuttle Radar Topography Mission (SRTM)(2013). Shuttle Radar Topography Mission (SRTM) Global. Distributed by OpenTopography. https://doi.org/10.5069/G9445JDF Accessed: 2020-05-03

Siza, Julius E., et al. "Prevalence of schistosomes and soil-transmitted helminths and morbidity associated with schistosomiasis among adult population in Lake Victoria Basin, Tanzania." The Korean journal of parasitology 53.5 (2015): 525.