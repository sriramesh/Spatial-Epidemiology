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
STH_kenya <- "Data/Kenya_STH.csv" %>% read.csv(header = T) %>%
  select(full_name_paper, lat, long, tri_prev)

STH_kenya <- STH_kenya %>%
  mutate(tri_prev = tri_prev/100,
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
roads <- readOGR("Data/KEN_rds", "KEN_roads") # roads
railroads <- readOGR("Data/KEN_rrd", "KEN_rails") # railroads
lakes <- readOGR("Data/KEN_wat", "KEN_water_areas_dcw") # inland lakes
canals <- readOGR("Data/KEN_wat", "KEN_water_lines_dcw") #inland canals and rivers

# set crs ----
wgs84 <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '
proj4string(STH_kenya) <- wgs84 %>% CRS()
proj4string(roads) <- wgs84 %>% CRS()
proj4string(railroads) <- wgs84 %>% CRS()
proj4string(lakes) <- wgs84 %>% CRS()
proj4string(canals) <- wgs84 %>% CRS()

roads_map <- roads %>% crop(STH_kenya %>% extent())
railroads_map <- railroads %>% crop(STH_kenya %>% extent())
lakes_map <- lakes %>% crop(STH_kenya %>% extent())
canals_map <- canals %>% crop(STH_kenya %>% extent())

STH_kenya_df <- read.csv("STH_kenya_covariates.csv")
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

#### Covariate transformations

I applied the following transformations to my covariates data:

* resampling the altitude raster, and
* calculating distances between prevalence points to lakes, rivers/canals, roads, and railroads.

##### 1) Resampling altitude raster

An initial look at my two raster layers - Bio12 and altitude - indicated that the Bio12 raster was at a higher resolution than the altitude, and that both had different extents:

```{r, echo=T}
alt <- raster::getData('alt', country = 'KEN', mask = F)
alt
```

```{r, echo=T}
bioclim <- raster::getData('worldclim', var='bio', res=2.5) # Bio12, higher resolution
bioclim[[12]]
```

```{r, include=F}
KEN_Adm_0 <- raster::getData("GADM", country="KEN", level=0)
KEN_bioclim <- bioclim[[12]] %>% crop(KEN_Adm_0 %>% extent())
KEN_alt <- alt %>% crop(KEN_bioclim %>% extent())
# resample higher res raster using lower res raster as reference
KEN_bioclim_resamp <- KEN_bioclim %>% resample(KEN_alt, method = "ngb")
```

I resampled the Bio12 raster to match the lower resolution of the altitude raster to ensure that the rasters were comparable in my forthcoming regression analysis. Because the annual precipitation data is continuous, I interpolated the data using a nearest-neighbor approach, which assigns the value of the nearest cell to cells in the prediction grid. The resampling procedure worked. We can see the altitude data's extent and resolution:

```{r, echo=T}
KEN_alt
```

And the resampled annual precipitation raster's extent and resolution:

```{r, echo=T}
KEN_bioclim_resamp
```

The resolution and extents now match up and can be used as covariates in regressions.

##### 2) Calculating distances

I calculated the distance in meters from each whipworm prevalence point to the nearest of each of the following 4 features: lakes, rivers/canals, roads, and railroads. I primarily used the `dist2line` function for this task. Each minimum distance was saved as a column in my dataset, and represents a unique distance-based covariate in the spatial and non-spatial regression analysis.

I also rasterized distance-based covariates which were included in my spatial regression to predict with them, and generate a risk map in the final step.

### Spatial and non-spatial regression analyses

For my regression analysis, I first specified non-spatial generalized linear models (GLMs) with various combinations of the 6 covariates (4 distance-based covariates, annual precipitation, and elevation). All models specified included binomial outputs. Some of the model specifications resulted in intercepts that had a negative coefficient, and I excluded these models from consideration. The primary candidates in terms of model specification were:

1. A non-spatial GLM with 2 covariates: altitude, and annual precipitation.
2. A non-spatial GLM with 4 covariates: altitude, distance to canals, distance to lakes, and distance to roads

Using a step-wise process of feature selection, I compared residual deviance values and cross-validated mean squared errors to identify the model with the lowest residual deviance.

I ultimately chose non-spatial GLM with 4 covariates: altitude, distance to canals, distance to lakes, and distance to roads. This model:

* had a positive coefficient on the intercept,
* the lowest residual deviance values, and
* the lowest cross-validated mean squared error.

Non-spatial GLM with 4 covariates:

```{r, echo = F}
options(scipen = 999)
glm_mod_1 <- glm(cbind(tri_np, round(tri_nneg)) ~
                   alt +
                   dist_to_canals + dist_to_lakes +
                   dist_to_roads,
                 data=STH_kenya_df, family=binomial())
glm_mod_1 %>% summary()
```

After choosing the 4-covariate GLM, I ran a few diagnostics to determine the extent of residual spatial autocorrelation. I first examined a residual plot for any linear patterns:

```{r, echo = F}
residual_plot <- ggplot() + geom_point(aes(glm_mod_1$fitted, STH_kenya_df$tri_prev_logit)) +
  xlab("Expected log prevalence values") +
  ylab("Observed log prevalence values") +
  ggtitle("Observed vs. expected log prevalence values, GLM")
residual_plot
```

The plot illustrates a lack of linearity. This indicates that these covariates, while statistically significant, still do not explain much of the spatial pattern underlying the whipworm prevalence data.

I then used a correlogram to examine clustering effects. I took note of the highest distance class at which statistically significant residual clustering can be detected. The 4-covariate, non-spatial GLM exhibits residual clustering up to 0.15 decimal degrees as shown in the correlogram below.

```{r, echo=F}
options(scipen = -1)
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

Finally, I specified a spatial regression to correct for the residual spatial autocorrelation. I used a Matern model, which includes a spatially-correlated fixed effect, as well as binomial outputs as shown below.

```{r, echo=F}
options(scipen = 999)
matern <- spaMM::fitme(cbind(tri_np, round(tri_nneg)) ~
                         alt +
                         dist_to_canals + dist_to_lakes +
                         dist_to_roads +
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

row <- row.names(coefs) %in% c('dist_to_canals')
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, 'Cond. SE']
or_ci_dist2canals <- c(lower, upper) %>% exp()

row <- row.names(coefs) %in% c('dist_to_lakes')
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, 'Cond. SE']
or_ci_dist2lakes <- c(lower, upper) %>% exp()

row <- row.names(coefs) %in% c('dist_to_roads')
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, 'Cond. SE']
or_ci_dist2roads <- c(lower, upper) %>% exp()

```

* Altitude: **`r or_ci_alt`** 
* Distance to rivers/canals: **`r or_ci_dist2canals`** 
* Distance to lakes: **`r or_ci_dist2lakes`** 
* Distance to roads: **`r or_ci_dist2roads`** 

I checked the residual plot and correlogram for the spatial regression to determine if the inclusion of the spatially-correlated fixed effect removed the residual spatial autocorrelation. Specifically, I looked for linearity in the residual plot, as shown below.

```{r, echo=F}
plot <- ggplot() + geom_point(aes(matern$fv, STH_kenya_df$tri_prev_logit)) +
  xlab("Expected log prevalence values") +
  ylab("Observed log prevalence values") +
  ggtitle("Observed vs. expected log prevalence values, Matern")
plot
```

I also looked for a lack of clustering across all distance classes of the correlogram corresponding to the spatial regression, as shown below.
 
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

The spatial regression residual plot indicates linearity between observed and expected values. Similarly, our correlogram indicates that at no distance classes do we have a statistically significant Moran's I statistic. This means that the inclusion of the spatially correlated fixed effect removed all residual clustering successfully. Our spatial model has independent and identically distributed residuals (iid), satisfying a key assumption of linear regression modeling.

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
                           alt +
                           dist_to_canals + dist_to_lakes +
                           dist_to_roads +
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

```{r, include=F}
latitude_raster <- longitude_raster <- raster(nrows = KEN_alt %>% nrow(),
                                              ncols = KEN_alt %>% ncol(),
                                              ext = KEN_alt %>% extent())

longitude_raster[] <- coordinates(longitude_raster)[,1]
latitude_raster[] <- coordinates(latitude_raster)[,2]

pred_stack <- stack(KEN_alt, KEN_bioclim_resamp, longitude_raster, latitude_raster)
names(pred_stack) <- c('alt', 'bio12', 'long', 'lat')

matern2 <- spaMM::fitme(cbind(tri_np, round(tri_nneg)) ~
               alt + bio12 +
               Matern(1|lat+long),
             data=STH_kenya_df, family=binomial())

predicted_risk <- raster::predict(pred_stack, matern2, type='response')
```

Finally, I used the spatial regression to generate a risk map of whipworm prevalence across Kenya, as shown below. To generate the risk map, I rasterized model covariates, created a rasterstack with column names that matched model covariate names, and used the Matern model specification to predict across the raster stack. The risk map is shown below, masked to Kenya Admin 1 bounds.

```{r, echo=F}
predicted_risk %>% mask(KEN_Adm_0) %>% plot(main="Whipworm prevalence risk map, Kenya",
                        xlab="Latitude", ylab="Longitude")
KEN_Adm_0 %>% lines()
```

## Results

Two key findings can be delineated from this analysis:

1. With the inclusion of a spatially-correlated fixed effect, **none** of our four covariates of interest (altitude, distance to lakes, distance to canals/rivers, and distance to roads) have statistically significant associations with *T.Trichiura* prevalence.

The 95% confidence intervals expressed as odds ratios of our spatial regression lead us to this conclusion. The coefficients indicate that the inclusion of the spatially correlated fixed effect effectively removed the statistically significant association between *T.Trichiura* prevalence and each of our four covariates of interest seen on the non-spatial GLM. This ultimately indicates that space accounts for much of that variation in the spatial pattern of whipworm prevalence.

We can evidence this using an anlysis of the odds ratios on each covariate coefficient in the Matern model as well:

```{r, include=F}
or_alt = ((1 - exp(-0.00382032))*100) %>% round(4)
or_dist2canals = ((1 - exp(-0.00009361))*100) %>% round(4)
or_dist2lakes = ((1 - exp(-0.00001498))*100) %>% round(4)
or_dist2roads = (1+0.000276 - 1) %>% round(4)
```

* *Altitude*: A one-unit increase in altitude is associated with only a `r or_alt` percent decrease in testing positive for whipworm.
* *Distance to rivers/canals*: A one-unit increase in distance to the nearest river or canal is associated with only a `r or_dist2canals` percent decrease in testing positive for whipworm.  
* *Distance to lakes*: A one-unit increase in distance to the nearest river or canal is associated with only a `r or_dist2lakes` percent decrease in testing positive for whipworm.  
* *Distance to roads*: A one-unit increase in distance to the nearest road is associated with only a `r or_dist2roads` percent increase in testing positive for whipworm.  

2. The risk map generated with our spatial regression indicates lower risk of whipworm prevalence around western Kenya. However, there are significant limitations to this conclusion: this is where our sample was taken from. Thus, our risk estimation suffers somewhat from the lack of availability of whipworm prevalence data from regions outside of western Kenya.

#### Results from non-spatial regression analysis

Our non-spatial regression analysis neither evidenced nor contradicted the literature which suggests that proximity to waterbodies, is associated with higher infection prevalence of *T.Trichiura*. It suggested that altitude, distance to waterbodies (lakes, canals, and rivers), and distance to roads were statistically significant predictors of whipworm prevalence, with p-values of 8.81e-09, 1.2348e-07, 0.00303, and <2e-16 respectively. Further, a relatively tiny coefficient on the altitude covariate of the non-spatial GLM also supported conclusions in prior literature that altitude has mixed impacts on *T.Trichiura* prevalence.

However, the odds ratio interpretation of the coefficients shows us that for a one-unit increase in minimum distance to a canal or river, we expect to see only a **`r or_dist2canals`** percent decrease in the odds of testing positive for whipworm. Similarly, a one-unit increase in minimum distance to a lake is correlated with only a **`r or_dist2lakes`** percent decrease in the odds of testing positive for whipworm.

In terms of the robustness of the non-spatial GLM, we can use the deviance residuals, null deviance, and residual deviance to conclude that while our non-spatial GLM specified with these 4 covariates (distance to rivers/canals, distance to lakes, altitude, and distance to roads) might be a relatively good fit for our data, much of the variance in infection prevalence is *not* explained by our covariates. Taking the non-spatial GLM deviance residuals, we can see that the third quartile (3Q) of the deviance residuals has a value of 1.7988, which is less than our benchmark of 2. This means that the first 75 percent of our model observations are within an acceptable distance of our predicted values, with "acceptable" being defined as a value less than 2. However, the model's null deviance totals 1287.54, and its  residual deviance totals 835.89. These are relatively high values and indicate our model could benefit from the inclusion of other, better fit covariates.

## Conclusion

All in all, the findings from this study neither evidenced nor contradicted the primary and secondary hypothesis. It neither evidenced nor contradicted the literature which suggests that proximity to waterbodies is associated with higher infection prevalence of *T.Trichiura*, as well as transit accessibility. The study supported conclusions in prior literature that altitude has mixed impacts on *T.Trichiura* prevalence, and it also found no statistically significant association between annual precipitation and *T.Trichiura* prevalence.

Potential reasons for these insignificant results include the two factors that affected our cluster analysis:

1. The data is not normally distributed even with a logit transformation, which is problematic for conducting a hypothesis test using the Moran's I.
2. There does not seem to be symmetry in the neighbor structure, which is problematic for defining spatial weights.

The first point limits our conclusions around the existence of clustering because the null hypothesis of the Moran's I is a normal distribution. Our outcomes data should be normally distributed for the hypothesis test to be conclusive.

The second point is problematic because it points to a logical fallacy; it means that if k is a neighbor of j, then j is *not* a neighbor of k. Constructing neighbor structures in a cluster analysis of the whipworm prevalence data yields asymmetrical results and this should be investigated.

Other potential reasons include lack of data availability in regions outside of western Kenya, and small sample size within western Kenya with regards to prevalence points (n = 68). These reasons also limit the conclusions we can draw around infection risk outside of western Kenya. 

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
