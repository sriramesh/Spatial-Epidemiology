---
title: "Spatial Regressions with 2016 Uganda BCG Vaccination Data"
author: "Sri Ramesh"
date: "4/19/2020"
output: html_document
---

The following provides a multivariate regression analysis with spatially-correlated fixed effects using Professor Hugh Sturrock's BCG vaccine prevalence data in Uganda in 2016. First, I was interested in examining how distance to infrastructure and distance to water bodies affected BCG vaccination prevalence. For this, I pulled infrastructure vector data from the Digital Chart of the World (DCW) database on Uganda's water bodies (polygons), roads (lines), and railroads (lines). The Digital Chart of the World was compiled by the United States Defense Mapping Agency around 1992; while it is quite outdated, it is freely available online. Second, I was interested in examining the impact of urbanization trends on BCG vaccination prevalence. To that end, I used Professor Sturrock et. al's remotely sensed night-time lights imagery across Africa for the year 2013.

* Professor Sturrock's BCG vaccination prevalence data: https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week6/assignment/BCG_vaccination_UGA.csv
* Digital Chart of the World database for Uganda water bodies, roads, and railroads: https://www.diva-gis.org/gdata
* Intercalibration and Gaussian Process Modeling of Nighttime Lights
Imagery for Measuring Urbanization Trends in Africa 2000–2013: https://geodata.globalhealthapp.net/

Using these data sets, I examined 4 model specifications with BCG vaccination prevalence as the dependent variable:

1. **Model 1:** Generalized Linear Model (GLM) of BCG vaccination prevalence in 2016 Uganda with 4 covariates: distance to water bodies, distance to roads, distance to railroads, and night-time lights across Africa in 2013.
2. **Model 2:** GLM of BCG vaccination prevalence in 2016 Uganda with 3 covariates: distance to water bodies, distance to railroads, and night-time lights across Africa in 2013.
3. **Model 3:** GLM of BCG vaccination prevalence in 2016 Uganda with 2 covariates: distance to water bodies and distance to railroads.
4. **Model 4:** GLM of BCG vaccination prevalence in 2016 Uganda with 1 covariate: night-time lights across Africa in 2013.

Because I only know how to predict with raster data, I used Model 2 to fit a spatial covariance model and conduct cross-validation (Questions 1-2 and 4-8 on the rubric), and Model 4 to predict over a raster (Questions 3 and 7 on the rubric).

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
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

# global options
options(stringsAsFactors = F) # import quality
par(mar=c(1,1,1,1)) # ensure plot margin sizes are always large enough

# load vaccination data and calculate vaccine prevalence
BCG_vaccination_UGA <- "BCG_vaccination_UGA.txt" %>% 
  read.table(header = T, sep = ",", dec = ".") %>%
  mutate(prevalence = number_positive / numer_examined)

# load covariates and set crs ----
ntl_2013_raster <- 'gp2013africa.tif' %>% raster() # Annual mighttime lights (NTL) images for Africa, 2013
roads <- readOGR("UGA_rds", "UGA_roads") #roads in Uganda
railroads <- readOGR("UGA_rrd", "UGA_rails") #railroads in Uganda
waterbodies <- readOGR("UGA_wat", "UGA_water_areas_dcw") #inland waterbodies in Uganda

# convert to spdf
BCG_vaccination_UGA_spdf <- SpatialPointsDataFrame(coords = BCG_vaccination_UGA %>% select(lng,lat),
                                                   data = BCG_vaccination_UGA %>% select(ID, numer_examined, number_positive, prevalence))

```

Here is a web map of the BCG vaccination prevalence points alongside the data on Uganda's water bodies, roads, and railroads:

```{r, echo = T, warning=F, message=F}
pal = colorNumeric("Oranges", BCG_vaccination_UGA$prevalence)

leaflet() %>% addProviderTiles("Stamen.TonerLite") %>%
  addCircleMarkers(data = BCG_vaccination_UGA,
                   fillOpacity=1,
                   fillColor= ~pal(prevalence),
                   radius=~prevalence*2, stroke=TRUE, weight=.1) %>%
  addPolylines(data = roads, weight=1.5, color = "dimgrey") %>%
  addPolylines(data = railroads, weight=2, color = "grey") %>%
  addPolygons(data = waterbodies, weight=1.5, color="dodgerblue", opacity = .5) %>%
  addLegend(pal = pal, values = BCG_vaccination_UGA$prevalence, title = "BCG Vaccine Prevalence (%)")

```

The light blue water bodies indicate water bodies, the thicker grey lines indicate railroads, and the lighter grey lines indicate roads in Uganda. Upon first look, we see potential clustering of higher prevalence in and around Kampala and other regions of Uganda. We also see higher densities of railroads and roads in some areas, indicated where all Uganda is more urbanized. I now investigate how well these covariates, alongside the 2013 night-time lights data (not pictured here) explains the spatial pattern of the BCG prevalence values.

# Analysis of GLM with Spatially-Correlated Fixed Effects (Q1-Q2 and Q4-Q8)

First, I created the dataset of covariates for the regression modeling by calculating the distance between each prevalence point and the closest water body, road and railroad and appending them as 3 separate columns. I also extracted the night-time lights data as values and appended them as 1 additonal column.

It took quite a while to calculate the distances to roads, water bodies, and railroads (especially the roads data), so instead of re-running the commands here, I saved the finalized dataset and have pre-loaded it here:

```{r, echo = T, warning=F, message=F}
BCG_vaccination_UGA <- read.csv("BCG_vaccination_UGA.csv")
BCG_vaccination_UGA %>% colnames()
```

## 1 Fit a GLM to the data

Now I fit the 4 model specifications described above to the data. All of the models are specified as binomial generalized linear models (GLM). The 1st term is the number of positive cases, and the 2nd term is the number of negative cases.

```{r, echo = T, warning=F, message=F}
# Model 1: Nighttime Lights in 2013, Distance to water bodies, distance to roads, and distance to railroads
glm_mod_1 <- glm(cbind(number_positive, numer_examined - number_positive) ~
                   ntl_2013 + dist_to_water + dist_to_roads + dist_to_railroads,
                 data=BCG_vaccination_UGA, family=binomial())
glm_mod_1 %>% summary()
```

Models 1 indicates that distance to roads is a statistically insignificant predictor of BCG vaccination prevalence, and night-time lights is is only slightly important in terms of statistical significance. We run Models 2 and 3, which exclude distance to roads and both distance to roads and night-time lights respectively: 

```{r, echo = T, warning=F, message=F}
# Model 2: all the above, without distance to roads
glm_mod_2 <- glm(cbind(number_positive, numer_examined - number_positive) ~
                   dist_to_water + dist_to_railroads + ntl_2013,
                 data=BCG_vaccination_UGA, family=binomial())
glm_mod_2 %>% summary()

# Model 3: only distance to water and distance to railroads
glm_mod_3 <- glm(cbind(number_positive, numer_examined - number_positive) ~
                   dist_to_water + dist_to_railroads,
                 data=BCG_vaccination_UGA, family=binomial())
glm_mod_3 %>% summary()
```

The results show that Models 2 and 3 have comparable residual deviances, but Model 2's is slightly lower. So we will go with Model 2's covariates to fit a spatial regression.

## 2 Check fitted values

I now plot values predicted by Model 2 against observed values of BCG vaccination prevalence:

```{r, echo = T, warning=F, message=F}
plot_2 <- ggplot() + geom_point(aes(glm_mod_2$fitted, BCG_vaccination_UGA$prevalence)) +
  xlab("Fitted BCG Vaccination Prevalence") +
  ylab("Observed BCG Vaccination Prevalence") +
  ggtitle("BCG Vaccination Prevalence Fitted v. Observed, Model 2, non-spatial")
plot_2
```

The plot indicates potential clustering at higher and lower prevalence points. There is no linear trend of any kind. This indicates poor model fit, upon first look.

## 4 Compute and plot a correlogram based on Moran's coefficient for the residuals

I now confirm the presence of clustering by generating a correlogram based on the global Moran's I test for spatial autocorrelation:
 
```{r, echo = T, warning=F, message=F}
nbc <- 10
cor_r_2 <- pgirmess::correlog(coords=BCG_vaccination_UGA %>% select(lng, lat),
                            z=glm_mod_2$residuals,
                            method="Moran", nbclass=nbc)
# view correlogram
cor_r_2 # definite clustering up to 1.7 decimal degrees

# plot correlogram
correlogram_2 <- cor_r_2 %>% as.data.frame
correlogram_2$variable <- "residuals_glm"

correlogram_2 %>% subset(variable %in% "residuals_glm") %>% ggplot(aes(dist.class, coef)) + 
  geom_hline(yintercept = 0, col="grey") +
  geom_line(col="steelblue") + 
  geom_point(col="steelblue") +
  xlab("distance") + 
  ylab("Moran's coefficient")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Correlogram for Model 2 (no spatial effect)")

```

The correlogram/ Moran's I statistics at each distance class confirm the presence of global spatial autocorrelation up to 1.7 decimal degrees, which is distance class 3 in this correlogram. The presence of clustering is problematic, because one of the key assumptions of linear regressions is independent and identically distributed (iid) residuals.

## 5 Add a spatial effect into the GLM and re-test for residual spatial autocorrelation

I now re-specify the model using a Matern covariance model, which relaxes the iid assumption through the inclusion of a spatially-correlated fixed effect. We can see if the Matern specification removes the clustering:

```{r, echo = T, warning=F, message=F}
# fit matern covariance model (i.e. spatial model)
glm_mod_2_spatial <- spaMM::fitme(cbind(number_positive, numer_examined - number_positive) ~
                                    dist_to_water + dist_to_railroads + ntl_2013 +
                                    Matern(1|lat+lng),
                                  data=BCG_vaccination_UGA, family=binomial())
glm_mod_2_spatial %>% summary()

# view spatial model fitted vs. observed values
plot_2 <- ggplot() + geom_point(aes(glm_mod_2_spatial$fv, BCG_vaccination_UGA$prevalence)) +
  xlab("Fitted BCG Vaccination Prevalence") +
  ylab("Observed BCG Vaccination Prevalence") +
  ggtitle("BCG Vaccination Prevalence Fitted v. Observed, Model 2 with spatial effect")
plot_2
```

The plot of the fitted values for the the spatial regression is encouraging because it indicates a much stronger fit between predicted and observed values than the non-spatial model with the same covariates. Looking at the correlograms after inclusion of the spatial effect will confirm whether all clustering was removed:

```{r}
# view and plot correlogram with spatially correlated fixed effect
nbc <- 10
cor_r_2_spatial <- pgirmess::correlog(coords = BCG_vaccination_UGA %>% select(lng, lat),
                            z = residuals(glm_mod_2_spatial),
                            method="Moran", nbclass=nbc)
cor_r_2_spatial

correlogram_2_spatial <- cor_r_2_spatial %>% as.data.frame
correlogram_2_spatial$variable <- "residuals_glm" 

correlogram_2_spatial %>% subset(variable %in% "residuals_glm") %>% ggplot(aes(dist.class, coef)) + 
  geom_hline(yintercept = 0, col="grey") +
  geom_line(col="steelblue") + 
  geom_point(col="steelblue") +
  xlab("distance") + 
  ylab("Moran's coefficient")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Correlogram, after spatial effect included")

```

The correlogram shows high p-values across all distance classes, meaning we must fail to reject the null at the 0.05 alpha level and conclude no global clustering at any of those distance classes. Our spatial model has removed the clustering effect seen in the non-spatial version of Model 2.

## 8 Validate the model and compute mean squared error

Finally, I use a v-fold cross-validation approach to compute the mean squared error for values predicted by my spatial model (i.e. Model 2 with spatially-correlated fixed effect), and determine which has the smallest generalization error.

```{r, echo = T, warning=F, message=F}
# take 20% to act as validation set
set.seed(1)
validation_rows <- 1:nrow(BCG_vaccination_UGA) %>% sample(40)
BCG_vaccination_UGA_train <- BCG_vaccination_UGA[-validation_rows,] # training sets
BCG_vaccination_UGA_valid <- BCG_vaccination_UGA[validation_rows,] # validation, will predict these using training sets

# Fit spatial GLM on only the training data
glm_mod_2_spatial_validation <- spaMM::fitme(cbind(number_positive, numer_examined - number_positive) ~
                                               dist_to_water + dist_to_railroads + ntl_2013 + Matern(1|lat+lng), data=BCG_vaccination_UGA_train,
                                             family=binomial())

# Use spatial GLM to predict the 20% validation data and compare with actual data to see how well the model did
predictions_validation <- glm_mod_2_spatial_validation %>% predict(BCG_vaccination_UGA_valid)

ggplot() + geom_point(aes(predictions_validation %>% as.vector(), BCG_vaccination_UGA_valid$prevalence)) +
  xlab("Validated predicted values of prevalence") +
  ylab("Validated observed values of prevalence") +
  ggtitle("Predicted vs. Observed BCG Vaccination Prevalence, Spatial GLM Model 2, Validated")

# Calculate mse
predictions_validation %>% mse(BCG_vaccination_UGA_valid$prevalence)
```

# Predicting over rasters (Q3 and Q7)

I now turn to predict using Model 4, which is simply a univariate GLM of BCG vaccination prevalence using Professor Sturrock et. al's remotely sensed nighttime lights imagery across Africa subsetted to the year 2013:

```{r, echo = T, warning=F, message=F}
glm_mod_4 <- glm(cbind(number_positive, numer_examined - number_positive) ~
                   ntl_2013,
                 data=BCG_vaccination_UGA, family=binomial())
glm_mod_4 %>% summary()
```

I chose the **Kampala** administrative region of Uganda as my window to predict over. I subsetted the 2013 nighttime raster layer to the Kampala bounding box:

```{r, echo = T, warning=F, message=F}
UGA_Adm_1 <- raster::getData("GADM", country="UGA", level=1) # Uganda boundaries
proj4string(UGA_Adm_1) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ')
Kampala <- UGA_Adm_1 %>% subset(NAME_1 %in% "Kampala")
ntl_2013_kampala <- ntl_2013_raster %>% crop(Kampala) %>% mask(Kampala)
ntl_2013_kampala %>% plot(main="Prediction raster, 2013 Nighttime lights in Kampala region",
                          xlab="Latitude", ylab="Longitude")
Kampala %>% lines()
```

I will now predict over this raster using model specifications with and without spatially-correlated fixed effects.

## 3 Predict prevalence over a grid using the raster package 

My non-spatial model is specified as a univariate GLM, and produced the following results:

```{r}
pred_raster <- ntl_2013_kampala
names(pred_raster) <- c("ntl_2013")
predicted_risk_masked <- raster::predict(pred_raster, glm_mod_4, type='response') %>% mask(Kampala)
predicted_risk_masked %>% plot(xlab="Latitude", ylab="Longitude",
                               main="Predicted prevalence, univariate GLM")
Kampala %>% lines()
```

We can now compare these results with the spatial model specification.

## 7 Predict prevalence using the stack of covariates and the model

I will now account for how space contributes to the observed spatial pattern in vaccination prevalence by specifying a Matern covariance model:

```{r, echo = T, warning=F, message=F}
latitude_raster <- longitude_raster <- raster(nrows = ntl_2013_kampala %>% nrow(),
                                              ncols = ntl_2013_kampala %>% ncol(),
                                              ext = ntl_2013_kampala %>% extent())

longitude_raster[] <- coordinates(longitude_raster)[,1]
latitude_raster[] <- coordinates(latitude_raster)[,2]

pred_stack <- stack(ntl_2013_kampala, longitude_raster, latitude_raster)
names(pred_stack) <- c("ntl_2013", "lng", "lat")

glm_mod_4_spatial <- spaMM::fitme(cbind(number_positive, numer_examined - number_positive) ~
                                    ntl_2013 + Matern(1|lat+lng),
                                  data=BCG_vaccination_UGA, family=binomial())

predicted_prevalence_raster <- pred_stack %>% predict(glm_mod_4_spatial)
predicted_prevalence_raster_kampala <- predicted_prevalence_raster %>% mask(Kampala)
predicted_prevalence_raster_kampala %>% plot(main="Predicted prevalence, univariate Matern model",
                                             xlab="Latitude", ylab="Longitude")
Kampala %>% lines()
```

We can see that the Matern model's predictions differed from the predictions generated by the GLM.

# Citations

Savory DJ, Gething PW, Bennett, Andrade-Pacheco R, Midekisa A, Sturrock HJW. **Intercalibration and Gaussian Process Modeling of Nighttime Lights Imagery for Measuring Urbanization Trends in Africa 2000–2013.** Remote Sens. 2017, 9(7), 713; doi:10.3390/rs9070713
