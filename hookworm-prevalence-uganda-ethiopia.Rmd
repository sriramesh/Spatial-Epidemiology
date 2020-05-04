---
title: 'Spatial Variation in Risk in Hookworm Prevalence Data from Uganda and Ethiopia'
author: "Sri Ramesh"
date: "3/23/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

# load libraries ----
library(Metrics)
library(spatstat)
library(sp)
library(raster)
library(geoR)
library(gtools)
library(lme4)
library(leaflet)
library(oro.nifti)
library(tidyverse)

# global options ----
par(mar=c(1,1,1,1)) # ensures figure margins are not too large

```

The following analyzes spatial variation in point process and point referenced data. The point process data is Professor Sturrock's (UCSF) case-control data representing Nairobi, Kenya. The point referenced data is Professor Sturrock's hookworm prevalence data from Uganda and Tanzania. The data sets can be found here:

* (https://github.com/HughSt/HughSt.github.io/tree/master/course_materials/week3/Assignment). 

This document specifically explores 3 methods, the first for analyzing the spatial variation in risk and the latter two for interpolating the data:

1. kernel density estimation (KDE) to calculate relative risk in the Nairobi case-control data, and 
2. interpolation techniques (kriging and inverse-distance weighting (IDW)) of the hookworm prevalence data.

## 1 Kernel Density Estimation (KDE)

KDE applies to point process data. With KDE, we are interested in generating a continuous surface that illustrates the density of cases, and specific regions where clustering may be occurring (i.e. where we see clusters not simply due to spatial randomness). We can also use KDE to determine the risk of being a case relative to the risk of being a control, known as a *relative risk density estimation.*

I used the following steps to produce a relative risk density estimation with [Professor Sturrock's Nairobi case-control data found here:](https://github.com/HughSt/HughSt.github.io/tree/master/course_materials/week3/Assignment)

1. Load Nairobi case data file (*cases_nairobi.csv*)
2. Create kernel density estimates at different bandwidths
3. Calculate risk ratios
4. Map risk ratio rasters on top of leaflet base map

### 1. Load Nairobi case data file and boundaries

```{r q1loaddata, echo = T, warning=F, message=F}
nairobi_cases <- read.csv("cases_nairobi.csv")
nairobi_cases_spdf <- SpatialPointsDataFrame(coords = nairobi_cases[,c("lng", "lat")], 
                                             data = nairobi_cases[,c("X", "case")])

KEN_Adm0 <- raster::getData('GADM',country='KEN',level=0)

```

### 2. Plot kernel density estimate of cases at different bandwidths

```{r q1kde, echo = T, warning=F, message=F}
# define window
ken_owin <- owin(xrange=range(nairobi_cases$lng), yrange=range(nairobi_cases$lat))
# make ppp object
nairobi_cases_ppp <- ppp(nairobi_cases$lng, nairobi_cases$lat, window = ken_owin, 
                         marks=as.factor(nairobi_cases$case))
```

The Nairobi case data contains events, or cases of interest. Over each of these cases, we draw what are called *kernels.* Kernels have their own distribution. The case lies at the midpoint of its kernel distribution.

The individual kernel distributions are determined by their *bandwidth*, or how wide/narrow they are. Another name for a wider kernel is a larger bandwidth, and another name for a narrower kernel is a smaller bandwidth. Most importantly, the wider or **larger the bandwidth**, the **smoother the KDE.** We can plot different bandwidths and use a cross-validation method to select the best bandwidth for our Nairobi case data.

```{r, echo = T, warning=F, message=F}
# plots
par(mfrow=c(2,2))
nairobi_cases_ppp %>% density(0.02) %>% plot(main = "Bandwidth of 0.02")
nairobi_cases_ppp %>% density(0.1) %>% plot(main = "Bandwidth of 0.1")
nairobi_cases_ppp %>% density(0.5) %>% plot(main = "Bandwidth of 0.5")
nairobi_cases_ppp %>% density(bw.ppl) %>% plot(main = "CV-based bandwidth selection")

```

At each point on the x-axis, we sum up the heights of the individual kernel distributions. This sum produces a *kernel density estimation* (KDE). The KDE represents the density of cases at a given point. **As you can see, the KDE is entirely a function of those individual kernels' distributions.** 

### 3. Calculate risk ratios

Each point process has its own KDE. In our case we have two point processes:

* the Nairobi case data, and
* the Nairobi control data.

We can first calculate the KDE for each point process, choosing the most appropriate bandwidth for each one. Then, we can divide the case data's KDE by the control data's KDE to get a *relative risk estimate.* The relative risk estimate ranges from 0 to 1.00 (i.e. a percentage) and illustrates the chance that a given event is a case versus a control.

```{r q1riskratios, echo = T, warning=F, message=F}
nairobi_cases_ppp <- ppp(nairobi_cases$lng, nairobi_cases$lat, window = ken_owin, 
                         marks=as.factor(nairobi_cases$case))
rel_risk_est <- nairobi_cases_ppp %>% relrisk(relative = T)
par(mfrow=c(1,1))
rel_risk_est %>% plot(main="Relative Risk Estimate")

```

As you can see here, the relative risk of being a case versus being a control in our Nairobi case-control data ranges from 15% to 45%.

### 4. Map risk ratio rasters ontop of leaflet base map

Finally, we can map our case-control data alongside our relative risk ratios on a leaflet web map:

* Red circles: cases
* Grey circles: controls

```{r q1map, echo = T, warning=F, message=F}
nairobi_cases <- nairobi_cases %>% mutate(case = as.integer(case))
nairobicases_spdf <- SpatialPointsDataFrame(coords = nairobi_cases %>% select(lng, lat),
                                            data = nairobi_cases %>% select(X, case))
case_color_scheme <- colorNumeric(c("grey", "red"), nairobicases_spdf$case)

rel_risk_raster <- rel_risk_est %>% raster(crs = crs(KEN_Adm0))
pal <- colorNumeric(palette=tim.colors(64), domain=values(rel_risk_raster), na.color = NA)
basemap <- leaflet() %>% addProviderTiles("CartoDB.Positron")
basemap %>%
  addRasterImage(rel_risk_raster, opacity=0.6, colors = pal) %>%
  addCircleMarkers(data=nairobicases_spdf,
                   color = case_color_scheme(nairobicases_spdf$case),
                   radius=1) %>%
  addLegend(pal = pal, values = values(rel_risk_raster), title = "Relative Risk Ratios")

```

The map makes sense. The relative risk ratios are higher where there seem to be more cases *relative to controls*, which is on the right-hand side of our box window here. Similarly, the relative risk ratios are lower where there seem to be less cases relative to controls, which is on the left-hand side of our window.

In sum, we used kernel density estimation to confirm statistically the spatial pattern of risk (i.e. higher on the right side of the window, lower on the left side of the window).
 
## 2 Interpolation

The component of interest in point process data, like the Nairobi case-control data above, is determining the underlying spatial pattern. This is because the point process is determined *by* the underlying spatial pattern, which makes it possible for a concept like "risk" to be relevant and make sense.

However, risk as a concept does not really apply to point referenced data, such as infection prevalence. This is because the underlying spatial pattern of interest in point referenced data is the *values* (of prevalence, of air quality, etc.) taken at those points. Therefore, we are more interested in what those values, are where all we could not do surveys to collect them.

This is the purpose of interpolation. Interpolation is where we predict values where all we did not run the prevalence/air quality survey, using the values taken at points where all we *did* do the prevalence/air quality surveys. We explore two methods of interpolation in the following:

* Inverse Distance Weighting (IDW), and
* Kriging.

I used the following steps to interpolate [Professor Sturrock's hookworm prevalence data from Uganda and Tanzania found here:](https://github.com/HughSt/HughSt.github.io/tree/master/course_materials/week3/Assignment)

1. Load hookworm data from Tanzania/Uganda (*tanzania_uganda_hkprev.csv*)
2. Calculate the best IDW scenario by testing different powers
3. Map the "best" IDW scheme raster files on top of leaflet basemap
4. Comparison map: map kriged raster for the chosen window (Mwanza region in TZA) and the prevalence points
5. Cross-validation to compare interpolation methods
6. Visualize where predictions from IDW differ to kriging
7. Include a trend surface to see if that improves kriging estimates

### 1. Load hookworm data from Tanzania/Uganda

```{r q2load, echo=T, warning=F, message=F}
HK <- read.csv("tanzania_uganda_hkprev.csv")
TZA_Adm_1 <- raster::getData("GADM", country="TZA", level=1) # used tanzania instead of uganda

```

Here is what our hookworm prevalence data looks like right now. We are interested in predicting values of hookworm prevalence where all we don't have data.

```{r, echo=T, warning=F, message=F}
pal = colorNumeric("Oranges", HK$Hookworm_prev_perc)

basemap_hk <- HK %>% leaflet() %>% addTiles() %>%
  addCircleMarkers(~x, ~y, fillOpacity=1,
                   stroke=T,
                   radius=3,
                   weight=0.5,
                   fillColor= ~pal(Hookworm_prev_perc)) %>% 
  addLegend(pal = pal, values = ~Hookworm_prev_perc, title = "Known Hookworm prevalence (%)")
basemap_hk

```

As you can see, we have data on hookworm prevalence across many regions in East Africa. So to focus our interpolation efforts, we will choose a box or window to which we will limit our predictions/interpolations. **The window I chose was the *Mwanza* Adm 1 level region in Tanzania. This window is close to the Tanzania-Uganda border.**

```{r q2bestidw, echo=T, warning=F, message=F}
# set window
mwanza <- TZA_Adm_1[TZA_Adm_1$NAME_1 %in% "Mwanza", ]
mwanza_window <- owin(mwanza@bbox[1,], mwanza@bbox[2,])
```

### 2. Calculate the best IDW scenario by testing different powers

The IDW method essentially takes each of our points, looks at all of the neighboring points, and averages the values of malaria prevalence taken at those neighboring points. Finally, it applies a weight to each of those averages, and comes up with a final weighted average. Those final weighted averages represent the interpolated values.

The weight applied to each of those averages is determined by what's called a *power function.* There can be a wide range of power functions. Most importantly, the higher the power, the more weight assigned to closer neighbors of a given point. I now examine different power functions and choose an optimal power to run the IDW function:

```{r, echo=T, warning=F, message=F}
# define ppp object
TZA_hookworm_data_ppp <- ppp(HK$x, HK$y, marks = HK$Hookworm_prev_perc, window = mwanza_window)
# set parameters for 4 panel display
par(mfrow=c(2,2))
TZA_hookworm_data_ppp %>% idw(power=0.2, at="pixels") %>% plot(col=heat.colors(20), main="Power = 0.2")
TZA_hookworm_data_ppp %>% idw(power=0.5, at="pixels") %>% plot(col=heat.colors(20), main="Power = 0.5")
TZA_hookworm_data_ppp %>% idw(power=1, at="pixels") %>% plot(col=heat.colors(20), main="Power = 1")
TZA_hookworm_data_ppp %>% idw(power=2, at="pixels") %>% plot(col=heat.colors(20), main="Power = 2")
# determine optimal power function using for loop
powers <- seq(0.05, 2, 0.05)
mse_result <- NULL
for (power in powers) {
  CV_idw <- TZA_hookworm_data_ppp %>% idw(power=power, at="points")
  mse_result <- c(mse_result, mse(TZA_hookworm_data_ppp$marks, CV_idw))
}
optimal_power <- powers[which.min(mse_result)]
optimal_power

```

The optimal power was 2. This was also the highest power tested. 

### 3. Map the "best" IDW scheme raster files on top of leaflet basemap

Now, we will use this power as an input to our IDW function, which will predict/interpolate values. The "pixels" option means generating interpolated estimates across a grid of pixels (i.e. prediction grid). Alternatively, we could use the "points" option to generate interpolated estimates at every point, using leave-one-out-cross validation.

```{r q2mapidw, echo=T, warning=F, message=F}
# perform idw
TZA_hookworm_data_ppp_idw <- TZA_hookworm_data_ppp %>% idw(power = optimal_power, at="pixels")
```

Finally, we can map our interpolated data on leaflet to see how it looks:

```{r mapit, echo=T, warning=F, message=F}
# Convert to a raster
TZA_hookworm_data_idw_raster <- TZA_hookworm_data_ppp_idw %>% raster(crs = crs(TZA_Adm_1))

# define color ramp
colPal <- colorNumeric(tim.colors(), TZA_hookworm_data_idw_raster[], na.color = NA)
# layer map
basemap_hk %>% addRasterImage(TZA_hookworm_data_idw_raster, col = colPal, opacity = 0.7) %>%
  addLegend(pal = colPal, values = TZA_hookworm_data_idw_raster[], title = "Predicted Infection Prevalence (%)")

```

If you zoom into the Mwanza window on the map above, you will see where all points were interpolated on the grid per the IDW method. For example, the IDW method has predicted infection prevalence to be around 30-50% on the upper regions of our window, and 70-80% on the bottom left region of our window.

As you can see, these predictions were based directly on the circles which indicate where we actually have inflection prevalence values taken (by way of surveys, etc.). For example, it is clear that the bottom left region of our window has 70-80% predicted infection prevalence primarily because we have 5-6 values of infection prevalence in this region at 80-90%. Similarly, we have 30-40% predicted infection prevalence on the middle-right region of our window, where we see a cluster of actual infection prevalence values ranging from 30-40%.

### 4. Comparison map: map kriged raster for the chosen window (Mwanza region in TZA) and the prevalence points

Now we can turn to kriging to see how interpolated values differ.

In this step, first I created a *distance matrix* which is essentially the distance between every location where the hookworm prevalence value was collected, and every other location where the hookworm prevalence value was collected. This gives us a matrix of distances. These distance values go on the x-axis of my graph.

Then, I measured how similar each of the hookworm prevalence values were to each other. I did this by calculating the variance between each hookworm prevalence value (which is the difference in a pair of values squared), then dividing that by 2. This means that similarity of two given hookworm prevalence values is measured by half of the variance, or the *semi-variance*. The semi-variance values go on the y-axis of my graph.

```{r q2kriging, echo=T, warning=F, message=F}
# create geodata object
TZA_hookworm_data_geo <- HK %>% select(x, y, Hookworm_prev_perc) %>% as.geodata()
# generate distance matrix
MaxDist <- (HK %>% select(x,y) %>% dist() %>% max()) /2

# generate variocloud
VarioCloud <- TZA_hookworm_data_geo %>% variog(option="cloud", max.dist= MaxDist)
VarioCloud %>% plot(main = "Vario Cloud, all pairwise computations")

```

The resulting graph is what's called a variogram displayed as a cloud, or a "vario cloud". This vario cloud essentially illustrates clustering. Values that are more similar (i.e. those that have lower semi variance) are also closer together (i.e. are at smaller distance classes). Values that are more different (i.e. have higher semi variance) are also farther apart (i.e. are at higher distance classes).

As you can see, the Vario cloud version of the variogram is overwhelming to look at; it shows all values on the distance matrix and all values of semi-variance. However, all we are interested in is discerning any linear relationship between lower semi-variance values and lower distance classes and vice versa. Such a linear trend would indicate clustering or spatial autocorrelation.

We can simplify the vario cloud by binning distances into "distance classes" using arbitrary cutoff points for distances. The bins I chose simply ensure that the distance classes are big/small enough such that each distance class has a minimum amount of distances included in it to show on a graph.

```{r, echo=T, warning=F, message=F}
# bin by distance
Vario <- TZA_hookworm_data_geo %>% variog(max.dist = MaxDist, uvec = seq(0.01, MaxDist, 0.2)) 
Vario %>% plot(main = "Variogram, simplified")
```

The simplified variogram above clearly illustrates a linear trend, meaning we are ready to fit regression models through this trend to see which model best captures the relationship between semi-variance, and distance. I examine spherical and exponential co-variance models, and choose the one that has the lowest sum of squares. It turns out that the spherical model has the lowest sum of squares and thus best captures that relationship.

```{r, echo=T, warning=F, message=F}
# fit variogram to each model
VarioMod_sph <- Vario %>% variofit(cov.model = "sph")
VarioMod_exp <- Vario %>% variofit(cov.model = "exp")
# plot each model
Vario %>% plot(pch=16)
VarioMod_sph %>% lines(col="blue",lwd=2)
VarioMod_exp %>% lines(col="red",lwd=2)
# test if spherical model has lower sum of squares than exponential model
sum_sph <- summary(VarioMod_sph)
sum_exp <- summary(VarioMod_exp)
sum_sph$sum.of.squares < sum_exp$sum.of.squares
```

Now that I decided on my co-variance model, I can create a prediction grid, and krige to those points on the prediction grid.

```{r, echo=T, warning=F, message=F}
# create prediction grid
IDW <- TZA_hookworm_data_ppp %>% idw(power = optimal_power, at="pixels")
par(mfrow=c(1,1))
pred_grid_x <- rep(IDW$xcol,length(IDW$yrow))
pred_grid_y <- rep(IDW$yrow,length(IDW$xcol)) %>% sort()
pred_grid <- cbind(pred_grid_x, pred_grid_y)
# krig to those points
KrigPred <- TZA_hookworm_data_geo %>% krige.conv(loc=pred_grid, krige=krige.control(obj.model=VarioMod_sph))

```

```{r, echo=T, warning=F, message=F}
# plot kriged points
KrigPred_raster <- data.frame(x=pred_grid_x, y=pred_grid_y, z=KrigPred$predict) %>% rasterFromXYZ(crs = crs(TZA_Adm_1))
KrigPred_raster %>% plot(main = "Kriged raster for Mwanza window and prevalence points")
```

We can see the final kriged predictions a leaflet map as well.

```{r}
basemap_hk %>% addRasterImage(KrigPred_raster, col = colPal, opacity = 0.7) %>%
  addLegend(pal = colPal, values = TZA_hookworm_data_idw_raster[], title = "Predicted Infection Prevalence (%)")
```

As you can see, the IDW method and kriging method produced different predictions on our Mwanza window.

### 5. Cross validation to compare interpolation methods

I cross validated the values and compared mean squared errors to find that kriging produces a lower mean squared error. This means that kriging is a better method for interpolating this data. 

```{r q2compare, echo=T, warning=F, message=F}
CV_idw_opt <- TZA_hookworm_data_ppp %>% idw(power = optimal_power, at= "points")
CV_kriging_opt <- TZA_hookworm_data_geo %>% xvalid(model = VarioMod_sph)
CV_idw_opt %>% mse(HK$Hookworm_prev_perc)
CV_kriging_opt$predicted %>% mse(HK$Hookworm_prev_perc)

```

### 6. Visualize where predictions from IDW differ to kriging

Because these rasters are in the same resolution and extent, and all values of hookworm prevalence are on the same scale (0 to 100) we can compare them using subtraction. Here, I subtracted kriged estimates from IDW estimates of hookworm prevalence values to see where all the two methods similar and un-similar predictions.

```{r differences, echo=T, warning=F, message=F}
IDW_raster <- IDW %>% raster()
plot(IDW_raster - KrigPred_raster, main = "Difference matrix comparing Kriged estimates to IDW estimates",
     xlab = "Longitude", ylab = "Latitude")

```

This plot illustrates the difference between using IDW versus kriging to interpolate our hookworm prevalence data. Regions of greatest difference are indicated in green. It seems like the methodologies differed the most at the southeastern and southwestern extents of our window.

### 7. Include a trend surface to see if that improves kriging estimates

I included a 1st order trend in the spherical model used to fit the variogram to see if that would improve kriged estimates. 

```{r trend, echo=T, warning=F, message=F}

# Fit (spherical) variogram with 1st order trend
Vario_trend <- TZA_hookworm_data_geo %>% variog(max.dist = MaxDist, trend="1st")
Vario_trend_sph <- Vario_trend %>% variofit(cov.model = "sph")

# Get kriged values with 1st order trend
CV_kriging_opt_trend <- TZA_hookworm_data_geo %>% xvalid(model = Vario_trend_sph)

# Now compare
CV_kriging_opt_trend$predicted %>% inv.logit() %>% mse(HK$Hookworm_prev_perc)  
CV_kriging_opt$predicted %>% inv.logit() %>% mse(HK$Hookworm_prev_perc)

```

As you can see, kriging with or without a 1st order trend produces the same mean squared error. Thus, I conclude that the inclusion of a trend surface does not improve kriging estimates.
