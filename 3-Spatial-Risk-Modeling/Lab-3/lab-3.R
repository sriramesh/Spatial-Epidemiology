## Objective: explore methods to understand and predict risk across space from point data
# point PROCESS or point PATTERN data (e.g. cases, controls)
# point REFERENCED or point LEVEL data (e.g. air quality values)

# load libraries for this week ----
library(Metrics)
library(spatstat)
library(sp)
library(raster)
library(geoR) # requires downloading XQuartz: https://www.xquartz.org/
library(gtools)
library(lme4)
library(leaflet)
library(oro.nifti)
library(tidyverse)


# utility data ----
url_casecontrol <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week3/Lab_files/CaseControl.csv"
url_malaria <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week1/Lab_files/Data/mal_data_eth_2009_no_dups.csv"


# global options ----
par(mfrow=c(2,2)) # figures on 2 rows
par(mar=c(1,1,1,1)) # ensures figure margins are not too large

# load data -----
# obfuscated case and control malaria data from namibia, lats and longs
CaseControl <- read.csv(url_casecontrol)

# boundaries for namibia, admin 0 level
NAM_Adm0 <- raster::getData('GADM',country='NAM',level=0)

# manipulate cc data ----
head(CaseControl) # looks like it's at the household level, and there case/control is a boolean value

# subset cases and controls
Cases <- CaseControl %>% filter(case %in% 1)
Controls <- CaseControl %>% filter(case %in% 0)

# also make a spatial points df
CaseControl_SPDF <- SpatialPointsDataFrame(coords = CaseControl[,c("long", "lat")],
                                           data = CaseControl[,c("household_id", "case")])

# plot cases and controls on leaflet map ----

# first create case v. control color scheme for legend
case_color_scheme <- colorNumeric(c("firebrick", "lawngreen"), CaseControl_SPDF$case)
# the way colornumeric works is cool.
# you first provide a vector of colors, then you provide the variable to apply the color vector to

# add basemap and point pattern layer
basemap <- leaflet() %>% addProviderTiles("CartoDB.Positron")

basemap %>% addCircleMarkers(data=CaseControl_SPDF,
                             color = case_color_scheme(CaseControl_SPDF$case),
                             radius = 2)


# 1 risk mapping using KDE ----

# define window for point PATTERN or point LEVEL data ----
Nam_Owin <- spatstat::owin(xrange=range(CaseControl$long),
                 yrange=range(CaseControl$lat))
# owin creates an object of class "owin" representing an observation window in the two-dimensional plane

# define object for the ppp (Point Pattern object of Points) using that window ----
Cases_ppp <- spatstat::ppp(Cases$long,
                 Cases$lat,
                 window = Nam_Owin)
# ppp function creates an object of class "ppp" representing a point pattern dataset in the two-dimensional plane.
plot(Cases_ppp)

# now you can plot a KDE for the cases (remember we need a different KDE for the controls) ----
case_density <- stats::density(Cases_ppp) # object of class list
plot(case_density) # Units are intensity of points per unit square

# so this is not the kernels for individual cases. this is the final KDE.

# we can change bandwidths of the individual kernel estimates, to make them larger or smaller.

# these are essentially heat maps around the point process cases:
Cases_ppp %>% density(0.02) %>% plot() #bandwidth of 0.02, quite small
Cases_ppp %>% density(0.1) %>% plot() #bandwidth of 0.1, quite large
Cases_ppp %>% density(bw.ppl) %>% plot() # automatic bandwidth selection based on cross-validation

# we can also put the heatmap/KDE on leaflet basemaps like so:

# convert the DENSITY object to a raster layer + apply a CRS
density_raster <- Cases_ppp %>% density(bw.ppl) %>% raster(crs = crs(NAM_Adm0)) # here, we're using the CRS of the namibia adm 1 boundaries data

# Plot
basemap %>% addRasterImage(density_raster, opacity=0.6) # opacity of the raster image file

# But this is just a density of CASES
# it doesn’t account for the denominator - the CONTROLS. We need another KDE for the controls tho

# To do this, we can use the kelsall & diggle et. al. method (it's in an academic paper)
# it calculates the RATIO of the density estimate of cases to controls
# cross-validation maybe?

# add marks to the points ----
# Marks are just VALUES associated with each point
# Case: 1, Control: 0

# just revising the original ppp object to include option for marks:
CaseControl_ppp <- ppp(CaseControl$long,
                       CaseControl$lat, 
                       window = Nam_Owin, # same window as cases KDE
                       marks=as.factor(CaseControl$case)) #marks is a factorized vector, 1 or 0

# use relrisk function from spatstat package
# will show you the "risk" (probability?) of being a case, relative to background pop
risk_est <- relrisk(CaseControl_ppp) #outputs an image!
plot(risk_est)

rel_risk_est <- CaseControl_ppp %>% relrisk(relative = T) # this is RELATIVE risk, if option is left unspecified, it's not relative
plot(rel_risk_est)

# maybe we can put these risk estimations on our basemap now? ----
risk_raster <- raster(risk_est, crs = crs(NAM_Adm0)) #convert to raster + add CRS
rel_risk_raster <- raster(rel_risk_est, crs = crs(NAM_Adm0)) #convert to raster + add CRS

# define a color ramp for the raster image
pal <- colorNumeric(palette=tim.colors(64), domain=values(risk_raster), na.color = NA)

basemap %>% addRasterImage(risk_raster, opacity=0.6, colors = pal)
basemap %>% addRasterImage(rel_risk_raster, opacity=0.6, colors = pal)


# 2 Interpolation of point referenced/level data ----

# load data ----
# ethiopia malaria prevalence data
ETH_malaria_data <- read.csv(url_malaria, header=T)

# ethiopia admin level 1
ETH_Adm_1 <- raster::getData("GADM", country="ETH", level=1)

# IDW (inverse distance-weighting) ----
## do this using spatstat package to create a PPP object, set observation window, marks, etc.
# bounding box around oromia state in ethiopia = our window for this data

oromia <- ETH_Adm_1[ETH_Adm_1$NAME_1 %in% "Oromia", ] # subsetting the spatial polygon
oromia_window <- owin(oromia@bbox[1,], oromia@bbox[2,]) # taking that data and wrapping it in class owin

# then define a ppp of the prevalence data
ETH_malaria_data_ppp <- ppp(ETH_malaria_data$longitude,
                            ETH_malaria_data$latitude,
                            marks = ETH_malaria_data$pf_pr, # prevalence? the OUTCOME we want to interpolate
                            window = oromia_window)

# marks = outcome variable we want to interpolate (predict, or fill in, basically)

# now plot different IDW results ----
idw(ETH_malaria_data_ppp, power=0.2, at="pixels") %>%
  plot(col=heat.colors(20), main="power = 0.2") 

# power represents the power function we want to use
# 'at' can be 'pixels' (where it generates estimates across a grid of pixels)
# 'at' can also be 'points' (where it interpolates values at every point using leave-one-out-cross validation)

idw(ETH_malaria_data_ppp, power=0.5, at="pixels") %>%
  plot(col=heat.colors(20), main="power = 0.5")

idw(ETH_malaria_data_ppp, power=1, at="pixels") %>%
  plot(col=heat.colors(20), main="power = 1")

idw(ETH_malaria_data_ppp, power=2, at="pixels") %>%
  plot(col=heat.colors(20), main="power = 2")

# 'power' is the power function. it determines which weight to assign per IDW method.
# 'at' is pixels where it generates estimates across a grid of pixels
# use CROSS VALIDATION to determine which power to use

powers <- seq(0.05, 2, 0.05)
mse_result <- NULL

# cross validation IDW
for (power in powers) {

  # the part that makes it cross validation is at="points" i think...  
  CV_idw <- idw(ETH_malaria_data_ppp, power=power, at="points") #CV_idw = cross validation
  # we want to minimize MSE (a measure of error)
  mse_result <- c(mse_result, mse(ETH_malaria_data_ppp$marks, CV_idw))
}

# see which produced the lowest error
optimal_power <- powers[which.min(mse_result)]
optimal_power
powers %>% plot(mse_result)

# plot observed vs. expected with the optimal power value
CV_idw_opt <- idw(ETH_malaria_data_ppp, power = optimal_power, at= "points")
ETH_malaria_data_ppp$marks %>% plot(CV_idw_opt)

# Convert to a raster
ETH_malaria_data_idw_raster <- idw(ETH_malaria_data_ppp, power = optimal_power, at="pixels") %>%
  raster(crs = crs(ETH_Adm_1))

colPal <- colorNumeric(tim.colors(), ETH_malaria_data_idw_raster[], na.color = NA)

# plot 
basemap %>% addRasterImage(ETH_malaria_data_idw_raster,
                           col = colPal,
                           opacity = 0.7) %>%
  addLegend(pal = colPal, values = ETH_malaria_data_idw_raster[])

# Kriging ----
## to do kriging in R, we need the geoR package
## remember, kriging is another form of interpolation for point referenced/level data (eg air quality)
## it is an alternative to inverse distance weighting (IDW)

# create GEODATA object using geoR package, import only select columns
ETH_malaria_data_geo <- ETH_malaria_data %>%
  select(longitude, latitude, pf_pr) %>% as.geodata()

# plot a summary plot using ther Lowes parameter
ETH_malaria_data_geo %>% plot(lowes=T)
# the Lowes option on 'plot' gives us lowes curves for the relationship between x and y

# should assess whether there is a first-order trend BEFORE kriging

# these particular plots show no trends

# you can add trend = "1st" or trend = "2nd" to plot command if there was evidence of trends in these plots

# generate and plot variogram (part of actually doing kriging) ---- 

MaxDist <- (ETH_malaria_data %>% select(longitude,latitude) %>% dist() %>% max()) /2

# dist = distance matrix computation, another part of kriging

VarioCloud <- ETH_malaria_data_geo %>% variog(option="cloud", max.dist= MaxDist)

VarioCloud %>% plot() # all pairwise comparisons

## rule of thumb: limit variogram estimation to half of the maximum interpoint distance
# can include trend surface using trend = "1st" or trend = "2nd"

# To make it easier to interpret, we can bin points by distance
Vario <- ETH_malaria_data_geo %>%
  variog(max.dist = MaxDist)

Vario %>% plot() # all pairwise comparisons
# x axis = distance classes, y axis = semivariance?

# can also change the way variograms are constructed
Vario <- ETH_malaria_data_geo %>%
  variog(max.dist = MaxDist, uvec = seq(0.01, MaxDist, 0.2)) 

# Just be careful not to have too few pairs of points in any distance class.

Vario$n # number in each bin (distance classes?)

# What is the minimum? A rule of thumb is 30 in each bin
min(Vario$n)

# Plot
Vario %>% plot(pch=16)
# y axis: semivariance
# x axis: distance

# fitting variogram model using different models ----
# We can now fit variogram model by minimized least squares using different covariance models.
# In this case we are just going to use a ‘spherical’ and ‘exponential’ model.

VarioMod_sph <- Vario %>% variofit(cov.model = "sph") # spherical model to fit variogram model

VarioMod_exp <- Vario %>% variofit(cov.model = "exp") # exponential model to fit variogram model

# plot results ----
Vario %>% plot(pch=16)
VarioMod_sph %>% lines(col="blue",lwd=2) # spherical model
VarioMod_exp %>% lines(col="red",lwd=2) # exponential model

# Get summaries of the fits ----
summary(VarioMod_sph) #spherical model
summary(VarioMod_exp) #exponential model

# Create prediction grid ----
IDW <- ETH_malaria_data_ppp %>% idw(power=0.2, at="pixels")
pred_grid_x <- rep(IDW$xcol,length(IDW$yrow))
pred_grid_y <- rep(IDW$yrow,length(IDW$xcol)) %>% sort()
pred_grid <- cbind(pred_grid_x, pred_grid_y)

# Now krig to those points ----
KrigPred <- ETH_malaria_data_geo %>% 
  krige.conv(loc=pred_grid,
             krige=krige.control(obj.model=VarioMod_sph))

# Visualize predictions ----
KrigPred %>% image(col=heat.colors(50))

# as raster
KrigPred_raster <- data.frame(x=pred_grid_x,
                              y=pred_grid_y,
                              z=KrigPred$predict) %>% rasterFromXYZ()

KrigPred_raster %>% plot()
ETH_malaria_data %>% select(longitude, latitude) %>% points(cex = ETH_malaria_data$pf_pr * 10)

# generate cross-validated predictions ----
xvalid_result <- ETH_malaria_data_geo %>% xvalid(model = VarioMod_sph)

# By default it xvalidates point by point....

# Plot observed versus expected ----
plot(xvalid_result$data,xvalid_result$predicted, asp=1)
abline(0,1)

# notice that some of the kriging values are less than 0.
# this can't be true since we're modeling probabilities
# so we should apply a transformation to our data before kriging, then back-transform our results

# we can use a logit regression

ETH_malaria_data <- ETH_malaria_data %>%
  mutate(pf_pr_adj = pf_pr + 0.001,
         pf_pr_logit = logit(pf_pr_adj))

ETH_malaria_data_geo_logit <- ETH_malaria_data %>% select(longitude, latitude, pf_pr_logit) %>%
  as.geodata()

# Fit (spherical) variogram
Vario_logit <- ETH_malaria_data_geo_logit %>% variog(max.dist = MaxDist)
VarioMod_sph_logit <- Vario_logit %>% variofit(cov.model = "sph")

# Get cross-validated (CV) kriged predictions
xvalid_result_logit <- ETH_malaria_data_geo_logit %>% xvalid(model = VarioMod_sph_logit)
xvalid_result_inv_logit <- xvalid_result_logit$predicted %>% inv.logit()
