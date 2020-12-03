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

# load data
nairobi_cases <- read.csv("cases_nairobi.csv")

# boundaries for kenya, admin 0 level
KEN_Adm0 <- raster::getData('GADM',country='KEN',level=0)


# Q1 ----

# 1 generate a relative risk density estimate ----

nairobi_cases <- read.csv("cases_nairobi.csv")
nairobi_cases_spdf <- SpatialPointsDataFrame(coords = nairobi_cases[,c("lng", "lat")], data = nairobi_cases[,c("X", "case")])

# 2 plot kernel density estimates at different bandwidths ----

# define window
ken_owin <- owin(xrange=range(nairobi_cases$lng), yrange=range(nairobi_cases$lat))

# make ppp object
nairobi_cases_ppp <- ppp(nairobi_cases$lng, nairobi_cases$lat, window = ken_owin, 
                         marks=as.factor(nairobi_cases$case))
case_density <- density(nairobi_cases_ppp)
plot(density(nairobi_cases_ppp,0.02), main = "Bandwidth of 0.02")
plot(density(nairobi_cases_ppp,0.1), main = "Bandwidth of 0.1")
plot(density(nairobi_cases_ppp,bw.ppl), main = "CV-based bandwidth selection")

# 3 calculate risk ratios ----
nairobi_cases_ppp <- ppp(nairobi_cases$lng, nairobi_cases$lat, 
                         window = ken_owin, 
                         marks=as.factor(nairobi_cases$case))

rel_risk_est <- nairobi_cases_ppp %>% relrisk(relative = T)
plot(rel_risk_est, main="Relative Risk Estimate")

# 4 map risk ratio rasters ontop of leaflet base map ----
rel_risk_raster <- raster(rel_risk_est, crs = crs(KEN_Adm0))
pal <- colorNumeric(palette=tim.colors(64), domain=values(rel_risk_raster), na.color = NA)

basemap <- leaflet() %>% addProviderTiles("CartoDB.Positron")
basemap %>% addRasterImage(rel_risk_raster, opacity=0.6, colors = pal) %>%
  addLegend(pal = pal, values = values(rel_risk_raster), 
            title = "Relative Risk Ratios")


# Q2 ----

# 1 load data ----
HK <- read.csv("tanzania_uganda_hkprev.csv")
TZA_Adm_1 <- raster::getData("GADM", country="TZA", level=1) # used tanzania instead of uganda


# 2 Calculate the best IDW scenario by testing different powers ----

# set window
mwanza <- TZA_Adm_1[TZA_Adm_1$NAME_1 %in% "Mwanza", ]
mwanza_window <- owin(mwanza@bbox[1,], mwanza@bbox[2,])

# define ppp object
TZA_hookworm_data_ppp <- ppp(HK$x, HK$y, marks = HK$Hookworm_prev_perc, window = mwanza_window)

# set parameters for 4 panel display
par(mfrow=c(2,2))
idw(TZA_hookworm_data_ppp, power=0.2, at="pixels") %>% plot(col=heat.colors(20), main="Power = 0.2")
idw(TZA_hookworm_data_ppp, power=0.5, at="pixels") %>% plot(col=heat.colors(20), main="Power = 0.5")
idw(TZA_hookworm_data_ppp, power=1, at="pixels") %>% plot(col=heat.colors(20), main="Power = 1")
idw(TZA_hookworm_data_ppp, power=2, at="pixels") %>% plot(col=heat.colors(20), main="Power = 2")

# determine optimal power function using for loop
powers <- seq(0.05, 2, 0.05)
mse_result <- NULL

for (power in powers) {
  
  CV_idw <- TZA_hookworm_data_ppp %>% idw(power=power, at="points")
  mse_result <- c(mse_result, mse(TZA_hookworm_data_ppp$marks, CV_idw))
}

optimal_power <- powers[which.min(mse_result)]
optimal_power


# 3 Map the "best" IDW scheme raster files on top of leaflet basemap ----

# Convert to a raster
TZA_hookworm_data_idw_raster <- TZA_hookworm_data_ppp %>%
  idw(power = optimal_power, at="pixels") %>% raster(crs = crs(TZA_Adm_1))
# define color ramp
colPal <- colorNumeric(tim.colors(), TZA_hookworm_data_idw_raster[], na.color = NA)
# layer map
basemap %>% addRasterImage(TZA_hookworm_data_idw_raster, col = colPal, opacity = 0.7) %>%
  addLegend(pal = colPal, values = TZA_hookworm_data_idw_raster[],
            title = "IDW Values")



# 4 plot kriged raster for the chosen window (Mwanza region in TZA) and the prevalence points ----
# create geodata object
TZA_hookworm_data_geo <- HK %>% select(x, y, Hookworm_prev_perc) %>% as.geodata()
# generate distance matrix
MaxDist <- (HK %>% select(x,y) %>% dist() %>% max()) /2
VarioCloud <- TZA_hookworm_data_geo %>% variog(option="cloud", max.dist= MaxDist)
# bin by distance
Vario <- TZA_hookworm_data_geo %>% variog(max.dist = MaxDist, uvec = seq(0.01, MaxDist, 0.2)) 
# fit variogram to each model
VarioMod_sph <- Vario %>% variofit(cov.model = "sph")
VarioMod_exp <- Vario %>% variofit(cov.model = "exp")
# plot each model
Vario %>% plot(pch=16)
VarioMod_sph %>% lines(col="blue",lwd=2)
VarioMod_exp %>% lines(col="red",lwd=2)
# summarize each model
summary(VarioMod_sph)
summary(VarioMod_exp)
# test if spherical model has lower sum of squares than exponential model
(summary(VarioMod_sph)$sum.of.squares < summary(VarioMod_exp)$sum.of.squares)[[1]]
# create prediction grid
IDW <- TZA_hookworm_data_ppp %>% idw(power = optimal_power, at="pixels")
par(mfrow=c(1,1))
pred_grid_x <- rep(IDW$xcol,length(IDW$yrow))
pred_grid_y <- rep(IDW$yrow,length(IDW$xcol)) %>% sort()
pred_grid <- cbind(pred_grid_x, pred_grid_y)
# krig to those points
KrigPred <- TZA_hookworm_data_geo %>%
  krige.conv(loc=pred_grid, krige=krige.control(obj.model=VarioMod_sph))
# plot kriged points
KrigPred_raster <- data.frame(x=pred_grid_x, y=pred_grid_y, z=KrigPred$predict) %>% rasterFromXYZ()
KrigPred_raster %>% plot(main = "Kriged raster for Mwanza window and prevalence points")
HK %>% select(x,y) %>% points(cex = .5)

# 5 Compare methods using mean squared error ----
# cross validate idw values
CV_idw_opt <- TZA_hookworm_data_ppp %>% idw(power = optimal_power, at= "points")
CV_kriging_opt <- TZA_hookworm_data_geo %>% xvalid(model = VarioMod_sph)

mse(CV_idw_opt, HK$Hookworm_prev_perc)
mse(kriged_result$predicted, HK$Hookworm_prev_perc)

# 6 Visualize where predictions from IDW differ to kriging ----
IDW_raster <- IDW %>% raster()
plot(IDW_raster - KrigPred_raster)

# 7 Include a trend surface to see if that improve kriging estimates -----

# Fit (spherical) variogram with 1st order trend
Vario_trend <- TZA_hookworm_data_geo %>% variog(max.dist = MaxDist, trend="1st")
Vario_trend_sph <- Vario_trend %>% variofit(cov.model = "sph")

# Get kriged values with 1st order trend
kriged_result_trend <- TZA_hookworm_data_geo %>% xvalid(model = Vario_trend_sph)
kriged_result_inv_logit_trend <- kriged_result_trend$predicted %>% inv.logit() 

# Now compare
CV_kriging_opt$predicted %>% inv.logit() %>% mse(HK$Hookworm_prev_perc)
kriged_result_inv_logit_trend %>% mse(HK$Hookworm_prev_perc)  
