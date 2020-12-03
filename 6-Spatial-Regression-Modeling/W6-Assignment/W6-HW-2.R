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

BCG_vaccination_UGA <- read.csv("BCG_vaccination_UGA.csv")

# 2 run glm specifications ----
glm_mod_1 <- glm(cbind(number_positive, numer_examined - number_positive) ~
                   ntl_2013 + dist_to_water + dist_to_roads + dist_to_railroads,
                 data=BCG_vaccination_UGA, family=binomial())
glm_mod_1 %>% summary()

glm_mod_2 <- glm(cbind(number_positive, numer_examined - number_positive) ~
                   dist_to_water + dist_to_railroads + ntl_2013,
                 data=BCG_vaccination_UGA, family=binomial())
glm_mod_2 %>% summary()

glm_mod_3 <- glm(cbind(number_positive, numer_examined - number_positive) ~
                   ntl_2013,
                 data=BCG_vaccination_UGA, family=binomial())
glm_mod_3 %>% summary()

# 2 check fitted values for last 2 models -----
ggplot() + geom_point(aes(glm_mod_2$fitted, BCG_vaccination_UGA$prevalence))

ggplot() + geom_point(aes(glm_mod_3$fitted, BCG_vaccination_UGA$prevalence))

# plot correlograms for last 2 models ----
nbc <- 10
cor_r_2 <- pgirmess::correlog(coords=BCG_vaccination_UGA %>% select(lng, lat),
                            z=glm_mod_2$residuals,
                            method="Moran", nbclass=nbc)
cor_r_2 # definite clustering up to 1.7 decimal degrees
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
  ggtitle("Correlogram, before spatial effect included")

cor_r_3 <- pgirmess::correlog(coords=BCG_vaccination_UGA %>% select(lng, lat),
                            z=glm_mod_3$residuals,
                            method="Moran", nbclass=nbc)
cor_r_3 # definite clustering up to 1.7 decimal degrees
correlogram_3 <- cor_r_3 %>% as.data.frame
correlogram_3$variable <- "residuals_glm" 

correlogram_3 %>% subset(variable %in% "residuals_glm") %>% ggplot(aes(dist.class, coef)) + 
  geom_hline(yintercept = 0, col="grey") +
  geom_line(col="steelblue") + 
  geom_point(col="steelblue") +
  xlab("distance") + 
  ylab("Moran's coefficient")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Correlogram, before spatial effect included")

# adding a spatial component to our model ----
glm_mod_spatial <- spaMM::fitme(cbind(number_positive, numer_examined - number_positive) ~
                                    dist_to_water + dist_to_railroads + ntl_2013 +
                                    Matern(1|lat+lng),
                                  data=BCG_vaccination_UGA, family=binomial())
glm_mod_spatial %>% summary()

# re-testing for residual spatial autocorrelation ----
nbc <- 10
cor_r <- pgirmess::correlog(coords = BCG_vaccination_UGA %>% select(lng, lat),
                            z = residuals(glm_mod_2_spatial),
                            method="Moran", nbclass=nbc)
cor_r

correlograms <- cor_r %>% as.data.frame
correlograms$variable <- "residuals_glm" 

correlograms %>% subset(variable %in% "residuals_glm") %>% ggplot(aes(dist.class, coef)) + 
  geom_hline(yintercept = 0, col="grey") +
  geom_line(col="steelblue") + 
  geom_point(col="steelblue") +
  xlab("distance") + 
  ylab("Moran's coefficient")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Correlogram, after spatial effect included")

# 3 plot prevalence over a grid ----

# choose window
UGA_Adm_1 <- raster::getData("GADM", country="UGA", level=1) # Uganda boundaries
proj4string(UGA_Adm_1) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ')
Kampala <- UGA_Adm_1 %>% subset(NAME_1 %in% "Kampala")

# crop ntl raster to window
ntl_2013_kampala <- ntl_2013_raster %>% crop(Kampala)
ntl_2013_kampala %>% plot()
Kampala %>% lines()

# predict over ntl raster
pred_raster <- ntl_2013_kampala
names(pred_raster) <- c("ntl_2013")

par(mar=c(1,1,1,1))

predicted_risk <- raster::predict(pred_raster, glm_mod_3, type='response')
predicted_risk %>% plot()

predicted_risk_masked <- predicted_risk %>% mask(Kampala)
predicted_risk_masked %>% plot()

# prediction using spatial GLM onto a stack of covariates ----

# Create an empty raster with the same extent and resolution as the bioclimatic layers
latitude_raster <- longitude_raster <- raster(nrows = ntl_2013_kampala %>% nrow(),
                                              ncols = ntl_2013_kampala %>% ncol(),
                                              ext = ntl_2013_kampala %>% extent())


# Change the values to be latitude and longitude respectively
longitude_raster[] <- coordinates(longitude_raster)[,1]
latitude_raster[] <- coordinates(latitude_raster)[,2]

# Now create a final prediction stack of the 4 (3?) variables we need
pred_stack <- stack(ntl_2013_kampala, longitude_raster, latitude_raster)

# Rename to ensure the names of the raster layers in the stack match those used in the model
names(pred_stack) <- c("ntl_2013", "lng", "lat")
pred_stack %>% plot()

# Make predictions using the stack of covariates and the spatial GLM, then clip and mask to Oromia state
predicted_prevalence_raster <- pred_stack %>% predict(glm_mod_2_spatial)
predicted_prevalence_raster %>% plot()
Kampala %>% lines()
predicted_prevalence_raster_kampala <- predicted_prevalence_raster %>% mask(Kampala)
predicted_prevalence_raster_kampala %>% plot()

# CV model validation -----
# take 20% to act as validation set
set.seed(1)
validation_rows <- 1:nrow(BCG_vaccination_UGA) %>% sample(40)
BCG_vaccination_UGA_train <- BCG_vaccination_UGA[-validation_rows,] # training sets
BCG_vaccination_UGA_valid <- BCG_vaccination_UGA[validation_rows,] # validation, will predict these using training sets

# Fit spatial GLM on only the training data
glm_mod_2_spatial_validation <- spaMM::fitme(cbind(number_positive, numer_examined - number_positive) ~
                                               dist_to_water + dist_to_railroads + Matern(1|lat+lng),
                                             data=BCG_vaccination_UGA_train, family=binomial())

# Use spatial GLM to predict the 20% validation data and compare with actual data to see how well the model did
predictions_validation <- glm_mod_2_spatial_validation %>% predict(BCG_vaccination_UGA_valid)

ggplot() + geom_point(aes(predictions_validation %>% as.vector(), BCG_vaccination_UGA_valid$prevalence))

# Calculate mse
predictions_validation %>% mse(BCG_vaccination_UGA_valid$prevalence)

