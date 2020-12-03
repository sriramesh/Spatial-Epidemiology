# How could you compare how well the best fitting IDW performs versus kriging? ----

## cross validation measures could help determine fit with IDW v. fit with kriging
## remember, both are types of interpolation that we use when are working with point prevalence data

# Which appears to be more accurate? -----

# 1 perform IDW by...getting window using Oromia boundaries, making ppp object with window, identifying the optimal power function

# oromia window
oromia <- ETH_Adm_1[ETH_Adm_1$NAME_1 %in% "Oromia", ]
oromia_window <- owin(oromia@bbox[1,], oromia@bbox[2,]) #make window using owin from geoR package

# make ppp object
ETH_malaria_data_ppp <- ppp(ETH_malaria_data$longitude,ETH_malaria_data$latitude,
                        marks=ETH_malaria_data$pf_pr,window=oromia_window)

# Get cross-validated values using the 'optimal' power
# cross-val idw with optimal power (CV_idw_opt)
CV_idw_opt <- ETH_malaria_data_ppp %>% idw(power=1.4, at="points")
# at="points" is for cross-validated!
# remember how we found the optimal power was 1.4 using that for loop? that's why it's 1.4

# 2 Now Krig
# First fit variogram ----
ETH_malaria_data <- ETH_malaria_data %>%
  mutate(pf_pr_adj = pf_pr + 0.001,
         pf_pr_logit = logit(pf_pr_adj))

ETH_malaria_data_geo_logit <- ETH_malaria_data %>%
  select(longitude, latitude, pf_pr_logit) %>% as.geodata()

# Fit (spherical) variogram
Vario_logit <- ETH_malaria_data_geo_logit %>% variog(max.dist = MaxDist)
VarioMod_sph_logit <- Vario_logit %>% variofit(cov.model = "sph")

# Get CV kriged predictions
xvalid_result_logit <- ETH_malaria_data_geo_logit %>% xvalid(model = VarioMod_sph_logit)
xvalid_result_inv_logit <- xvalid_result_logit$predicted %>% inv.logit() 

# Now compare
CV_idw_opt %>% mse(ETH_malaria_data$pf_pr)
xvalid_result_inv_logit %>% mse(ETH_malaria_data$pf_pr)

# the one with the lower mse = better fit

# Can you visualize where predictions from IDW differ to kriging? ----
# Here we can generate rasters of predictions from both methods and visualize the difference

IDW <- ETH_malaria_data_ppp %>% idw(power=0.2, at="pixels") # load IDW with optimal power
pred_grid_x <- rep(IDW$xcol,length(IDW$yrow))
pred_grid_y <- rep(IDW$yrow,length(IDW$xcol)) %>% sort()
pred_grid <- cbind(pred_grid_x,pred_grid_y) # create prediction grid

# Now krig to those points on the grid
KrigPred_logit <- ETH_malaria_data_geo_logit %>% krige.conv(loc=pred_grid,
                             krige=krige.control(obj.model=VarioMod_sph_logit))

# Create a raster of inv.logit values
KrigPred_logit_raster <- data.frame(x=pred_grid_x,
                                    y=pred_grid_y,
                                    z=inv.logit(KrigPred_logit$predict)) %>%
  rasterFromXYZ()

IDW_raster <- IDW %>% raster()
plot(IDW_raster - KrigPred_logit_raster)

# Does inclusion of a trend surface improve kriging estimates? ----

# Fit (spherical) variogram
# this is where we include the trend surface: trend = "1st"
Vario_logit_trend <- ETH_malaria_data_geo_logit %>%
  variog(max.dist = MaxDist, trend="1st")

VarioMod_sph_logit_trend <- Vario_logit_trend %>%
  variofit(cov.model = "sph")

# Get CV kriged predictions (using spherical variogram)
xvalid_result_trend_logit <- ETH_malaria_data_geo_logit %>%
  xvalid(model = VarioMod_sph_logit_trend)

xvalid_result_trend_inv_logit <- xvalid_result_trend_logit$predicted %>%
  inv.logit() 

# Now compare
xvalid_result_inv_logit %>% mse(ETH_malaria_data$pf_pr)  
xvalid_result_trend_inv_logit %>% mse(ETH_malaria_data$pf_pr)  
# Doesn't seem to improve things in this case since the mse is essentially the same

  