# load libraries ----
library(raster)
library(spaMM)
library(caret)
library(ggplot2)
library(ModelMetrics)
library(tidyverse)

# Get Ethiopia malaria data
malaria_url <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week1/Lab_files/Data/mal_data_eth_2009_no_dups.csv"
ETH_malaria_data <- malaria_url %>% read.csv(header=T)

# Get Biolclim layers for Ethiopia
bioclim_layers <- raster::getData('worldclim', var='bio', res=0.5, lon=38.7578, lat=8.9806) # lng/lat for Addis Ababa

# crop layers to make them easier to handle
ETH_Adm_1 <- raster::getData("GADM", country="ETH", level=1) # Admin boundaries
Oromia <- ETH_Adm_1 %>% subset(NAME_1 %in% "Oromia")
bioclim_layers_oromia <- bioclim_layers %>% raster::crop(Oromia)

# extract bioclim data as variables in dataset
ETH_malaria_data <- ETH_malaria_data %>%
  mutate(bioclim1 = raster::extract(bioclim_layers_oromia[[1]], ETH_malaria_data[,c("longitude", "latitude")]),
         bioclim2 = raster::extract(bioclim_layers_oromia[[2]], ETH_malaria_data[,c("longitude", "latitude")]))

# split data into 10 folds (default) ----
folds_list <- caret::createFolds(ETH_malaria_data$pf_pr)
folds_list[[1]]

# split data into 20 folds
# folds_list <- createFolds(ETH_malaria_data, k = 20)

# loop through folds ----
# before we loop through folds, lets create empty objects of things we want to keep 

cross_validated_prediction <- NULL
observed <- NULL

# Now we need to 1) loop through each fold, 2) fit a model to the training data of each fold, 2) predict to validation data 3) store validation stats

# Now loop through each fold 
for(fold in 1:length(folds_list)){
  
  training_data <- ETH_malaria_data[-folds_list[[fold]], ]
  validation_data <- ETH_malaria_data[folds_list[[fold]], ]

  # fitting the model to each training dataset of each fold
  fold_mod_spatial <- spaMM::fitme(cbind(pf_pos, examined - pf_pos) ~
                                     bioclim2 + Matern(1|latitude+longitude),
                                   data=training_data, family=binomial())
  
  # predict to the validation data
  x_valid_pred <- fold_mod_spatial %>% predict(validation_data)
  
  # Add to cross_validated_prediction
  cross_validated_prediction <- c(cross_validated_prediction, x_valid_pred)
  observed <- c(observed, validation_data$pf_pr)
}

# We now have a vector of cross validated predictions
# and their corresponding observed prevalence
# We can use this to create a scatter plot and/or estimate MSE

# plot
ggplot() + geom_point(aes(cross_validated_prediction, observed))

# calculate MSE
cross_validated_prediction %>% mse(observed)

# If you want to compare models in terms of cross-validated prediction error,
# you can repeat this process but swap the model/covariates.



