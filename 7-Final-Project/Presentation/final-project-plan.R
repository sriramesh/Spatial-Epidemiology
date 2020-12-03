# (A) Intro to dataset - Kenya hookworm/whipworm/ringworm data

# (B) Exploratory analysis
# DO: - Histogram of data (normal or not? logit transform or not?)
# DO: - Plot prevalence points onto leaflet, make comments on clustering or not
# - plan to use global Moran's I test and define weights using k-nearest neighbors
# - plan to use monte carlo simulation to confirm the value of the Moran's I
# - plan to use local moran's I and LISAs to test for local clustering

# (C) External Data
# 1) climatic data (altitude, mean temperature)
# 2) proximity to transport infrastructure (distance to railroads, distance to roads)

# geostatistical studies out of Bolivia and Thailand that prevalence of soil-transmitted helthminths is correlated with
# climatic factors, namely altitude and mean temperature, as well as land use, land cover, and land surface temperature
# one study out kenya finds that a school's proximity to lake victoria was correlated with distribution of soil-transmitted healthminths

# (D) Statistical tests you are going to use
# hypothesis: In Kenya, higher mean temperature, higher altitude, and closer proximity to both water bodies and transport infrastructre is associated with higher prevalence of soil-transmitted helminths
# generalized linear model with and without spatially correlated fixed effects to see if space can explain away any of the spatial pattern in the point-prevalence data
# global moran's I test

# (E) validation methods
# - v folds cross validation and MSE
