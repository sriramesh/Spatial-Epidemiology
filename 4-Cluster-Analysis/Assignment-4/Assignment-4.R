# load libraries ----

# visualization
library(rgdal)
library(raster)
library(ggplot2)
library(spatstat)
library(plotrix)
library(fields)
library(leaflet)
library(maptools)
library(RColorBrewer)
library(lattice)
library(geoR)
library(plotrix) 

# spatial data management and point process analysis
library(sp)
library(tidyverse)
library(gtools)

# Moran's I and spatial dependencies
library(spdep) # Spatial Dependence: Weighting Schemes, Statistics and Models
library(ape) # Analyses of Phylogenetics and Evolution
library(pgirmess) # Data Analysis in Ecology

# point processes
library(spatstat)
library(splancs) # K-function
library(smacpod) # Spatial scanning statistic
library(car) # contains a function for logistic transformation (log odds ratio) to make more normal


# load data -----
url_ethmalaria <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week1/Lab_files/Data/mal_data_eth_2009_no_dups.csv"
ETH_malaria_data <- read.csv(url_ethmalaria, header=T)

# checking out what it looks like - normal or not
ETH_malaria_data <- ETH_malaria_data %>%
  mutate(log_odds = car::logit(pf_pr))

# ETH_malaria_data %>% ggplot(aes(x=pf_pr)) + geom_histogram()

# ETH_malaria_data %>% ggplot(aes(x=log_odds)) + geom_histogram()


# Q1 create inverse distance matrix ----
ETH.dists <- ETH_malaria_data %>%
  select(latitude, longitude) %>% 
  bind_cols() %>% dist() %>% as.matrix()

ETH.dists %>% dim() # 203 x 203 matrix of distance between all sets of points

# take the inverse of the distance matrix values so that CLOSER values = LARGER weight
ETH.dists.inv <- 1/ETH.dists
diag(ETH.dists.inv) <- 0

# Q1 Calculate the maximum distance between points ----
maxDist <- ETH_malaria_data %>%
  select(latitude, longitude) %>% bind_cols() %>% dist() %>% max()
maxDist

# Q1 Plot data in correlogram. Include a relevant title and clear axis labels ----
xy <- ETH_malaria_data %>%
  select(longitude, latitude) %>% bind_cols() %>% as.matrix()

pgi.cor <- correlog(coords=xy, z=ETH_malaria_data$log_odds,
                    method="Moran", nbclass=10) # the "distance class" columns

pgi.cor %>% plot(xlab = "Distance Classes", ylab = "Moran's I Statistic",
                 main = "Correlogram: Distance Classes v. Moran's I")
pgi.cor # distclass = midpoint for the bin

# Q1 Comment on what the correlogram reveals. Is there evidence of spatial autocorrelation? ----
# If so, over what spatial lags? Is the clustering positive or negative?

# Q2 Create different neighbor structures and assign weights to them ----

coords <- xy %>% coordinates() # set spatial coordinates to create a spatial object
IDs <- coords %>% as.data.frame() %>% row.names() # set row names as IDs

# k-nearest neighbors:
Neigh_nb <- knearneigh(coords, k=1, longlat = TRUE) %>% knn2nb(row.names=IDs)

dsts <- Neigh_nb %>% nbdists(coords) %>% unlist()
dsts %>% summary()

# maximum distance to provide at least one neighbor to each point (biggest nearest neighbor distance in the dataset)
max_1nn <- dsts %>% max()
max_1nn

# neighbors within maximum distance:
Neigh_kd1 <- coords %>% dnearneigh(d1=0, d2=max_1nn, row.names=IDs)
# neighbors within 2X maximum distance:
Neigh_kd2 <- coords %>% dnearneigh(d1=0, d2=2*max_1nn, row.names=IDs)

# list of neighbor structures
nb_1 <- list(d1=Neigh_kd1, d2=Neigh_kd2)
nb_1 %>% sapply(function(x) is.symmetric.nb(x, verbose=F, force=T))

# Does not always hold for k-nearest neighbours
nb_1 %>% sapply(function(x) n.comp.nb(x)$nc)

# plot neighbors comparing the two distances
par(mfrow=c(2,1), mar= c(1, 0, 1, 0))
xy %>% plot(pch=16)
Neigh_kd1 %>% plot(coords, col="green",add=T)
xy %>% plot(pch=16)
Neigh_kd2 %>% plot(coords,col="green", add=T)

weights <- Neigh_kd1 %>% nb2listw(style="W")
weights

# Q2 Run the Moran's I test ----

# row standardized binary weights, using minimum distance for one neighbor
moran.test(ETH_malaria_data$log_odds, listw=weights)  #using row standardised weights

# Q3 Compare your results from each approach (correlogram vs. distance-based neighbor function) ----


set.seed(1234)
bperm <- moran.mc(ETH_malaria_data$log_odds, listw=weights,nsim=999)
bperm

# Plot simulated test statistics
par(mfrow=c(1,1), mar= c(5, 4, 4, 2))
hist(bperm$res, freq=T, breaks=20, xlab="Moran's I statistics",
     main = "Distribution of Moran's I statistics over 999 simulations")
abline(v=0.27, col="red")

# Local Spatial Autocorrelation test: Local Moran's I test ----

# First calculate the local Moran's I around each point (binary weights, not row-standardized weights)
I <- ETH_malaria_data$log_odds %>% localmoran(weights) # "spdep" package

# Print 'LISA' for each point
Coef <- data.frame(I[IDs,], row.names=row.names(coords), check.names=FALSE) %>%
  printCoefmat()

# Plot the spatial data against its spatially lagged values (the weighted mean of its neighbors)
nci <- ETH_malaria_data$log_odds %>%
  moran.plot(listw=weights,
             xlab="Log prevalence",
             ylab="Spatially lagged log prev",
             labels=T, pch=16, col="grey")
text(x=c(3,3,-5,-5),
     y=c(0.9,-1.9,0.9,-1.9),
     labels=c("High-High", "High-Low", "Low-High", "Low-Low"),
     cex=0.8)

# Map points that are local outliers in the plot
# find which points are statistically significant outliers
infl <- apply(nci$is.inf, 1, any)
sum(infl %in% T) # 22 true (more than would expect by chance)

# !! what is happening here -----

x <- ETH_malaria_data$log_odds
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=T)
wx <- lag(weights, ETH_malaria_data$log_odds)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)), labels=c("L", "H"), include.lowest=T)
lhlh <- interaction(lhx, lhwx, infl, drop=T)

names<-rep("none", length(lhlh))
names[lhlh=="L.L.TRUE"]<-"LL"
names[lhlh=="H.L.TRUE"]<-"HL"
names[lhlh=="L.H.TRUE"]<-"LH"
names[lhlh=="H.H.TRUE"]<-"HH"

# map local clusters ----
ETH_malaria_localM <- cbind(xy,names) %>% as.data.frame()
colnames(ETH_malaria_localM) <- c("longitude", "latitude", "names")

ETH_malaria_localM <- ETH_malaria_localM %>%
  mutate(latitude = lapply(latitude, function(x) as.numeric(as.character(x))),
         longitude = lapply(longitude, function(x) as.numeric(as.character(x))))

Pal <- c( "cyan4","coral4","coral","cyan","lightgrey")
factpal <- Pal %>% colorFactor(names)

ETH_malaria_localM %>% leaflet() %>% addTiles() %>%
  addCircleMarkers(~longitude, ~latitude, fillOpacity=1,
                   color= ~factpal(names), radius=4, stroke=TRUE, weight=1) %>% 
  addLegend(pal = factpal, values = ~names, title="Class")

