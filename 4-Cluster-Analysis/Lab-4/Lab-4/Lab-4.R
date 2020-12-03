## OBJECTIVE: intro to methods for identifying CLUSTERS in R ----

# load libraries for visualization ----
library(rgdal)
library(raster)
library(ggplot2)
library(spatstat)
library(plotrix)
library(fields)

library(leaflet)
#library(plotGoogleMaps)
library(maptools)
library(RColorBrewer)
library(lattice)
library(geoR)
library(plotrix) 

# load libraries for spatial data management and point process analysis
library(sp)
library(tidyverse)
library(gtools)

# Moran's I and spatial dependencies
library(spdep) # Spatial Dependence: Weighting Schemes, Statistics and Models
library(ape) # Analyses of Phylogenetics and Evolution
library(pgirmess) # Data Analysis in Ecology

# Attach libraries for point processes
library(spatstat)
library(splancs) # K-function
library(smacpod) # Spatial scanning statistic

# load malaria prevalence data from Burkina Faso and visualize it using leaflet
# take a look. is there evidence of spatial clustering?

# Open BF malaria data
malariaprev_url <- "https://www.dropbox.com/s/bfs3pinxe1lvvxr/data_bf2_binomial.csv?dl=1"
BF_malaria_data <- malariaprev_url %>% read.csv(header=T)

# load admin 1 levels for Burkina Faso
BF_Adm_1 <- raster::getData("GADM", country="BFA", level=1)
# set crs
proj4string(BF_Adm_1) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Calculate prevalence
BF_malaria_data <- BF_malaria_data %>% mutate(prevalence = positives/examined)

# put it on leaflet
pal = colorNumeric("Oranges", BF_malaria_data$prevalence)

BF_malaria_data %>%
  leaflet() %>% addTiles() %>% addCircleMarkers(~longitude, ~latitude, fillOpacity=1,
                                                fillColor= ~pal(prevalence), radius=~prevalence*10,
                                                stroke=TRUE, weight=1) %>%
  addLegend(pal = pal, values = ~prevalence)

# yes, this looks like it contains global spatial autocorrelation. we should test for it

## GLOBAL Spatial Autocorrelation
# (1) Calculate Moran's I using a distance-based matrix
# process: transform data if needed -> dist matrix -> inverse of dist matrix -> weight matrix -> moran's I
## 1 look at the distribution of the prevalence column
# are they normally distributed? (if so, complete spatial randomness. if not, clustering, maybe)
BF_malaria_data %>% ggplot(aes(x=prevalence)) + geom_histogram()
## honestly the histogram kinda looks normally distributed...?
## i guess not? professor says it looks skewed.
### so if the data is super skewed, we need to transform them before running Moran's I
### in this case we'll run a logit
### Moran's I produces a comparison that is close to the normal distribution

library(car) # contains a function for logistic transformation (log odds ratio) to make more normal
BF_malaria_data <- BF_malaria_data %>%
  mutate(log_odds = logit(prevalence))
BF_malaria_data %>% ggplot(aes(x=log_odds)) + geom_histogram() ## i guess this is better!

# now generate a distance matrix. calc dist between all points ----
BF.dists <- BF_malaria_data %>% select(latitude, longitude) %>% 
  bind_cols() %>% dist() %>% as.matrix()
BF.dists %>% dim() # 109 x 109 matrix of distance between all sets of points

# take the inverse of the distance matrix values so that CLOSER values = LARGER weight ----
BF.dists.inv <- 1/BF.dists
diag(BF.dists.inv) <- 0

# Computes Moran's I autocorrelation coefficient of x, giving a matrix of WEIGHTS
# (here based on distance! can be based on other things...)
Moran.I(BF_malaria_data$log_odds, BF.dists.inv) # Moran.I comes from the "ape" package

# (2) Create a correlogram to explore Moran's I over different spatial lags ----

## calculate the max distance between points
maxDist <- BF_malaria_data %>% select(latitude, longitude) %>% 
  bind_cols() %>% dist() %>% max()
maxDist

xy <- BF_malaria_data %>% select(longitude, latitude) %>% 
  bind_cols() %>% as.matrix()

# calculate correlogram
pgi.cor <- correlog(coords=xy, z=BF_malaria_data$log_odds, method="Moran", nbclass=10) # "pgirmess" package
pgi.cor %>% plot()
pgi.cor # distclass = midpoint for the bin

# based on the correlogram, over what SPATIAL LAGS are there evidence of spatial autocorrelation?
# is this positive or negative spatial autocorrelation/clustering?

# compare correlogram to a semivariogram

# variogram
# make geodata object
BF_malaria_data_geo <- BF_malaria_data %>% 
  select(longitude, latitude, log_odds) %>% as.geodata()
# bin and plot variogram (10 bins)
Vario <- BF_malaria_data_geo %>% 
  variog(max.dist=7.53,uvec=seq(0.4121237,7.1595572,l=10))

# plot
par(mfrow = c(2,1))
Vario %>% plot()
pgi.cor %>% plot()

# (3) Calculate Moran's I using a BINARY distance matrix
## for this approach we will create sets of neighbors based on their proximity to each other
## this can be used with point data but is especially useful for areal data

coords <- xy %>% coordinates() # set spatial coordinates to create a spatial object
IDs <- coords %>% as.data.frame() %>% row.names() # set row names as IDs

# in this approach we chose a distance d such that:
# pairs of points with distances LESS THAN d = neighbors (neighbors = 1)
# pairs of points with distances GREATER THAN d = not neighbors (neighbors = 0)

# using the "spdep" package...
Neigh_nb <- knearneigh(coords, k=1, longlat = TRUE) %>%
  knn2nb(row.names=IDs)

# assigns at least one neighbor to each, and calculates the distances between
# returns the distance between nearest neighbors for each point
dsts <- Neigh_nb %>% nbdists(coords) %>% unlist()
dsts %>% summary()


# We create different neighbor structures based upon distance
# maximum distance to provide at least one neighbor to each point
max_1nn <- dsts %>% max()
max_1nn

# neighbors within maximum distance:
Neigh_kd1 <- coords %>% dnearneigh(d1=0, d2=max_1nn, row.names=IDs)
# neighbors within 2X maximum distance:
Neigh_kd2 <- coords %>% dnearneigh(d1=0, d2=2*max_1nn, row.names=IDs)

# list of neighbor structures
nb_1 <- list(d1=Neigh_kd1, d2=Neigh_kd2)
nb_1 %>% sapply(function(x) is.symmetric.nb(x, verbose=F, force=T))

# Checks for symmetry (i.e. if i is a neighbor of j, then j is a neighbor of i).
# Does not always hold for k-nearest neighbours
nb_1 %>% sapply(function(x) n.comp.nb(x)$nc)

# number of disjoint connected subgraphs
# plot neighbors comparing the two distances
par(mfrow=c(2,1), mar= c(1, 0, 1, 0))
xy %>% plot(pch=16)
Neigh_kd1 %>% plot(coords, col="green",add=T)
xy %>% plot(pch=16)
Neigh_kd2 %>% plot(coords,col="green", add=T)

# assign weights to neighbor list
# we will use the neighbor structure with all neighbors within the maximum neighbor distance between any two points
# row standardized binary weights, using minimum distance for one neighbor
weights <- Neigh_kd1 %>% nb2listw(style="W")
# "B" is simplest binary weights
weights

# using this weights matrix, we can now run the Moran's I test on the logit transformed prevalence using the neighborhood matrix ----
moran.test(BF_malaria_data$log_odds, listw=weights)  #using row standardised weights

# sampling for monte carlo simulation
# monte carlo simulation to get the p-value
set.seed(1234)
bperm <- moran.mc(BF_malaria_data$log_odds, listw=weights,nsim=999)
bperm

#statistic = 0.15, observed rank = 1000, p-value = 0.001

# Plot simulated test statistics
par(mfrow=c(1,1), mar= c(5, 4, 4, 2))
hist(bperm$res, freq=T, breaks=20, xlab="Simulated Moran's I")
abline(v=0.15, col="red")

# Moran's I for polygons/areal data ----
# we can now also take a look at running Moran's I for areal data (polygons)
# can use a dataset on leukemia from New York (Turnbull et al 1990)

# geo json file:
nydata_url <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week7/Lab_files/nydata.geojson"
nydata <- rgdal::readOGR(nydata_url)

# let's take a look at the data
head(nydata@data)

# for now, we are only interested in seeing if there is global clustering in the area-level case incidence
## Cases variable gives us the estimated number of cases per area
## We also need to consider the population in each area
## areas with higher populations are more likely to have more cases, just due to pop size 
## so we will first normalize this by creating an incidence var per 1000 people

nydata$inc_per_1000 <- (nydata$Cases / nydata$POP8) * 1000

# as these are areas not points, we will NOT use distance to define the neighbors
# instead, we will figure out which polygons are directly touching each other along a boundary or boundary point

# contiguity neighbors - all that share a boundary point
nydata_nb <- nydata %>% poly2nb()  #queen contiguity

nydata_nbr <- nydata %>% poly2nb(queen=F)  #rook contiguity

# coordinates
coords_ny <- nydata %>% coordinates()

# view and compare the neighbors
par(mfrow=c(1,2))
par(mar=c(1,1,1,1)) # ensures figure margins are not too large
nydata %>% plot()
nydata_nb %>% plot(coords_ny,col="blue",add=T) # queen contiguity
nydata %>% plot()
nydata_nbr %>% plot(coords_ny,col="green",add=T) # rook contiguity

## set weights - contiguity
# weights style W - row standardized
nydata_w <- nydata_nb %>% nb2listw()
nydata_w

# weights style B - binary
nydata_wB <- nydata_nb %>% nb2listw(style="B")
nydata_wB

## moran's tests of global spatial autocorrelation -- for polygon/areal data!
moran.test(nydata$inc_per_1000, listw=nydata_w)  # using row standardized

# LOCAL Spatial Autocorrelation ----

# first calculate local moran's I around each point based on spatial weights object
# this is BINARY, based on at least 1 neighbor
# use spdep package
# First calculate the local Moran's I around each point based on the spatial weights object (binary based on at least one neighbor)
I <- localmoran(BF_malaria_data$log_odds, weights)                         # "spdep" package

# Print 'LISA' for each point
Coef <- I[IDs,] %>% data.frame() %>% printCoefmat(row.names=row.names(coords), check.names=F)

# Plot the spatial data against its spatially lagged values (the weighted mean of its neighbors)                         
nci <- moran.plot(BF_malaria_data$log_odds, listw=weights, 
                xlab="Log prevalence",
                ylab="Spatially lagged log prev",
                labels=T, pch=16, col="grey")
text(c(3,3, -5,-5), c(0.9, -1.9,0.9,-1.9),
     c("High-High", "High-Low", "Low-High", "Low-Low"),
     cex=0.8)

# map points that are local outliers in the plot
# find which points are statistically significant outliers
# Map points that are local outliers in the plot
infl <- apply(nci$is.inf, 1, any) # find which points are statistically significant outliers
sum(infl %in% T) #13 true (12% - more than would expect by chance)

x <- BF_malaria_data$log_odds
lhx <- x %>% cut(breaks=c(min(x), mean(x), max(x)),
                 labels=c("L", "H"), include.lowest=T)

wx <- weights %>% stats::lag(BF_malaria_data$log_odds)

lhwx <- wx %>%
  cut(breaks=c(min(wx), mean(wx), max(wx)),
      labels=c("L", "H"), include.lowest=T)

lhlh <- interaction(lhx,lhwx,infl,drop=T)

names <- rep("none", length(lhlh))

names[lhlh %in% "L.L.TRUE"] <- "LL"
names[lhlh %in% "H.L.TRUE"] <- "HL"
names[lhlh %in% "L.H.TRUE"] <- "LH"
names[lhlh %in% "H.H.TRUE"] <- "HH"

BF_malaria_localM <- cbind(xy, names) %>% as.data.frame()
colnames(BF_malaria_localM) <- c("longitude", "latitude", "names")

BF_malaria_localM[c("longitude", "latitude")] <- lapply( BF_malaria_localM[c("longitude", "latitude")], function(x) as.numeric(as.character(x)) )

factpal <- colorFactor(c( "cyan4","coral4","coral","cyan","lightgrey"), names)

## leaflet of the LOCAL clusters
BF_malaria_localM %>% leaflet() %>%
  addTiles() %>% addCircleMarkers(~longitude, ~latitude, fillOpacity=1,
                                  color= ~factpal(names), radius=4, stroke=TRUE, weight=1) %>% 
  addLegend(pal = factpal, values = ~names, title="Class")




