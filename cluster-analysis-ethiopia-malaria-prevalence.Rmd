---
title: "Global and local incidences of clustering in Ethiopia malaria prevalence data"
author: "Sri Ramesh"
date: "4/4/2020"
output: html_document
---

```{r setup, include=F}
knitr::opts_chunk$set(echo = TRUE)
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

```

In the following, I will discern whether Professor Sturrock's (UCSF) malaria prevalence data out of Ethiopia (*mal_data_eth_2009_no_dups.csv*) contains instances of global and local clustering. The data can be found here: *https://github.com/HughSt/HughSt.github.io/tree/master/course_materials/week1/Lab_files/Data*.

I will focus on the Moran's I test for global spatial autocorrelation (primarily used on point-referenced data like this) using four different methods:

1. Creating an inverse-distance weight matrix to calculate Moran's I
2. Creating a correlogram to plot Moran's I statistics at different distance classes
3. Using k-nearest neighbors and a distance-based approach in tandem to calculate Moran's I
4. Using a Monte Carlo simulation to examine Moran's I

I will then turn to tests for local spatial autocorrelation using LISAs.

## Components of the Moran's I test

The Moran's I test for global spatial autocorrelation is a key tool for global cluster identiifcation in point referenced or continous data. In this case, "clustering" refers not to the distribution of points, which in and of themselves only refer to survey locations. **Clustering here refers to the values of malaria prevalence taken at those points.**

The Moran's I test requires:

1. A set of continuous data. Let's say n = 203 where each point has a continuous value.
2. A defined community of neighbors for each of the 203 points. This requires us to define what's considered a "neighbor," and each methodology has its own definition.
3. A spatial weight for each neighbor linkage.

## Map of Ethiopia Prevalence Data

Here is a map of what the point-process data on malaria prevalence in Ethiopia currently looks like.

```{r}
# load data
url_ethmalaria <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week1/Lab_files/Data/mal_data_eth_2009_no_dups.csv"

ETH_malaria_data <- read.csv(url_ethmalaria, header=T)

# layer map
pal = colorNumeric("Oranges", ETH_malaria_data$pf_pr)

ETH_malaria_data %>%
  leaflet() %>% addProviderTiles("CartoDB.DarkMatter") %>%
  addCircleMarkers(~longitude, ~latitude, fillOpacity=1,
                   fillColor= ~pal(pf_pr),
                   radius=3,
                   weight=0.1,
                   stroke=TRUE) %>%
  addLegend(pal = pal, values = ~pf_pr, title = "Malaria Prevalence")
```

Because we care about clustering in our values, what we are looking for on this map is if there are clusters of *colors* (that is to say, we don't care if we see clusters of points). There do seem to be clusters of lower prevalence values here and there, and one or two clusters of higher prevalence values. We can now confirm this observation statistically using the Global Moran's I test.

Before we jump into the Moran's I test, we first examine our data to see if it is normally distributed.

```{r, warning=F, message=F}

# see if malaria prevalence variable is normally distributed
ETH_malaria_data %>% ggplot(aes(x=pf_pr)) + geom_histogram() +
  xlab("Ethiopia malaria prevalence") +
  ylab("Frequency") +
  ggtitle("Histogram of Ethiopia malaria prevalence values")
```

I saw that the malaria prevalence variable was not normally distributed, so I used a logit transformation and added a minuscule value to it to see if that would correct it. 

```{r, warning=F, message=F}
# see if malaria prevalence variable is normally distributed
ETH_malaria_data <- ETH_malaria_data %>% mutate(log_odds = logit(pf_pr) + 0.0001)

ETH_malaria_data %>% ggplot(aes(x=log_odds)) + geom_histogram() +
  xlab("Log odds of Ethiopia malaria prevalence") +
  ylab("Frequency") +
  ggtitle("Histogram of Ethiopia malaria prevalence values")
```

This did not happen even with a logit transformation. This is not ideal because the null hypothesis of the Moran’s I test is a normal distribution. The lack of normally distributed values will limit the robustness of the Moran's I test, as noted in the conclusion.

## Q1 Create inverse distance matrix

An inverse distance approach defines "neighbors" as essentially every point being a neighbor to every other point in the dataset. This is why we are left with a matrix of n by n values.

In this case, I first take every point of malaria prevalence (n = 203) and calculated the distance to every other point of malaria prevalence. We are left with a matrix of 203 x 203 values where the values represent the distance between each pair of points. These pairs of points are our neighbor linkages.

I then take the inverse of every distance in the matrix. **This essentially weights shorter neighbor linkages higher** as 1 divided by a small number results in a large number, and vice-versa.

Invariably, the original 203x203 distance matrix includes 203 neighbor linkages of value 0, representing the 203 points' distances to themselves. When take the inverse of 0, we are left with Infinity, thus leaving with 203 infinity values on our inverse distance matrix. They appear as a diagonal on our inverse distance matrix. We replace these 203 infinity values with 0 to finish our inverse distance-based weighting process, simply as a data cleaning step.

```{r}
# calculate inverse distance matrix
ETH.dists <- ETH_malaria_data %>% select(latitude, longitude) %>% bind_cols() %>% dist() %>% as.matrix()

# Take the inverse of the matrix values so that closer values have a larger weight and vice versa
ETH.dists.inv <- 1/ETH.dists

# replace the diagonal values with zero
diag(ETH.dists.inv) <- 0
```


```{r}
# From the ape package, computes Moran's I where weights = inverse distance matrix
ETH_malaria_data$log_odds %>% Moran.I(weight=ETH.dists.inv)

```

The key component of the Moran's I is defining the spatial weights. Using an inverse distance based approach to define the spatial weights, we derive a Moran's I statistic of 0.1847263. The associated p-value is lower than 0.01, indicating statistically significant global spatial autocorrelation overall.

## Q1 Calculate the maximum distance between points

We can take a quick look at the maximum distance between any two points on our matrix.

```{r, warning=F, message=F}
maxDist <- ETH_malaria_data %>% select(latitude, longitude) %>% bind_cols() %>% dist() %>% max()
maxDist

```

The max distance is 7.957449 decimal degrees.

## Q1 Plot data in correlogram

Now, we can turn to correlograms which illustrate correlation between distance classes and global spatial autocorrelation. The correlograms essentially tell us what is the maximum distance at which we can compute statistically significant Moran's I values.

```{r, warning=F, message=F}
xy <- ETH_malaria_data %>% select(longitude, latitude) %>% bind_cols() %>% as.matrix()

pgi.cor <- correlog(coords=xy, z=ETH_malaria_data$log_odds, method="Moran", nbclass=10)

pgi.cor %>% plot(xlab = "Distance Classes", ylab = "Moran's I Statistic",
                 main = "Correlogram: 10 Distance Classes v. Moran's I")

pgi.cor

pgi.cor <- correlog(coords=xy, z=ETH_malaria_data$log_odds, method="Moran", nbclass=15)

pgi.cor %>% plot(xlab = "Distance Classes", ylab = "Moran's I Statistic",
                 main = "Correlogram: 15 Distance Classes v. Moran's I")

pgi.cor

pgi.cor <- correlog(coords=xy, z=ETH_malaria_data$log_odds, method="Moran", nbclass=20)

pgi.cor %>% plot(xlab = "Distance Classes", ylab = "Moran's I Statistic",
                 main = "Correlogram: 20 Distance Classes v. Moran's I")

pgi.cor
```

## Q1 What the correlogram reveals

The correlograms reveal global spatial autocorrelation and instances of positive and negative clustering at different distance classes (where all there is a red dot on the correlograms). I compared correlation when 10 vs. 15 vs. 20 categories were used. With 10 categories, we have statistically significant positive clustering at a maximum distance of 0.3979675 decimal degrees. With 15, we have statistically significant positive clustering at a max distance of 0.7958354 decimal degrees, and with 20, we have it at a max distance of 0.9947692 decimal degrees. This is all at a statistical significance level of 0.05. At further distances, we seem to have primarily negative spatial autocorrelation but none of it is statistically significant.

## Q2 Create different neighbor structures and assign weights to them

An inverse-distance based approach to calculating Moran's I involves defining "neighbors" as every point to every other point, no matter how far apart they are. This means each of our 203 points of malaria prevalence has 202 neighbors, and thus 202 neighbor linkages. These linkages are then each assigned a weight using an inverse distance approach, which greatly decreases the importance of many of the weights. However, it remains to be the case that the inverse distance based approach still uses the greatest number of linkages that are possible.

We can set a limit to this, where not every point is a neighbor to every other point, by defining a threshold distance. This is the crux of the concept of *k-nearest neighbors* and *distance-based neighbors.* k-nearest neighbors has to do with limiting the **number of neighbor linkages.** When k = 1, that means only immediately adjacent points are called "neighbors." We will use k=1 or immediate neighbors to set neighbor limits:

```{r, warning=F, message=F}
# define k-nearest neighbors
coords <- xy %>% coordinates()
IDs <- coords %>% as.data.frame() %>% row.names()
Neigh_nb <- coords %>% knearneigh(k=1, longlat = TRUE) %>% knn2nb(row.names=IDs)

dsts <- Neigh_nb %>% nbdists(coords) %>% unlist()
dsts %>% summary()
```

The five-number summary above illustrates the distribution of distances used to provide *at least one neighbor (k = 1)* to each malaria prevalence point. Now is where the **distance-based approach** comes in. We take the max distance used to provide at least one neighbor and use that as an input to a distance-based neighbor approach.

```{r, warning=F, message=F}
# find maximum distance used to provide at least one neighbor to each point
max_1nn <- dsts %>% max()
max_1nn
```

The max distance to provide at least one neighbor to each point is 0.6581074 units. We can input this into our distance-based approach for defining neighbors to see how many neighbor linkages it gives us. For comparison sake, we can also increase the limit of what counts as a "neighbor" by computing twice the max distance. We can see how many more neighbor linkages twice the max distance gives us too. 

```{r, warning=F, message=F}
# define distance-based neighbors within maximum distance:
Neigh_kd1 <- coords %>% dnearneigh(d1=0, d2=max_1nn, row.names=IDs) # when k=1 (immediate neighbors)
# define distance-based neighbors within 2 times the maximum distance:
Neigh_kd2 <- coords %>% dnearneigh(d1=0, d2=2*max_1nn, row.names=IDs) # when k=1

# list of neighbor structures
nb_1 <- list(d1=Neigh_kd1, d2=Neigh_kd2) #when k=1

# Check for symmetry
nb_1 %>% sapply(function(x) is.symmetric.nb(x, verbose=F, force=T))
nb_1 %>% sapply(function(x) n.comp.nb(x)$nc)
```

As an aside here, the outputs of "5" and "1" indicate that there does not appear to be symmetry between neighbors (i.e. if *i* is a neighbor of *j*, then *j* is a neighbor of *i*). This should be investigated.

Now to see how many neighbor linkages each distance parameter leaves us with:

```{r}
# plot neighbor links with different definitions of "neighbors" (max distance vs. 2 times max distance)
par(mfrow=c(2,1), mar= c(1, 0, 1, 0))
xy %>% plot(pch=16)
Neigh_kd1 %>% plot(coords, col="green",add=T, main = "Neighbor Definition: Maximum Distance")
Neigh_kd1

xy %>% plot(pch=16)
Neigh_kd2 %>% plot(coords,col="green", add=T, main = "Neighbor Definition: 2 times Maximum Distance")
Neigh_kd2
```

**Clearly, looking at the "number of nonzero links" we see that twice the max distance gives us more than double the number of neighbor linkages than the max distance.** The number of neighbor linkages would only increase as we increase the distance parameter used. This might tap out at n-squared (in our case 203 * 203 = 41,209) in which case every point is seen as a  "neighbor" to every other point.

Moving forward, let's use the number of neighbor linkages as set by the **maximum distance** parameter, or 2956 neighbor linkages. Now that we have our number of neighbor linkages, we need to assign weights to these linkages.

```{r}
weights <- Neigh_kd1 %>% nb2listw(style="W")
weights
```

Here, n = 203 indicates that we have 203 malaria prevalence points. By definition per the k-nearest neighbor method, every point has at least one neighbor. Now the number of neighbors that each point has can range quite widely. Here is the total number of neighbors, and corresponding spatial weights, for our 1st malaria prevalence point:

```{r}
# list of neighbors (each integer represents a unique neighbor/neighbor linkage to point 1) 
weights[["neighbours"]][[1]]

# list of weights (one weight per neighbor linkage)
weights[["weights"]][[1]]
```

The above output shows that point 1 has 8 neighbors (ie. there are 8 integers in that vector), and each of those 8 linkages were assigned a weight of 0.125. This is what's called **row standardized weights;** the value of each weight is the output of simply dividing 1 over the number of neighbors (1/8 = 0.125). Row standardized weights are the default option in R when assigning spatial weights.

## Q2 Run the Moran's I test

Now that we have our defined list of neighbor linkages and corresponding row-standardized weights matrix, we can calculate Moran's I.

```{r, warning=F, message=F}
# spdep package, Moran's I using k=1 and row-standardized weights
ETH_malaria_data$log_odds %>% moran.test(listw=weights)
```

Given the parameters of k=1 and row-standardized weights, we find a statistically significant Moran's I statistic of 0.2709168871, indicating that our data has positive global clustering.

## Q3 Compare results from inverse distance matrix and distance-based neighbor approaches, and draw conclusions

Both the inverse distance based and k-nearest/distance-based neighbor approaches show evidence of statistically significant global spatial autocorrelation (ie positive global clustering). The p-value associated with the Moran's I statistic using the k-nearest/distance based approach (where k = 1, and row-standardized weights are used) is less than 2.2e-16. This is much less than even the 0.01 threshold for p-values. The Moran's I statistic with this approach takes on a value of 0.2709168871, which means positive clustering.

The inverse distance matrix and correlograms tell a slightly more granular story of global spatial autocorrelation than the distance-based approach. The correlograms that I generated for 10, 15, and 20 bins indicate that we have statistically significant positive clustering when the max distance is 0.3979675 to 0.9947692 decimal degrees.

Overall, the distance-based neighbor approach shows higher spatial autocorrelation (0.2709168871) than the inverse-distance based approach (0.1847263). This difference makes sense because we approached the neighbor definition and spatial weights questions in different ways with each approach. In the former, we used all neighbor linkages possible and simply took the inverse of the distance as weights. In the latter, we restricted the definition of what counts as a "neighbor" to the max distance needed to give every point at least one neighbor (k=1), and also, we used row-standardized weights instead of weights based on inverse distances.

Perhaps because the k-nearest neighbor method is more popular/common, we can stick with that method and use its parameters for running a Monte Carlo simulation. This might confirm what the value of the Moran's I statistic is (again, using row-standardized weights and k=1):

```{r}
set.seed(1234)
bperm <- ETH_malaria_data$log_odds %>% moran.mc(listw=weights,nsim=999)
bperm

# Plot simulated test statistics
par(mfrow=c(1,1), mar= c(5, 4, 4, 2))
hist(bperm$res, freq=T, breaks=20, xlab="Moran's I statistics",
     main = "Distribution of Moran's I statistics over 999 simulations")
abline(v=0.27092, col="red")

```

The Monte Carlo simulation results indicate that with a high degree of statistical significance, the Moran's I statistic is 0.27092.

In conclusion, we could tentatively say there is global spatial autocorrelation - specifically positive clustering, as the Moran's I statistic is positive. This is a tentative conclusion only because the data was not normally distributed to begin with. This reduces the robustness of the Moran's I test because the null hypothesis of the Moran's I test is that the data is normally distributed.

## Local Spatial Autocorrelation using LISAs

Clustering in the underlying spatial process (in this case being *values taken at each point*) can vary from location to location. We can decompose the Global Moran's I into its localized components at each point because the Moran's I is ultimately a sum of its local parts.

With LISA, our "localized areas" are simply 203 windows, one window drawn around each prevalence point in our dataset.

```{r lisa, warning=F, echo=T, error=F}
# First calculate the local Moran's I around each point based on row standardized weights where k = 1
I <- ETH_malaria_data$log_odds %>% localmoran(weights) # "spdep" package

# Print 'LISA' for each point
Coef <- I %>% as.data.frame(row.names=row.names(coords), check.names=FALSE) %>% printCoefmat()

```

The LISA output provided here includes:

* Observed local Moran's I *(Ii)*
* Expected local Moran's I *(E.Ii)* (!! based on the Global Moran's I value? or assumption of stationarity?)
* Variance of the local Moran's I *(V.Ii)*
* Z-value which can be used to determine statistical significance of the clustering *(Z.Ii)*

We can map the observed vs. expected local moran's I values in a 4-quadrant grid to see where deviations are the greatest in terms of the observed vs. expected Moran's I. Areas of higher deviation would indicate that the underlying spatial process is *not* stationary, as it deviates at those prevalence points.

```{r maplisas, echo=T, warning=F, error=F}
# Plot the spatial data against its spatially lagged values (the weighted mean of its neighbors)                         
nci <- ETH_malaria_data$log_odds %>% moran.plot(listw=weights,
                                                main = "Observed vs. Expected Log prevalence",
                                               xlab="Log prevalence (log_odds variable)",
                                               ylab="Spatially lagged log prev",
                                               labels=T, pch=16, col="grey")
```

This plot illustrates all 203 of our malaria point prevalence points in terms of observed vs. expected log prevalence. "Expected" prevalence refers to **the weighted mean of each point's neighbors.** The gray spots are used to represent those expected values.

Reading the quadrants can be done as follows:

* Bottom right: High-Low (HL) meaning **abnormally high: high observed value of log prevalence, but low expected value, where expectation is based on the mean in the local region**. The area around the prevalence point has low values, so we expect the observed prevalence point we have to be low too. However, the observed prevalence point we have is actually a high value.
* Top left: Low-High (LH) meaning **abnormally low: high expected value, but low observed value**. The area *around the point* have high values, but the observed point itself has a low value.
* Bottom left: Low-Low (LL): observed value matches the mean in the localized region, which is low
* Top right: High-High (HH): observed value matches the mean in the localized region, which is high

In the case of this malaria point prevalence data, most of our data sits in the LH and HH quadrants. This means that most of our deviances from the stationarity assumption are in the LH region: we have abnormally low observed prevalence in regions where mean prevalence is higher.

Now we can identify the number of outlier points in our data:

```{r outliers}
# find which points are statistically significant outliers
infl <- apply(nci$is.inf,1,any) 
sum(infl %in% T)

```

It seems that 22 of our 203 points are outliers, which is more than one might expect is due to chance.

```{r lisamapping, echo=F, include=F}
# x <- ETH_malaria_data$log_odds
# lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=T)
# wx <- lag(weights, ETH_malaria_data$log_odds)
# 
# lhwx <- wx %>% cut(breaks=c(min(wx), mean(wx), max(wx)),
#                    labels=c("L", "H"), include.lowest=T)
# 
# lhlh <- interaction(lhx, lhwx, infl, drop=T)
# 
# names <- rep("none", length(lhlh))
# 
# names[lhlh %in% "L.L.TRUE"] <- "LL"
# names[lhlh %in% "H.L.TRUE"] <- "HL"
# names[lhlh %in% "L.H.TRUE"] <- "LH"
# names[lhlh %in% "H.H.TRUE"] <- "HH"
# 
# ETH_malaria_localM <- cbind(xy, names) %>% as.data.frame()
# 
# ETH_malaria_localM <- ETH_malaria_localM %>% 
#   mutate(longitude = as.numeric(as.character(longitude)),
#          latitude = as.numeric(as.character(longitude)))
# 
# factpal <- c( "cyan4","coral4","coral","cyan","lightgrey") %>% colorFactor(names)
#  
# ETH_malaria_localM %>% leaflet() %>% addTiles() %>%
#   addCircleMarkers(~longitude, ~latitude, fillOpacity=1,
#                    color= ~factpal(names), radius=4, stroke=TRUE, weight=1) %>%
#   addLegend(pal = factpal, values = ~names, title="Class")

```
