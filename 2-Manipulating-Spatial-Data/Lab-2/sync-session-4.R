# Semivariogram v Correlogram ----
## Semivariogram: how semivariance CHANGES as distance CHANGES
### it gives you of what the INVERSE of the correlogram is
### semivariogram = for understanding the shape of the spatial process
## how variance between points changes over a distance between points

## as you move away from points (increase distance), how does the variance change? = semivariogram
## inverse of a correlogram

## when you are very close, you should have LOW variance
## that shape = spherical, elliptical, etc.
## that shape = modeling that spatial process

## once we accurately estimate that shape (ie elliptial model, spherical model, etc)
## that model is what allows us to then predict points where we don't actually have them
## putting potential values where we don't actually have data

# Correlogram ----
### semivariogram maps over DISTANCE, how variance changes)
## increased distance means HIGH variance

## the opposite with correlogram, which has to do with CORRELOGRAM
## increased distance means low correlation (no spatial autocorrelation), and 
## reduced distance means high correlation (more spatial autocorrelation)

## semivariogram shows the inverse of the correlogram

## clustering over different distances = what correlogram
## measure of global clustering = Moran's I

## in Lab 4, the first 3 distance classes = statistically significant Moran's

## units of correlogram: the units are whatever you have defined in your CRS
## in this case, use decimal degrees!
## decimal degrees are one reference system
## things are closer together at the poles than they are at the equator
## at the equator, one decimal degree = 110 or 111 kilometers

## how many distance classes to set? there's a rule of thumb...
# depends on how finely you want to look at the clustering


## semivariogram = just trying to give you an estimate of the SHAPE of that decay (based on distance)
## semivariogram should slope up because variance (ie semivariance) increases as distance increases

## calculating moran's I using a binary distance matrix ----

# this is different from k-nearest neighbors

# checking for symmetry: should always have symmetry basically (NOT for k-nearest neighbors though)
# symmetry = i.e. if i is a neighbor of j, then j is a neighbor of i
# d1 and d2 are TRUE, and d1 = 1 and d2 = 1
# Anything other than 1 means you don't have symmetry going on...
# THIS MEANS if i is a neighbor of j (d1?), j is NOT a neighbor of i (d2?)
# that can't be true obviously, so something's wrong.

# we need to assign weights to the neighbor list
# We will use function where neighbors are defined using maximum distance


# increasing distance within which one is considered a neighbor = BIG increase in neighbor linkages

# k-nearest neighbors: number of neighbors
# other ways to define neighbors: by distance


# if you use more distance classes, you will have smaller classes
# use the distance instead (not the bins)

# correlogram: if we increase the window over which we calc this statistic, different stat. signif.
# over what SPATIAL SCALE that's occuring
# that is basically the distance classes, what distances (over 1 decimal degree)

# distance based: global spatial autocorrelation

# qualifier: these data were highly skewed so that influences the robustness of the Moran's I.





