# number of ways rasters can be resampled (changed for their extent and resolution) and it depends on
# if your data is continuous
# if your data is categorical
# etc...

library(raster)
library(tidyverse)

bk_url <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week2/Lab_files/BF_land_use.tif"
pop_url <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week2/Assignment/BF_pop.tif"

# load burkina faso data and boundaries
BF_elev <- raster::getData("alt", country="BF")
BF_land_use <- raster(bk_url)
pop <- raster(pop_url)

crs(BF_land_use)
crs(BF_elev)

# change crs ----
# reproject - can time a little time
BF_land_use <- BF_land_use %>% 
  projectRaster(crs=crs(BF_elev), method="ngb") 

# crop rasters to the same bounding box ----
# For demo purposes, we are going to crop the rasters across an arbitray bounding boox
BF_elev_crop <- BF_elev %>% crop(extent(BF_elev,499,507,301,309)) 
pop_crop <- pop %>% crop(extent(BF_elev_crop)) #crop population layer to the extent of BF_elev layer
BF_land_use_crop <- BF_land_use %>% crop(extent(BF_elev_crop)) #crop land use layer to extent of BF_elev layer

# make ref raster ----
# create a reference raster that is low-res: we will resample using this reference raster
new_raster <- BF_elev_crop %>% aggregate(fact=3) # aggregate BF elev layer by a factor of 3, new extent
new_raster <- new_raster %>% raster::shift(dx=0.004, dy=0.002) #new extent
new_raster_outline <- new_raster %>% rasterToPolygons(dissolve=TRUE) #new extent

# plot the cropped elevation raster and overlay the grid of the new raster we are resampling to
BF_elev_crop %>% plot()

# plot cell outlines
new_raster_outline %>% lines()

# resample: change res and extent to new raster ----
# bilinear interpolation since it is continuous data
# old raster: BF_elev_crop
# new raster: BF_elev_crop aggregated by a function of 3 AND shifted to dx=0.004, dy=0.002
BF_elev_crop_resampled_bilin <- BF_elev_crop %>% resample(new_raster, method="bilinear")
BF_elev_crop_resampled_bilin %>% plot()
new_raster_outline %>% lines()

# repeat process with landuse raster ----
BF_land_use_crop %>% plot() #plot original raster
table(BF_land_use_crop[]) # look at land use classifications -- it is categorical data

# aggregate land use raster by factor of 9 and plot and resample that aggregated raster (not the ref raster)
# we aggregate to get it to match the reference raster?
BF_land_use_crop_aggregated <- BF_land_use_crop %>% aggregate(fun='modal', fact = 9)
BF_land_use_crop_aggregated %>% plot()

new_raster_outline %>% lines()

# resample aggregated 9 by 9 raster to reference raster from old problem
# when you resample, make sure you use the new RASTER itself, not the new raster's outline
# ngb = categorical data (land use classifications)
# binomial = continuous data (elevation data)
BF_land_use_crop_aggregated_resamp <- BF_land_use_crop %>% resample(new_raster, method = "ngb")
BF_land_use_crop_aggregated_resamp %>% plot()
new_raster_outline %>% lines()

# repeat resampling process on population ----
# The raw data in the pop_crop shows the number of individuals per cell
# If we want to resample to a coarser resolution,
# we are probably most interested in maintaining these as counts and
# therefore we want to SUM the numbers as opposed to use interpolation or nearest neighbour

# aggregation prior to resampling ----
pop_crop[]
pop_crop %>% plot(col=topo.colors(64))
pop_crop %>% cellStats(sum) # we are summing per cell
pop_crop_aggregated <- pop_crop %>% aggregate(fact = 3.012048, fun = sum) #use sum because its population per capita
pop_crop_aggregated %>% cellStats(sum) # Check the total population is the same as the original raster

# outline the new extent ----
pop_crop_aggregated %>% plot(col=topo.colors(64))
new_raster_outline %>% lines()

pop_crop_aggregated_resamp <- pop_crop_aggregated %>% resample(new_raster)
pop_crop_aggregated_resamp %>% plot(col=topo.colors(64))
new_raster_outline %>% lines()

# The total should more or less be equal to the original
# its not exact as we have to use interpolation to predict populations in the new cell.
pop_crop_aggregated_resamp %>% cellStats(sum) # Not bad..
