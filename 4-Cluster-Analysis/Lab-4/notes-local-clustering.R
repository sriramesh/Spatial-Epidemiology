# case control data
casecontrol_url <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week3/Lab_files/CaseControl.csv"
casecontrol <- read.csv(casecontrol_url)

# convert to an spdf and subset cases and controls
casecontrol_spdf <- SpatialPointsDataFrame(coords = casecontrol %>% select(long, lat),
                                           data = casecontrol %>% select(household_id, case))
cases <- casecontrol_spdf[casecontrol$case %in% 1, ]
controls <- casecontrol_spdf[casecontrol$case %in% 0, ]

# namibia adm 1 boundaries
NAM_Adm0 <- raster::getData('GADM',country='NAM',level=0)

# Let's plot and see what we have
pal <- c("blue", "red")
case_color_scheme <- colorNumeric(pal, casecontrol_spdf$case)

leaflet() %>% addTiles() %>%
  addCircleMarkers(data=casecontrol_spdf,
                   color = case_color_scheme(casecontrol_spdf$case),
                   radius=3)

# In the previous lecture, you already generated first order KDEs
# and calculated the ratio of the density estimate of cases:controls
# Now you will look at 2nd order functions, summarizing the spatial dependence between events

# convert object class of cases subset into ppp
cases_ppp <- cases %>% as("ppp")

# apply ripley's k
# summarize spatial dependence between events at a wide range of spatial scales

ripleys_k <- cases_ppp %>% Kest(correction=c("isotropic", "Ripley")) #uses the "spatstat" package

# Plot the estimate of K(r); note different border-corrected estimates ('iso', 'border' and 'trans')
par(mfrow=c(1,1))

# Red dashed line is expected K value computed for a CRS process
ripleys_k %>% graphics::plot(xlab="d (dd)", ylab="K(dd)", main="Ripley's K")

# Plot confidence envelope using monte carlo simulation
E <- cases_ppp %>% envelope(Kest, nsim=999)
E %>% graphics::plot()


# The K-function computed for cases assumes that H0 is complete spatial randomness.
# What are the limitations of this assumption?

# * global measures assume the underlying spatial process is STATIONARY (consistent?) across the entire spatial frame
# * if this false, the global spatial autocorrelation tests can leave out localized differences in the underlying spatial process
# * this is particularly true if we had a large dataset (lets say, 1 million malaria prevalence points instead of 203)
# 
# * lisas scan the whole dataset, but only MEASURE clustering at local boxes/windows?
#   
#   * what is "spatial dependence"? is that the same as "clustering" or "spatial autocorrelation"?
#   
#   * LISA = Local Indicators of Spatial Association/Autocorrelation
# * Local Moran's I = decomposed version of global Moran's I
# * looks for deviations in localized points

# Most point processes follow a Poisson distribution (i.e. we can model the distribution of current and future cases using a Poisson random variable). Complete spatial randomness (CSR) is equivalent to a *homogenous Poisson process.* CSR means that events are distributed independently over the entire region. To detect clustering, we are looking for a departure from CSR. In other words, we are testing whether the underlying point process is best modeled using a homogenous poisson distribution, or something else. The former would indicate no clustering, and the latter would indicate clustering.

