par(mar=c(1,1,1,1))

library(sp)
library(ggplot2)
library(rgdal)
library(spdep)
library(leaflet)
library(spaMM)
library(viridis)

# Let's look at the data
nydata <- rgdal::readOGR("https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week7/Lab_files/nydata.geojson")

head(nydata@data)

# Let's create an incidence column (count per 1000 people)
nydata$inc_per_1000 <- (nydata$Cases / nydata$POP8) * 1000

# plot incidence data
cases_pal <- colorBin(viridis(4), bins = c(0, 0.1, 0.5, 1, 8), nydata$inc_per_1000)

plot(nydata, col = cases_pal(nydata$inc_per_1000), asp = 1)
legend('bottomright', legend = c('0 - 0.1', '0.1 - 0.5', '0.5 - 1', '1 - 8'),
       fill = cases_pal(c(0, 0.1, 0.5, 1, 8)),
       title = 'Cases / 1000')

# For more info on the dataset type ?spData::nydata
?spData::nydata
data <- spData::nydata

# use poisson regression which is suitable for count outcomes

# First round the case numbers
nydata$CASES <- round(nydata$TRACTCAS) #TRACTCAS is the count data, not a typo

# identify covariates that could also predict count (exposure to toxic gas sites, etc)
# remember these covariates are at the AREAL level (eg census tracts),
# so we can use things like median family income per census tract or whatever as covariates

nyc_glm_mod <- glm(CASES ~ PEXPOSURE + PCTOWNHOME + PCTAGE65P, offset = log(POP8), 
                   data = nydata, family = 'poisson')
summary(nyc_glm_mod)
# PEXPOSURE is statistically significant, so is percent age over 65


# plot fitted vs. observed values of the model
# Scatter plot
ggplot() + geom_point(aes(nyc_glm_mod$fitted.values, nydata$CASES))

# Create maps
nydata$fitted <- nyc_glm_mod$fitted.values

col_pal <- colorNumeric(topo.colors(64), c(0,9))
par(mfrow=c(2,1), mar=c(rep(0.8,4)))
plot(nydata, col = col_pal(nydata$fitted), asp=1, main = 'Fitted')
plot(nydata, col = col_pal(nydata$CASES), asp=1, main = 'Observed')
legend("bottomright", inset = 0.2,
       legend=0:9, fill = col_pal(0:9),
       title = 'Counts')

# Contiguity neighbors - all that share a boundary point
nydata_nb <- poly2nb(nydata)  #queen contiguity
nydata_nb

# coordinates
coords<-coordinates(nydata)

#view the neighbors
plot(nydata, asp = 1)
plot(nydata_nb,coords,col="blue",add=T)

# convert neighborhood list to a neighborhood MATRIX (aka adjacency matrix)
adj_matrix <- nb2mat(nydata_nb, style="B")
# remember to include style=B. for CAR, we have to use binary weights

# Annoyingly we have to remove the row.names, otherwise R complains later
row.names(adj_matrix) <- NULL

# fit the CAR model
nyc_car_mod <- fitme(CASES ~ PEXPOSURE + PCTOWNHOME + PCTAGE65P + adjacency(1|AREAKEY), 
                     offset = log(nydata$POP8),
                     adjMatrix = adj_matrix,
                     data = nydata@data, family = 'poisson')

summary(nyc_car_mod)

# so now, we do no longer have statistically significant covariates
# the inclusion of the spatial component removed those associations

# remember to generate 95% CIs
terms <- c('PEXPOSURE', 'PCTOWNHOME', 'PCTAGE65P')
coefs <- summary(nyc_car_mod)$beta_table %>% as.data.frame()

row <- row.names(coefs) %in% terms
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, 'Cond. SE']
data.frame(terms = terms,
           IRR = exp(coefs[row,'Estimate']),
           lower = exp(lower),
           upper = exp(upper))

# Scatter plot
ggplot() + geom_point(aes(fitted(nyc_car_mod), nydata$CASES))

# Create maps
nydata$fitted_car <- fitted(nyc_car_mod)

col_pal <- colorNumeric(topo.colors(64), c(0,9))
par(mfrow=c(2,1), mar=rep(2,4))
plot(nydata, col = col_pal(nydata$fitted_car), asp=1, main = 'Fitted - CAR')
legend("bottomright", inset = 0.1, cex = 0.8,
       legend=0:9, fill = col_pal(0:9),
       title = 'Counts')
plot(nydata, col = col_pal(nydata$CASES), asp=1, main = 'Observed')



