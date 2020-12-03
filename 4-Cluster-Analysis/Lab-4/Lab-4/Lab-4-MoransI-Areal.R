# load NY data
nydata_url <- "https://raw.githubusercontent.com/HughSt/HughSt.github.io/master/course_materials/week7/Lab_files/nydata.geojson"
nydata <- rgdal::readOGR(nydata_url)

# lets take a look at the data
head(nydata@data)

# append column corresponding to incidence per capita
nydata$inc_per_1000 <- (nydata$Cases / nydata$POP8) * 1000 #can't use dplyr since its an object of class spdf, not a df

# neighbors are defined as all polygons that share a BOUNDARY
nydata_queen <- nydata %>% poly2nb() #queen contiguity, from spdep package
nydata_rook <- nydata %>% poly2nb(queen = F) #queen contiguity, only slightly diff from queen contiguity

# store coordinates
coords_ny <- nydata %>% coordinates() #this is the lat and long column taken together

# plot and compare neighbors
par(mfrow=c(1,2))
nydata %>% plot(main = "Queen Contiguity")
nydata_queen %>% plot(coords_ny,col="blue",add=T)
nydata %>% plot(main = "Rook Contiguity")
nydata_rook %>% plot(coords_ny,col="green",add=T)

# 3 assign weights ----
nydata_queen_w <- nydata_queen %>% nb2listw() # row-standardized weights, where EACH ROW SUMS TO ONE
nydata_queen_wB <- nydata_queen %>% nb2listw(style="B") # binary style weights

nydata_rook_w <- nydata_rook %>% nb2listw()
nydata_rook_wB <- nydata_rook %>% nb2listw(style="B")

# moran's I test for global spatial autocorrelation in areal data ----
nydata$inc_per_1000 %>% moran.test(listw=nydata_queen_w)  #using row standardized


