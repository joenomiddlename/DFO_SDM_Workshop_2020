# Exploratory Analysis Script
library(rgeos)
library(rgdal)
library(sp)
library(maptools)
library(spatstat)
library(spdep)
library(INLA)
library(inlabru)
library(readxl)
library(lubridate)
library(ggmap)
library(raster)

# load in data
list2env(readRDS('./Data/Compiled_Data.rds'),globalenv())

# Always check the coordinate reference system
Sightings_DRWW_sp@proj4string # latlon
Sightings_survey@proj4string # latlon
Effort_survey@proj4string # latlon

# Try plotting the data. The gg and gm functions works great for this
# If the data are in lat/lon format then gmap allows google maps to be downloaded.
source('utility_functions.R')
# The bespoke gg.spatiallines_mod function plots transect lines too!

gmap(Sightings_survey) +
  gg(Domain_latlon) +
  gg.spatiallines_mod(Effort_survey) +
  gg(Sightings_survey, colour='blue') +
  gg(Sightings_DRWW_sp, colour='purple') +
  gg(WW_ports, colour='red')

# Use statistics canada Lambert projection CRS, but change to km
#Can_proj <- CRS('+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +units=km +no_defs')
# Can_proj <- CRS(SRS_string = 'PROJCS["NAD27 / UTM zone 20N",
#     GEOGCS["NAD27",
#         DATUM["North_American_Datum_1927",
#             SPHEROID["Clarke 1866",6378206.4,294.9786982139006,
#                 AUTHORITY["EPSG","7008"]],
#             AUTHORITY["EPSG","6267"]],
#         PRIMEM["Greenwich",0,
#             AUTHORITY["EPSG","8901"]],
#         UNIT["degree",0.0174532925199433,
#             AUTHORITY["EPSG","9122"]],
#         AUTHORITY["EPSG","4267"]],
#     PROJECTION["Transverse_Mercator"],
#     PARAMETER["latitude_of_origin",0],
#     PARAMETER["central_meridian",-63],
#     PARAMETER["scale_factor",0.9996],
#     PARAMETER["false_easting",500000],
#     PARAMETER["false_northing",0],
#     UNIT["metre",1,
#         AUTHORITY["EPSG","9001"]],
#     AXIS["Easting",EAST],
#     AXIS["Northing",NORTH],
#     AUTHORITY["EPSG","26720"]]')
 Can_proj <- CRS(SRS_string='PROJCS["NAD83(CSRS) / UTM zone 20N",
     GEOGCS["NAD83(CSRS)",
         DATUM["NAD83_Canadian_Spatial_Reference_System",
             SPHEROID["GRS 1980",6378137,298.257222101,
                 AUTHORITY["EPSG","7019"]],
             AUTHORITY["EPSG","6140"]],
         PRIMEM["Greenwich",0,
             AUTHORITY["EPSG","8901"]],
         UNIT["degree",0.01745329251994328,
             AUTHORITY["EPSG","9122"]],
         AUTHORITY["EPSG","4617"]],
     UNIT["metre",1,
         AUTHORITY["EPSG","9001"]],
     PROJECTION["Transverse_Mercator"],
     PARAMETER["latitude_of_origin",0],
     PARAMETER["central_meridian",-63],
     PARAMETER["scale_factor",0.9996],
     PARAMETER["false_easting",500000],
     PARAMETER["false_northing",0],
     AUTHORITY["EPSG","2961"],
     AXIS["Easting",EAST],
     AXIS["Northing",NORTH]]')
Can_proj <- CRS("+init=EPSG:2961")
Can_proj <- fm_crs_set_lengthunit(Can_proj, unit='km')

# Transform all the spatial objects to the new format
Sightings_DRWW_sp <- spTransform(Sightings_DRWW_sp, Can_proj)
Sightings_Opp_sp <- spTransform(Sightings_Opp_sp, Can_proj)
Sightings_survey <- spTransform(Sightings_survey, Can_proj)
Effort_survey <- spTransform(Effort_survey, Can_proj)
WW_ports <- spTransform(WW_ports, Can_proj)
Slope_tmp <- spTransform(Slope, Can_proj)
# Note the warning message. This is due to the curvature of the spatial transform. 
# We can bypass this by using bilinear interpolation in the raster package
rm(Slope_tmp)
Slope <- raster(Slope)
Slope <- projectRaster(Slope, crs=Can_proj)
Slope <- as(Slope, 'SpatialPixelsDataFrame')
# repeat for Bathym
Bathym <- as(projectRaster(raster(Bathym), crs=Can_proj), 'SpatialPixelsDataFrame')
Domain <- spTransform(Domain, Can_proj)
Dist_Brier <- as(projectRaster(raster(Dist_Brier), crs=Can_proj), 'SpatialPixelsDataFrame')
Dist_Quoddy <- as(projectRaster(raster(Dist_Quoddy), crs=Can_proj), 'SpatialPixelsDataFrame')

# Plot to see the new CRS - note the gmap function no longer works
ggplot() +
  gg(Domain) +
  gg.spatiallines_mod(Effort_survey) +
  gg(Sightings_survey, colour='blue') +
  gg(Sightings_DRWW_sp, colour='purple') +
  gg(WW_ports, colour='red')

ggplot() +
  gg(Domain) +
  gg(Bathym)

ggplot() +
  gg(Domain) +
  gg(Dist_Brier)

ggplot() +
  gg(Domain) +
  gg(Dist_Quoddy)

# Is there autocorrelation between sighting locations?
View(Sightings_DRWW_sp@data[2:4,])
# There is evidence of repeated sightings being made within a day of the same individuals
# How far away are these sightings in space (in meters)?
gDistance(Sightings_DRWW_sp[2:4,], byid = T)
# How far away are these sightings in time (in hours)?
Sightings_DRWW_sp[3:4,]$WS_TIME-Sightings_DRWW_sp[2:3,]$WS_TIME
# Could this be the same animal? How fast is this implied movement?
gDistance(Sightings_DRWW_sp[2:4,], byid = T)[cbind(c(1,2),c(2:3))] / 
  (c(18, 9)/60)
# Speeds of ~ 4km and 16 km Could be same animal
Sightings_DRWW_sp$PLATFORM <- as.factor(Sightings_DRWW_sp$PLATFORM)

# There is no overlap between the sightings from the two WW companies
ggplot() + gg(Domain) + gg(Sightings_DRWW_sp, colour=Sightings_DRWW_sp$PLATFORM_CODE)

# We need to subset data to be safe
Sightings_Brier_nodup <- Sightings_DRWW_sp[Sightings_DRWW_sp$PLATFORM=='BRIER ISLAND WHALEWATCH',]
Sightings_Brier_nodup <- Sightings_Brier_nodup[!duplicated(Sightings_Brier_nodup$WS_DATE),]

Sightings_Quoddy_nodup <- Sightings_DRWW_sp[Sightings_DRWW_sp$PLATFORM=='QUODDY LINK',]
Sightings_Quoddy_nodup <- Sightings_Quoddy_nodup[!duplicated(Sightings_Quoddy_nodup$WS_DATE),]

dim(Sightings_DRWW_sp)[1]; dim(Sightings_Brier_nodup)[1]; dim(Sightings_Quoddy_nodup)[1];

# Test the hypothesis of autocorrelation between sightings/encounters locations
# Compute the distance in time and the euclidean distance in space between sightings
# Again
difftime_Brier<- dist(as.numeric(ymd(Sightings_Brier_nodup$WS_DATE)))
diffspace_Brier <- dist(Sightings_Brier_nodup@coords)

difftime_Quoddy<- dist(as.numeric(ymd(Sightings_Quoddy_nodup$WS_DATE)))
diffspace_Quoddy <- dist(Sightings_Quoddy_nodup@coords)

ggplot(data=data.frame(difftime = as.numeric(difftime_Brier)[as.numeric(difftime_Brier)<85],
                       diffspace = as.numeric(diffspace_Brier)[as.numeric(difftime_Brier)<85]),
       aes(x=difftime, y=diffspace)) +
  geom_point() + geom_smooth(method='loess')

ggplot(data=data.frame(difftime = as.numeric(difftime_Quoddy)[as.numeric(difftime_Quoddy)<85],
                       diffspace = as.numeric(diffspace_Quoddy)[as.numeric(difftime_Quoddy)<85]),
       aes(x=difftime, y=diffspace)) +
  geom_point() + geom_smooth(method='loess')
# Caused by 2 discrete clusters/regions visited
ggplot() + gg(Domain) + gg(Sightings_DRWW_sp, colour=Sightings_DRWW_sp$PLATFORM_CODE)

# Focus on autocorrelations within clusters
ggplot(data=data.frame(difftime = as.numeric(difftime_Quoddy)[as.numeric(difftime_Quoddy)<85 & as.numeric(diffspace_Quoddy)<50],
                       diffspace = as.numeric(diffspace_Quoddy)[as.numeric(difftime_Quoddy)<85 & as.numeric(diffspace_Quoddy)<50]),
       aes(x=difftime, y=diffspace)) +
  geom_point() + geom_smooth(method='loess')

# Note that the autocorrelations could be due to a seasonal change in density
# It may not be problematic at all!

# Split the data into training and test
unique(Sightings_survey$YEAR)
unique(Effort_survey$YEAR)

Sightings_Brier_nodup_test <- Sightings_Brier_nodup[Sightings_Brier_nodup$YEAR==2011,]
Sightings_Brier_nodup <- Sightings_Brier_nodup[Sightings_Brier_nodup$YEAR!=2011,]
Sightings_Quoddy_nodup_test <- Sightings_Quoddy_nodup[Sightings_Quoddy_nodup$YEAR==2011,]
Sightings_Quoddy_nodup <- Sightings_Quoddy_nodup[Sightings_Quoddy_nodup$YEAR!=2011,]
Sightings_survey_test <- Sightings_survey[Sightings_survey$YEAR==2011,]
Sightings_survey <- Sightings_survey[Sightings_survey$YEAR!=2011,]
Effort_survey_test <- Effort_survey[Effort_survey$YEAR=='2011',]
Effort_survey <- Effort_survey[Effort_survey$YEAR!='2011',]

# How many survey sightings by year
table(Sightings_survey$YEAR)
# How much survey effort (in trackline length) is there by year
by(gLength(Effort_survey,byid = T), Effort_survey$YEAR, sum)
# Crude estimate of relative `CPUE'
(table(Sightings_survey$YEAR) / 
       by(gLength(Effort_survey,byid = T), Effort_survey$YEAR, sum))/
  min(table(Sightings_survey$YEAR) / 
        by(gLength(Effort_survey,byid = T), Effort_survey$YEAR, sum))
# Highly variable - but survey locations differ each year! This is an advantage 
# of taking a model-based approach
# 2008 was focused in region with high density

# How many WW sightings by year
table(Sightings_Quoddy_nodup$YEAR)
table(Sightings_Brier_nodup$YEAR)

source('utility_functions.R')

ggplot() +
  gg(Domain) +
  gg.spatiallines_mod(Effort_survey[Effort_survey$YEAR=='2007',], colour='red') +
  gg.spatiallines_mod(Effort_survey[Effort_survey$YEAR=='2008',], colour='black') +
  gg.spatiallines_mod(Effort_survey[Effort_survey$YEAR=='2009',], colour='blue') +
  gg.spatiallines_mod(Effort_survey_test, colour='purple')

# Note that the 2011 data (in grey) overlap with all the previous years and WW data
# Perfect test dataset for comparing models

# The 3 years of surveys have actually been conducted under 3 different surveys
# under similar (but different) protocols
unique(Effort_survey$DATASET)
# We will model differences with a unique intercept, but assume same det fun

# Are there any obvious trends with covariates bathymetry or Slope?
# Evaluate the empirical densities of each covariate at encounter locations
# and across domain
# Convert covariates to im format
Slope_im <- as(as(Slope, 'SpatialGridDataFrame'),'im')
Bathym_im <- as(as(Bathym, 'SpatialGridDataFrame'),'im')
# Convert Observation locations to ppp object
All_obs_ppp <- ppp(x=c(Sightings_survey@coords[,1],
                       Sightings_Brier_nodup@coords[,1],
                       Sightings_Quoddy_nodup@coords[,1]),
                   y=c(Sightings_survey@coords[,2],
                       Sightings_Brier_nodup@coords[,2],
                       Sightings_Quoddy_nodup@coords[,2]),
                   window = as(Domain, 'owin'))

Slope_obs <- Slope_im[All_obs_ppp]
ggplot(data=data.frame(Slope_obs = Slope_obs),
       aes(x=Slope_obs)) +
  geom_density() +
  geom_density(data=data.frame(Slope_dom = as.numeric(Slope_im$v)),
            aes(x=Slope_dom), colour='red')
# Perhaps a log transform is desired here
ggplot(data=data.frame(Slope_obs = log(Slope_obs)),
       aes(x=Slope_obs)) +
  geom_density() +
  geom_density(data=data.frame(Slope_dom = log(as.numeric(Slope_im$v))),
               aes(x=Slope_dom), colour='red')
# Perhaps there is evidence of whales preferring areas of higher slope?

Bathym_obs <- Bathym_im[All_obs_ppp]
ggplot(data=data.frame(Bathym_obs = Bathym_obs),
       aes(x=Bathym_obs)) +
  geom_density() +
  geom_density(data=data.frame(Bathym_dom = as.numeric(Bathym_im$v[Bathym_im$v<0])),
               aes(x=Bathym_dom), colour='red')
# Perhaps a log transform is desired here. Max bathym is 313
# Lets only consider water
# Flip the sign to make it depth (larger is deeper)
ggplot(data=data.frame(Bathym_obs = log(1-Bathym_obs)),
       aes(x=Bathym_obs)) +
  geom_density() +
  geom_density(data=data.frame(Bathym_dom = log(as.numeric(1-Bathym_im$v[Bathym_im$v<0]))),
               aes(x=Bathym_dom), colour='red') +
  geom_vline(xintercept = 0, colour='blue')
# Perhaps there is evidence of whales preferring areas of shallower waters?

# Ignoring effort for the time being what would a Maxent-style model tell us?
log_Slope_im <- Slope_im
log_Slope_im$v <- log(log_Slope_im$v)
log_Bathym_im <- Bathym_im
log_Bathym_im[log_Bathym_im >= 0] <- -1
log_Bathym_im$v <- log(1-log_Bathym_im$v)
Maxent_mod <- ppm(All_obs_ppp~log_Slope_im+I(log_Slope_im^2)+log_Bathym_im+I(log_Bathym_im^2),
                  data=list(log_Slope_im=log_Slope_im, log_Bathym_im=log_Bathym_im,
                            Domain=as(Domain, 'owin')),
                  subset = Domain)
summary(Maxent_mod)
dev.off()
plot(predict(Maxent_mod, ngrid=c(300,300), se=T)$estimate,superimpose=F,
     hcl.colors(60, "YlOrRd", rev = TRUE))
plot(predict(Maxent_mod, ngrid=c(300,300), se=T)$se,superimpose=F,
     hcl.colors(60, "YlOrRd", rev = TRUE))

plot(log(predict(Maxent_mod, ngrid=c(300,300), se=T)$estimate),superimpose=F,
     hcl.colors(60, "YlOrRd", rev = TRUE))
plot(log(predict(Maxent_mod, ngrid=c(300,300), se=T)$se),superimpose=F,
     hcl.colors(60, "YlOrRd", rev = TRUE))

# Save the data required for modelling
saveRDS(list(Sightings_survey=Sightings_survey,
             Sightings_survey_test=Sightings_survey_test,
             Sightings_Brier_nodup=Sightings_Brier_nodup,
             Sightings_Brier_nodup_test=Sightings_Brier_nodup_test,
             Sightings_Quoddy_nodup=Sightings_Quoddy_nodup,
             Sightings_Quoddy_nodup_test=Sightings_Quoddy_nodup_test,
             Sightings_Opp_sp=Sightings_Opp_sp,
             Effort_survey=Effort_survey,
             Effort_survey_test=Effort_survey_test,
             WW_ports=WW_ports,
             Slope=Slope,
             Bathym=Bathym,
             Dist_Brier=Dist_Brier,
             Dist_Quoddy=Dist_Quoddy,
             Can_proj=Can_proj,
             Domain=Domain,
             Domain_latlon=Domain_latlon,
             Maxent_mod=Maxent_mod),
        './Data/Modelling_Data.rds')

