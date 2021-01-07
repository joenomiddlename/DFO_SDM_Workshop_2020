# Install inlabru
#remotes::install_github("fbachl/inlabru", ref="stable")
# Load the packages
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
library(tweedie)
library(measurements)
library(raster)
#library(sf)

# read in whale watch and directed research sightings
# Two of the platforms have data for only one year each 
# (Grand Manan Whale and Seabird Research Station & Quoddy Link Marine and Wildlife Cruises) 
# – the other two have data over multiple years (9 and 10 years I believe)
Sightings_DRWW <- read_excel("Data/1_MarWSDB_BpSightings_DRWW.xlsx")
# add a zero before the times before 10am
Sightings_DRWW$WS_TIME[!(substr(Sightings_DRWW$WS_TIME,1,1) %in% c('1','N'))] <- 
  paste0('0',Sightings_DRWW$WS_TIME[!(substr(Sightings_DRWW$WS_TIME,1,1) %in% c('1','N'))])
Sightings_DRWW$WS_TIME[!(substr(Sightings_DRWW$WS_TIME,1,1) %in% c('N'))] <- 
  paste0(substr(Sightings_DRWW$WS_TIME[!(substr(Sightings_DRWW$WS_TIME,1,1) %in% c('N'))],1,2),
         ':',
         substr(Sightings_DRWW$WS_TIME[!(substr(Sightings_DRWW$WS_TIME,1,1) %in% c('N'))],3,4))
Sightings_DRWW$WS_TIME <- hm(Sightings_DRWW$WS_TIME)
# are there any duplicate events
max(as.numeric(table(Sightings_DRWW$WS_EVENT_ID))) # no
View(Sightings_DRWW) 
# read in opportunistic sightings from vessels from two sources – 
# our regional database and the Ocean Biodiversity Information System
Sightings_Opp <- read_excel("Data/2_MarWSDB_OBIS_BpSightings_OV.xlsx")
# add a zero before the times before 10am
Sightings_Opp$WS_TIME[!(substr(Sightings_Opp$WS_TIME,1,1) %in% c('1','N'))] <- 
  paste0('0',Sightings_Opp$WS_TIME[!(substr(Sightings_Opp$WS_TIME,1,1) %in% c('1','N'))])
Sightings_Opp$WS_TIME[!(substr(Sightings_Opp$WS_TIME,1,1) %in% c('N'))] <- 
  paste0(substr(Sightings_Opp$WS_TIME[!(substr(Sightings_Opp$WS_TIME,1,1) %in% c('N'))],1,2),
         ':',
         substr(Sightings_Opp$WS_TIME[!(substr(Sightings_Opp$WS_TIME,1,1) %in% c('N'))],3,4))
Sightings_Opp$WS_TIME <- hm(Sightings_Opp$WS_TIME)
# look at the observer types
unique(Sightings_Opp$PLATFORM)
# change the names to 'Opportunistic'
Sightings_Opp$PLATFORM <- 'Opportunistic'
Sightings_Opp$PLATFORM_CODE <- 999
# Add arbitrary PO_LATITUDE and PO_LONGITUDE
Sightings_Opp$PO_LATITUDE <- NA
Sightings_Opp$PO_LONGITUDE <- NA
# Add 1 to N_VESSELS
Sightings_Opp$N_VESSELS <- 1
Sightings_Opp$WS_DATE <- ymd(paste0(Sightings_Opp$YEAR,'-',
                                    Sightings_Opp$MONTH,'-',
                                    Sightings_Opp$DAY))
View(Sightings_Opp) 
# read in survey data from NOAA
Sightings_NOAA <- readOGR('./Data/3_NOAA_NARW_Surveys/NOAA_NARWS_FIWHsightings.shp')
# check the coordinate reference system (CRS)
Sightings_NOAA@proj4string # latlon
Sightings_NOAA$DATE_LO <- dmy(Sightings_NOAA$DATE_LO)
# read in effort data from NOAA
Effort_NOAA <- readOGR('./Data/3_NOAA_NARW_Surveys/NOAA_NARWS_Effort.shp')
# check the CRS
Effort_NOAA@proj4string
# convert the earlier sightings data to SpatialPointsDataFrame Objects with 
# correct/consistent CRS
Sightings_DRWW_sp <- SpatialPointsDataFrame(coords = cbind(Sightings_DRWW$LONGITUDE,
                                                           Sightings_DRWW$LATITUDE),
                                            data = Sightings_DRWW,
                                            proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'))
Sightings_Opp_sp <- SpatialPointsDataFrame(coords = cbind(Sightings_Opp$LONGITUDE,
                                                          Sightings_Opp$LATITUDE),
                                           data = Sightings_Opp,
                                           proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'))

# bind the first two sightings datasets together
Sightings_nonsurvey <- 
  rbind(Sightings_DRWW_sp[,c('DAY','MONTH','YEAR','WS_TIME','COUNT',
                             'DATATYPE','DATASOURCE','PLATFORM','PO_LATITUDE',
                             'PO_LONGITUDE','N_VESSELS')], 
        Sightings_Opp_sp[,c('DAY','MONTH','YEAR','WS_TIME','COUNT',
                            'DATATYPE','DATASOURCE','PLATFORM','PO_LATITUDE',
                            'PO_LONGITUDE','N_VESSELS')])
Sightings_nonsurvey$WS_DATE <- ymd(paste0(Sightings_nonsurvey$YEAR,'-',
                                              Sightings_nonsurvey$MONTH,'-',
                                              Sightings_nonsurvey$DAY))

# Create DRWW Port SpatialPointsDataFrame
DRWW_Ports_sp <- SpatialPointsDataFrame(coords = cbind(c(-66.3488, -66.7488, -66.7488, -67.0543),
                                                       c(44.2637, 44.7619, 44.7619, 45.0714)),
                                        data = data.frame(PLATFORM_CODE=c(15,11,73,69),
                                                          COMPANY=c('BRIER ISLAND WHALEWATCH',
                                                                    'GRAND MANAN WHALE AND SEABIRD',
                                                                    'WHALES N SAILS',
                                                                    'QUODDY LINK'),
                                                          N_VESSELS=c(3,1,1,1)),
                                        proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs'))

# read in shapefile of coastline
DOMAIN <- readOGR('./Data/0_Study_Area_and_Coastline/FIWH_MAR_StudyArea.shp')
DOMAIN_latlon <- spTransform(DOMAIN, CRS('+proj=longlat +datum=WGS84 +no_defs'))
# COAST <- readOGR('./Data/0_Study_Area_And_Coastline/FIWH_MAR_Coastline.shp')
# COAST <- st_as_sf(COAST) 
# COAST <- st_polygonize(COAST)
# COAST <- as(COAST, "Spatial") # If you want sp
# class(shp_airports)
#ggmap::ggmap(get_map(apply(coordinates(Sightings_Opp_sp), MARGIN = 2, range))) + gg(DOMAIN_latlon)
gmap(Sightings_Opp_sp) + gg(DOMAIN_latlon) 

# read in the DFO survey data 
# read in survey data from NOAA
Sightings_DFO <- readOGR('./Data/4_DFO_2007TNASS/DFO_TNASS2007_FIWHsightings.shp')
# check the coordinate reference system (CRS)
Sightings_DFO@proj4string # latlon
# read in effort data from NOAA
Effort_DFO <- readOGR('./Data/4_DFO_2007TNASS/DFO_TNASS2007_Effort.shp')
# check the CRS
Effort_DFO@proj4string
Effort_DFO$Year

Sightings_DFO$Dt_tm_l[!(substr(Sightings_DFO$Dt_tm_l,12,12) %in% c('1'))] <- 
  paste0(substr(Sightings_DFO$Dt_tm_l[!(substr(Sightings_DFO$Dt_tm_l,12,12) %in% c('1'))],1,11),
         '0', 
         substr(Sightings_DFO$Dt_tm_l[!(substr(Sightings_DFO$Dt_tm_l,12,12) %in% c('1'))],12,100))
Sightings_DFO$DATE_LO <- dmy_hm(Sightings_DFO$Dt_tm_l)
Sightings_DFO$TIME_LO <- time(Sightings_DFO$DATE_LO)
Sightings_DFO$DATE_LO <- date(Sightings_DFO$DATE_LO)
Sightings_DFO$DAY <- day(Sightings_DFO$DATE_LO)
Sightings_DFO$MONTH <- month(Sightings_DFO$DATE_LO)
Sightings_DFO$YEAR <- year(Sightings_DFO$DATE_LO)
Sightings_DFO$LINE_ID <- Sightings_DFO$Trns_ID
Sightings_DFO$ALTITUD <- Sightings_DFO$Alti
Sightings_DFO$SPEED <- Sightings_DFO$Speed
Sightings_DFO$VISIBIL <- Sightings_DFO$Visblty
Sightings_DFO$BEAUFOR <- Sightings_DFO$Beaufrt
Sightings_DFO$SPCODE <- Sightings_DFO$Sp_code
Sightings_DFO$COUNT <- Sightings_DFO$Clustsz
Sightings_DFO$ANGLE <- Sightings_DFO$Rel_ang
Sightings_DFO$RID <- Sightings_DFO$Observr

Sightings_DFO_sp <- Sightings_DFO[,colnames(Sightings_DFO@data) %in% c(colnames(Sightings_NOAA@data),'Prp_dst')]

# Read in NOAA Cetacean Survey data
Sightings_NOAA_Cet <- readOGR('./Data/5_NOAA_Cetacean_Surveys/NNOAA20062017_CanFIWHsightings.shp')
# check the coordinate reference system (CRS)
Sightings_NOAA_Cet@proj4string # latlon
# read in effort data from NOAA
Effort_NOAA_Cet <- readOGR('./Data/5_NOAA_Cetacean_Surveys/NOAA20062017_CanEffort.shp')
# check the CRS
Effort_NOAA_Cet@proj4string

Sightings_NOAA_Cet$DATE_LO <- date(dmy_hms(Sightings_NOAA_Cet$DATETIM))
Sightings_NOAA_Cet$TIME_LO <- time(dmy_hms(Sightings_NOAA_Cet$DATETIM))
Sightings_NOAA_Cet$COUNT <- Sightings_NOAA_Cet$SIZEBES
Sightings_NOAA_Cet$ANGLE <- Sightings_NOAA_Cet$DECANGL

Sightings_NOAA_Cet_sp <- Sightings_NOAA_Cet[,colnames(Sightings_NOAA_Cet@data) %in% c(colnames(Sightings_NOAA@data),'PERPDIS')]

# Find a good selection of 3 years of data
# NOAA effort and sightings by year
by(Effort_NOAA$km, Effort_NOAA$YEAR, sum)
by(Sightings_NOAA$YEAR, Sightings_NOAA$YEAR, length)

# nonsurvey sightings by year
by(Sightings_nonsurvey$YEAR, Sightings_nonsurvey$YEAR, length)

# dfo sightings by year
by(Sightings_DFO$YEAR, Sightings_DFO$YEAR, length) # 2007 only!
sum(as.numeric(gLength(spTransform(Effort_DFO,CRS('EPSG:3005')), byid = T)))/1000
# use DFO data as validation

# NOAA cetaceans survey data sightings and effort by year
by(Effort_NOAA_Cet$km, Effort_NOAA_Cet$YEAR, sum)
by(Sightings_NOAA_Cet$YEAR, Sightings_NOAA_Cet$YEAR, length)

# 2007, 2008 and 2009? How about using 2011 data for validation?
Sightings_nonsurvey <- Sightings_nonsurvey[Sightings_nonsurvey$YEAR %in% c(2007, 2008, 2009, 2011),]
Sightings_DRWW_sp <- Sightings_DRWW_sp[Sightings_DRWW_sp$YEAR %in% c(2007, 2008, 2009, 2011),]
Sightings_Opp_sp <- Sightings_Opp_sp[Sightings_Opp_sp$YEAR %in% c(2007, 2008, 2009, 2011),]
Sightings_DFO_sp <- Sightings_DFO_sp[Sightings_DFO_sp$YEAR %in% c(2007, 2008, 2009, 2011),]
Sightings_NOAA <- Sightings_NOAA[Sightings_NOAA$YEAR %in% c(2007, 2008, 2009, 2011),]
Sightings_NOAA_Cet_sp <- Sightings_NOAA_Cet_sp[Sightings_NOAA_Cet_sp$YEAR %in% c(2007, 2008, 2009, 2011),]

# Note that we only have sightings from BRIER ISLAND AND QUODDY LINK WW COMPANIES

Effort_DFO$YEAR <- 2007
Effort_NOAA <- Effort_NOAA[Effort_NOAA$YEAR %in% as.character(c(2007, 2008, 2009, 2011)),]
Effort_NOAA_Cet <- Effort_NOAA_Cet[Effort_NOAA_Cet$YEAR %in% as.character(c(2007, 2008, 2009, 2011)),]
DRWW_Ports_sp <- DRWW_Ports_sp[c(1,4),]

# Add a distance variable for each survey dataset in the same units (km)
# Sightings_NOAA$DISTANCE <- (gDistance(spTransform(Effort_NOAA, CRS('EPSG:3005')),
#                                       spTransform(Sightings_NOAA, CRS('EPSG:3005')),
#                                      byid = T) / 1000)[cbind(1:dim(Sightings_NOAA)[1], 
#                                                              pmatch(Sightings_NOAA$LINE_ID,
#                                                                     Effort_NOAA$LINE_ID,
#                                                                     duplicates.ok = T))]
Sightings_NOAA$DISTANCE <- tan(conv_unit(Sightings_NOAA$ANGLE, 'degree','radian')) * Sightings_NOAA$ALTITUD # In feet
# Convert to meters
Sightings_NOAA$DISTANCE <- conv_unit(Sightings_NOAA$DISTANCE, 'ft','m')
# max sighting distance 10371m https://www.mass.gov/files/documents/2016/08/uw/rwhale00.pdf
# avg sighting distance reported in similar surveys 2482m.
Sightings_NOAA$W <- 10371

Sightings_NOAA_Cet_sp$DISTANCE <- Sightings_NOAA_Cet_sp$PERPDIS
# These are in meters, they truncated to 2.1km 
# They flew at height of 183m at 200kph

Sightings_DFO_sp$DISTANCE <- Sightings_DFO_sp$Prp_dst
# In meters. Flew at 198m at 185kmh

hist(Sightings_DFO_sp$DISTANCE)
hist(Sightings_NOAA$DISTANCE)
hist(Sightings_NOAA_Cet_sp$DISTANCE)

# Merge the survey datasets together, keeping the survey as a factor variable
Sightings_DFO_sp$DATASET <- 'DFO'
Sightings_NOAA$DATASET <- 'NOAA_1'
Sightings_NOAA_Cet_sp$DATASET <- 'NOAA_2'

surv_names <- intersect(names(Sightings_NOAA_Cet_sp),
  intersect(names(Sightings_DFO_sp),
            names(Sightings_NOAA)))
Sightings_survey <- rbind(Sightings_DFO_sp[,surv_names],
                          Sightings_NOAA[,surv_names], 
                          Sightings_NOAA_Cet_sp[,surv_names])
Sightings_survey$DATASET <- factor(Sightings_survey$DATASET)

hist(Sightings_survey$DISTANCE)

Effort_DFO$DATASET <- 'DFO'
Effort_NOAA$DATASET <- 'NOAA_1'
Effort_NOAA_Cet$DATASET <- 'NOAA_2'
Effort_survey <- rbind(Effort_DFO[,c('DATASET','YEAR')],
                       Effort_NOAA[,c('DATASET','YEAR')],
                       Effort_NOAA_Cet[,c('DATASET','YEAR')])
Effort_survey$DATASET <- factor(Effort_survey$DATASET)

# load the covariate data
# reduce the size by a factor of 4
Bathym <- aggregate(raster('./Data/Covariates/FIWH_MAR_Bathymetry.asc'), fact=4)
Bathym <- as(Bathym,'SpatialPixelsDataFrame')
Bathym@proj4string <- CRS('+proj=longlat +datum=WGS84 +no_defs')

Slope <- aggregate(raster('./Data/Covariates/FIWH_MAR_Slope.asc'), fact=4)
Slope <- as(Slope,'SpatialPixelsDataFrame')
Slope@proj4string <- CRS('+proj=longlat +datum=WGS84 +no_defs')

# Create the inla mesh
mesh <- inla.mesh.2d(boundary = gSimplify(DOMAIN_latlon,tol=0.5),
                     offset = c(0.2,0.2),
                     max.edge = c(0.2,0.4),
                     cutoff = 0.05,
                     min.angle = 21,
                     crs = CRS('+proj=longlat +datum=WGS84 +no_defs'))
mesh$n
plot(mesh)

# Create the SpatialPixelsDataFrame of distance from port (around land)
Dist_BRIER <- pixels(mesh, mask = gBuffer(DOMAIN_latlon,width=0.06))
Dist_QUODDY <- pixels(mesh, mask = gBuffer(DOMAIN_latlon,width=0.06))

DOMAIN_R <- raster(extent(gBuffer(DOMAIN_latlon,width=0.06)), nrow=100, ncol=100) ## 2nd to RasterLayer ...
DOMAIN_RR <- rasterize(DOMAIN_latlon, DOMAIN_R)                       ## ...
DOMAIN_RR[is.na(DOMAIN_RR)]<-0                               ## Set cells on land to "0"

library(gdistance)
## gdistance requires that you 1st prepare a sparse "transition matrix"
## whose values give the "conductance" of movement between pairs of
## adjacent and next-to-adjacent cells (when using directions=16)
tr1 <- transition(DOMAIN_RR, transitionFunction=mean, directions=16)
tr1 <- geoCorrection(tr1,type="c")

## Compute a matrix of pairwise distances between points
Dist_BRIER$Dist_Brier <- costDistance(tr1, Dist_BRIER@coords,DRWW_Ports_sp@coords[1,]) / 1000
Dist_QUODDY$Dist_Quoddy <- costDistance(tr1, Dist_QUODDY@coords,DRWW_Ports_sp@coords[2,]) / 1000
ggplot() + gg(Dist_BRIER) + gg(DRWW_Ports_sp[1,], color='red') + gg(DOMAIN_latlon)
ggplot() + gg(Dist_QUODDY) + gg(DRWW_Ports_sp[2,], color='red') + gg(DOMAIN_latlon)

# plot the data

source('utility_functions.R')

gmap(Sightings_nonsurvey) +
  gg(DOMAIN_latlon) +
  gg.spatiallines_mod(Effort_NOAA, colour='orange') +
  gg.spatiallines_mod(Effort_NOAA_Cet, colour='brown') +
  gg.spatiallines_mod(Effort_DFO, colour='grey') +
  gg(Sightings_nonsurvey, colour='blue') +
  gg(Sightings_NOAA, colour='black') +
  gg(DRWW_Ports_sp, colour='red') +
  gg(Sightings_NOAA_Cet_sp, colour='green') +
  gg(Sightings_DFO_sp, colour='purple') 

gmap(Sightings_nonsurvey) +
  gg(DOMAIN_latlon) +
  gg.spatiallines_mod(Effort_survey) +
  gg(Sightings_nonsurvey, colour='blue') +
  gg(Sightings_survey, colour='purple') +
  gg(DRWW_Ports_sp, colour='red')
  
ggplot() +
  gg(DOMAIN_latlon) +
  gm(Bathym) +
  gg(Sightings_nonsurvey, colour='yellow') +
  gg(Sightings_survey, colour='pink') +
  gg.spatiallines_mod(Effort_survey) +
  gg(DRWW_Ports_sp, colour='red') 

ggplot() +
  gg(DOMAIN_latlon) +
  gm(Bathym) +
  gg(Sightings_nonsurvey[Sightings_nonsurvey$DATATYPE=='WhaleWatch',], colour='yellow') +
  gg(Sightings_survey, colour='pink') +
  gg.spatiallines_mod(Effort_survey) +
  gg(DRWW_Ports_sp, colour='red')
# Save the results for later use

Compiled_Data <- saveRDS(
  list(Sightings_survey = Sightings_survey,
       #Sightings_nonsurvey = Sightings_nonsurvey,
       Sightings_Opp_sp = Sightings_Opp_sp,
       Sightings_DRWW_sp = Sightings_DRWW_sp,
       Effort_survey = Effort_survey,
       WW_ports = DRWW_Ports_sp,
       Bathym=Bathym,
       Slope=Slope,
       Dist_Quoddy = Dist_QUODDY,
       Dist_Brier = Dist_BRIER,
       Domain=DOMAIN,
       Domain_latlon=DOMAIN_latlon,
       Mesh = mesh),
  './Data/Compiled_Data.rds'
)

save.image('Read_Data.RData')
