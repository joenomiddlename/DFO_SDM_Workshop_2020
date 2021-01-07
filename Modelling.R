# Modelling using the general framework and inlabru
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

rgdal::set_rgdal_show_exportToProj4_warnings(FALSE)
rgdal::set_thin_PROJ6_warnings(TRUE)
options("rgdal_show_exportToProj4_warnings"="none")

bru_options_set(
  bru_verbose = TRUE,
  verbose = TRUE,
  control.inla=list(int.strategy='eb')
)

# Load data
list2env(readRDS('./Data/Modelling_Data.rds'),globalenv())

# Step 1 build the mesh
mesh_land <- inla.mesh.2d(boundary = Domain,
                               min.angle = 25,
                               max.edge = 45,
                               cutoff = 18)
ggplot() + gg(mesh_land) + gg.spatiallines_mod(Effort_survey,colour='red')
# This builds a mesh with minimum angle 21 degrees, maximum triangle edge 15km,
# minimum triangle edge 10km. Note that spatial correlations <15km will not be detectable
mesh_land$n # lots of triangle vertices can lead to slow computation time
plot(mesh_land)

# We have a potential issue - the coastline is very jagged - this could lead to
# numerical instabilities. Solution - buffer the coastline
source('utility_functions.R')
# Original
ggplot() + gg(Domain) + gg.spatiallines_mod(Effort_survey,colour='red')
# Use gSimplify to smooth the edges
ggplot() + gg(gSimplify(Domain,tol=8)) + gg.spatiallines_mod(Effort_survey,colour='red')

# Alternative is buffering
# Buffered by 5km
ggplot() + gg(gBuffer(Domain,width=5)) + gg.spatiallines_mod(Effort_survey,colour='red')
# It still doesn't contain all the survey lines
# Buffered by 10km
ggplot() + gg(gBuffer(Domain,width=10)) + gg.spatiallines_mod(Effort_survey,colour='red')
gContains(gBuffer(Domain,width=10), Effort_survey) # success - all effort lines are in domain
# Buffered by 10km and then smoothed
ggplot() + gg(gSimplify(gBuffer(Domain,width=15),tol=5)) + gg.spatiallines_mod(Effort_survey,colour='red')
gContains(gSimplify(gBuffer(Domain,width=15),tol=5), Effort_survey) # success - all effort lines are in domain

Domain_smoothed <- gSimplify(gBuffer(Domain,width=15),tol=5)

mesh_land <- inla.mesh.2d(boundary = Domain_smoothed,
                         min.angle = 25,
                         max.edge = 45,
                         cutoff = 18)
mesh_land$crs <- Can_proj
plot(mesh_land)
mesh_land$n
ggplot() + #gg(gSimplify(gBuffer(Domain,width=15),tol=3)) +
  gg(mesh_land) + 
  gg.spatiallines_mod(Effort_survey,colour='red')

# remove all edges of the survey lines that fall outside of the boundary
# Effort_survey_clip <- gIntersection(Effort_survey, Domain_smoothed, byid = T)
# Effort_survey_clip <- sp::spChFIDs(Effort_survey_clip, rownames(Effort_survey@data))
# Effort_survey_clip <- SpatialLinesDataFrame(Effort_survey_clip,
#                                             Effort_survey@data)
# plot(Domain_smoothed); plot(Effort_survey_clip,col='red',add=T)
# # ggplot no longer works


# Compare with a mesh that ignores the land barrier
mesh_noland <- inla.mesh.2d(boundary = Domain_smoothed,
                          min.angle = 25,
                          max.edge = c(45,60),
                          cutoff = 30,
                          offset = c(10,20))
mesh_noland$crs <- Can_proj
plot(mesh_noland)
mesh_noland$n

# Create the log depth and log slope covariates
log_Depth <- Bathym
log_Depth$FIWH_MAR_Bathymetry[log_Depth$FIWH_MAR_Bathymetry >= 0] <- 0
log_Depth$FIWH_MAR_Bathymetry <- log(1-log_Depth$FIWH_MAR_Bathymetry)
names(log_Depth) <- 'log_Depth' 
# Need to buffer the covariate a little to help with computation
pixels_Domain <- pixels( mesh_land,mask=gSimplify(gBuffer(Domain,width=19),tol=3),
                         nx=400,ny=400)
pixels_Domain$log_Depth <- over(pixels_Domain,log_Depth)$log_Depth
# impute missing values with the nearest neighbour value
pixels_Domain$log_Depth[is.na(pixels_Domain$log_Depth)] <- 
  log_Depth$log_Depth[nncross(as(SpatialPoints(pixels_Domain@coords[which(is.na(pixels_Domain$log_Depth)),]),'ppp'),
                          as(SpatialPoints(log_Depth@coords),'ppp'),
                          what = 'which')]
log_Depth <- pixels_Domain
log_Depth_sq <- log_Depth
names(log_Depth_sq) <- 'log_Depth_sq'
log_Depth_sq$log_Depth_sq <- log_Depth_sq$log_Depth_sq^2
# Repeat for slope
log_Slope <- Slope
log_Slope$FIWH_MAR_Slope <- log(1+log_Slope$FIWH_MAR_Slope)
names(log_Slope) <- 'log_Slope'
pixels_Domain <- pixels( mesh_land,mask=gSimplify(gBuffer(Domain,width=19),tol=3),
                         nx=400,ny=400)
pixels_Domain$log_Slope <- over(pixels_Domain,log_Slope)$log_Slope
pixels_Domain$log_Slope[is.na(pixels_Domain$log_Slope)] <- 
  log_Slope$log_Slope[nncross(as(SpatialPoints(pixels_Domain@coords[which(is.na(pixels_Domain$log_Slope)),]),'ppp'),
                              as(SpatialPoints(log_Slope@coords),'ppp'),
                              what = 'which')]
log_Slope <- pixels_Domain
log_Slope_sq <- log_Slope
names(log_Slope_sq) <- 'log_Slope_sq'
log_Slope_sq$log_Slope_sq <- log_Slope_sq$log_Slope_sq^2

# Plot the covariates with mesh overlayed
ggplot() + gg(log_Depth) + gg(mesh_land)
ggplot() + gg(log_Slope) + gg(mesh_land)
# Problem - we need to give values of the covariates on boundary
# Note that this won't affect model fit as integration only performed on tracklines

# # Step 1 create SpatialPixelsDataFrame across expanded domain
# pixels_domain <- pixels(mesh_land, nx=400, ny=400)
# pixels_domain$log_Depth <- 0
# gD
# pixels_domain[gWithin(as(pixels_domain,'SpatialPoints'), as(log_Depth,'SpatialPolygonsDataFrame'), byid = T),] 

### Fitting a model to the survey data

# First define the SPDE model
matern_land <- inla.spde2.pcmatern(mesh_land, 
                              prior.sigma = c(2, 0.01), 
                              prior.range = c(35, 0.01))
matern_noland <- inla.spde2.pcmatern(mesh_noland, 
                                   prior.sigma = c(2, 0.01), 
                                   prior.range = c(35, 0.01))
# Define a half-normal detection probability function (code from inlabru tutorials)
log_hn = function(distance, lsig){ 
  -0.5*(distance/exp(lsig))^2
}

# (taken from inlabru) We need to now separately define the components of the model 
# (the SPDE, the Intercept and the detection function parameter `lsig`)
cmp_land1 <- ~ mySpatial(main = coordinates, model = matern_land) + 
  lsig(1) + Intercept(1)

cmp_noland1 <- ~ mySpatial(main = coordinates, model = matern_noland) + 
  lsig(1) + Intercept(1)

# and the formula, which describes how these components are combined to 
# form the linear predictor (remembering that we need an offset due to the unknown
# direction of the detections!):

formula_1 = coordinates + distance ~ mySpatial +
  log_hn(distance, lsig) + 
  Intercept + log(2)

# Then fit the model, passing both the components and the formula
# and specify integration domains for the spatial and distance dimensions:
# Note that we are going to threshold our upper detection distances at 2km
hist(Sightings_survey$DISTANCE)
Sightings_survey$DISTANCE[Sightings_survey$DISTANCE>2000 & !is.na(Sightings_survey$DISTANCE)] <- 
  rnorm(sum(Sightings_survey$DISTANCE>2000, na.rm=T),mean=1999,sd=0.1)
# needs renaming
Sightings_survey$distance <- Sightings_survey$DISTANCE
Sightings_survey$distance[is.na(Sightings_survey$distance)] <- mean(Sightings_survey$distance,na.rm=T)
Sightings_survey <- Sightings_survey[,-c(8)]
Sightings_survey$distance <- Sightings_survey$distance / 1000

# The survey effort lines need splitting into length 2 lines
# This will let inlabru easily crop tracklines to fit mesh
Effort_survey_split <- spatial_lines_splitter(Effort_survey)
ggplot() + gg(Domain) + gg(Effort_survey_split)
ggplot() + gg(Domain) + gg.spatiallines_mod(Effort_survey)

bru_options_set(
  bru_verbose = TRUE,
  verbose = TRUE,
  control.inla=list(int.strategy='eb')
)
# inlabru removes all edges of the survey lines that fall outside mesh

fit = like(data = Sightings_survey,
           samplers = Effort_survey_split,
           domain = list(
             coordinates = mesh_land,
             distance = INLA::inla.mesh.1d(seq(0, 2, length.out = 30))),
           formula = formula_1,
           family = 'cp',
           options = list(bru_initial = list(lsig = -1)))
fit <- bru(fit, components = cmp_land1) # can be slow
fit <- readRDS('./Model_Files/fit.rds')
fit$dic$dic #DIC is 840
# Alternatively call lgcp function
# fit = lgcp(components = cmp_land1, 
#            data = Sightings_survey,
#            samplers = Effort_survey,
#            domain = list(
#              coordinates = mesh_land,
#              distance = INLA::inla.mesh.1d(seq(0, 2, length.out = 30))),
#            formula = formula_1,
#            options = list(bru_initial = list(lsig = -1),
#                           bru_verbose = 3))

# predict the intensity from this model
# Use only 20 samples for demonstration
# Increase in real applications!
plot_pixels <- pixels(mesh_land, mask = Domain)
pred_int <- predict(fit, plot_pixels, ~ exp(mySpatial + Intercept), n.samples = 20)
colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(...))
}

ggplot() + gg(pred_int[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
  colsc(pred_int@data$mean)
ggplot() + gg(pred_int[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
  colsc(pred_int@data$sd)

# Look at parameters of spde model
spde.range <- spde.posterior(fit, "mySpatial", what = "range")
plot(spde.range)
spde.var <- spde.posterior(fit, "mySpatial", what = "variance")
plot(spde.var)

# Plot detection probability function
distdf <- data.frame(distance = seq(0, 2, length = 100))
dfun <- predict(fit, distdf, ~ exp(log_hn(distance, lsig)), n.samples = 20)
plot(dfun)

# We can look at the posterior for expected number of whales:
predpts <- ipoints(Domain, mesh_land)
Lambda <- predict(fit, predpts, ~ sum(weight * exp(mySpatial + Intercept)), n.samples = 20)
Lambda
#around 700 individuals with a large variance. We want more MC samples here

# Note that we are extrapolating far away from our observations! This could prove
# disastrous for models with covariates. Shrink our prediction area to lie within
# 30km of the nearest trackline
predpts_restricted <- predpts[
  which(apply(gWithinDistance(predpts,Effort_survey,dist = 30,byid = T),2,sum)>0),]
Lambda_restricted <- predict(fit, predpts_restricted, ~ sum(weight * exp(mySpatial + Intercept)), n.samples = 20)
Lambda_restricted
### Adding covariates

# Define the new components and object
# Add the log_Depth variable 
cmp_land2 <- ~ mySpatial(main = coordinates, model = matern_land) + 
  Intercept(1) + log_Depth(main = log_Depth, model='linear') +
  lsig(1)

formula_2 = coordinates + distance ~ mySpatial +
  Intercept + log_Depth +
  log_hn(distance, lsig) + 
  log(2)

fit2 = like(data = Sightings_survey,
            samplers = Effort_survey_split,
            domain = list(
              coordinates = mesh_land,
              distance = INLA::inla.mesh.1d(seq(0, 2, length.out = 30))),
            formula = formula_2,
            family = 'cp')
#fit2 <- bru(fit2, components = cmp_land2)
fit2 <- readRDS('./Model_Files/fit2.rds')
fit2$dic$dic # DIC is 840 - no improvement
summary(fit2)

# add log_Slope variable 
## Let them do this
cmp_land3 <- ~ mySpatial(main = coordinates, model = matern_land) + 
  Intercept(1) + log_Slope(main = log_Slope, model='linear') +
  lsig(1)

formula_3 = coordinates + distance ~ mySpatial +
  Intercept + log_Slope +
  log_hn(distance, lsig) + 
  log(2)

fit3 = like(data = Sightings_survey,
            samplers = Effort_survey_split,
            domain = list(
              coordinates = mesh_land,
              distance = INLA::inla.mesh.1d(seq(0, 2, length.out = 30))),
            formula = formula_3,
            family = 'cp')
#fit3 <- bru(fit3, components = cmp_land3)
fit3 <- readRDS('./Model_Files/fit3.rds')
fit3$dic$dic # DIC is 841.6 - no improvement
summary(fit3)

# What happens if we add both variables?
cmp_land4 <- ~ mySpatial(main = coordinates, model = matern_land) + 
  Intercept(1) + log_Depth(main = log_Depth, model='linear') +
  log_Slope(main = log_Slope, model = 'linear') + 
  #log_Depth_sq(main = log_Depth_sq, model='linear') +
  #log_Slope_sq(main = log_Slope_sq, model = 'linear') +
  lsig(1)

formula_4 = coordinates + distance ~ mySpatial +
  Intercept + log_Depth + #log_Depth_sq +
  log_Slope + #log_Slope_sq +
  log_hn(distance, lsig) + 
  log(2)

fit4 = like(data = Sightings_survey,
            samplers = Effort_survey_split,
            domain = list(
              coordinates = mesh_land,
              distance = INLA::inla.mesh.1d(seq(0, 2, length.out = 30))),
            formula = formula_4,
            family = 'cp')
#fit4 <- bru(fit4, components = cmp_land4)
fit4 <- readRDS('./Model_Files/fit4.rds')
fit4$dic$dic # DIC is 840 - no improvement
summary(fit4)

# Extra work for the keeners
# Notice that the log_Depth_Sq variable leads to instability!
# Fortunately this is detected in a higher DIC
# Notice how wild the predicted population size becomes!

# predict the intensity from the best covariate model (fit2)
## Let them do this
plot_pixels2 <- pixels(mesh_land, mask = Domain)

pred_int2 <- predict(fit2, plot_pixels2, 
                     ~ exp(mySpatial + Intercept + log_Depth), n.samples = 20)
pred_field2 <- predict(fit2, plot_pixels2, 
                     ~ exp(mySpatial), n.samples = 20)
colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(...))
}
multiplot(ggplot() + gg(pred_int2[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int2@data$mean,pred_int@data$mean)) + ggtitle('Covariate Model Mean'),
          ggplot() + gg(pred_int2[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int2@data$sd,pred_int@data$sd)) + ggtitle('Covariate Model SD'),
          ggplot() + gg(pred_int[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int2@data$mean,pred_int@data$mean)) + ggtitle('Covariate-Free Model Mean'),
          ggplot() + gg(pred_int[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int2@data$sd,pred_int@data$sd)) + ggtitle('Covariate-Free Model SD'),
          layout = matrix(c(1:4),nrow=2,ncol=2,byrow = F))

# Look at parameters of spde model
spde.range2 <- spde.posterior(fit2, "mySpatial", what = "range")
spde.var2 <- spde.posterior(fit2, "mySpatial", what = "variance")
multiplot(plot(spde.range2)+ylim(range(spde.range$pdf)),plot(spde.range)+xlim(range(spde.range2$range)),
          plot(spde.var2)+ylim(range(spde.var$pdf)),plot(spde.var)+xlim(range(spde.var2$variance)),
          layout = matrix(1:4,2,2,byrow = F))
ggplot() + gg(pred_field2[1])

# plot the covariate effects
depthdf <- data.frame(log_Depth = seq(from=min(log_Depth@data$log_Depth),
                                      to=max(log_Depth@data$log_Depth),
                                      length.out = 100))
depthpred <- predict(fit2, depthdf, ~ exp(log_Depth))
plot(depthpred)
# Manually if the function breaks
samples <- generate(fit2)
# depthpred <- sapply(samples,FUN = function(x){return(depthdf$log_Depth*x$log_Depth)})
# ggplot(data.frame(mean=apply(depthpred,1, mean),
#                   LCL=apply(depthpred,1, quantile, probs=c(0.025)),
#                   UCL=apply(depthpred,1, quantile, probs=c(0.975)),
#                   logdepth=seq(from=min(log_Depth@data$log_Depth),
#                             to=max(log_Depth@data$log_Depth),
#                             length.out = 100)),
#        aes(y=mean,x=logdepth,ymax=UCL,ymin=LCL)) +
#   geom_line() + geom_ribbon(alpha=0.4)

# Plot detection probability function
distdf <- data.frame(distance = seq(0, 2, length = 100))
dfun2 <- predict(fit2, distdf, ~exp(log_hn(distance,lsig)),n.samples = 20)
# If this calls an error run below

dfun2 <- sapply(samples,FUN = function(x){return(exp(log_hn(distdf$distance,x$lsig)))})
ggplot(data.frame(mean=apply(dfun2,1, mean),
                  LCL=apply(dfun2,1, quantile, probs=c(0.025)),
                  UCL=apply(dfun2,1, quantile, probs=c(0.975)),
                  distance=distdf$distance),
       aes(y=mean,x=distance,ymax=UCL,ymin=LCL)) +
  geom_line() + geom_ribbon(alpha=0.4)

# We can look at the posterior for expected number of whales:
predpts <- ipoints(Domain, mesh_land)
# It looks like we do not need to do this!

Lambda2 <- predict(fit2, predpts, ~ sum(weight * exp(mySpatial + Intercept +
                                                       log_Depth)))
Lambda2
Lambda
# 625 (421 - 920) individuals

# Avoid extrapolation
Lambda2_restricted <- predict(fit2, predpts_restricted, 
                              ~ sum(weight * exp(mySpatial + Intercept +
                                                 log_Depth)))
Lambda2_restricted
Lambda_restricted

# The covariate model estimates fewer whales. Credible intervals overlap.

## Now add the whale watch sightings

# Assume that effort decreases with distance from port
ggplot() + gg(Dist_Brier) + gg(mesh_land)
# Buffer
pixels_Domain <- pixels( mesh_land,mask=gSimplify(gBuffer(Domain,width=19),tol=3),
                         nx=400,ny=400)
pixels_Domain$Dist_Brier <- over(pixels_Domain,Dist_Brier)$Dist_Brier
pixels_Domain$Dist_Brier[is.na(pixels_Domain$Dist_Brier)] <- 
  Dist_Brier$Dist_Brier[nncross(as(SpatialPoints(pixels_Domain@coords[which(is.na(pixels_Domain$Dist_Brier)),]),'ppp'),
                              as(SpatialPoints(Dist_Brier@coords),'ppp'),
                              what = 'which')]
Dist_Brier <- pixels_Domain
ggplot() + gg(Dist_Brier) + gg(mesh_land)
Dist_Brier$Dist_Brier[is.infinite(Dist_Brier$Dist_Brier)] <- 0
max(Dist_Brier$Dist_Brier)
Dist_Brier$Dist_Brier <- Dist_Brier$Dist_Brier / 980.7996
# repeat for distance from Quoddy port
pixels_Domain <- pixels( mesh_land,mask=gSimplify(gBuffer(Domain,width=19),tol=3),
                         nx=400,ny=400)
pixels_Domain$Dist_Quoddy <- over(pixels_Domain,Dist_Quoddy)$Dist_Quoddy
pixels_Domain$Dist_Quoddy[is.na(pixels_Domain$Dist_Quoddy)] <- 
  Dist_Quoddy$Dist_Quoddy[nncross(as(SpatialPoints(pixels_Domain@coords[which(is.na(pixels_Domain$Dist_Quoddy)),]),'ppp'),
                                as(SpatialPoints(Dist_Quoddy@coords),'ppp'),
                                what = 'which')]
Dist_Quoddy <- pixels_Domain
Dist_Quoddy$Dist_Quoddy[is.infinite(Dist_Quoddy$Dist_Quoddy)] <- 0
Dist_Quoddy$Dist_Quoddy <- Dist_Quoddy$Dist_Quoddy / 980.7996
ggplot() + gg(Dist_Quoddy) + gg(mesh_land)
Dist_Quoddy_sq <- pixels( mesh_land,mask=gSimplify(gBuffer(Domain,width=19),tol=3),
                          nx=400,ny=400)
Dist_Brier_sq <- Dist_Quoddy_sq
Dist_Quoddy_sq$Dist_Quoddy_sq <- Dist_Quoddy$Dist_Quoddy^2
Dist_Brier_sq$Dist_Brier_sq <- Dist_Brier$Dist_Brier^2
# Area increases with squared radius. Create a quadratically decreasing covariate
# Use the half-norm function again

# Plot the sightings with distance from the port
hist(over(Sightings_Brier_nodup,Dist_Brier)$Dist_Brier)
hist(over(Sightings_Quoddy_nodup,Dist_Quoddy)$Dist_Quoddy)

# Define the components and formulae
# First model assumption lets the data inform the effort quantity AND
# the shape of the detection probability surface

formula_WW_B1 = coordinates ~ mySpatial_Brier +
  Dist_Brier_sq + 
  Intercept_WW_Q +
  Diff_WW_Q_B +
  log_Depth

fit_WW_B1 = like(data = Sightings_Brier_nodup,
           domain = list(
             coordinates = mesh_land),
           formula = formula_WW_B1,
           family = 'cp')

# and the formula, which describes how these components are combined to 
# form the linear predictor (remembering that we need an offset due to the unknown
# direction of the detections!):

formula_WW_Q1 = coordinates ~ mySpatial_Quoddy +
  Dist_Quoddy_sq + 
  Intercept_WW_Q +
  log_Depth

fit_WW_Q1 = like(data = Sightings_Quoddy_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_Q1,
                 family = 'cp')

# Use brew to combine both like objects! Note that we need to define a joint 
# set of components for this
# Repeat for Quoddy link
cmp_WW_Combined1 <- ~ 
  mySpatial(main = coordinates, model = matern_land) + 
  mySpatial_Brier(main = coordinates, copy='mySpatial', fixed=T) + 
  mySpatial_Quoddy(main = coordinates, copy='mySpatial', fixed=T) + 
  Dist_Quoddy_sq(main=Dist_Quoddy_sq, model='linear') +
  Dist_Brier_sq(main=Dist_Brier_sq, model='linear') +
  Diff_WW_Q_B(1) +
  Intercept_WW_Q(1) +
  Intercept(1) +
  lsig(1) +
  log_Depth(main = log_Depth, model='linear')
# How many boats do each company have in their fleet?
unique(Sightings_Brier_nodup$N_VESSELS) # 3
unique(Sightings_Quoddy_nodup$N_VESSELS) # 1

# Next model consider fixing the total effort of Brier equal to 3.22 times that of Quoddy
# Achieve this by adding log(3.22) offset
# We need to define an integration function to integrate the distance function
# with the current guess of Dist_Brier_sq 
int_dist_Brier <- function(par){
  return(-log(sum(exp(Dist_Brier_sq@data[,1]*par))/2000))
}
int_dist_Quoddy <- function(par){
  return(-log(sum(exp(Dist_Quoddy_sq@data[,1]*par))/2000))
}
# Divide by 2000 to put numbers closer to 1

formula_WW_B2 = coordinates ~ mySpatial_Brier +
  Dist_Brier_sq*Brier_par + 
  Intercept_WW_Q + log(3.22) +
  log_Depth +
  int_dist_Brier(Brier_par)

formula_WW_Q2 = coordinates ~ mySpatial_Quoddy +
  Dist_Quoddy_sq*Quoddy_par + 
  Intercept_WW_Q  +
  log_Depth +
  int_dist_Quoddy(Quoddy_par)

fit_WW_B2 = like(data = Sightings_Brier_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_B2,
                 family = 'cp',
                 allow_latent = T)

fit_WW_Q2 = like(data = Sightings_Quoddy_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_Q2,
                 family = 'cp',
                 allow_latent = T)
# Define the components
cmp_WW_Combined2 <- ~ mySpatial_Brier(main = coordinates, copy='mySpatial', fixed=T) + 
  mySpatial_Quoddy(main = coordinates, copy='mySpatial', fixed=T) + 
  Dist_Quoddy_sq(main=Dist_Quoddy_sq, model='offset') +
  Dist_Brier_sq(main=Dist_Brier_sq, model='offset') +
  Intercept_WW_Q(1) +
  mySpatial(main = coordinates, model = matern_land) + 
  Intercept(1) +
  lsig(1) +
  Quoddy_par(1) +
  Brier_par(1) +
  log_Depth(main = log_Depth, model='linear')

# Third model - put an informative prior on the scale factor of effort between
# Brier and Quoddy, centered at 3
formula_WW_B3 = coordinates ~ mySpatial_Brier +
  Dist_Brier_sq*Brier_par + 
  Intercept_WW_Q + 
  Diff_QQ_B_Q +
  log_Depth +
  int_dist_Brier(Brier_par)

formula_WW_Q3 = coordinates ~ mySpatial_Quoddy +
  Dist_Quoddy_sq*Quoddy_par + 
  Intercept_WW_Q +
  log_Depth +
  int_dist_Quoddy(Quoddy_par)

fit_WW_B3 = like(data = Sightings_Brier_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_B3,
                 family = 'cp')

fit_WW_Q3 = like(data = Sightings_Quoddy_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_Q3,
                 family = 'cp')
# Define the components
cmp_WW_Combined3 <- ~ mySpatial_Brier(main = coordinates, copy='mySpatial', fixed=T) + 
  mySpatial_Quoddy(main = coordinates, copy='mySpatial', fixed=T) + 
  Dist_Quoddy_sq(main=Dist_Quoddy_sq, model='offset') +
  Dist_Brier_sq(main=Dist_Brier_sq, model='offset') +
  Intercept_WW_Q(1) +
  Diff_QQ_B_Q(main=1,model='linear',mean.linear=1.17, prec.linear=100) +
  mySpatial(main = coordinates, model = matern_land) + 
  Intercept(1) +
  lsig(1) +
  Quoddy_par(1) +
  Brier_par(1) +
  log_Depth(main = log_Depth, model='linear')

## Now fit these joint models!
fit2_like = like(data = Sightings_survey,
            samplers = Effort_survey_split,
            domain = list(
              coordinates = mesh_land,
              distance = INLA::inla.mesh.1d(seq(0, 2, length.out = 30))),
            formula = formula_2,
            family = 'cp')

fit_combined1 <- bru(fit2_like, fit_WW_B1, fit_WW_Q1, components = cmp_WW_Combined1)
fit_combined1 <- readRDS('./Model_Files/fit_combined1.rds')
fit_combined1$dic$dic #DIC is 2259 - can't compare with previous DIC's
summary(fit_combined1)
# Plot the distance function for the two ports
# Note that there are 3 vessels from Brier vs 1 from Quoddy
# We have data June-August. 
# Quoddy does 2 trips per day @ 2.75hours and 3 trips from mid July (10th)
# Brier has 2 x "Mega Nova" 4 hour trips and 3 from July (11th)-August. 
# Also has 5 x Zodiac 2.25 hour trips and 6 in mid July (11th) - August
# Quoddy has 2.75 *( 2 * (30 + 9) + 3 * (22 + 31) ) = 651.75 boat hours per year
# Brier has 4 * (2 * (30 + 10) + 3 * (21 + 31)) + 2.25 * (5 * (30 + 10) + 6 * (21 + 31) )
# = 2096 boat hours
# Note that 2096 / 651.75 = 3.22 times more effort
# Expect Diff_WW_Q_B to equal ~log(3.22). It equals ~ -log(3)!
Dist_Brier_sq_plot <- pixels(mesh=mesh_land,mask=Domain)
Dist_Brier_sq_plot <- predict(fit_combined1,Dist_Brier_sq_plot,~exp(Dist_Brier_sq),n.samples = 20)
ggplot() + gg(Dist_Brier_sq_plot[1])
Dist_Quoddy_sq_plot <- pixels(mesh=mesh_land,mask=Domain)
Dist_Quoddy_sq_plot <- predict(fit_combined1,Dist_Quoddy_sq_plot,~exp(Dist_Quoddy_sq),n.samples = 20)
ggplot() + gg(Dist_Quoddy_sq_plot[1])

plot_pixels_WW1 <- pixels(mesh_land, mask = Domain)

pred_int_WW1 <- predict(fit_combined1, plot_pixels_WW1, 
                     ~ exp(mySpatial + Intercept + log_Depth), n.samples = 20)
pred_field_WW1 <- predict(fit_combined1, plot_pixels_WW1, 
                       ~ exp(mySpatial), n.samples = 20)
colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(...))
}
multiplot(ggplot() + gg(pred_int2[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int2@data$mean,pred_int_WW1@data$mean)) + ggtitle('Covariate Model Mean'),
          ggplot() + gg(pred_int2[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int2@data$sd,pred_int_WW1@data$sd)) + ggtitle('Covariate Model SD'),
          ggplot() + gg(pred_int_WW1[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int2@data$mean,pred_int_WW1@data$mean)) + ggtitle('Joint Model Mean'),
          ggplot() + gg(pred_int_WW1[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int2@data$sd,pred_int_WW1@data$sd)) + ggtitle('Joint Model SD'),
          layout = matrix(c(1:4),nrow=2,ncol=2,byrow = F))

# We may be extrapolating too far here. Let's restrict our plot to those pixels within
# 30 km of nearest trackline
# compute the indices of the pixels that lie < 30km of nearest trackline
restricted_ind <- which(apply(gWithinDistance(as(pred_int2,'SpatialPoints'),
                                              Effort_survey,dist = 30,byid = T),2,sum)>0)

pred_int2_restricted <- pred_int2[restricted_ind,]
pred_int_WW1_restricted <- pred_int_WW1[restricted_ind,]

multiplot(ggplot() + gg(pred_int2_restricted[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int2_restricted@data$mean,pred_int_WW1_restricted@data$mean),probs = c(0.000,1))) + ggtitle('Covariate Model Mean'),
          ggplot() + gg(pred_int2_restricted[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int2_restricted@data$sd,pred_int_WW1_restricted@data$sd),probs = c(0.001,1))) + ggtitle('Covariate Model SD'),
          ggplot() + gg(pred_int_WW1_restricted[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int2_restricted@data$mean,pred_int_WW1_restricted@data$mean),probs = c(0.001,1))) + ggtitle('Joint Model Mean'),
          ggplot() + gg(pred_int_WW1_restricted[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int2_restricted@data$sd,pred_int_WW1_restricted@data$sd),probs = c(0.001,1))) + ggtitle('Joint Model SD'),
          layout = matrix(c(1:4),nrow=2,ncol=2,byrow = F))

# What is the estimated population size?
Lambda_WW1 <- predict(fit_combined1, predpts, ~ sum(weight * exp(mySpatial + Intercept +
                                                                 log_Depth)), n.samples = 20)
Lambda_WW1
Lambda_df <- rbind(Lambda_WW1,Lambda,Lambda2)
Lambda_df$Model <- c('WW1','Covariate-Free','Covariate')
ggplot(Lambda_df,aes(x=Model, y=mean, ymax=q0.975, ymin=q0.025)) +
  geom_point() + geom_errorbar()
# and in the restricted region?
Lambda_WW1_restricted <- predict(fit_combined1, predpts_restricted, ~ sum(weight * exp(mySpatial + Intercept +
                                                                              log_Depth)), n.samples = 20)
Lambda_WW1_restricted
Lambda_df_restricted <- rbind(Lambda_WW1_restricted,Lambda_restricted,Lambda2_restricted)
Lambda_df_restricted$Model <- c('WW1','Covariate-Free','Covariate')
ggplot(Lambda_df_restricted,aes(x=Model, y=mean, ymax=q0.975, ymin=q0.025)) +
  geom_point() + geom_errorbar()

### FIT MODEL 2

fit_combined2 <- bru(fit2_like, fit_WW_B2, fit_WW_Q2, components = cmp_WW_Combined2,
                     options = list(bru_initial = list(Brier_par = -100,
                                                       Quoddy_par = -100)))
fit_combined2 <- readRDS('./Model_Files/fit_combined2.rds')
fit_combined2$dic$dic #DIC is 2326.411 - substantially worse fit.
summary(fit_combined2)

Dist_Brier_sq_plot2 <- pixels(mesh=mesh_land,mask=Domain)
Dist_Brier_sq_plot2 <- predict(fit_combined2,Dist_Brier_sq_plot2,~exp(Dist_Brier_sq*Brier_par),n.samples = 20)
ggplot() + gg(Dist_Brier_sq_plot2[1]) + colsc(Dist_Brier_sq_plot2[1]$mean)
Dist_Quoddy_sq_plot2 <- pixels(mesh=mesh_land,mask=Domain)
Dist_Quoddy_sq_plot2 <- predict(fit_combined2,Dist_Quoddy_sq_plot2,~exp(Dist_Quoddy_sq*Quoddy_par),n.samples = 20)
ggplot() + gg(Dist_Quoddy_sq_plot2[1]) + colsc(Dist_Quoddy_sq_plot2[1]$mean)

plot_pixels_WW2 <- pixels(mesh_land, mask = Domain)

pred_int_WW2 <- predict(fit_combined2, plot_pixels_WW2, 
                        ~ exp(mySpatial + Intercept + log_Depth), n.samples = 20)
pred_field_WW2 <- predict(fit_combined2, plot_pixels_WW2, 
                          ~ exp(mySpatial), n.samples = 20)

multiplot(ggplot() + gg(pred_int_WW1[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int_WW1@data$mean,pred_int_WW2@data$mean),probs = c(0.005,0.995))) + ggtitle('Joint Model 1 Mean'),
          ggplot() + gg(pred_int_WW1[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int_WW1@data$sd,pred_int_WW2@data$sd),probs = c(0.005,0.995))) + ggtitle('Joint Model 1 SD'),
          ggplot() + gg(pred_int_WW2[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int_WW1@data$mean,pred_int_WW2@data$mean),probs = c(0.005,0.995))) + ggtitle('Joint Model 2 Mean'),
          ggplot() + gg(pred_int_WW2[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int_WW1@data$sd,pred_int_WW2@data$sd),probs = c(0.005,0.995))) + ggtitle('Joint Model 2 SD'),
          layout = matrix(c(1:4),nrow=2,ncol=2,byrow = F))

# We may be extrapolating too far here. Let's restrict our plot to those pixels within
# 30 km of nearest trackline

pred_int_WW2_restricted <- pred_int_WW2[restricted_ind,]
multiplot(ggplot() + gg(pred_int_WW1_restricted[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int_WW1_restricted@data$mean,pred_int_WW2_restricted@data$mean)) + ggtitle('Joint Model 1 Mean'),
          ggplot() + gg(pred_int_WW1_restricted[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int_WW1_restricted@data$sd,pred_int_WW2_restricted@data$sd)) + ggtitle('Joint Model 1 SD'),
          ggplot() + gg(pred_int_WW2_restricted[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int_WW1_restricted@data$mean,pred_int_WW2_restricted@data$mean)) + ggtitle('Joint Model 2 Mean'),
          ggplot() + gg(pred_int_WW2_restricted[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int_WW1_restricted@data$sd,pred_int_WW2_restricted@data$sd)) + ggtitle('Joint Model 2 SD'),
          layout = matrix(c(1:4),nrow=2,ncol=2,byrow = F))

# What is the estimated popualation size?
Lambda_WW2 <- predict(fit_combined2, predpts, ~ sum(weight * exp(mySpatial + Intercept +
                                                                   log_Depth)), n.samples = 20)
Lambda_WW2
Lambda_df <- rbind(Lambda_WW1,Lambda_WW2,Lambda,Lambda2)
Lambda_df$Model <- c('WW1','WW2','Covariate-Free','Covariate')
ggplot(Lambda_df,aes(x=Model, y=mean, ymax=q0.975, ymin=q0.025)) +
  geom_point() + geom_errorbar()
# and in the restricted region?
Lambda_WW2_restricted <- predict(fit_combined2, predpts_restricted, ~ sum(weight * exp(mySpatial + Intercept +
                                                                                         log_Depth)), n.samples = 20)
Lambda_WW2_restricted
Lambda_df_restricted <- rbind(Lambda_WW1_restricted,Lambda_WW2_restricted,Lambda_restricted,Lambda2_restricted)
Lambda_df_restricted$Model <- c('WW1','WW2','Covariate-Free','Covariate')
ggplot(Lambda_df_restricted,aes(x=Model, y=mean, ymax=q0.975, ymin=q0.025)) +
  geom_point() + geom_errorbar()

# Finally test the model with the informative prior
fit_combined3 <- bru(fit2_like, fit_WW_B3, fit_WW_Q3, components = cmp_WW_Combined3)
fit_combined3 <- readRDS('./Model_Files/fit_combined3.rds')
fit_combined3$dic$dic #DIC is 2308 - performance in-between the other 2 models
summary(fit_combined3)

Dist_Brier_sq_plot3 <- pixels(mesh=mesh_land,mask=Domain)
Dist_Brier_sq_plot3 <- predict(fit_combined3,Dist_Brier_sq_plot3,~exp(Dist_Brier_sq*Brier_par),n.samples = 20)
ggplot() + gg(Dist_Brier_sq_plot3[1]) + colsc(Dist_Brier_sq_plot3[1]$mean) + 
  gg(Sightings_Brier_nodup,colour='green') + gg(Domain)
Dist_Quoddy_sq_plot3 <- pixels(mesh=mesh_land,mask=Domain)
Dist_Quoddy_sq_plot3 <- predict(fit_combined3,Dist_Quoddy_sq_plot3,~exp(Dist_Quoddy_sq*Quoddy_par),n.samples = 20)
ggplot() + gg(Dist_Quoddy_sq_plot3[1]) + colsc(Dist_Quoddy_sq_plot3[1]$mean) + 
  gg(Sightings_Quoddy_nodup,colour='green') + gg(Domain)

plot_pixels_WW3 <- pixels(mesh_land, mask = Domain)

pred_int_WW3 <- predict(fit_combined3, plot_pixels_WW3, 
                        ~ exp(mySpatial + Intercept + log_Depth), n.samples = 20)
pred_field_WW3 <- predict(fit_combined3, plot_pixels_WW3, 
                          ~ exp(mySpatial), n.samples = 20)

multiplot(ggplot() + gg(pred_int_WW1[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int_WW1@data$mean,pred_int_WW3@data$mean),probs = c(0.005,0.995))) + ggtitle('Joint Model 1 Mean'),
          ggplot() + gg(pred_int_WW1[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int_WW1@data$sd,pred_int_WW3@data$sd),probs = c(0.005,0.995))) + ggtitle('Joint Model 1 SD'),
          ggplot() + gg(pred_int_WW3[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int_WW1@data$mean,pred_int_WW3@data$mean),probs = c(0.005,0.995))) + ggtitle('Joint Model 3 Mean'),
          ggplot() + gg(pred_int_WW3[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(quantile(c(pred_int_WW1@data$sd,pred_int_WW3@data$sd),probs = c(0.005,0.995))) + ggtitle('Joint Model 3 SD'),
          layout = matrix(c(1:4),nrow=2,ncol=2,byrow = F))
# Again, this third model seems to be predicting unrealistically high whale densities away
# from the tracklines!

# We may be extrapolating too far here. Let's restrict our plot to those pixels within
# 30 km of nearest trackline

pred_int_WW3_restricted <- pred_int_WW3[restricted_ind,]

multiplot(ggplot() + gg(pred_int_WW1_restricted[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int_WW1_restricted@data$mean,pred_int_WW3_restricted@data$mean)) + ggtitle('Joint Model 1 Mean'),
          ggplot() + gg(pred_int_WW1_restricted[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int_WW1_restricted@data$sd,pred_int_WW3_restricted@data$sd)) + ggtitle('Joint Model 1 SD'),
          ggplot() + gg(pred_int_WW3_restricted[1]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int_WW1_restricted@data$mean,pred_int_WW3_restricted@data$mean)) + ggtitle('Joint Model 3 Mean'),
          ggplot() + gg(pred_int_WW3_restricted[2]) + gg(Domain) + gg.spatiallines_mod(Effort_survey, colour='red') +
            colsc(c(pred_int_WW1_restricted@data$sd,pred_int_WW3_restricted@data$sd)) + ggtitle('Joint Model 3 SD'),
          layout = matrix(c(1:4),nrow=2,ncol=2,byrow = F))

# What is the estimated popualation size?
Lambda_WW3 <- predict(fit_combined3, predpts, ~ sum(weight * exp(mySpatial + Intercept +
                                                                   log_Depth)), n.samples = 20)
Lambda_WW3
Lambda_df <- rbind(Lambda_WW1,Lambda_WW2,Lambda_WW3,Lambda,Lambda2)
Lambda_df$Model <- c('WW1','WW2','WW3','Covariate-Free','Covariate')
ggplot(Lambda_df,aes(x=Model, y=mean, ymax=q0.975, ymin=q0.025)) +
  geom_point() + geom_errorbar()
# and in the restricted region?
Lambda_WW3_restricted <- predict(fit_combined3, predpts_restricted, ~ sum(weight * exp(mySpatial + Intercept +
                                                                                         log_Depth)), n.samples = 20)
Lambda_WW3_restricted
Lambda_df_restricted <- rbind(Lambda_WW1_restricted,Lambda_WW2_restricted,Lambda_WW3_restricted,Lambda_restricted,Lambda2_restricted)
Lambda_df_restricted$Model <- c('WW1','WW2','WW3','Covariate-Free','Covariate')
ggplot(Lambda_df_restricted,aes(x=Model, y=mean, ymax=q0.975, ymin=q0.025)) +
  geom_point() + geom_errorbar()

# What are the effects of depth?
samples_WW1 <- generate(fit_combined1,n.samples = 20)
samples_WW2 <- generate(fit_combined2,n.samples = 20)
samples_WW3 <- generate(fit_combined3,n.samples = 20)

# What is the maximum log depth seen by transect line?
max_depth_obs <- max(over(Effort_survey,log_Depth))

depthpred_WW1 <- sapply(samples_WW1,FUN = function(x){return(exp(depthdf$log_Depth*x$log_Depth))})
ggplot(data.frame(mean=apply(depthpred_WW1,1, mean),
                  LCL=apply(depthpred_WW1,1, quantile, probs=c(0.025)),
                  UCL=apply(depthpred_WW1,1, quantile, probs=c(0.975)),
                  logdepth=seq(from=min(log_Depth@data$log_Depth),
                            to=max(log_Depth@data$log_Depth),
                            length.out = 100)),
       aes(y=mean,x=logdepth,ymax=UCL,ymin=LCL)) +
  geom_line() + geom_ribbon(alpha=0.4) +
  geom_vline(xintercept = max_depth_obs)
depthpred_WW2 <- sapply(samples_WW2,FUN = function(x){return(exp(depthdf$log_Depth*x$log_Depth))})
ggplot(data.frame(mean=apply(depthpred_WW2,1, mean),
                  LCL=apply(depthpred_WW2,1, quantile, probs=c(0.025)),
                  UCL=apply(depthpred_WW2,1, quantile, probs=c(0.975)),
                  logdepth=seq(from=min(log_Depth@data$log_Depth),
                               to=max(log_Depth@data$log_Depth),
                               length.out = 100)),
       aes(y=mean,x=logdepth,ymax=UCL,ymin=LCL)) +
  geom_line() + geom_ribbon(alpha=0.4) +
  geom_vline(xintercept = max_depth_obs)
depthpred_WW3 <- sapply(samples_WW1,FUN = function(x){return(exp(depthdf$log_Depth*x$log_Depth))})
ggplot(data.frame(mean=apply(depthpred_WW3,1, mean),
                  LCL=apply(depthpred_WW3,1, quantile, probs=c(0.025)),
                  UCL=apply(depthpred_WW3,1, quantile, probs=c(0.975)),
                  logdepth=seq(from=min(log_Depth@data$log_Depth),
                               to=max(log_Depth@data$log_Depth),
                               length.out = 100)),
       aes(y=mean,x=logdepth,ymax=UCL,ymin=LCL)) +
  geom_line() + geom_ribbon(alpha=0.4) +
  geom_vline(xintercept = max_depth_obs)


# Consider a quadratic effect of slope? Use the prior model 3
# Let them do this
formula_WW_B4 = coordinates ~ mySpatial_Brier +
  Dist_Brier_sq*Brier_par + 
  Intercept_WW_Q + 
  Diff_QQ_B_Q +
  log_Depth +
  log_Depth_sq +
  int_dist_Brier(Brier_par)

formula_WW_Q4 = coordinates ~ mySpatial_Quoddy +
  Dist_Quoddy_sq*Quoddy_par + 
  Intercept_WW_Q +
  log_Depth +
  log_Depth_sq +
  int_dist_Quoddy(Quoddy_par)

fit_WW_B4 = like(data = Sightings_Brier_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_B4,
                 family = 'cp')

fit_WW_Q4 = like(data = Sightings_Quoddy_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_Q4,
                 family = 'cp')

formula_survey_4 <- 
  coordinates + distance ~ mySpatial + Intercept + log_Depth + 
  log_hn(distance, lsig) + log(2) + log_Depth_sq

# We need to redefine the survey likelihood!
fit4_survey = like(data = Sightings_survey,
                 samplers = Effort_survey_split,
                 domain = list(
                   coordinates = mesh_land,
                   distance = INLA::inla.mesh.1d(seq(0, 2, length.out = 30))),
                 formula = formula_survey_4,
                 family = 'cp')

# Define the components
cmp_WW_Combined4 <- ~ mySpatial_Brier(main = coordinates, copy='mySpatial', fixed=T) + 
  mySpatial_Quoddy(main = coordinates, copy='mySpatial', fixed=T) + 
  Dist_Quoddy_sq(main=Dist_Quoddy_sq, model='offset') +
  Dist_Brier_sq(main=Dist_Brier_sq, model='offset') +
  Intercept_WW_Q(1) +
  Diff_QQ_B_Q(main=1,model='linear',mean.linear=1.17, prec.linear=100) +
  mySpatial(main = coordinates, model = matern_land) + 
  Intercept(1) +
  lsig(1) +
  Quoddy_par(1) +
  Brier_par(1) +
  log_Depth(main = log_Depth, model='linear') +
  log_Depth_sq(main = log_Depth_sq, model='linear')

fit_combined4 <- bru(fit4_survey, fit_WW_B4, fit_WW_Q4, components = cmp_WW_Combined4)
# crashes again. log_Depth_sq bad once again

# What about adding slope?
# Let them do this
formula_WW_B5 = coordinates ~ mySpatial_Brier +
  Dist_Brier_sq*Brier_par + 
  Intercept_WW_Q + 
  Diff_QQ_B_Q +
  log_Depth +
  log_Slope +
  int_dist_Brier(Brier_par)

formula_WW_Q5 = coordinates ~ mySpatial_Quoddy +
  Dist_Quoddy_sq*Quoddy_par + 
  Intercept_WW_Q +
  log_Depth +
  log_Slope +
  int_dist_Quoddy(Quoddy_par)

formula_survey_5 <- 
  coordinates + distance ~ mySpatial + Intercept + log_Depth + 
  log_hn(distance, lsig) + log(2) +
  log_Slope

# We need to redefine the survey likelihood!
fit5_survey = like(data = Sightings_survey,
                   samplers = Effort_survey_split,
                   domain = list(
                     coordinates = mesh_land,
                     distance = INLA::inla.mesh.1d(seq(0, 2, length.out = 30))),
                   formula = formula_survey_5,
                   family = 'cp')

fit_WW_B5 = like(data = Sightings_Brier_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_B5,
                 family = 'cp')

fit_WW_Q5 = like(data = Sightings_Quoddy_nodup,
                 domain = list(
                   coordinates = mesh_land),
                 formula = formula_WW_Q5,
                 family = 'cp')
# Define the components
cmp_WW_Combined5 <- ~ mySpatial_Brier(main = coordinates, copy='mySpatial', fixed=T) + 
  mySpatial_Quoddy(main = coordinates, copy='mySpatial', fixed=T) + 
  Dist_Quoddy_sq(main=Dist_Quoddy_sq, model='offset') +
  Dist_Brier_sq(main=Dist_Brier_sq, model='offset') +
  Intercept_WW_Q(1) +
  Diff_QQ_B_Q(main=1,model='linear',mean.linear=1.17, prec.linear=100) +
  mySpatial(main = coordinates, model = matern_land) + 
  Intercept(1) +
  lsig(1) +
  Quoddy_par(1) +
  Brier_par(1) +
  log_Depth(main = log_Depth, model='linear') +
  log_Slope(main = log_Slope, model='linear')

fit_combined5 <- bru(fit5_survey, fit_WW_B5, fit_WW_Q5, components = cmp_WW_Combined5)

# Crashes again!

# How would we fit a spatio-temporal model?
# we would need to first define integration points manually
ips <- ipoints(Effort_survey_split, mesh_land, group = "YEAR")
# Next we need to redefine the components
cmp <- LHSstuff ~ RHSstuff + mySpatial(main = coordinates, model = matern_land,
                                       group = YEAR, ngroup = 3,
                                       control.group(model='?'))
# where ? can be iid for independent spatial fields each year
# ar1 for an autoregressive separable spatio-temporal field 
# see ?control.group for possible models

# We cannot fit a spatio-temporal model to our data due to 
# insufficient data and due to the survey tracklines visiting 
# non-overlapping regions (non-randomly) each year!

saveRDS(fit,'./Model_Files/fit.rds')
saveRDS(fit2,'./Model_Files/fit2.rds')
saveRDS(fit3,'./Model_Files/fit3.rds')
saveRDS(fit4,'./Model_Files/fit4.rds')
saveRDS(fit_combined1,'./Model_Files/fit_combined1.rds')
saveRDS(fit_combined2,'./Model_Files/fit_combined2.rds')
saveRDS(fit_combined3,'./Model_Files/fit_combined3.rds')




