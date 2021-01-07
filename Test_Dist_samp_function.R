## Test distance sampling
data(mexdolphin, package = "inlabru")
mesh <- mexdolphin$mesh
hn <- function(distance, lsig) {
  exp(-0.5 * (distance / exp(lsig))^2)
}
matern <- inla.spde2.pcmatern(mexdolphin$mesh,
                              prior.sigma = c(2, 0.01),
                              prior.range = c(50, 0.01)
)
cmp <- ~ mySPDE(main = coordinates, model = matern) +
  lsig(1) + Intercept(1)
form <- coordinates + distance ~ mySPDE +
  log(hn(distance, lsig)) +
  Intercept + log(2)
samplers <- mexdolphin$samplers
samplers@data <- data.frame(depth=samplers@data[,5])
points <- mexdolphin$points
points@data <- data.frame(distance=points@data$distance)
fit <- lgcp(
  components = cmp,
  points,
  samplers = samplers,
  domain = list(
    coordinates = mesh,
    distance = INLA::inla.mesh.1d(seq(0, 8, length.out = 30))
  ),
  formula = form
)
pxl <- pixels(mexdolphin$mesh, nx = 100, ny = 50)
pr.int <- predict(fit, pxl, ~ exp(mySPDE + Intercept))
ggplot() + gg(pr.int)

# Test again with different CRS
proj <- CRS('PROJCS["NAD27 / UTM zone 15N",
    GEOGCS["NAD27",
        DATUM["North_American_Datum_1927",
            SPHEROID["Clarke 1866",6378206.4,294.9786982139006,
                AUTHORITY["EPSG","7008"]],
            AUTHORITY["EPSG","6267"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.0174532925199433,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4267"]],
    PROJECTION["Transverse_Mercator"],
    PARAMETER["latitude_of_origin",0],
    PARAMETER["central_meridian",-93],
    PARAMETER["scale_factor",0.9996],
    PARAMETER["false_easting",500000],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AXIS["Easting",EAST],
    AXIS["Northing",NORTH],
    AUTHORITY["EPSG","26715"]]')
points2 <- spTransform(points, proj)
samplers2 <- spTransform(samplers, proj)

mesh2 <- inla.mesh.2d(boundary = spTransform(mexdolphin$ppoly,
                                             proj),
                      max.edge = 150000,
                      min.angle = 21,
                      cutoff = 120000,
                      offset = c(300000,300000))
plot(mesh2)
ggplot() + gg(mesh2) + gg(samplers2)
matern2 <- inla.spde2.pcmatern(mexdolphin$mesh,
                              prior.sigma = c(2, 0.01),
                              prior.range = c(50000, 0.01)
)
cmp2 <- ~ mySPDE(main = coordinates, model = matern2) +
  lsig(1) + Intercept(1)
fit <- lgcp(
  components = cmp2,
  points2,
  samplers = samplers2,
  domain = list(
    coordinates = mesh2,
    distance = INLA::inla.mesh.1d(seq(0, 8, length.out = 30))
  ),
  formula = form
)
pxl <- pixels(mesh2, nx = 100, ny = 50)
pr.int <- predict(fit, pxl, ~ exp(mySPDE + Intercept))
ggplot() + gg(pr.int)
