# utility functions
gg.spatiallines_mod <- function (data, mapping = NULL, crs = NULL, ...) 
{
  if (!is.null(crs)) {
    data = spTransform(data, crs)
  }
  qq = coordinates(data)
  cnames = coordnames(data)
  if (is.null(cnames)) {
    cnames = c("x", "y")
  }
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, 
                                                     lapply(k, function(x) x[1:(nrow(x) - 1), ]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, 
                                                     lapply(k, function(x) x[2:(nrow(x)), ]))))
  colnames(sp) = cnames
  colnames(ep) = paste0("end.", cnames)
  df = data.frame(cbind(sp, ep), data@data[rep(1:dim(data)[1], 
                                                      times=sapply(qq, FUN=function(x){dim(x[[1]])[1]-1})),])
  dmap = aes_string(x = cnames[1], y = cnames[2], xend = paste0("end.", 
                                                                cnames[1]), yend = paste0("end.", cnames[2]))
  if (!is.null(mapping)) {
    dmap = modifyList(dmap, mapping)
  }
  geom_segment(data = df, mapping = dmap, ...)
}

spatial_lines_splitter <- function (data_tmp) 
{
  qq = coordinates(data_tmp)
  cnames = coordnames(data_tmp)
  if (is.null(cnames)) {
    cnames = c("x", "y")
  }
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, 
                                                     lapply(k, function(x) x[1:(nrow(x) - 1), ]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, 
                                                     lapply(k, function(x) x[2:(nrow(x)), ]))))
  colnames(sp) = cnames
  colnames(ep) = paste0("end.", cnames)
  df = data.frame(cbind(sp, ep), data_tmp@data[rep(1:dim(data_tmp)[1], 
                                               times=sapply(qq, FUN=function(x){dim(x[[1]])[1]-1})),])
  lines <- SpatialLinesDataFrame(SpatialLines(apply(df,1,function(x){Lines(list(Line(coords=rbind(as.numeric(cbind(x[1], x[2])),
                                                    as.numeric(cbind(x[3], x[4]))))),ID=sample(1:1e12, size=1))} )),
                 data=df[,-c(1:4)], match.ID = FALSE)
  lines@proj4string <- data_tmp@proj4string
  return(lines)
}

