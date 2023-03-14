#Load libraries
library(sp)
library(raster)
library(lidR)
library(rgdal)
library(downloader)
library(httr)
library(jsonlite)

# FUNCTIONS TO CALCULATE METRICS


#Given a distributed plot spatial dataframe and the LAS object of a lidar point cloud file covering the plot,
#  this function returns an LAS object containing the plot and a surrounding buffer zone

make_buffered_las <- function(plot_spdf, tile_las, buffer = 50){

  #Convert plot SPdf to same CRS as tile_las
  plot_converted <- spTransform(plot_spdf, as.character(tile_las@proj4string))
  plot_extent <- extent(plot_converted)

  
  #Define Plot center and buffer
  x <- plot_spdf$easting
  y <- plot_spdf$northing

  
  in_radius = 20
  out_radius = in_radius + buffer

  
  #Clip out plot and buffer
  data.buffered <- lidR::lasclipRectangle(tile_las, xleft = (x - out_radius), xright = (x + out_radius),
                                      ytop = (y + out_radius), ybottom = (y - out_radius))
  

  
  #Normalize buffered plot LAS
  dtm <- grid_terrain(data.buffered, 1, kriging(k = 10L))
  data.buffered <- lasnormalize(data.buffered, dtm)
  

  
  data.plot <- lidR::lasclipRectangle(data.buffered, xleft = (x - in_radius), xright = (x + in_radius),
                                      ytop = (y + in_radius), ybottom = (y - in_radius))
  

  #Remove unreliable values
  #data.plot@data$Z[data.plot@data$Z <= .5] <- NA
  
  return(data.plot)



}

#Given am normalized LAS object for a plot, this function calculates and returns diversity indices
structural_diversity_metrics <- function(data.40m, plot_spatPolyDF) {
  x <- plot_spatPolyDF$easting
  y <- plot_spatPolyDF$northing
  chm <- grid_canopy(data.40m, res = 1, dsmtin()) 
  mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 
  max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 
  rumple <- rumple_index(chm) 
  top.rugosity <- sd(chm@data@values, na.rm = TRUE) 
  cells <- length(chm@data@values) 
  chm.0 <- chm
  chm.0[is.na(chm.0)] <- 0 
  zeros <- which(chm.0@data@values == 0) 
  deepgaps <- length(zeros) 
  deepgap.fraction <- deepgaps/cells 
  cover.fraction <- 1 - deepgap.fraction 
  vert.sd <- cloud_metrics(data.40m, sd(Z, na.rm = TRUE)) 
  sd.1m2 <- grid_metrics(data.40m, sd(Z), 1) 
  sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
  Zs <- data.40m@data$Z
  Zs <- Zs[!is.na(Zs)]
  entro <- entropy(Zs, by = 1) 
  gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3)
  GFP.AOP <- mean(gap_frac$gf) 
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
  VCI.AOP <- VCI(Zs, by = 1, zmax=100) 
  out.plot <- list(x, y, mean.max.canopy.ht,max.canopy.ht, 
             rumple,deepgaps, deepgap.fraction, 
             cover.fraction, top.rugosity, vert.sd, 
             sd.sd, entro, GFP.AOP, VAI.AOP,VCI.AOP) 
  names(out.plot) <- 
    c("easting", "northing", "mean.max.canopy.ht.aop",
      "max.canopy.ht.aop", "rumple.aop", "deepgaps.aop",
      "deepgap.fraction.aop", "cover.fraction.aop",
      "top.rugosity.aop","vert.sd.aop","sd.sd.aop", 
      "entropy.aop", "GFP.AOP.aop",
      "VAI.AOP.aop", "VCI.AOP.aop")
  return(out.plot)
}


#Given a distributed plot spatial dataframe and the LAS object of a lidar point cloud covering tile with plot, 
#  this function should calculate diversity metrics for the area around the plot.


plot_diversity_metrics <- function(plot_spdf, in_LAS){
  data.distributed_plot <- make_buffered_las(plot_spdf, in_LAS)
  return(structural_diversity_metrics(data.distributed_plot, plot_spdf))
}

