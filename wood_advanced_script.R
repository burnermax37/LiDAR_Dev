##############
library(sp)
library(raster)
library(neonUtilities)
library(geoNEON)



##Variant of get veg, used when other functions call veglist
###Function takes vst_mappingandtagging dataframe and site code, returns veg dataframe
make_veg <- function(maptag_df, apind_df, SITECODE){
  #Get vegmap with byLocTOS using veglist$mappingandtagging
  vegmap <- getLocTOS(maptag_df, 'vst_mappingandtagging')
  
  
  #Combine veglist and vegmap
  out_veg <- merge(apind_df, vegmap, by = c('individualID', 'namedLocation', 'domainID', 'siteID', 'plotID'))
  
  
  #Filter out rows with missing heighgt measurements, missing northing values, or missing easting values
  out_veg <- out_veg[which(!is.na(out_veg$height)),]
  out_veg <- out_veg[which(!is.na(out_veg$adjEasting)),]
  out_veg <- out_veg[which(!is.na(out_veg$adjNorthing)),]
  
  
  
  #Add column to each row indicating coordinates of tile in which the plot falls
  out_veg$tile_coordinates <- paste0(as.character(round1000(out_veg$adjEasting)),'_',as.character(round1000(out_veg$adjNorthing)))
  
  return(out_veg)
}






#Function takes row of plot polygons spdf and cooresponding LAS object, calculates CHM for plot
make_plot_chm <- function(plot_polygon, tile_LAS){
  plot_data <- make_buffered_las(plot_polygon, tile_LAS)
  x <- plot_polygon$easting
  y <- plot_polygon$northing
  chm <- grid_canopy(plot_data, res = 1, dsmtin()) 
  return(chm)
}







#Function takes plot polygons spdf, coords df, and tile row of files df, and makes a chm for each plot in the tile
make_plot_chm_ldf <- function(in_file, in_spdf, in_coords){
  #Get plot ids and polygons of plots falling in current tile
  plot_ids <- in_coords[in_coords$coord_String == in_file$coords,]$plotID
  coord_polys <- in_spdf[which(in_spdf$plotID %in% plot_ids),]
  
  #Read LAS
  tile_LAS<- readLAS(paste0(getwd(),'/',in_file$name))
  
  #Check number of plots falling in tile
  coord_N <- length(plot_ids)

  
  chm_list <- list()
  
  #If tile fully encloses any plots
  if(coord_N != 0){
    #For each plot
    for(i in 1:coord_N){
      poly <- coord_polys[i,]
      chm_list <- append(chm_list, list(plotID = plot_ids[i], chm = make_plot_chm(poly, tile_LAS), coords = in_file$coords))
    }
    
  }
  return(chm_list)
}


##Function takesveg, file dataframe, plot polygons, and coord dataframe to make nested list of plot chms by tile
make_tile_chms <- function(in_veg, in_files, in_coords, in_spdf){
  
  
  chm_list <- list()
  for(i in 1:nrow(in_files)){
    chm_list <- append(chm_list, list(make_plot_chm_ldf(in_files[i,],in_coords = in_coords, in_spdf = in_spdf)))
  }
  return(chm_list)
}










###Function to create boundary CHM
###Function takes file_df, plot polygon (spdf row), and entry from boundary_ldf, creates a CHM for each set of coordinates in that entry,
  ##and combines them into one CHM
make_boundary_chm <- function(boundary_entry, in_files, plot_poly){
  #Get required files
  boundary_files <- in_files[in_files$coords %in% boundary_entry$coords,]
  
  #Create plot LAS objects
  
  my_las <- readLAS(paste0(getwd(),'/',boundary_files[1,'name']))

  
  for(i in 2:nrow(boundary_files)){
    my_las <- rbind(my_las, readLAS(paste0(getwd(),'/',boundary_files[i,'name'])))
  }
  
  #Create CHM
  return(make_plot_chm(plot_polygon = plot_poly, tile_LAS = my_las))

}





















###Function takes vegetation dataframe and CHM, plots relationship between CHM heights and measured heights of trees
height_discrepency <- function(in_veg, in_chm){
  #For each chm/plot in list, extract tree heights, and compare them to measured tree heights
  #Remove rows from in_veg that lack height values

  
  #Extract veg data that falls inside CHM
  vegsub <- in_veg[which(in_veg$adjEasting >= extent(in_chm)[1] &
                        in_veg$adjEasting <= extent(in_chm)[2] &
                        in_veg$adjNorthing >= extent(in_chm)[3] & 
                        in_veg$adjNorthing <= extent(in_chm)[4]),]
  #Extract CHM data for each tree, with buffer matching coordinate uncertainty
  bufferCHM <- extract(in_chm, cbind(vegsub$adjEasting, 
                                  vegsub$adjNorthing),
                       buffer=vegsub$adjCoordinateUncertainty, 
                       fun=max)
  
  
  
  #Comparison Method Two: Tree Based
  #  Sort veg structure data by height
  vegsub <- vegsub[order(vegsub$height, decreasing = T),]
  
  #  For each tree x, we want to estimate which nearby trees may lay beneath the canopy of x. Iterate over all trees:
  #Caluclate the distance of each tree from target tree
  #Choose an estimate for canopy size, and remove trees within that distance of main tree
  
  vegfil <- vegsub
  
  for(i in 1:nrow(vegsub)){
    if(is.na(vegfil$height[i])){
      next
    }
    #Calculate distance of each tree from target tree
    dist <- sqrt((vegsub$adjEasting[i] - vegsub$adjEasting)^2+
                   (vegsub$adjNorthing[i] - vegsub$adjNorthing)^2)
    
    #Choose estimate for canopy size as fraction of height, and remove trees within that distance
    vegfil$height[which(dist<0.3*vegsub$height[i] & 
                          vegsub$height<vegsub$height[i])] <- NA
  }
  
  #Remove empty entries
  vegfil <- vegfil[which(!is.na(vegfil$height)),]
  
  
  #Exclude dead or broken trees
  vegfil <- vegfil[which(vegfil$plantStatus=="Live"),]
  
  if(length(vegfil$height) > 0){
    #Extract raster values based on filtered trees
    filterCHM <- extract(in_chm, cbind(vegfil$adjEasting, vegfil$adjNorthing),
                         buffer=vegfil$adjCoordinateUncertainty+1, fun=max)
    
    plot(filterCHM~vegfil$height, pch=20, 
         xlab="Height", ylab="Canopy height model")
    lines(c(0,50), c(0,50), col="grey")
    out <- mean(vegfil$height - filterCHM)
  }else{
    print('After filtering, there were no Non-NA values')
  }
  return(out)
}



















#####MAIN


SITECODE = 'SOAP'


#Make request
veglist <- loadByProduct(dpID='DP1.10098.001', site=SITECODE, package="basic", check.size = F)

#Make veg
veg <- make_veg(veglist$vst_mappingandtagging, veglist$vst_apparentindividual, SITECODE = SITECODE)



#Subset to get veg_coords dataframe, veg_boundary list dataframe, veg_unique, and veg_files



  #Subset coord_df to get plots for which there is vegetation structure data
veg_coords <- coord_df[coord_df$plotID %in% unique(veg$plotID),]

  #Subset boundary pseudo-dataframe, getting only entries with vegetation structure data
veg_boundary_ldf <- list()
for(entry in boundary_list){
  if(entry$plotID %in% unique(veg$plotID)){
    veg_boundary_ldf <- append(veg_boundary_ldf, list(entry))
  }
}

rm(entry)

#Make list of unique coordinates needed for all plots with vegetation structure, and subset files df for files 
veg_unique <- get_unique_coordinates(veg_coords,veg_boundary_ldf)
veg_files <- files_df[which(files_df$coords %in% veg_unique),]


#Make non_boundary plots CHM list
tile_chm_list <- make_tile_chms(veg, veg_files, veg_coords, base_plots_SPDF)


DF$height_discrepency <- NA

#For every plot
for(i in 1:nrow(DF)){
  #if plot appears in veg
  if(DF[i,'plotID'] %in% veg_coords$plotID){
    print(DF[i,'plotID'])
    #Retrieve required CHM
    chm <- NA
    for(entry in tile_chm_list){
      if(DF[i, 'plotID'] %in% entry$plotID){
        chm <- entry[['chm']]
      }
    }
    
    #Add height discrepency value of plot to plots dataframe
    DF[i,'height_discrepency'] <- height_discrepency(veg, in_chm = chm)
  }
}
 
rm(tile_chm_list)


#Create boundary CHM list


  
#For every entry in veg_boundary_ldf
for(entry in veg_boundary_ldf){
  #Create joined CHM
  b_poly <- base_plots_SPDF[base_plots_SPDF$plotID == entry$plotID,]
  b_chm <- make_boundary_chm(entry, veg_files, b_poly)
  
  #Use joined CHM to calculate height discrepancy
  DF[DF$plotID == entry$plotID,'height_discrepency'] == height_discrepency(veg, b_chm)
}

