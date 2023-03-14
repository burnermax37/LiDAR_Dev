tile_na_counts <- function(file_tile, coorddf, plot_sdf){
  
  plot_ids <- coorddf[(coorddf$coord_String == file_tile$coords),]$plotID
  
  plots <- plot_sdf[(plot_sdf$plotID %in% plot_ids),]
  N <- length(plots)
  tile_LAS<- readLAS(paste0(getwd(),'/',file_tile$name))
  
  DF <- data.frame(plotID <- rep('',N),
                   NAcount <- rep(NA,N))

  
  if(N != 0){
    for(i in seq(1,N)){
      plo_it <- make_buffered_las(plots[i,], tile_LAS)
      count <- length(plo_it@data$Z[is.na(plo_it@data$Z)])
      
      DF[i,1] <- plot_sdf[i,]$plotID
      DF[i,2] <- count
    }
  }
  
  
  return(DF)
  
}


##Function takes row from boundary plots nested list, plot polygons spDF, and file df
## returns number of NA height values in plot poitn cloud
calculate_boundary_na <- function(boundary_plot, plots_spdf, file_df){
  
  #Extract rows for required files, record number of extracted rows
  files_sub <- file_df[file_df$coords %in% boundary_plot$coords,]
  N <- length(files_sub$coords)
  
  plot_data <- plots_spdf[plots_spdf$plotID == boundary_plot$plotID,]
  
  #Get first LAS object
  merged_LAS <- readLAS(paste0(getwd(),'/',files_sub[1,]$name))
  
  #Get additional LAS objects and combine with first
  for(i in seq(2,N)){
    new_LAS <- readLAS(paste0(getwd(),'/',files_sub[i,]$name))
    merged_LAS <- rbind(merged_LAS, new_LAS)
  }
  
  #calculate las
  las <- make_buffered_las(plot_data, merged_LAS)
  count <- length(las@data$Z[is.na(las@data$Z)])
  
  out <- list(plotID = plot_data$plotID, NAcount = count)
  
  return(out)
}



build_na_df <- function(this_plots_spdf, this_coords_df, this_files_df, this_boundaries_df){
  #Initialize empty dataframe
  bound_N <- length(this_boundaries_df)
  
  DF <- data.frame(plotID = rep('',bound_N),
                   NAcount = rep(NA, bound_N)) 
  
  print('DF Made')
  #Add rows for boundary plots
  if(bound_N > 0){
    for(i in seq(1, bound_N)){
      DF[i,] <- calculate_boundary_na(this_boundaries_df[[i]], this_plots_spdf, this_files_df)
    }
  }
  
  names(DF) <- c('1','2')
  
  print('DF Filled')
  #Create metrics data for non-boundary plots by file, and merge with DF
  file_N <- length(this_files_df)
  
  for(i in seq(1,file_N)){
    new_file_DF <- tile_na_counts(this_files_df[i,],this_coords_df,this_plots_spdf)
    names(new_file_DF) <- c('1','2')
    DF <- rbind(DF, new_file_DF)
  }
  
  names(DF) <- c('plotID','NAcount')
  
  DF <- arrange(DF, plotID)
  
  return(DF)
}

base <- 'NEON_D17_SOAP_DP1_XX_classified_point_cloud_colorized.laz'

for(i in seq(1,length(coord_unique))){
  files_df[i,]$name <- gsub('XX', coord_unique[i],base)
  files_df[i,]$coords <- coord_unique[i]
}

#Function takes row of file df, coords df, and plots sp_df; calculates structural diversity metrics for each
#  plot located in tile and returns results as dataframe
assemble_tile_metrics <- function(file_tile, coorddf, plot_sdf){
  
  plot_ids <- coorddf[(coorddf$coord_String == file_tile$coords),]$plotID
  
  plots <- plot_sdf[(plot_sdf$plotID %in% plot_ids),]
  tile_LAS<- readLAS(paste0(getwd(),'/',file_tile$name))
  
  
  #Prepare empty data frame
  N <- length(plots)
  DF <- data.frame(plotID = rep('',N),
                   x = rep(NA,N),
                   y = rep(NA,N),
                   mean.max.canopy.ht = rep(NA,N),
                   max.canopy.ht = rep(NA,N), 
                   rumple = rep(NA,N),
                   deepgaps = rep(NA,N),
                   deepgap.fraction = rep(NA,N),
                   cover.fraction = rep(NA,N),
                   top.rugosity = rep(NA,N),
                   vert.sd = rep(NA,N), 
                   sd.sd = rep(NA,N),
                   entro = rep(NA,N),
                   GFP.AOP = rep(NA,N),
                   VAI.AOP = rep(NA,N),
                   VCI.AOP = rep(NA,N)) 
  
  if(N != 0){
    for(i in seq(1,N)){
      DF[i,1] <- as.character(plots[i,]$plotID)
      DF[i,(2:16)] <- plot_diversity_metrics(plots[i,], tile_LAS)
    }
  }
  
  
  return(DF)
}






NEON_all_plots <- readOGR('All_NEON_TOS_Plots_V7/All_NEON_TOS_Plot_Polygons_V7.shp')





###Function for getting number of baseplots and number of files for a NEON site
site_description <- function(sitecode, in_spdf){
  #Plot base plots
  base_plots_SPDF <- in_spdf[(in_spdf$siteID == sitecode)&(in_spdf$subtype == 'basePlot'),]

  
  #Get coordinates for all plots not on a tile boundary, then split by whether plot is on a boundary
  coord_df <- build_plot_frame(base_plots_SPDF)
  boundary_df <- coord_df[coord_df$coord_String == 'Plot crosses a tile boundary',]
  coord_df <- coord_df[coord_df$coord_String != 'Plot crosses a tile boundary',]
  
  #Get coordinates for plots on boundary tiles
  boundary_list <- build_plot_pseudo_frame(base_plots_SPDF, boundary_df)
  
  #Get unique coordinates, coordinates for each necessary tile
  coord_unique <- get_unique_coordinates(coord_df, boundary_list)
  
  #Get most recent date for which data is available
  date <- get_most_recent_date(sitecode, 'DP1.30003.001')
  
  #If data is available
  if(!is.nan(date)){
    #Get file names and urls for all necessary tiles, idetnifying tile coordinates for which files are not available
    files_df <- build_file_df(coord_unique, date, sitecode)
    
    bad_coords <- files_df[files_df$name == 'bad tile',]$coords
    files_df <- files_df[files_df$name != 'bad_tile',]
    
    
    #Make list of bad coordinates, and remove bad tiles from file df
    bad_coords <- files_df[files_df$name == 'bad tile',]$coords
    files_df <- files_df[files_df$name != 'bad tile',]
    fileNum <- nrow(files_df)
    
    #Remove any non-boundary plots with bad coordinates
    if(length(bad_coords) > 0){
      bad_plot_num <- length(coord_df[coord_df$coord_String %in% bad_coords,'coord_string'])
      if(bad_plot_num != 0){
        print(paste0('Removed ',bad_plot_num,' non-boundary plots with missing tiles')) 
      }
      
      coord_df <- coord_df[!(coord_df$coord_String %in% bad_coords),]
      
      
      #Remove any boundary plots with bad coordinates
      bad_boundaries <- length(boundary_list)
      boundary_list <- remove_bad_boundary_plots(boundary_list, bad_coords)
      
      bad_boundaries <- bad_boundaries - length(boundary_list)
      bad_plot_num <- bad_plot_num + bad_boundaries
      
      

      bad_tilesNum <- length(bad_coords)
      
    }else{
      
      bad_tilesNum <- 0
      bad_plot_num <- 0
      
    }
    
    
  }else{
    fileNum <- NA
    bad_tilesNum <- NA
    bad_plot_num <- NA
  }
  
  
  
  out <- list(siteID = sitecode, baseplotNum = length(base_plots_SPDF$plotID), fileNum = fileNum, 
              bad_tilesNum = bad_tilesNum, bad_plotsNum = bad_plot_num)
  
  return(out)

}



site_df <- data.frame(siteID = unique(NEON_all_plots$siteID),
                      basePlotNum = rep(NA, 47),
                      fileNum = rep(NA,47),
                      bad_tilesNum = rep(NA, 47),
                      bad_plotsNum = rep(NA, 47))



main_files_quick <- function(sitecode, in_spdf){
  #Plot base plots
  base_plots_SPDF <- in_spdf[(in_spdf$siteID == sitecode)&(in_spdf$subtype == 'basePlot'),]
  
  
  #Get coordinates for all plots not on a tile boundary, then split by whether plot is on a boundary
  coord_df <- build_plot_frame(base_plots_SPDF)
  boundary_df <- coord_df[coord_df$coord_String == 'Plot crosses a tile boundary',]
  coord_df <- coord_df[coord_df$coord_String != 'Plot crosses a tile boundary',]
  
  #Get coordinates for plots on boundary tiles
  boundary_list <- build_plot_pseudo_frame(base_plots_SPDF, boundary_df)
  
  #Get unique coordinates, coordinates for each necessary tile
  coord_unique <- get_unique_coordinates(coord_df, boundary_list)
  
  #Get most recent date for which data is available
  date <- get_most_recent_date(sitecode, 'DP1.30003.001')
  if(is.nan(date)){
    print('No available data')
    return(0)
  }
  
  #Get file names and urls for all necessary tiles, idetnifying tile coordinates for which files are not available
  files_df <- build_file_df(coord_unique, date, sitecode)
  
  bad_coords <- files_df[files_df$name == 'bad tile',]$coords
  files_df <- files_df[files_df$name != 'bad_tile',]
  
  
  return(files_df)
  
}






for(i in 1:47){
  attempt <- try(site_description(site_df[i,'siteID'], NEON_all_plots))
  if(typeof(attempt) != 'try-error'){
    site_df[i,] <- attempt
  }
}











#Function to determin if a baseplot crosses tile boundaries, v2
is_on_boundary_2 <- function(x,y){
  the_names <- c('xMin', 'xMax', 'yMin', 'yMax')
  
  tile_diameter <- 1000
  plot_radius <- 20
  
  tile_x <- round1000(x)
  tile_y <- round1000(y)
  
  plot_extent <-   c((x - plot_radius),(x + plot_radius),(y - plot_radius),(y + plot_radius))
  names(plot_extent) <- the_names
  
  tile_extent <- c((tile_x),(tile_x + tile_diameter),(tile_y),(tile_y + tile_diameter))
  names(tile_extent) <- the_names
  
  
  out <- c(FALSE, FALSE, FALSE, FALSE)
  names(out) <- the_names
  
  #If plot left is less than tile left boundary, or plot lower boundary is less than tile lower boundary
  for(name in c('xMin','yMin')){
    if(plot_extent[name] <= tile_extent[name]){
      out[name] <- TRUE
    }
  }
  for(name in c('xMax', 'yMax')){
    if(plot_extent[name] >= tile_extent[name]){
      out[name] <- TRUE
    }
  }
 
  return(out) 
}





#Function takes spatialPolygons data frame for one plot that crosses a boundary, returns coordinates for 
#  each tile the plot crosses. Output is vector of coordinate strings
get_boundary_plot_coords_2 <- function(plot_spdf_row){
  options(scipen = 99999)
  x <- as.numeric(plot_spdf_row[19])
  y <- as.numeric(plot_spdf_row[20])
  
  east <- round1000(x)
  north <- round1000(y)
  
  
  boundary_vec <- is_on_boundary_2(x,y)
  
  #If tile is on left boundary, get easting coordinate of tile to left. If plot crosses right boundary, get
  # easting coordinate of tile to right
  if(boundary_vec['xMin']){
    east <- c(east, east - 1000)
  } else if(boundary_vec['xMax']){
    east <- c(east, east + 1000)
  }
  
  
  
  #If plot crosses lower tile boundary, get northing of tile below. Otherwise, if plot crosses upper tile boundary, get 
  # northing of tile above
  if(boundary_vec['yMin']){
    north <- c(north, north - 1000)
  } else if(boundary_vec['yMax']){
    north <- c(north, north + 1000)
  }
  
  
  
  #If tile crosses both boundaries, make vector of four empty strings. Otherwise, make vector of two empty strings
  out <- ifelse(boundary_vec[1]&boundary_vec[2],c('','','',''),c('',''))
  
  #Build vector of coordinate strings
  i <- 1
  for(n_coord in north){
    for(e_coord in east){
      out[i] <- paste0(as.character(e_coord),'_',as.character(n_coord))
      i = i + 1
    }
  }
  
  return(out)
}




