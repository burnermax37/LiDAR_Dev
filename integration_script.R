




#FUNCTIONS THAT GENERATE METRICS DATA FRAME



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

#Function takes row of boundary coordinates nested list, plot polygons spatial df, and file df
#Gets two point cloud files needed for plot and combines them in one object, calculates diversity metrics
calculate_boundary_metrics <- function(boundary_plot, plots_spdf, file_df){
  
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
  
  #calculate metrics
  metrics <- plot_diversity_metrics(plot_data, merged_LAS)
  return(metrics)
}



#Function takes plots sp_df, coords df, and file df. Calculates and compiles diversity metrics for each plot.
  #For every file, install file, calculate diveristy metrics of matching plots, remove file

build_metrics_df <- function(this_plots_spdf, this_coords_df, this_files_df, this_boundaries_df){
  #Initialize empty dataframe
  bound_N <- length(this_boundaries_df)
  
  DF <- data.frame(plotID = rep('',bound_N),
                   x = rep(NA,bound_N),
                   y = rep(NA,bound_N),
                   mean.max.canopy.ht = rep(NA,bound_N),
                   max.canopy.ht = rep(NA,bound_N), 
                   rumple = rep(NA,bound_N),
                   deepgaps = rep(NA,bound_N),
                   deepgap.fraction = rep(NA,bound_N),
                   cover.fraction = rep(NA,bound_N),
                   top.rugosity = rep(NA,bound_N),
                   vert.sd = rep(NA,bound_N), 
                   sd.sd = rep(NA,bound_N),
                   entro = rep(NA,bound_N),
                   GFP.AOP = rep(NA,bound_N),
                   VAI.AOP = rep(NA,bound_N),
                   VCI.AOP = rep(NA,bound_N)) 
  
  #Add rows for boundary plots
  if(bound_N > 0){
    for(i in seq(1, bound_N)){
      DF[i,1] <- as.character(this_boundaries_df[[i]]$plotID)
      DF[i,(2:16)] <- calculate_boundary_metrics(this_boundaries_df[[i]], this_plots_spdf, this_files_df)
    }
  }


  
  #Create metrics data for non-boundary plots by file, and merge with DF
  file_N <- length(this_files_df)
  
  for(i in seq(1,file_N)){
    new_file_DF <- assemble_tile_metrics(this_files_df[i,],this_coords_df,this_plots_spdf)
    DF <- rbind(DF, new_file_DF)
  }
  
  DF <- arrange(DF, plotID)
  
  return(DF)
}

