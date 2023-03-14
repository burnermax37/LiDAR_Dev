

main <- function(sitecode){
  #Load TOS Plot Shape Files into R
  base_plots_SPDF <- generate_baseplot_spdf('All_NEON_TOS_Plots_V7/All_NEON_TOS_Plot_Polygons_V7.shp', sitecode)
  
  #Plot base plots
  plot(base_plots_SPDF, border = 'blue')
  
  #Get coordinates for all plots not on a tile boundary, then split by whether plot is on a boundary
  coord_df <- build_plot_frame(base_plots_SPDF)
  
  boundary_df <- coord_df[coord_df$coord_String == 'Plot crosses a tile boundary',]
  coord_df <- coord_df[coord_df$coord_String != 'Plot crosses a tile boundary',]
  
  #Get coordinates for plots on boundary tiles
  boundary_list <- build_plot_pseudo_frame(base_plots_SPDF, boundary_df)
  
  #Get unique coordinates, coordinates for each necessary tile
  coord_unique <- get_unique_coordinates(coord_df, boundary_list)
  
  #Get file names and urls for all necessary tiles
  files_df <- ldply(coord_unique, find_point_cloud, SITECODE = sitecode, DATE = '2019-06')
  
  
  #build final dataframe of stuctural diversity metrics
  metric_df <- build_metrics_df(base_plots_SPDF, coord_df, files_df, boundary_list)
  
  return(metric_df)
}


DF <- main('SOAP')
