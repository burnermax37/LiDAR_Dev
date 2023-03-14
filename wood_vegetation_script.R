library(sp)
library(raster)
library(neonUtilities)
library(downloader)
library(httr)
library(jsonlite)
library(dplyr)





#Function to calculate number of unique genus in a plot
#Function takes a mapping and tagging woody vegetation structure data table as a dataframe, and a plot ID
#Function number of unique taxonomic genus in a plot
#Can also take a shrubgroup dataframe and calculate number of unique shurub genera
count_unique_genus <- function(maptag, plot_id){
  
  #Filter mappingandtagging data to rows matching plot ID
  maptag <- maptag[maptag$plotID == plot_id,]
  
  #Get all unique sicentific names at plot, and prepare empty vector
  scientific <- unique(maptag$scientificName)
  genus <- rep('', length(scientific))
  
  #Split scientific names
  scientific <- strsplit(scientific, ' ')
  
  for(i in 1:length(genus)){
    genus[i] <- scientific[[i]][[1]]
  }
  
  #Filter out repeated and unknown genera
  genus <- unique(genus)
  genus <- genus[genus != 'Unknown']
  
  return(length(genus))
}





#Function takes dataframe with plotID attribute and mappingandtagging datagframe, and returns former dataframe with added
#  column listing number of unique plant genera recorded in each plot
add_genus_column <- function(struct_df, mapandtag_df){
  
  #Add column for genus count
  struct_df$genus_count <- NA
  
  
  for(i in 1:nrow(struct_df)){
    plot <- struct_df[i,'plotID']
    
    #If the current plot appears in the mappingandtagging dataframe
    if((plot %in% mapandtag_df$plotID)){
      #Get unique genus counts
      struct_df[i,'genus_count'] <- count_unique_genus(mapandtag_df, struct_df[i,'plotID'])
    }
  }
  return(struct_df)
}





###Function takes veg structure apparentIndividual dataframe and plotID, returns named list of parameters calculated
# based on trees in the plot
get_tree_parameters <- function(in_plants, plot_id){
  
  #Filter apparentIndividual dataframe to rows matching plot ID
  in_plants <- in_plants[in_plants$plotID == plot_id,]
  
  #Remove dead individuals
  in_plants <- in_plants[in_plants$plantStatus == 'Live',]
  
  #Extract trees
  in_trees <- in_plants[which(in_plants$growthForm %in% c('small tree','sapling','single-bole tree', 'multi-bole tree')),]
  
  #Assemble metrics
  woody_plant_count <- nrow(in_plants)
  tree_count <- nrow(in_trees)
  tree_wood_plant_ratio <- round(tree_count/woody_plant_count, 2)
  mean_tree_height <- round(mean(in_trees$height, na.rm = T), 2)
  max_measured_tree_height <- max(in_trees$height, na.rm = T)

  


  
  out <- list(woody_plant_count, tree_count, tree_wood_plant_ratio, mean_tree_height, max_measured_tree_height)
  names(out) <- c('woody_plants.count', 'trees.count', 'tree_to_woody_plant.ratio', 'tree.mean_height',
                  'tree.max_measured_height')

  return(out)
}

##Function takes a dataframe with plotID attribute and a vegetation structure apparentIndividual dataframe
###Returns first dataframe but with added columns for tree metrics
add_tree_parameters <- function(in_plots, in_apparent){
  
  #Add empty columns
  in_plots$wood_plants.count <- NA
  in_plots$trees.count <- NA
  in_plots$tree_to_woody_plant.ratio <- NA
  in_plots$tree.mean_height <- NA
  in_plots$tree.max_measured_height <- NA
  
  #Create vector of column names
  columns <- c('wood_plants.count','trees.count','tree_to_woody_plant.ratio', 'tree.mean_height',
               'tree.max_measured_height')
  
  
  #For each row in the input dataframe
  for(i in 1:nrow(in_plots)){
    #If there are entries in the apparent Individual dataframe occuring in that plot
    if(in_plots[i,'plotID'] %in% unique(in_apparent$plotID)){
      in_plots[i,columns] <- get_tree_parameters(in_apparent, in_plots[i,'plotID'])
    }
  }
  
  return(in_plots)
}


###Function takes dataframe with plotID attribute, and site code, gets woody vegetation structure data 
##  attaches species counts by plot and tree parameters to input dataframe and returns result
add_woody_tree_data <- function(in_data, site){
  
  SERVER <- 'http://data.neonscience.org/api/v0/'
  PRODUCTCODE <- 'DP1.10098.001'
  
  

  DATE <- get_most_recent_date(SITECODE = site, PRODUCTCODE = PRODUCTCODE)
  if(!is.nan(DATE)){
    
    list2env(loadByProduct(dpID=PRODUCTCODE, site=site, package="basic", check.size = F), .GlobalEnv)
    
    #Add species counts
    in_data <- add_genus_column(in_data, vst_mappingandtagging)
    
    #Add tree parameters
    in_data <- add_tree_parameters(in_data, vst_apparentindividual)
  }


  
  return(in_data)
}






























###############
library(sp)
library(raster)
library(neonUtilities)
library(geoNEON)






#Get vegetation structure data
veglist <- loadByProduct(dpID="DP1.10098.001", site="STER", package="basic")

#Get geolocation data for vegetation structure
vegmap <- getLocTOS(veglist$vst_mappingandtagging, 'vst_mappingandtagging')

#Merge data on tree measurements with tree location data
veg <- merge(veglist$vst_apparentindividual, vegmap, by = c('individualID', 'namedLocation', 'domainID', 'siteID', 'plotID'))

#Plot data; circle size indicates stem diameter
plot_id <- 'SOAP_059'
#plot_id <- 'WREF_075'

symbols(veg$adjEasting[which(veg$plotID==plot_id)], 
        veg$adjNorthing[which(veg$plotID==plot_id)], 
        circles=veg$stemDiameter[which(veg$plotID==plot_id)]/100, 
        inches=F, xlab="Easting", ylab="Northing")

#Add estimated uncertainty to plot
symbols(veg$adjEasting[which(veg$plotID==plot_id)], 
        veg$adjNorthing[which(veg$plotID==plot_id)], 
        circles=veg$stemDiameter[which(veg$plotID==plot_id)]/100, 
        inches=F, xlab="Easting", ylab="Northing")
symbols(veg$adjEasting[which(veg$plotID==plot_id)], 
        veg$adjNorthing[which(veg$plotID==plot_id)], 
        circles=veg$adjCoordinateUncertainty[which(veg$plotID==plot_id)], 
        inches=F, add=T, fg="lightblue")

#Download CHM data for plot WREF_075
byTileAOP(dpID="DP3.30015.001", site="SOAP", year="2017", 
          easting=veg$adjEasting[which(veg$plotID==plot_id)], 
          northing=veg$adjNorthing[which(veg$plotID==plot_id)],
          savepath="C:/R/data/lidar_data")

chm <- raster("C:/R/data/lidar_data/DP3.30015.001/2017/FullSite/D17/2017_SOAP_2/L3/DiscreteLidar/CanopyHeightModelGtif/NEON_D17_SOAP_DP3_297000_4101000_CHM.tif")


plot(chm,col = topo.colors(5))




#Subset data to remove values outside of extent (speeds processing)
vegsub <- veg[which(veg$adjEasting >= extent(chm)[1] &
                      veg$adjEasting <= extent(chm)[2] &
                      veg$adjNorthing >= extent(chm)[3] & 
                      veg$adjNorthing <= extent(chm)[4]),]


#Extract CHM data for each tree location, with buffer size matching coordinate uncertainty,
#  and calculate max CHM value in each of the extracted areas
bufferCHM <- extract(chm, cbind(vegsub$adjEasting, 
                                vegsub$adjNorthing),
                     buffer=vegsub$adjCoordinateUncertainty, 
                     fun=max)

#Plot extracted data: tree height vs. max CHM value at tree location
plot(bufferCHM~vegsub$height, pch=20, xlab="Height", 
     ylab="Canopy height model")
lines(c(0,50), c(0,50), col="grey")

#Comparison Method One: Use Map-Centric approach to filter out values from understory vegetation
#  Aggregate vegetation structure data and CHM data to find tallest point in a grid. Grid size 10m cells

#  Round data down to nearest multiple of ten
vegsub$easting10 <- 10*floor(vegsub$adjEasting/10)
vegsub$northing10 <- 10*floor(vegsub$adjNorthing/10)

#  Aggregate vegetation structure data into 10m bins
vegbin <- stats::aggregate(vegsub, by = list(vegsub$easting10, vegsub$northing10), FUN = max)

#  Aggregate CHM data into 10m bins
CHM10 <- raster::aggregate(chm, fact = 10, fun = max)
plot(CHM10, col = topo.colors(5))

#  Adjust coordinates to match coordinates of pixel corners
vegbin$easting10 <- vegbin$easting10 + 5
vegbin$northing10 <- vegbin$northing10 + 5

#  Extract values from raster, extractig each pixel
binCHM <- extract(CHM10, cbind(vegbin$easting10, vegbin$northing10))

plot(binCHM~vegbin$height, pch = 20,
     xlab = 'Height', ylab = 'Canopy height model')
lines(c(0,50), c(0,50), col = 'grey')

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


#Extract raster values based on filtered trees
filterCHM <- extract(chm, cbind(vegfil$adjEasting, vegfil$adjNorthing),
                     buffer=vegfil$adjCoordinateUncertainty+1, fun=max)



plot(filterCHM~vegfil$height, pch=20, 
     xlab="Height", ylab="Canopy height model")
lines(c(0,50), c(0,50), col="grey")





