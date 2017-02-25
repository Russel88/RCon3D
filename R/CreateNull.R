#' Create RDS for all empty pixels
#'
#' Create arrays which represents all the empty pixels in an image
#' @param path The path of the _Array.R files
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the .tif files
#' @keywords array image
#' @return Creates arrays as RDS files in the specified path, and outputs the paths for these files
#' @export

CreateNulls <- function(path, channels){

  setwd(path)
  
  # Find the files
  files <- list.files(path, "_Array.R", full.names = T)
  
  # Remove channel names
  for(i in 1:length(channels)){
      files <- lapply(files, function(x) gsub(channels[i],"RCon3D.Null",x))
  }
  files.nc <- as.character(unique(files))
 
  for(i in unique(files.nc)) {
    
    files.sub <- list()
    arrays <- list()
    
    # Load the images
    for(j in 1:length(channels)){
      files.sub[[j]] <- gsub("RCon3D.Null",channels[j],i)
    }
    
    for(j in 1:length(channels)){
      arrays[[j]] <- readRDS(files.sub[[j]])
    }
    
    # Aggregate channels
    array.sum <- do.call("+",arrays)
    
    # Create null channel
    array.sum[array.sum == 0] <- NA
    array.sum[array.sum > 0] <- 0
    array.sum[is.na(array.sum)] <- 1

    saveRDS(array.sum, file = i)
    
  }
  
  # Find the files
  files <- list.files(path, "RCon3D.Null", full.names = T)
  
  return(files)
  
}