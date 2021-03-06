#' Create RDS for all empty pixels
#'
#' Create arrays which represents all the empty pixels in an image
#' @param files The paths of the _Array.R files
#' @param path The path of where to save the created RDS files
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the .tif files
#' @keywords array image
#' @return Creates arrays as RDS files in the specified path, and outputs the paths for these files
#' @export

create_nulls <- function(files, path, channels){

  workdir <- getwd()
  setwd(path)
  
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
  
  setwd(workdir)
  
  return(files)
  
}