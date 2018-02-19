#' Make Array from .tif file(s) 
#'
#' Create an array from one or more .tif files. Looks for .tif in the working directory
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the .tif files
#' @param split Logical indicating if all z-stacks are in one .tif (FALSE) or in a folder with a .tif for each z-stack (TRUE)
#' @param multi Integer. If the .tif contains multiple images this varible should indicate how to divide them up. Three images, each composed of 50, 100, and 150 layers: c(50,100,150). 10 images with each 50 layers. rep(50, 10). If NULL, input is assumed to be a single 3D image
#' @param multi.name If multi is not NULL, a vector with names for images.  
#' @keywords tif array image
#' @keywords tif array
#' @details If the .tif file has a color.space attribute saying "white is zero" (and split=FALSE) the binary coding of the image is reversed 
#' @return Arrays as RDS files in working directory
#' @importFrom tiff readTIFF
#' @export

tiff_to_array <- function(channels,split,multi,multi.name) {
  
  if(split){
    
    # Find subdirs
    subdirs <- list.dirs(getwd(), T, F)
    
    for(f in 1:length(subdirs)){
    
    message(paste("Loading image",f))
    
    # list the file paths
    files <- dir(subdirs[f], ".tif", full.names = T)
    
    cFiles <- list()
    
    for(j in 1:length(channels)){
      
      cFiles[[j]] <- files[grep(channels[j], files)]
      
      # load
      side <- dim(readTIFF(cFiles[[j]]))[1]
      cArray <- array(0, c(side, side, length(cFiles[[j]])))
      for(i in 1:length(cFiles[[j]])){
        cArray[,,i] <- readTIFF(cFiles[[j]][i])
      }
      saveRDS(cArray, file = paste0(subdirs[f], "_", channels[j],"_Array.R"))
    }
  }
    
  } else {
    
    files <- dir(getwd(), ".tif", full.names = T)
    
    for(f in 1:length(files)){
      
      message(paste("Loading image",f))
      
      if(!is.null(multi)){
        
        # Load
        tiflist <- readTIFF(files[f],all=TRUE)
        check <- readTIFF(files[f],info=TRUE)
        for(k in 1:length(multi)){
           sub.list <- tiflist[(sum(multi[1:k])-multi[k]+1):sum(multi[1:k])]
           side <- dim(tiflist[[1]])[1]
           cArray <- array(0, c(side, side, length(tiflist)))
           for(l in 1:length(sub.list)){
             cArray[,,l] <- sub.list[[l]]
           }
           # Fix reverse images
           if(attributes(check)$color.space == "white is zero") {
             cArray[cArray == 0] <- NA
             cArray[cArray > 0] <- 0
             cArray[is.na(cArray)] <- 1
           }
           saveRDS(cArray, file = gsub(".tif",paste0("_",multi.name[k],"_Array.R"),files[f]))
        }
        
      } else {
        # Load
        tiflist <- readTIFF(files[f],all=TRUE)
        side <- dim(tiflist[[1]])[1]
        cArray <- array(0, c(side, side, length(tiflist)))
        for(l in 1:length(tiflist)){
          cArray[,,l] <- tiflist[[l]]
        }
        
        # Fix reverse images
        check <- readTIFF(files[f],info=TRUE)
        if(attributes(check)$color.space == "white is zero") {
          cArray[cArray == 0] <- NA
          cArray[cArray > 0] <- 0
          cArray[is.na(cArray)] <- 1
        }
        saveRDS(cArray, file = gsub(".tif","_Array.R",files[f]))
      }
    }
  }
}
