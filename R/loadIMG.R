#' Load .tif images into R
#'
#' Create an array from one or more .tif files.
#' @param path The path of either a folder with all .tif files to be loaded (if split=FALSE), or a folder where each image is in a separate subfolder with all z-stacks and channels in that folder (if split=TRUE)
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the .tif files
#' @param split Logical indicating if all z-stacks are in one .tif (FALSE) or in a folder with a .tif for each z-stack (TRUE)
#' @param multi Integer. If the .tif contains multiple images this varible should indicate how to divide them up. Three images, each composed of 50, 100, and 150 layers: c(50,100,150). 10 images with each 50 layers. rep(50, 10). If NULL input is assumed to be a single 3D image
#' @param multi.name If multi is not NULL, a vector with names for images.  
#' @keywords tif array image
#' @details If the .tif file has a color.space attribute saying "white is zero" (and split=FALSE) the binary coding of the image is reversed 
#' @return Creates arrays as RDS files in the set path, and outputs the paths for these files
#' @export

loadIMG <- function(path, channels, split = FALSE, multi = NULL, multi.name = NULL){

  if(!is.null(multi)){
    stopifnot(!is.null(multi.name))
    stopifnot(length(multi) == length(multi.name))
  }
  
  # Save working dir
  workdir <- getwd()
  
  # Set working dir
  setwd(path)

  # Convert images to arrays
  tiff_to_array(channels,split=split,multi=multi,multi.name=multi.name)

  # Find the files
  files <- list.files(getwd(), "_Array.R", full.names = T)

  # Return to original working dir
  setwd(workdir)
  
  return(files)

}
