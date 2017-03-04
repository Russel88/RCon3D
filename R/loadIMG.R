#' Load .tif images into R
#'
#' Create an array from one or more .tif files. Looks for .tif in the working directory
#' @param path The path of either a folder with all .tif files to be loaded (if split=FALSE), or a folder where each image is in a separate subfolder with all z-stacks and channels in that folder (if split=TRUE)
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the .tif files
#' @param split Logical indicating if all z-stacks are in one .tif (FALSE) or in a folder with a .tif for each z-stack (TRUE)
#' @keywords tif array image
#' @return Creates arrays as RDS files in working directory, and outputs the paths for these files
#' @export

loadIMG <- function(path, channels, split=FALSE){

  # Set working dir
  setwd(path)

  # Convert images to arrays
  tiff_to_array(channels,split=split)

  # Find the files
  files <- list.files(getwd(), "_Array.R", full.names = T)

  return(files)

}
