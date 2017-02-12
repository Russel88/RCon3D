#' Find arrays created with loadIMG or tiffToArray
#'
#' Find arrays already created by the loadIMG or tiffToArray functions
#' @param path The path of the _Array.R files
#' @keywords tif array image
#' @return The paths for the array files
#' @export

findIMG <- function(path){
  
  # Find the files
  files <- list.files(path, "_Array.R", full.names = T)
  
  return(files)
  
}
