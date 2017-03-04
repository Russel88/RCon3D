#' Find arrays created with \code{loadIMG} or \code{tiff_to_array}
#'
#' Find arrays already created by the \code{loadIMG} or \code{tiff_to_array} functions
#' @param path The path of the _Array.R files
#' @keywords tif array image
#' @return The paths for the array files
#' @export

findIMG <- function(path){
  
  # Find the files
  files <- list.files(path, "_Array.R", full.names = T)
  
  return(files)
  
}
