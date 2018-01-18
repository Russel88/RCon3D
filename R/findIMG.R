#' Find arrays created with \code{loadIMG} or \code{tiff_to_array}
#'
#' Find arrays already created by the \code{loadIMG} or \code{tiff_to_array} functions
#' @param path The path of the _Array.R files
#' @param channels Optional. Channel names to look for in the image files
#' @keywords tif array image
#' @return The paths for the array files
#' @export

findIMG <- function(path, channels = NULL){
  
  # Find the files
  files <- list.files(path, "_Array.R", full.names = T)

  if(!is.null(channels)){
    files.new <- c()
    for(i in seq_along(channels)){
      files.new <- c(files.new,files[grep(channels[i],files)])
    }
    files <- files.new
  }
    
  return(files)
  
}
