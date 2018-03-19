#' Find Center of Mass in images
#'
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the array files
#' @keywords array image quantify
#' @return A list with the centers of mass in x, y, z format for each image
#' @export

center_of_mass <- function(imgs,channels) {

  # Loop for each channel
  centers <- list()

  ch_files <- lapply(channels, function(x) imgs[grep(x, imgs)])
  
  for(k in 1:length(ch_files[[1]])){
  
    img_sub <- lapply(1:length(ch_files), function(x) ch_files[[x]][[k]])
    imags <- lapply(img_sub, function(x) readRDS(x))
    img_fin <- do.call("+",imags)
    img_fin[img_fin != 0] <- 1
    
    where <- which(img_fin == 1, arr.ind = TRUE)
    pos <- round(apply(where, 2, median))
    names(pos) <- c("x","y","z")

    centers[[k]] <- pos  
  }
  return(centers)
}

