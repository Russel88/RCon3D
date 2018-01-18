#' Create RDS for a combination of channels
#'
#' @param path The path of where to create the _Array.R files
#' @param imgs The paths for the images to merge
#' @param channels Character vector with names of the channels. Channel names should be in the names of the .tif files
#' @param method How to merge. Either "union" (presence is coded if at least one channel is present), "intersect" (presence is coded only if all channels are present), or "subtract" (presence is coded as those in the first channel but not in the second channel).  
#' @keywords array image
#' @return Creates arrays as RDS files in the specified path, and outputs the paths for image files containing the name of the method and channels
#' @export

merge_channels <- function(path, imgs, channels, method){
  
  setwd(path)
  
  # Find the files
  files <- imgs
  files.l <- lapply(channels, function(x) files[grep(x, files)])

  # For each image
  for(k in 1:length(files.l[[1]])){
    
    imags <- lapply(files.l, function(x) readRDS(x[[k]]))
    
    if(method == "union"){
      img <- do.call("+",imags)
      img[img != 0] <- 1
      saveRDS(img, file = gsub(channels[1],paste0("RCon3D.",channels[1],".",channels[2],".union"),files.l[[1]][k]))
      files <- list.files(path, paste0("RCon3D.",channels[1],".",channels[2],".union"), full.names = T)
    }
    
    if(method == "intersect"){
      img <- do.call("+",imags)
      img[img != length(channels)] <- 0
      img[img == length(channels)] <- 1
      saveRDS(img, file = gsub(channels[1],paste0("RCon3D.",channels[1],".",channels[2],".intersect"),files.l[[1]][k]))
      files <- list.files(path, paste0("RCon3D.",channels[1],".",channels[2],".intersect"), full.names = T)
    }
    
    if(method == "subtract"){
      img <- imags[[1]] - imags[[2]]
      img[img == -1] <- 0
      saveRDS(img, file = gsub(channels[1],paste0("RCon3D.",channels[1],".",channels[2],".subtract"),files.l[[1]][k]))
      files <- list.files(path, paste0("RCon3D.",channels[1],".",channels[2],".subtract"), full.names = T)
    }
    
  }
  
  return(files)
  
}