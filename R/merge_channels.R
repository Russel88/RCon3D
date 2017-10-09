#' Create RDS for a combination of channels
#'
#' @param path The path of the _Array.R files
#' @param channels Character vector with names of the channels. Channel names should be in the names of the .tif files
#' @param method How to merge. Either "union" (presence is coded of at least one channel is present), "intersect" (presence is coded only if all channels are present), or "subtract" (presence is coded as those in the first channel but not in the second channel).  
#' @keywords array image
#' @return Creates arrays as RDS files in the specified path, and outputs the paths for these files
#' @export

merge_channels <- function(path, channels, method){

  setwd(path)
  
  # Find the files
  files <- list.files(path, "_Array.R", full.names = T)
  files.l <- lapply(channels, function(x) files[grep(x, files)])
  
  # For each image
  for(k in 1:length(files.l[[1]])){
    
    imgs <- lapply(files.l, function(x) readRDS(x[[k]]))
    
    if(method == "union"){
      img <- do.call("+",imgs)
      img[img != 0] <- 1
      saveRDS(img, file = gsub(channels[1],paste0("RCon3D.union.",channels[1],".",channels[2]),files.l[[1]][k]))
      files <- list.files(path, "RCon3D.union.", full.names = T)
    }
    
    if(method == "intersect"){
      img <- do.call("+",imgs)
      img[img != length(channels)] <- 0
      img[img == length(channels)] <- 1
      saveRDS(img, file = gsub(channels[1],paste0("RCon3D.intersect.",channels[1],".",channels[2]),files.l[[1]][k]))
      files <- list.files(path, "RCon3D.intersect.", full.names = T)
    }
    
    if(method == "subtract"){
      img <- imgs[[1]] - imgs[[2]]
      img[img == -1] <- 0
      saveRDS(img, file = gsub(channels[1],paste0("RCon3D.subtract.",channels[1],".",channels[2]),files.l[[1]][k]))
      files <- list.files(path, "RCon3D.subtract.", full.names = T)
    }

  }
  
  return(files)

}