#' Find 3D aggregates
#'
#' Function to group adjacent pixels in aggregates
#' @param imgs The paths of array files; i.e. output from loadIMG or findIMG functions. A dataframe from the quantify function
#' @param channel Name of the channel to find aggregates in. Should be in the names of the array files
#' @param kern.neighbour Numeric vector indicating range of neighbouring pixels to aggregate in the x,y,z directions. Has to be odd intergers. c(1,1,1) means no aggregating.
#' @param kern.smooth Numeric vector indicating range of median smoothing in the x,y,z directions. Has to be odd intergers. c(1,1,1) means no smoothing.
#' @param layers Optional. Should the function only look in a subset of layers. A list with lists of layers to use for each image. Can also be the output from ELayers 
#' @param pwidth Optional. Width of pixels in microns to calculate aggregate size in microns instead of pixels
#' @param zstep Optional. z-step in microns to calculate aggregate size in microns instead of pixels
#' @keywords array image aggregate
#' @return A list with two parts. First part is a dataframe with ID and size of aggregates and name of image, second part is a list of the arrays in which pixels are NA if empty or given a number indicating the aggregate ID
#' @import mmand
#' @export

Agg <- function(imgs,channel,kern.neighbour=c(3,3,3),kern.smooth=c(3,3,3),layers=NULL,pwidth=NULL,zstep=NULL) {
  
  # Load image
  ch_files <- imgs[grep(channel, imgs)]
  
  # Lists for results
  arrays <- list()
  results <- list()
  
  # For each replica
  for(k in 1:length(ch_files)) {
    
    message(paste("Running replica",k))
    
    # Load
    ch_t <- readRDS(ch_files[k])
    side <- dim(ch_t)[1]
    h <- dim(ch_t)[3]
    
    # Subset layers
    if(!is.null(layers[[k]])){
      ch_t <- ch_t[,,layers[[k]]]
    }
    
    # Extend parameter
    ep <- max(c(kern.neighbour,kern.smooth))
    
    # Extend
    ch_new <- array(NA, dim=c(side+(2*ep), side+(2*ep), h+(2*ep)))
    ch_new[(ep+1):(side+ep),(ep+1):(side+ep),(ep+1):(h+ep)] <- ch_t 
    
    # Smooth
    if(kern.smooth[1] %% 1 == 0 & kern.smooth[1] %% 2 != 0 &
       kern.smooth[2] %% 1 == 0 & kern.smooth[2] %% 2 != 0 &
       kern.smooth[3] %% 1 == 0 & kern.smooth[3] %% 2 != 0) {
      kern.s <- kernelArray(array(1,dim=kern.smooth))
      ch_new <- medianFilter(ch_new,kern.s)
      ch_new[ch_new > 0] <- 1 } else stop("Kernel smooth has to be an odd integers in all directions")
    
    # Find aggregates
    kern.n <- kernelArray(array(1,dim=kern.neighbour))
    ch_agg <- components(ch_new,kern.n)
    
    # Crop
    ch_fin <- ch_agg[(ep+1):(side+ep),(ep+1):(side+ep),(ep+1):(h+ep)]
    
    afr <- as.data.frame(table(ch_fin))
    colnames(afr) <- c("ID","Size")
    if(!is.null(pwidth) & !is.null(zstep)) afr$Size.micron <- afr$Size * pwidth^2 * zstep
    afr$Img <- sub(paste0("_.*"),"",sub(".*/", "", ch_files[k]))
    results[[k]] <- afr  
    arrays[[k]] <- ch_fin
  }
  
  afrx <- do.call(rbind,results)
  
  All <- list(Aggregates=afrx,Arrays=arrays)
  
  return(All)
}



