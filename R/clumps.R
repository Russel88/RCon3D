#' Find 3D aggregates/clumps
#'
#' Function to group adjacent pixels in aggregates/clumps
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param channel Name of the channel to find aggregates in. Should be in the names of the array files
#' @param kern.neighbour Numeric vector indicating range of neighbouring pixels to aggregate in the x,y,z directions. Has to be odd intergers. c(1,1,1) means no aggregating.
#' @param kern.smooth Optional. Numeric vector indicating range of median smoothing in the x,y,z directions. Has to be odd intergers. c(1,1,1) means no smoothing.
#' @param layers Optional. Should the function only look in a subset of layers. A list with lists of layers to use for each image. Can also be the output from \code{extract_layers} 
#' @param pwidth Optional. Width of pixels in microns to calculate aggregate size in microns instead of pixels
#' @param zstep Optional. z-step in microns to calculate aggregate size in microns instead of pixels
#' @param naming Optional. Add metadata to the output dataframe by looking through names of array files. Should be a list of character vectors, each list element will be added as a variable. Example: naming=list(Time=c("T0","T1","T2")). The function inserts a variable called Time, and then looks through the names of the array files and inserts characters mathcing either T0, T1 or T2
#' @keywords array image aggregate
#' @return A list with two parts. First part is a dataframe with ID and size of aggregates and name of image, second part is a list of the arrays in which pixels are NA if empty or given a number indicating the aggregate ID
#' @import mmand
#' @export

clumps <- function(imgs,channel,kern.neighbour=c(3,3,3),kern.smooth=NULL,layers=NULL,pwidth=NULL,zstep=NULL,naming=NULL) {
  
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
    if(!is.null(kern.smooth)) {
      if(kern.smooth[1] %% 1 == 0 & kern.smooth[1] %% 2 != 0 &
         kern.smooth[2] %% 1 == 0 & kern.smooth[2] %% 2 != 0 &
         kern.smooth[3] %% 1 == 0 & kern.smooth[3] %% 2 != 0) {
        kern.s <- kernelArray(array(1,dim=kern.smooth))
        ch_new <- medianFilter(ch_new,kern.s)
        ch_new[ch_new > 0] <- 1 } else stop("Kernel smooth has to be odd integers in all directions")
    }
    
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
  
  if(!is.null(naming)){
    afrxn <- afrx
    for(i in 1:length(naming)){
      name.temp <- naming[[i]]
      afrxn <- cbind(afrxn,NA)
      for(j in 1:length(name.temp)){
        afrxn[grep(name.temp[j],afrxn$Img),ncol(afrx)+i] <- name.temp[j]
      }
    }
    
    colnames(afrxn) <- c(colnames(afrx),names(naming))
    
  } else afrxn <- afrx
  
  All <- list(Aggregates=afrxn,Arrays=arrays)
  
  return(All)
}


