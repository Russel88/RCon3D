#' Find 3D aggregates/clumps
#'
#' Function to group adjacent pixels in aggregates/clumps
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param channel Name of the channel to find aggregates in. Should be in the names of the array files
#' @param kern.neighbour Numeric vector indicating range of neighbouring pixels to aggregate in the x,y,z directions. Has to be odd intergers. c(1,1,1) means no aggregating.
#' @param type.neighbour Type of kernel for neighbourhood. "box" includes diagonals, "diamond" is without diagonals
#' @param kern.smooth Optional. Numeric vector indicating range of median smoothing in the x,y,z directions. Has to be odd intergers. c(1,1,1) means no smoothing.
#' @param type.smooth Optional. Type of kernel for smooth "box" includes diagonals, "diamond" is without diagonals
#' @param layers Optional. Should the function only look in a subset of layers. A list with lists of layers to use for each image. Can also be the output from \code{extract_layers} 
#' @param pwidth Optional. Width of pixels in microns to calculate aggregate size in microns instead of pixels
#' @param zstep Optional. z-step in microns to calculate aggregate size in microns instead of pixels
#' @param naming Optional. Add metadata to the output dataframe by looking through names of array files. Should be a list of character vectors, each list element will be added as a variable. Example: naming=list(Time=c("T0","T1","T2")). The function inserts a variable called Time, and then looks through the names of the array files and inserts characters mathcing either T0, T1 or T2
#' @param coords Logical. Return coordinates of the centroids of each aggregate. This can be somewhat time-consuming if there are many aggregates
#' @keywords array image aggregate
#' @return A list with two parts. First part is a dataframe with ID, size of aggregates in pixels, size of aggregates in microns if pwidth and zstep are provided, coordinates if coords is TRUE, if coords is TRUE also a logical variable, Edge, indicating whether the aggregate touches the edge, and name of image. Second part is a list of the arrays in which pixels are NA if empty or given a number indicating the aggregate ID
#' @import mmand
#' @export

clumps <- function(imgs,channel,kern.neighbour=c(3,3,3),type.neighbour="box",kern.smooth=NULL,type.smooth="box",layers=NULL,pwidth=NULL,zstep=NULL,naming=NULL,coords=FALSE) {
  
  stopifnot(length(channel)==1)
  
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
    
    # Smooth
    if(!is.null(kern.smooth)) {
      if(kern.smooth[1] %% 1 == 0 & kern.smooth[1] %% 2 != 0 &
         kern.smooth[2] %% 1 == 0 & kern.smooth[2] %% 2 != 0 &
         kern.smooth[3] %% 1 == 0 & kern.smooth[3] %% 2 != 0) {
        kern.s <- shapeKernel(kern.smooth, type = type.smooth)
        ch_t <- medianFilter(ch_t,kern.s)
        ch_t[ch_t > 0] <- 1 } else stop("Kernel smooth has to be odd integers in all directions")
    }
    
    # Find aggregates
    kern.n <- shapeKernel(kern.neighbour, type = type.neighbour)
    ch_agg <- components(ch_t,kern.n)
    
    afr <- as.data.frame(table(ch_agg))
    colnames(afr) <- c("ID","Size")
    if(!is.null(pwidth) & !is.null(zstep)) afr$Size.micron <- afr$Size * pwidth^2 * zstep
    afr$Img <- sub(paste0("_.*"),"",sub(".*/", "", ch_files[k]))
    
    if(coords){
      coordsl <- lapply(as.numeric(as.character(afr$ID)),function(x) which(ch_agg == x, arr.ind = TRUE))
      centroids <- t(sapply(coordsl, function(x) apply(x, 2, median)))
      colnames(centroids) <- c("x","y","z")
      afr <- as.data.frame(cbind(afr, centroids))
      
      # Touch edge?
      edge.x <- c(1,dim(ch_agg)[1])
      edge.y <- c(1,dim(ch_agg)[2])
      edge.z <- c(1,dim(ch_agg)[3])
      
      touch <- sapply(coordsl, function(ac) any(any(sapply(ac[,1], function(x) x %in% edge.x)),
                                               any(sapply(ac[,2], function(x) x %in% edge.y)),
                                               any(sapply(ac[,3], function(x) x %in% edge.z))))
      
      afr$Edge <- touch
    }
    
    results[[k]] <- afr
    arrays[[k]] <- ch_agg
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





