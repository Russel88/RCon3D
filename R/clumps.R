#' Find 3D aggregates/clumps
#'
#' Function to group adjacent pixels in aggregates/clumps
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param channels Name of the channel to find aggregates in. Should be in the names of the array files
#' @param kern.neighbour Numeric vector indicating range of neighbouring pixels to aggregate in the x,y,z directions. Has to be odd intergers. c(1,1,1) means no aggregating.
#' @param type.neighbour Type of kernel for neighbourhood. "box" includes diagonals, "diamond" is without diagonals
#' @param kern.smooth Optional. Numeric vector indicating range of median smoothing in the x,y,z directions. Has to be odd integers. c(1,1,1) means no smoothing.
#' @param type.smooth Optional. Type of kernel for smooth "box" includes diagonals, "diamond" is without diagonals
#' @param layers Optional. Should the function only look in a subset of layers. A list with lists of layers to use for each image. Can also be the output from \code{extract_layers} 
#' @param pwidth Optional. Width of pixels in microns to calculate aggregate size in microns instead of pixels
#' @param zstep Optional. z-step in microns to calculate aggregate size in microns instead of pixels
#' @param naming Optional. Add metadata to the output dataframe by looking through names of array files. Should be a list of character vectors, each list element will be added as a variable. Example: naming=list(Time=c("T0","T1","T2")). The function inserts a variable called Time, and then looks through the names of the array files and inserts characters mathcing either T0, T1 or T2
#' @param coords Optional. Logical. Return coordinates of the centroids of each aggregate. This can be somewhat time-consuming if there are many aggregates
#' @param thresh Optional. Numeric. Clumps containing fewer than this number of pixels are discarded. Default 0
#' @param points Optional. List. 3D coordinates of points from which the distance to each aggregate should be calculated. The list should have an element for each image in imgs in the same order. Each element is a vector of length 3 with the coordinates. E.g. points = list(c(250,200,10),c(500,300,20)) 
#' @param cores Integer. Number of cores to use for parallel computing.
#' @keywords array image aggregate
#' @return A list with two parts. First part is a dataframe with ID, size of aggregates in pixels, size of aggregates in microns if pwidth and zstep are provided, coordinates if coords is TRUE, if coords is TRUE also variables, Edge.x, Edge.y and Edge.z, indicating how many pixels of each aggregate is on the edge of the image, and name of image. Second part is a list of the arrays in which pixels are NA if empty or given a number indicating the aggregate ID
#' @import mmand foreach doSNOW
#' @export

clumps <- function(imgs,channels,kern.neighbour=c(3,3,3),type.neighbour="box",kern.smooth=NULL,type.smooth="box",layers=NULL,pwidth=NULL,zstep=NULL,naming=NULL,coords=FALSE,thresh=0,points=NULL,cores=1) {
  
  # Find images
  if(length(channels) == 1){
    ch_files.x <-  imgs[grep(channels, imgs)]
  } else {
    ch_files <- lapply(channels, function(x) imgs[grep(x, imgs)])
    ch_files.x <- lapply(1:length(ch_files[[1]]), function(x) sapply(ch_files, function(y) y[[x]]))
  }
  
  # For each image
  if(cores > 1){
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }
  
  pb <- txtProgressBar(max = length(ch_files.x), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  res.list <- foreach(k = 1:length(ch_files.x), .options.snow = opts, .packages = "mmand") %dopar% {
    
    # Load
    if(length(channels) > 1){
      ch_t <- lapply(unlist(ch_files.x[k]), function(x) readRDS(x))
      ch_t <- do.call("+", ch_t)
      ch_t[ch_t != 0] <- 1
    } else {
      ch_t <- readRDS(ch_files.x[k])
    }

    # Subset layers
    if(!is.null(layers[[k]])){
      ch_t <- ch_t[,,layers[[k]]]
    }
    
    # Smooth
    if(!is.null(kern.smooth)) {
      if(kern.smooth[1] %% 2 == 1 &
         kern.smooth[2] %% 2 == 1 &
         kern.smooth[3] %% 2 == 1) {
        kern.s <- shapeKernel(kern.smooth, type = type.smooth)
        ch_t <- medianFilter(ch_t,kern.s)
        ch_t[ch_t > 0] <- 1 } else stop("Kernel smooth has to be odd integers in all directions")
    }
    
    # Find aggregates
    kern.n <- shapeKernel(kern.neighbour, type = type.neighbour)
    ch_agg <- components(ch_t,kern.n)

    # Extract aggregates
    afr <- as.data.frame(table(ch_agg))
    colnames(afr) <- c("ID","Size")
    if(!is.null(pwidth) & !is.null(zstep)) afr$Size.micron <- afr$Size * pwidth^2 * zstep
    afr$Img <- gsub(channels[1],"",sub(paste0("_Array.*"),"",sub(".*/", "", ch_files.x[[k]][1])))
    afr <- afr[afr$Size >= thresh,]
    
    # Extract if multiple channels
    if(length(channels) > 1){
      afr.list <- list()
      for(zz in 1:length(channels)){
        imag <- readRDS(ch_files.x[[k]][zz])
        imag <- ch_agg * imag
        afr.sub <- as.data.frame(table(imag))
        colnames(afr.sub) <- c("ID","Size")
        afr.sub <- afr.sub[afr.sub$ID != 0,] # All those removed by smoothing
        afr.sub$ID <- as.character(afr.sub$ID)
        if(!is.null(pwidth) & !is.null(zstep)) afr.sub$Size.micron <- afr.sub[,2] * pwidth^2 * zstep
        afr.sub$Img <- gsub(channels[1],"",sub(paste0("_Array.*"),"",sub(".*/", "", ch_files.x[[k]][1])))
        afr.sub <- afr.sub[afr.sub$ID %in% afr$ID,]
        afr.list[[zz]] <- afr.sub
      }
      names(afr.list) <- channels
    }
    
    if(coords){
      coordsl <- lapply(as.numeric(as.character(afr$ID)),function(x) which(ch_agg == x, arr.ind = TRUE))
      centroids <- t(sapply(coordsl, function(x) apply(x, 2, median)))
      colnames(centroids) <- c("x","y","z")
      afr <- as.data.frame(cbind(afr, centroids))
      
      # Touch edge?
      edge.x <- c(1,dim(ch_agg)[1])
      edge.y <- c(1,dim(ch_agg)[2])
      edge.z <- c(1,dim(ch_agg)[3])
      
      touch <- t(sapply(coordsl, function(ac) c(sum(sapply(ac[,1], function(x) x %in% edge.x)),
                                               sum(sapply(ac[,2], function(x) x %in% edge.y)),
                                               sum(sapply(ac[,3], function(x) x %in% edge.z)))))
      
      colnames(touch) <- c("Edge.x","Edge.y","Edge.z")
      
      afr <- cbind(afr, touch)
      
      # Distance
      if(!is.null(points) & !is.null(pwidth) & !is.null(zstep)){
        afr$dist <- sqrt(((afr$x-points[[k]][1])*pwidth)^2 + 
                           ((afr$y-points[[k]][2])*pwidth)^2 + 
                           ((afr$z-points[[k]][3])*zstep)^2)
      }

      if(length(channels) > 1){
        for(vv in 1:length(channels)){
          coordsl.sub <- lapply(as.numeric(as.character(afr.list[[vv]]$ID)),function(x) which(ch_agg == x, arr.ind = TRUE))
          centroids.sub <- t(sapply(coordsl.sub, function(x) apply(x, 2, median)))
          colnames(centroids.sub) <- c("x","y","z")
          afr.list[[vv]] <- as.data.frame(cbind(afr.list[[vv]], centroids.sub))
          
          # Touch edge?
          edge.x <- c(1,dim(ch_agg)[1])
          edge.y <- c(1,dim(ch_agg)[2])
          edge.z <- c(1,dim(ch_agg)[3])
          
          touch.sub <- t(sapply(coordsl.sub, function(ac) c(sum(sapply(ac[,1], function(x) x %in% edge.x)),
                                                    sum(sapply(ac[,2], function(x) x %in% edge.y)),
                                                    sum(sapply(ac[,3], function(x) x %in% edge.z)))))
          
          colnames(touch.sub) <- c("Edge.x","Edge.y","Edge.z")
          
          afr.list[[vv]] <- cbind(afr.list[[vv]], touch.sub)
          
          # Distance
          if(!is.null(points) & !is.null(pwidth) & !is.null(zstep)){
            afr.list[[vv]]$dist <- sqrt(((afr.list[[vv]]$x-points[[k]][1])*pwidth)^2 + 
                                          ((afr.list[[vv]]$y-points[[k]][2])*pwidth)^2 + 
                                          ((afr.list[[vv]]$z-points[[k]][3])*zstep)^2)
          }
        }
      }
      
    }
    
    if(length(channels) > 1){
      out <- list(list(afr,afr.list), ch_agg)
    } else {
      out <- list(afr, ch_agg)
    }
    return(out)
  }
  
  if(length(channels) > 1){
    results <- lapply(res.list, function(x) x[[1]][[1]])
    res.split <- lapply(res.list, function(x) x[[1]][[2]])
    res.split <- lapply(1:length(res.split[[1]]), function(x) do.call(rbind, lapply(res.split, function(y) y[[x]])))
    names(res.split) <- channels
  } else {
    results <- lapply(res.list, function(x) x[[1]])
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
  
  arrays <- lapply(res.list, function(x) x[[2]])
  
  if(length(channels) > 1){
    All <- list(Aggregates=afrxn,Aggregates_split=res.split,Arrays=arrays)
  } else {
    All <- list(Aggregates=afrxn,Arrays=arrays)
  }
  
  return(All)
}
