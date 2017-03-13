#' 3D Occupancy
#'
#' Function to calculate the occupancy of a target channel at certain distances from a focal channel
#' @param ... Arguments for the \code{occupancy.default} function
#' @param R Number of times to run the occupancy analysis
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param focal.channel Name of the channel from which the distance is calculated. Should be in the names of the array files
#' @param target.channel Name of the channel for calculate occupancy of. Should be in the names of the array files
#' @param size The maximum distance (microns) to examine. Has to be a multiple of both pwidth and zstep. Beware, increasing size will increase runtime exponetially!
#' @param npixel Number of random pixels to examine. Increasing this will increase precision (and runtime in a linear fashion)
#' @param dstep The interval between examined distances (microns). Increasing this decreases resolution but speeds up function linearly. Defaults to 1
#' @param pwidth Width of pixels in microns
#' @param zstep z-step in microns
#' @param freec The number of cores NOT to use. Defaults to 1
#' @param kern.smooth Optional. Numeric vector indicating range of median smoothing in the x,y,z directions. Has to be odd intergers. c(1,1,1) means no smoothing.
#' @param layers Optional. Should the function only look in a subset of layers. A list with lists of layers to use for each image. Can also be the output from \code{extract_layers} 
#' @param naming Optional. Add metadata to the output dataframe by looking through names of array files. Should be a list of character vectors, each list element will be added as a variable. Example: naming=list(Time=c("T0","T1","T2")). The function inserts a variable called Time, and then looks through the names of the array files and inserts characters mathcing either T0, T1 or T2
#' @keywords array image occupancy
#' @return A dataframe with the occupancy values for each distance
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export

occupancy <- function(...,R=NULL){
  
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  OCC.all <- list()
  for(r in 1:R){
    OCC <- occupancy.default(...)
    OCC$R <- r
    OCC.all[[r]] <- OCC
    # Update progress bar
    setTxtProgressBar(pb, r)
  }
  OCCall <- do.call("rbind", OCC.all)
  close(pb)
  return(OCCall)
}

#' @rdname occupancy
#' @export

occupancy.default <- function(imgs,focal.channel,target.channel,size,npixel,dstep=1,pwidth,zstep,freec=1,kern.smooth=NULL,layers=NULL,naming=NULL) {
  
  if(size%%zstep != 0) if(all.equal(size%%zstep,zstep)!=TRUE) stop("size not a multiple of zstep")
  if(size%%pwidth != 0) if(all.equal(size%%pwidth,pwidth)!=TRUE) stop("size not a multiple of pwidth")
  
  # Null box (pixels)
  null_box <- expand.grid(x = seq((-size/pwidth), (size/pwidth), by = 1), 
                        y = seq((-size/pwidth), (size/pwidth), by = 1), 
                          z = seq((-size/zstep), (size/zstep), by = 1))
  
  # Distances from focal pixel to each pixel in box (distances in microns, position in pixel space)
  d <- apply(null_box,1,function(nbox) sqrt((0 - (nbox[1]*pwidth))^2 +
                                              (0 - (nbox[2]*pwidth))^2 +
                                              (0 - (nbox[3]*zstep))^2))
  
  # Bind positions in box with distances
  boxd <- cbind(null_box,d)
  
  # Bins of distances (microns)
  ds <- seq(0, size, by = dstep)
  
  # Load images
  ch.f_files <- imgs[grep(focal.channel, imgs)]
  ch.t_files <- imgs[grep(target.channel, imgs)]
  
  # Start parallel for each replica
  cl <- makeCluster(detectCores()-freec)
  registerDoParallel(cl)  
  OCC_results <- foreach(i=1:length(ch.f_files),.combine=rbind) %dopar% {
    
    # Load
    ch.f <- readRDS(ch.f_files[i])
    ch.t <- readRDS(ch.t_files[i])
    
    # Subset
    if(!is.null(layers[[i]])){
      ch.f <- ch.f[,,layers[[i]]]
      ch.t <- ch.t[,,layers[[i]]]
    }
    
    # Smooth
    if(!is.null(kern.smooth)) {
      if(kern.smooth[1] %% 1 == 0 & kern.smooth[1] %% 2 != 0 &
         kern.smooth[2] %% 1 == 0 & kern.smooth[2] %% 2 != 0 &
         kern.smooth[3] %% 1 == 0 & kern.smooth[3] %% 2 != 0) {
        kern.s <- kernelArray(array(1,dim=kern.smooth))
        ch.f <- medianFilter(ch.f,kern.s)
        ch.f[ch.f > 0] <- 1
        ch.t <- medianFilter(ch.t,kern.s)
        ch.t[ch.t > 0] <- 1} else stop("Kernel smooth has to be odd integers in all directions")
    }
    
    # Density of channel
    dens <- length(which(ch.t == 1))/length(ch.t)
    
    # Addresses in array (pixels)
    side <- dim(ch.f)[1]
    address_array <- array(1:(side*side*dim(ch.f)[3]), 
                           c(side, side, dim(ch.f)[3]))
    
    # Coordinates of pixels in focal channel (pixels)
    chf_add <- data.frame(which(ch.f == 1, T))
    colnames(chf_add) <- c("x", "y", "z")
    
    # Randomly sample pixels (pixels)
    these <- sample(1:dim(chf_add)[1], size = npixel)
    
    # Get their addresses
    ch_pix <- chf_add[these,]
    
    # Matrix to collect results
    hits <- matrix(NA, length(ds), npixel)
    totals <- matrix(NA, length(ds), npixel)
    
    # Loop through pixels
    for(j in 1:npixel){
      
      # Focal pixel position
      p <- ch_pix[j,]
      
      # Coordinates of the box (pixels)
      xrange <- c(p$x-(size/pwidth), p$x+(size/pwidth))
      yrange <- c(p$y-(size/pwidth), p$y+(size/pwidth))
      zrange <- c(p$z-(size/zstep), p$z+(size/zstep))
      
      # Edge handling
      xrange[xrange < 1] <- 1
      xrange[xrange > dim(address_array)[1]] <- dim(address_array)[1]
      
      yrange[yrange < 1] <- 1
      yrange[yrange > dim(address_array)[2]] <- dim(address_array)[2]
      
      zrange[zrange < 1] <- 1
      zrange[zrange > dim(address_array)[3]] <- dim(address_array)[3]
      
      # Addresses of the array (pixels)
      box <- address_array[c(xrange[1]:xrange[2]), 
                           c(yrange[1]:yrange[2]), 
                           c(zrange[1]:zrange[2])]  
      
      # New box for specific pixel (changes depending on position relative to the edges)
      new_box <- expand.grid(x = (seq(xrange[1], xrange[2], by = 1)-p$x), 
                             y = (seq(yrange[1], yrange[2], by = 1)-p$y), 
                             z = (seq(zrange[1], zrange[2], by = 1)-p$z))
      
      # Find distances to positions in new box
      sub_box <- boxd[boxd$x %in% new_box$x & boxd$y %in% new_box$y & boxd$z %in% new_box$z,]
      
      new_d <- sub_box$d
      
      # Presence absence data for the box (pixels)
      id <- ch.t[box]
      
      # Postitions of pixels
      positions <- list()
      for(p in 1:(length(ds))){
        positions[[p]] <- which(new_d <= ds[p]+dstep/2 & new_d > ds[p]-dstep/2)
      }
      
      # Total counts
      for(l in 1:(length(ds))){
        totals[l,j] <- length(which(id[positions[[l]]]==0))+length(which(id[positions[[l]]]==1))
      }
      
      # Hits
      for(l in 1:(length(ds))){
        hits[l,j] <- length(which(id[positions[[l]]]==1))
      }
    }
    
    # Sum totals and hits for all pixels
    totals.sum <- apply(totals,1,sum,na.rm=T)
    hits.sum <- apply(hits,1,sum,na.rm=T)
    
    # Calculate occupancy
    occup <- hits.sum/totals.sum
    occup.norm <- occup/dens
    
    theseOCC <- cbind(gsub(focal.channel,"",sub(".*/", "", ch.f_files[i])),ds, occup, occup.norm)
    theseOCC <- as.data.frame(theseOCC)
    colnames(theseOCC) <- c("Img", "Distance", "Occupancy","Occupancy.Normalized")
    return(theseOCC)
  }
  stopCluster(cl)
  
  # Final polishing
  colnames(OCC_results) <- c("Img", "Distance", "Occupancy","Occupancy.Normalized")
  OCC_results$Occupancy <- as.numeric(as.character(OCC_results$Occupancy))
  OCC_results$Occupancy.Normalized <- as.numeric(as.character(OCC_results$Occupancy.Normalized))
  OCC_results$Distance <- as.numeric(as.character(OCC_results$Distance))
  
  OCC_results$Focal <- focal.channel
  OCC_results$Target <- target.channel
  
  if(!is.null(naming)){
    OCC_resultsx <- OCC_results
    for(i in 1:length(naming)){
      name.temp <- naming[[i]]
      OCC_resultsx <- cbind(OCC_resultsx,NA)
      for(j in 1:length(name.temp)){
        OCC_resultsx[grep(name.temp[j],OCC_resultsx$Img),ncol(OCC_results)+i] <- name.temp[j]
      }
    }
    
    colnames(OCC_resultsx) <- c(colnames(OCC_results),names(naming))
    
  } else OCC_resultsx <- OCC_results
  
  
  return(OCC_resultsx)
  
}
