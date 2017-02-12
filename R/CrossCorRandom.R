#' 3D random cross-correlation
#'
#' Function to calculate the pairwise 3D cross-correlation between two random channels
#' @param overlap Are the channels allowed to overlap, such that they can be at the same pixels
#' @param size The maximum distance (microns) to examine. Has to be a multiple of both pwidth and zstep. Beware, increasing size will increase runtime exponetially!
#' @param npixel Number of random pixels to examine. Increasing this will increase precision (and runtime in a linear fashion)
#' @param dstep The interval between examined distances (microns). Increasing this decreases resolution but speeds up function linearly. Defaults to 1
#' @param side Pixels in x,y dimension, a numeric 
#' @param pwidth Width of pixels in microns
#' @param zstep z-step in microns
#' @param R Number of runs
#' @param p1 For each pixel, the probability of being assigned to channel 1. 
#' @param p2 For each pixel, the probability of being assigned to channel 2. 
#' @param freec The number of cores NOT to use.Defaults to 1
#' @keywords array image cross-correlation
#' @return A dataframe with the cross-correlation vaules for each distance
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export

CrossCorRandom <- function(overlap,size,npixel,dstep,side,h,pwidth,zstep,R,p1=1/3,p2=1/3,freec=1) {
  
  stopifnot(size%%zstep == 0)
  stopifnot(size%%pwidth == 0)
  
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
  
  # Create empty array (images)
  ch0_t <- array(0, c(side, side, h))
  
  # Start parallel for each replica
  cl <- makeCluster(detectCores()-freec)  
  registerDoParallel(cl)  
  cc_results <- foreach(i=1:R,.combine=rbind) %dopar% {
    
    # Randomly place 1s in array 
    if(overlap==TRUE){  
      ch1_t <- apply(ch0_t, c(1,2,3), function(x) sample(c(0,1),1,prob=c((1-p1),p1)))
      ch2_t <- apply(ch0_t, c(1,2,3), function(x) sample(c(0,1),1,prob=c((1-p2),p2))) } else {
        ch <- apply(ch0_t, c(1,2,3), function(x) sample(c(0,1,2),1,prob=c((1-p1-p2),p1,p2)))
        ch1_t <- ch
        ch2_t <- ch
        ch1_t[ch1_t == 2] <- 0
        ch2_t[ch2_t == 1] <- 0
        ch2_t[ch2_t == 2] <- 1
      }
    
    # Densities of channels
    d1 <- length(which(ch1_t == 1))/length(ch1_t)
    d2 <- length(which(ch2_t == 1))/length(ch2_t)
    
    # Addresses in array (pixels)
    address_array <- array(1:(side*side*h), 
                           c(side, side, h))
    
    # Coordinates of pixels in channel1 (pixels)
    ch1_add <- data.frame(which(ch1_t == 1, T))
    colnames(ch1_add) <- c("x", "y", "z")
    
    # Coordinates of pixels in channel2 (pixels)
    ch2_add <- data.frame(which(ch2_t == 1, T))
    colnames(ch2_add) <- c("x", "y", "z")
    
    # Coordinates of all pixels (pixels)
    ch_add <- data.frame(which(ch1_t == 1 | ch1_t == 0, T))
    colnames(ch_add) <- c("x", "y", "z")
    
    # Randomly sample pixels (pixels)
    these <- sample(1:dim(ch_add)[1], size = npixel)
    
    # Get their addresses
    ch_pix <- ch_add[these,]
    
    # Matrix to collect results
    hits <- matrix(NA, length(ds), npixel)
    totals <- matrix(NA, length(ds), npixel)
    
    # Loop through pixels
    for(j in 1:npixel){
      
      # Focal pixel position
      p <- ch_pix[j,]
      
      # Does the focal pixel have a colour?
      test1 <- ch1_add[ch1_add$x == p$x & ch1_add$y == p$y & ch1_add$z == p$z,]
      test2 <- ch2_add[ch2_add$x == p$x & ch2_add$y == p$y & ch2_add$z == p$z,]
      
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
      id1 <- ch1_t[box]
      id2 <- ch2_t[box]
      
      # Postitions of pixels
      positions <- list()
      for(p in 1:(length(ds))){
        positions[[p]] <- which(new_d <= ds[p]+dstep/2 & new_d > ds[p]-dstep/2)
      }
      
      # Total counts
      for(l in 1:(length(ds))){
        totals[l,j] <- length(which(id1[positions[[l]]]==0))+length(which(id1[positions[[l]]]==1))
      }
      
      # Hits
      if((nrow(test1)+nrow(test2))==0) {hits[,j] <- 0} else {
        
        for(l in 1:(length(ds))){
          if(nrow(test1)>0 & nrow(test2)==0) hits[l,j] <- length(which(id2[positions[[l]]]==1))
          if(nrow(test2)>0 & nrow(test1)==0) hits[l,j] <- length(which(id1[positions[[l]]]==1))
          if(nrow(test1)>0 & nrow(test2)>0) hits[l,j] <- length(which(id1[positions[[l]]]==1))+length(which(id2[positions[[l]]]==1))      
        }
        
      }
      
    }
    
    # Sum totals and hits for all pixels
    totals.sum <- apply(totals,1,sum,na.rm=T)
    hits.sum <- apply(hits,1,sum,na.rm=T)
    
    # Calculate probability and cross-correlation
    Prop <- hits.sum/totals.sum
    CC <- Prop/(2*d1*d2)
    
    theseCC <- cbind(i,ds, CC)
    theseCC <- as.data.frame(theseCC)
    colnames(theseCC) <- c("R","Distance", "CC")
    
    return(theseCC)
    
  }
  
  return(cc_results)
  
}