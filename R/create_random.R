#' Create random image
#'
#' Create array of random image
#' @param path The path of where to save the RDS files
#' @param overlap Logical. Are the channels allowed to overlap, such that they can be at the same pixels
#' @param side Pixels in x,y dimension, a numeric 
#' @param h Height of image, i.e. number of layers
#' @param probs Numeric vector with an element for each channel, these elements gives the probability for each pixel of being assigned to each channel. 
#' @keywords array image
#' @return Creates arrays as RDS files in the specified path, and outputs the paths for these files
#' @export

create_random <- function(path, overlap, side, h, probs){

  # Create empty array
  ch0 <- array(0, c(side, side, h))
  
  chs <- list()
  
  # Randomly place 1s in array 
  if(overlap==TRUE){  
    
    for(c in 1:length(probs)){
      chs[[c]] <- apply(ch0, c(1,2,3), function(x) sample(c(0,1),1,prob=c((1-probs[c]),probs[c])))
    }  

    } else {
      ch <- apply(ch0, c(1,2,3), function(x) sample(seq(0,length(probs),1),1,prob=c((1-sum(probs)),probs)))
      
      for(i in 1:length(probs)){
        
        ch.x <- ch
        ch.x[ch.x != i] <- 0
        ch.x[ch.x == i] <- 1

        chs[[i]] <- ch.x
        
      }
      
    }
  
  # Save RDS
  for(k in 1:length(probs)){
    saveRDS(chs[[k]], file = paste0(path,"/RandomImage_Ch",k,".R"))
  }
  
  # Find the files
  files <- list.files(path, "RandomImage_Ch", full.names = T)
  
  return(files)
  
}