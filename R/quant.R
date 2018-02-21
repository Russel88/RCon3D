#' Quantifies the pixels in images
#'
#' Quantifies the number of pixels for each channel in each layer (z-stack)
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the array files
#' @param naming Optional. Add metadata to the output dataframe by looking through names of array files. Should be a list of character vectors, each list element will be added as a variable. Example: naming=list(Time=c("T0","T1","T2")). The function inserts a variable called Time, and then looks through the names of the array files and inserts characters mathcing either T0, T1 or T2
#' @param cores Number of cores for parallel computing
#' @keywords array image quantify
#' @return A dataframe with number of pixels for each layer
#' @import foreach doSNOW
#' @export

quant <- function(imgs,channels,naming = NULL,cores = 1) {
  
  # Loop for each channel
  qua.fin <- list()
  for(c in 1:length(channels)){
    
    message(paste("Channel",channels[c]))
    
    # Images
    ch_files <- imgs[grep(channels[c], imgs)]
    
    # Parallel
    if(cores > 1){
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
      on.exit(stopCluster(cl))
    } else {
      registerDoSEQ()
    }
    pb <- txtProgressBar(max = length(ch_files), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Loop foreach image
    qua <- foreach(k = 1:length(ch_files), .combine = rbind, .options.snow = opts) %dopar% {

      # Load RDS
      ch_t <- readRDS(ch_files[k])
      
      # Count for each layer
      q <- apply(ch_t, 3, sum)
      
      # Output
      temp <- cbind(sub(paste0("_Array.*"),"",sub(".*/", "", ch_files[k])),channels[c],q,1:length(q))

    }
    qua.fin[[c]] <- qua
  }
  
  qua.final <- as.data.frame(do.call(rbind, qua.fin))
  
  # Outputting results
  colnames(qua.final) <- c("Img","Channel","Count","Layer")  
  qua.final$Count <- as.numeric(as.character(qua.final$Count))
  
  if(!is.null(naming)){
    quan <- qua.final
    for(i in 1:length(naming)){
      name.temp <- naming[[i]]
      quan <- cbind(quan,NA)
      for(j in 1:length(name.temp)){
        quan[grep(name.temp[j],quan$Img),4+i] <- name.temp[j]
      }
    }
    
    colnames(quan) <- c(colnames(qua.final),names(naming))
    
  } else quan <- qua.final
  
  quan$Layer <- as.numeric(as.character(quan$Layer))
  
  return(quan)
  
}

