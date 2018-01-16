#' Quantifies the pixels in images
#'
#' Quantifies the number of pixels for each channel in each layer (z-stack)
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the array files
#' @param naming Optional. Add metadata to the output dataframe by looking through names of array files. Should be a list of character vectors, each list element will be added as a variable. Example: naming=list(Time=c("T0","T1","T2")). The function inserts a variable called Time, and then looks through the names of the array files and inserts characters mathcing either T0, T1 or T2
#' @keywords array image quantify
#' @return A dataframe with number of pixels for each layer
#' @export

quant <- function(imgs,channels,naming=NULL) {
  
  # DF for storing data
  qua <- data.frame(Img = NULL, Channel = NULL, Count = NULL, Layer = NULL)
  
  # Loop for each channel
  for(c in 1:length(channels)){
    
    # Load images
    ch_files <- imgs[grep(channels[c], imgs)]
    
    # Loop for each image
    for(k in 1:length(ch_files)){
      
      # Load RDS
      ch_t <- readRDS(ch_files[k])
      
      # Count for each layer
      for(i in 1:dim(ch_t)[3]){
        
        q <- length(which(ch_t[,,i] > 0))
        
        temp <- cbind(sub(paste0("_Array.*"),"",sub(".*/", "", ch_files[k])),channels[c],q,i)
        
        qua <- rbind(qua,temp)
      }
      
    }
    
  }
  
  # Outputting results
  colnames(qua) <- c("Img","Channel","Count","Layer")  
  qua$Count <- as.numeric(as.character(qua$Count))
  
  if(!is.null(naming)){
    quan <- qua
    for(i in 1:length(naming)){
      name.temp <- naming[[i]]
      quan <- cbind(quan,NA)
      for(j in 1:length(name.temp)){
        quan[grep(name.temp[j],quan$Img),4+i] <- name.temp[j]
      }
    }
    
    colnames(quan) <- c(colnames(qua),names(naming))
    
  } else quan <- qua
  
  quan$Layer <- as.numeric(as.character(quan$Layer))
  
  return(quan)
  
}

