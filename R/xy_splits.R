#' Analyses for each xy position
#'
#' Splits the image in each xy position and then run an analysis
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param channels Character vector with name(s) of channels. Channel names should be in the names of the array files
#' @param do What to do? \itemize{
#'  \item "max" and "min" find the max or min layer where the first 1 is found.
#'  \item "sum" sums the pixels (e.g. for biomass xy-distribution). 
#'  \item "which" returns z-positions of pixels for each xy dimension. 
#'  \item "section" first sections the image in an upper/lower part and then sums for each part. Sectioning is done for each xy-position and therefore takes a variable biomass distribution across the image into account.
#'  \item Alternatively a user-defined function that will be applied for each xy-position. E.g. function(x) median(which(x > 0)) will return the median position of each channel at each xy-position.
#'  } 
#' @param upper.part Only if do="section". Either the fraction to assign to the upper section (between 0 and 1), or a fixed number of layers that are assigned to the upper part (an integer, 1L, 2L, 3L etc.. Remember the L!). The upper fraction is inclusive (i.e. x >= upper.part)
#' @param layer.start Only if do="section". Are first layers the "Top" or the "Bottom" of the specimen?
#' @details For sectioning, NA is returned if a channel is absent at a specific xy-position, but a zero is returned if the channel is present at that xy-position but absent in the section.
#' @keywords array image
#' @return A recursive list. First layer are the images, within those are the channels. Each element is a matrix with a value for each xy position. If do="section" an extra layer is added denoting the upper and lower section.
#' @export

xy_splits <- function(imgs,channels,do,upper.part=0.5,layer.start=NULL) {
  
  # Check
  if(do == "section") {
    if(!layer.start %in% c("Bottom","Top")) stop("Set layer.start to either Top or Bottom")
  }
  
  # Results
  results <- list()

  # Options
  if(!is.function(do)){
    if(do == "section") opts <- "which" else opts <- do 
  } else opts <- FALSE
    
  # Loop for each channel
  for(ch in 1:length(channels)){
    
    # Load images
    ch_files <- imgs[grep(channels[ch], imgs)]
    
    # Results
    result.image <- list()
    
    # Loop for each image
    for(k in 1:length(ch_files)){
      
      # Load RDS
      ch_t <- readRDS(ch_files[k])
      
      # Analyse
      if(!is.function(do)){
        res <- switch(opts,
                      min = suppressWarnings(apply(ch_t,c(1,2),function(x) min(which(x > 0)))),
                      max = suppressWarnings(apply(ch_t,c(1,2),function(x) max(which(x > 0)))),
                      sum = apply(ch_t,c(1,2),sum),
                      which = apply(ch_t,c(1,2),function(x) which(x > 0)))
        
        if(opts != "which"){
          res[is.infinite(res)] <- NA
        }
        
      } else {
        res <- apply(ch_t,c(1,2),do)
      }
      
      # Collect results
      result.image[[k]] <- res
      
    }
    results[[ch]] <- result.image
  }

  # Sectioning
  if(is.function(do)) do <- FALSE
  if(do == "section"){
    
    section.results <- list()
    
    for(s in 1:length(ch_files)){
      
      # Extract
      chs <- lapply(1:length(channels),function(ch){
        results[[ch]][[s]]
      })
      
      # Lengths
      chs.l <- simplify2array(lapply(1:length(channels),function(ch){
        apply(chs[[ch]],c(1,2),function(x) length(unlist(x)))
      }))

      # Combine the channels
      chs.combined <- apply(simplify2array(chs),c(1,2),unlist)
      
      # Upper and lower parts
      if(is.integer(upper.part)){
        
        # Functions
        top <- function(x, n=upper.part){
          if(length(x) > 0) {
            result <- list()
            for(i in 1:n){
              j <- which(x == max(x, na.rm = TRUE))
              result[[i]] <- x[j]
              x[j] <- -Inf
            }
            return(unlist(result))
          } else {
            return(0)
          }
        }
        bottom <- function(x, n=upper.part){
          if(length(x) > 0) {
            result <- list()
            for(i in 1:n){
              j <- which(x == min(x, na.rm = TRUE))
              result[[i]] <- x[j]
              x[j] <- Inf
            }
            return(unlist(result))
          } else {
            return(0)
          }
        }
        
        if(layer.start == "Bottom") {

          chs.upper <- apply(chs.combined,c(1,2),function(x) {
            unlist(x) %in% top(unlist(x))})

          chs.lower <- apply(chs.combined,c(1,2),function(x) {
            !unlist(x) %in% top(unlist(x))})
         
        }
        
        if(layer.start == "Top") {
          
          chs.upper <- apply(chs.combined,c(1,2),function(x) {
            unlist(x) %in% bottom(unlist(x))})

          chs.lower <- apply(chs.combined,c(1,2),function(x) {
            !unlist(x) %in% bottom(unlist(x))})
          
        }
        
      } else {
        
        if(layer.start == "Bottom") {
          
          chs.upper <- apply(chs.combined,c(1,2),function(x) {
            unlist(x) >= quantile(unlist(x),(1-upper.part))})
          
          chs.lower <- apply(chs.combined,c(1,2),function(x) {
            unlist(x) < quantile(unlist(x),(1-upper.part))})
          
        } 
        
        if(layer.start == "Top") {

          chs.upper <- apply(chs.combined,c(1,2),function(x) {
            unlist(x) < quantile(unlist(x),upper.part)})
          
          chs.lower <- apply(chs.combined,c(1,2),function(x) {
            unlist(x) >= quantile(unlist(x),upper.part)})
          
        }

      }
      
      # Sums
      upper.sums <- lapply(1:length(channels),function(ch){
        sapply(1:dim(chs.upper)[1],function(x){
          sapply(1:dim(chs.upper)[2],function(y){
            sum(chs.upper[x,y][[1]][(sum(chs.l[x,y,0:(ch-1)])+1):sum(chs.l[x,y,(ch-1):ch])])
          })
        })
      })
      
      lower.sums <- lapply(1:length(channels),function(ch){
        sapply(1:dim(chs.lower)[1],function(x){
          sapply(1:dim(chs.lower)[2],function(y){
            sum(chs.lower[x,y][[1]][(sum(chs.l[x,y,0:(ch-1)])+1):sum(chs.l[x,y,(ch-1):ch])])
          })
        })
      })
    
      names(upper.sums) <- channels
      names(lower.sums) <- channels
      
      # Collect
      section.results[[s]] <- list(Upper=upper.sums,Lower=lower.sums)
      
    }
    
    final <- section.results
    
  } else {
    
    final <- list()
    for(i in 1:length(ch_files)){
      sublist <- list()
      for(ch in 1:length(results)){
        sublist[[ch]] <- results[[ch]][[i]]
      }
      names(sublist) <- channels
      final[[i]] <- sublist
    }
    
  }
  
  # Names
  naming <- unique(grep(paste(channels,collapse="|"),imgs, value=TRUE))
  for(i in 1:length(channels)){
    naming <- gsub(channels[i],"", naming)
  }
  names(final) <- unique(gsub("_Array.R","",gsub(".*/","",naming)))
  
  return(final)
  
}



