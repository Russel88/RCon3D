#' Find distances between 3D aggregates/clumps
#'
#' With two outputs from the \code{clumps} function compute the distances between pairs of aggregates. 
#' Computes the distances between aggregates in \code{clumps.x} with those in \code{clumps.y}.
#' If you want distances between aggregates from one output, just input the data.frame in \code{clumps.x}.
#' Provide pwidth and zstep to get the distance in microns. Alternatively it's the pixel distance
#' @param clumps.x Data.frame output from the \code{clumps} function (...$Aggregates) 
#' @param clumps.y Optional. Data.frame output from the \code{clumps} function (...$Aggregates) 
#' @param pwidth Width of pixels in microns
#' @param zstep z-step in microns
#' @param cores Number of cores to use for parallel computing
#' @keywords array image aggregate
#' @return A dataframe with the distances
#' @import foreach doSNOW
#' @export
clumps_dist <- function(clumps.x, clumps.y = clumps.x, pwidth = 1, zstep = 1, cores = 1){
  
  if(!all(unique(clumps.x$Img) == unique(clumps.y$Img))){
    stop("clumps.x and clumps.y has to be from the same images (`Img` columns should be similar)")
  }
  stopifnot(is.data.frame(clumps.x),
            is.data.frame(clumps.y))
  
  if(cores > 1){
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }
  
  pb <- txtProgressBar(max = length(unique(clumps.x$Img)), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  out <- foreach(i = 1:length(unique(clumps.x$Img)), .combine = rbind, .options.snow = opts) %dopar% {
    
    # Subset
    afr.x <- clumps.x[clumps.x$Img == unique(clumps.x$Img)[i],]
    afr.y <- clumps.y[clumps.y$Img == unique(clumps.y$Img)[i],]
    
    # Combinations
    agg.mat <- expand.grid(afr.x$ID, afr.y$ID, stringsAsFactors = FALSE)
    # agg.mat <- as.data.frame(apply(agg.mat, 2, as.numeric))
    
    # Distances
    dist.fun <- function(x,y){sqrt(((x[1]-y[1])*pwidth)^2 + ((x[2]-y[2])*pwidth)^2 + ((x[3]-y[3])*zstep)^2)}
    
    agg.mat$dist <- apply(agg.mat, 1, function(am) dist.fun(unlist(afr.x[afr.x$ID == am[1],c("x","y","z")]),
                                                                  unlist(afr.y[afr.y$ID == am[2],c("x","y","z")])))
    colnames(agg.mat) <- c("ID.x","ID.y","dist")
    agg.mat$Img <- unique(clumps.x$Img)[i]
    
    return(agg.mat)
  }

  return(out)  

}
