#' Median smooth image and save as RDS
#'
#' Smooth images as a preparation for subsequent analyses
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param kern.smooth Numeric vector indicating range of median smoothing in the x,y,z directions. Has to be odd intergers.
#' @param type.smooth Type of kernel for smooth "box" includes diagonals, "diamond" is without diagonals
#' @param cores Integer. Number of cores to use for parallel computing.
#' @keywords array image
#' @return Creates arrays as RDS files in the set path, and outputs the paths for these files
#' @import foreach doSNOW mmand
#' @export

smoothIMG <- function(imgs,kern.smooth=c(3,3,3),type.smooth="box",cores=1) {

  # For each image
  if(cores > 1){
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }
  
  pb <- txtProgressBar(max = length(imgs), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  res.list <- foreach(k = 1:length(imgs), .options.snow = opts, .packages = "mmand") %dopar% {
    
    # Load
    ch_t <- readRDS(imgs[k])

    # Smooth
    if(kern.smooth[1] %% 2 == 1 &
       kern.smooth[2] %% 2 == 1 &
       kern.smooth[3] %% 2 == 1) {
    kern.s <- shapeKernel(kern.smooth, type = type.smooth)
    ch_t <- medianFilter(ch_t,kern.s)
    ch_t[ch_t > 0] <- 1 } else stop("Kernel smooth has to be odd integers in all directions")

    # Save
    saveRDS(ch_t, file = gsub("_Array.R","_SmoothArray.R",imgs[k]))
    
    return(NULL)
  }
  
  files <- gsub("_Array.R","_SmoothArray.R",imgs)

  return(files)
}
