#' Apply mathematical morphology kernel on image and save as RDS
#'
#' Morph images as a preparation for subsequent analyses
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param morph One of "erode", "dilate", "opening", or "closing"
#' @param kern Numeric vector indicating range of morpinh in the x,y,z directions. Has to be odd integers.
#' @param type Type of kernel for morphing "box" includes diagonals, "diamond" is without diagonals
#' @param cores Integer. Number of cores to use for parallel computing.
#' @keywords array image
#' @return Creates arrays as RDS files and outputs the paths for these files
#' @import foreach doSNOW mmand
#' @export

morphIMG <- function(imgs,morph=NULL,kern=c(3,3,3),type="box",cores=1) {

  stopifnot(!is.null(type))
  
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

    # Morph
    if(kern[1] %% 2 == 1 &
       kern[2] %% 2 == 1 &
       kern[3] %% 2 == 1) {
    kern.s <- shapeKernel(kern, type = type)
    
    ch_t <- switch(morph,
      erode = erode(ch_t,kern.s),
      dilate = dilate(ch_t,kern.s),
      opening = opening(ch_t,kern.s),
      closing = closing(ch_t,kern.s)
    )
    ch_t[ch_t > 0] <- 1 } else stop("Kernel has to be odd integers in all directions")

    # Save
    saveRDS(ch_t, file = gsub("Array.R","MorphArray.R",imgs[k]))
    
    return(NULL)
  }
  
  files <- gsub("Array.R","MorphArray.R",imgs)

  return(files)
}
