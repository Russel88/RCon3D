#' Apply thresholding on image and save as RDS
#'
#' A single threshold is found for each image/z-stack. 
#' @param imgs The paths of array files; i.e. output from \code{loadIMG} or \code{findIMG} functions.
#' @param method Either Manual, Otsu, or BEM
#' @param cores Integer. Number of cores to use for parallel computing.
#' @param breaks Integer. Number of breaks for Otsu's threshold. Increase for increased precision
#' @param opt Numeric. Between 0 and 1. Parameter for BEM method. Increase to decrease threshold
#' @param thresh Numeric. Vector of thresholds for manual thresholding, one for each image.
#' @param BEM.opt List. A list of options for model fitting for the BEM method. Passed to the nlsLM function.
#' @keywords array image
#' @return Prints thresholds. saves thresholded RDS files and outputs the path for the thresholded images.
#' @import foreach doSNOW minpack.lm stats utils
#' @importFrom graphics hist.default
#' @export

threshIMG <- function(imgs, method = "Manual", cores = 1, breaks = 100, opt = 0.1, thresh = rep(0,length(imgs)), BEM.opt = list(start = list(a=-1e6, b=0.1, c=1e6), control = nls.lm.control(maxiter=1e3))) {
  
  if(is.null(method)) stop("Method is missing")
  if(!method %in% c("Otsu","BEM","Manual")) stop("Method has to be either Otsu, BEM or Manual")
  if(method == "Manual" & length(imgs) != length(thresh)) stop("Length of thresh has to equal the number of images")
  
  # Progress bar
  pb <- txtProgressBar(max = length(imgs), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Start parallel
  if(cores == 1) {
    registerDoSEQ() 
  } else {
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  }
  
  i <- NULL
  res <- foreach(i=1:length(imgs), .options.snow = opts, .packages = "minpack.lm") %dopar% {
   
    # Define functions
    BEM <- function(x, opt = 0.1, BEM.opt){
      bit8 <- TRUE
      if(max(x) <= 1){
        bit8 <- FALSE
        x <- round(x*255)
      }
      if(max(x) > 255 | min(x) < 0) stop("Intensity values has to either be between 0 and 255 or between 0 and 1")
      tt <- 0:max(x)
      xs <- sapply(tt, function(i) sum(x > i))
      z <- do.call(nlsLM, c(list(xs ~ a*tt^b+c),BEM.opt))
      zd <- coef(z)[1] * (tt^(coef(z)[2] - 1) * coef(z)[2])
      thr.adj <- (c(0,zd[1:(length(zd)-1)]) - zd)/zd
      if(bit8){
        return(tt[which(thr.adj < opt)[1]])
      } else {
        return(tt[which(thr.adj < opt)[1]]/255)
      }
    }
    
    otsu <- function(y, breaks) {
      # Borrowed from the EBImage package
      h = hist.default(y, breaks = breaks, plot = FALSE)
      counts = as.double(h$counts)
      mids = as.double(h$mids)
      len = length(counts)
      w1 = cumsum(counts)
      w2 = w1[len] + counts - w1
      cm = counts * mids
      m1 = cumsum(cm)
      m2 = m1[len] + cm - m1
      var = w1 * w2 * (m2/w2 - m1/w1)^2
      maxi = which(var == max(var, na.rm = TRUE))
      return((mids[maxi[1]] + mids[maxi[length(maxi)]] ) / 2)
    }
    
    # Find thresholds
    ch <- readRDS(imgs[i])
    
    if(method == "Otsu") thr <- otsu(as.vector(ch), breaks)
    if(method == "BEM") thr <- BEM(as.vector(ch), opt, BEM.opt)
    if(method == "Manual") thr <- thresh[i]
    
    ch[ch >= thr] <- 1
    ch[ch < thr] <- 0
    
    fnam <- gsub("Array.R","ThreshArray.R",imgs[i])
    
    saveRDS(ch, file = fnam)
    
    return(list(fnam,thr))
    
  }
  if(cores > 1) stopCluster(cl)
  
  thresh <- sapply(res, function(x) x[[2]])
  files <- sapply(res, function(x) x[[1]])
  
  print(data.frame(Image = gsub(".*\\/","",files), Threshold = thresh),row.names = FALSE)
  
  return(files)
}