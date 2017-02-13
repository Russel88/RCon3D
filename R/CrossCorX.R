#' 3D cross-correlation
#'
#' Wrapper function to run CrossCor several times
#' @param ... Arguments to the CrossCor function
#' @param R Number of times to run the cross-correlation analysis
#' @keywords array image cross-correlation
#' @return A dataframe with the cross-correlation vaules for each distance
#' @export
#' 
CrossCorX <- function(...,R=NULL){
  
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  CC.all <- list()
  for(r in 1:R){
    CC <- CrossCor(...)
    CC$R <- r
    CC.all[[r]] <- CC
    # Update progress bar
    setTxtProgressBar(pb, r)
  }
  CCall <- do.call("rbind", CC.all)
  close(pb)
  return(CCall)
}
