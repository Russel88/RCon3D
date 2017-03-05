#' Extract layers 
#'
#' Function to extract layers from the \code{layer_split} function for use in the \code{co_agg}, \code{occupancy} or \code{clumps} analysis
#' @param splitdf A dataframe from the \code{layer_split} function
#' @param part What part should be extracted? Example: part="Top"
#' @param ... Arguments for subsetting the splitdf dataframe first
#' @keywords array image split
#' @return A list with lists of the extracted layers, for use in the \code{co_agg}, \code{occupancy} or \code{clumps} functions
#' @export

extract_layers <- function(splitdf,part,...){

  Split <- subset(splitdf,splitdf$Split == part,...)
  Part <- foreach(i = unique(Split$Exp)) %do% {
    Temp <- Split[Split$Exp == i,]
    Llist <- c(min(as.numeric(Temp$Layer)):max(as.numeric(Temp$Layer)))
    return(Llist)
  }
  return(Part)
}

