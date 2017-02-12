#' Extract layers 
#'
#' Function to extract layers from the LayerSplit function for use in the CrossCor or Agg analysis
#' @param splitdf A dataframe from the LayerSplit function
#' @param part What part should be extracted? Example: part="Top"
#' @param ... Arguments for subsetting the splitdf dataframe first
#' @keywords array image split
#' @return A list with lists of the extracted layers, for use in the Agg and CrossCor function
#' @export

ELayers <- function(splitdf,part,...){

  Split <- subset(splitdf,splitdf$Split == part,...)
  Part <- foreach(i = unique(Split$Exp)) %do% {
    Temp <- Split[Split$Exp == i,]
    Llist <- c(min(as.numeric(Temp$Layer)):max(as.numeric(Temp$Layer)))
    return(Llist)
  }
  return(Part)
}

