#' Split quantification in Top, Middle and Bottom
#'
#' Function to split quantification in Top, Middle and Bottom layer by proportion occupied
#' @param qua A dataframe from the \code{quant} function
#' @param side xy dimension of image in pixels
#' @param pt The Top is defined from layer 1 (or max if layer.start="Bottom") until, but not including, the layer where the image is filled by the proportion pt (0.9 = 90 percent). If pt is set higher than found in any layer in the image then add.b specifies the Bottom and the rest is the Top. Note that if pt, add.t and add.b is set such that there is overlap between Bottom, Middle and Top, then Top will override the others 
#' @param add.t How many layers to add to the top, after pt calculation (e.g. pt=0 and add.t=10 will specify the top as the top 10 layers)
#' @param add.b How many layers are specified as the bottom
#' @param channel Name of channel to focus on. If NULL will aggregate all channels before splitting
#' @param trim Should the bottom be trimmed? Specifies the bottom as the first layer (from the bottom) with max counts. Defaults to FALSE
#' @param layer.start Are first layers the "Top" or the "Bottom" of the specimen?
#' @keywords array image split
#' @return The original qua dataframe including an Exp variable, with a unique name for each image (based on the naming given in the \code{quant} function). And a Split varaible denoting the split
#' @export

layer_split <- function(qua,side,pt,add.t,add.b,channel=NULL,trim=FALSE,layer.start=NULL){
  
  # Make Exp variable in qua dataframe
  if(length(colnames(qua)) > 4){
  cols <- colnames(qua)[5:length(colnames(qua))]
  qua$Exp <- do.call(paste, c(qua[cols], sep="_"))
  } else qua$Exp <- 1
  
  LS <- foreach(i=unique(qua$Exp),.combine=rbind) %do% {
    
    Temp <- qua[qua$Exp == i,]
    
    # Subset relevant channel  
    if(!is.null(channel)) {
      Sub <- Temp[Temp$Channel == channel,]
      Layer <- aggregate(Count ~ Layer,data=Sub,FUN=sum)
    } else {
      Layer <- aggregate(Count ~ Layer,data=Temp,FUN=sum)
    }
    
    # If layer zero is the top of the specimen
    if(layer.start == "Top"){
      
      # Trimming
      if(trim==TRUE) Layerx <- Layer[1:Position(function(x) max(Layer$Count) == x,Layer$Count,right=TRUE),] else Layerx <- Layer
      
      # Calculate fill
      Layerx$Fill <- Layerx$Count / side^2
      
      # Split in layers
      if(max(Layerx$Fill)>=pt){
        Bottom <- Layerx[(nrow(Layerx)-add.b+1):nrow(Layerx),]
        Middle <- Layerx[(Position(function(x) x >= pt,Layerx$Fill)+add.t):(nrow(Layerx)-add.b),]
        Top <- Layerx[1:(Position(function(x) x >= pt,Layerx$Fill)-1+add.t),]} else{
          Top <- Layerx[1:nrow(Layerx),]
          Middle <- NULL
          Bottom <- NULL
        }
      
      
      Temp$Split <- NA
      Temp[Temp$Layer %in% Middle$Layer,"Split"] <- "Middle"
      Temp[Temp$Layer %in% Bottom$Layer,"Split"] <- "Bottom"
      Temp[Temp$Layer %in% Top$Layer,"Split"] <- "Top"
      Temp[is.na(Temp$Split),"Split"] <- "Trimmed"
    }
    
    # If layer zero is the bottom of the specimen
    if(layer.start == "Bottom"){
      
      # Trimming
      if(trim==TRUE) Layerx <- Layer[Position(function(x) max(Layer$Count) == x,Layer$Count,right=TRUE):nrow(Layer),] else Layerx <- Layer
      
      # Calculate fill
      Layerx$Fill <- Layerx$Count / side^2
      
      # Split in layers
      if(max(Layerx$Fill)>=pt){
        Bottom <- Layerx[1:add.b,]
        Middle <- Layerx[(add.b+1):(Position(function(x) x >= pt,Layerx$Fill,right=TRUE)-add.t),]
        Top <- Layerx[(Position(function(x) x >= pt,Layerx$Fill,right=TRUE)+1-add.t):nrow(Layerx),]} else{
          Top <- Layerx[(add.b+1):nrow(Layerx),]
          Middle <- NULL
          Bottom <- Layerx[1:add.b,]
        }
      
      Temp$Split <- NA
      Temp[Temp$Layer %in% Middle$Layer,"Split"] <- "Middle"
      Temp[Temp$Layer %in% Bottom$Layer,"Split"] <- "Bottom"
      Temp[Temp$Layer %in% Top$Layer,"Split"] <- "Top"
      Temp[is.na(Temp$Split),"Split"] <- "Trimmed"
    }
    
    return(Temp)
  }
  return(LS)
}

