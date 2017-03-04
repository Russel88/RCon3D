#' Trim and standardize quantification 
#'
#' Function to trim and standardize quantification
#' @param qua A dataframe from the \code{quant} function
#' @param trim Should the bottom be trimmed? Specifies the bottom as the first layer (from the bottom) with max counts. Defaults to FALSE
#' @param layer.start Are first layers the "Top" or the "Bottom" of the specimen? If equal to "Top", the layers will be reversed so that low numbers denote the bottom of the specimen
#' @keywords array image split
#' @return The original qua dataframe including an Exp variable, with a unique name for each image (based on the naming given in the quantify function). And a NewLayer varaible denoting the layer from the bottom of the specimen after trimming 
#' @export

layer_stand <- function(qua,trim=FALSE,layer.start=NULL){
  
  # Make Exp variable in qua dataframe
  if(length(colnames(qua)) > 4){
    cols <- colnames(qua)[5:length(colnames(qua))]
    qua$Exp <- do.call(paste, c(qua[cols], sep="_"))
  } else qua$Exp <- 1
  
  # Layer 0 is Top
  if(layer.start=="Top"){
    # Trimming
    if(trim==TRUE){
      qua.layers <- foreach(i=unique(qua$Exp),.combine=rbind) %do% {
        qua.temp <- qua[qua$Exp == i,]
        if(length(unique(qua.temp$Channel))>1){
          Layer.sum <- aggregate(Count~Layer,FUN=sum,data=qua.temp)
          Layer.bottom <- as.numeric(Layer.sum[which.max(Layer.sum$Count),"Layer"])
          qua.temp2 <- qua.temp[as.numeric(qua.temp$Layer) <= Layer.bottom,]
        } else qua.temp2 <- qua.temp[1:which.max(qua.temp$Count),]
        qua.temp2$NewLayer <- rev(qua.temp2$Layer)
        return(qua.temp2)
      }} else {
        # No trimming
        qua.layers <- foreach(i=unique(qua$Exp),.combine=rbind) %do% {
          qua.temp <- qua[qua$Exp == i,]
          qua.temp$NewLayer <- rev(qua.temp$Layer)
          return(qua.temp)
        }
      }
  }
  
  # Layer 0 is Bottom
  if(layer.start=="Bottom"){
    # Trimming
    if(trim==TRUE){
      qua.layers <- foreach(i=unique(qua$Exp),.combine=rbind) %do% {
        qua.temp <- qua[qua$Exp == i,]
        if(length(unique(qua.temp$Channel))>1){
          Layer.sum <- aggregate(Count~Layer,FUN=sum,data=qua.temp)
          Layer.bottom <- as.numeric(Layer.sum[which.max(Layer.sum$Count),"Layer"])
          qua.temp2 <- qua.temp[as.numeric(qua.temp$Layer) >= Layer.bottom,]
        } else qua.temp2 <- qua.temp[which.max(qua.temp$Count):nrow(qua.temp),]
        qua.temp2$NewLayer <- qua.temp2$Layer-(min(qua.temp2$Layer)-1)
        return(qua.temp2)
      }} else {
        # No trimming
        qua.layers <- qua
      }
  }
  return(qua.layers)
}
