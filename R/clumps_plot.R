#' Plot 3D aggregates
#'
#' Provide pwidth and zstep to get true aspect ratio. Provide center and radius to draw a circular outline
#' @param clumps.out Output from \code{clumps} 
#' @param replica Image replica to plot
#' @param col Colour by aggregate ID ("agg") or by distance to a point ("dist")
#' @param pwidth Pixel width in microns
#' @param zstep Microns between the layers
#' @param center 3D coordinates (c(x,y,z)) of where to center circular outline
#' @param radius Radius of the circular outline, in microns if pwidth and zstep are provided, else in pixels
#' @param thresh Clumps smaller than this are not plotted (in pixels)
#' @param thresm.m If TRUE will apply threshold to sizes in microns instead of pixels
#' @keywords array image aggregate
#' @export
clumps_plot <- function(clumps.out, replica = 1, col = "agg", pwidth = NULL, zstep = NULL, center = NULL, radius = NULL, thresh = 0, thresh.m = FALSE){
  
  # Find positions
  M <- reshape2::melt(clumps.out[[length(clumps.out)]][[replica]])
  M <- M[!is.na(M$value),]

  # Remove by threshold
  subsize <- clumps.out[[1]][clumps.out[[1]]$Img == unique(clumps.out[[1]]$Img)[1],]
  if(thresh.m){
    if("Size.micron" %in% colnames(subsize)){
      rem <- as.character(subsize[subsize$Size.micron < thresh,"ID"])
    } else {
      stop("Clump size in microns not found. Either set thresh.m = FALSE or give 'clumps' both the zstep and pwidth argument")
    }
  } else {
    rem <- as.character(subsize[subsize$Size < thresh,"ID"])
  }
  M <- M[!M$value %in% rem,]
  
  img.dims <- dim(clumps.out[[length(clumps.out)]][[replica]])
  
  # Colour by aggregate
  if(col == "agg"){
    if(!is.null(pwidth) & !is.null(zstep)){
      rgl::plot3d(x = M$Var1*pwidth, y = M$Var2*pwidth, z = M$Var3*zstep, col = M$value, 
             aspect = c(img.dims[1]*pwidth,img.dims[2]*pwidth,img.dims[3]*zstep),
             xlab = expression(paste("x (",mu,m,")")), 
             ylab = expression(paste("y (",mu,m,")")), 
             zlab = expression(paste("z (",mu,m,")")))
    } else {
      rgl::plot3d(x = M$Var1, y = M$Var2, z = M$Var3, col = M$value,
             xlab = "x", ylab = "y", zlab = "z",
             aspect = c(img.dims[1],img.dims[2],img.dims[3]))
    }
  }
  
  # Colour by distance
  if(col == "dist"){
    
    # Colour function
    myColorRamp <- function(colors, values) {
      v <- (values - min(values))/diff(range(values))
      x <- colorRamp(colors)(v)
      rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
    }
    
    # Merge pixel positions with distance data.frame
    dist.df <- merge(M, clumps.out$Aggregates[clumps.out$Aggregates$Img == unique(clumps.out$Aggregates$Img)[replica],], 
                     by.x = "value", by.y = "ID")
    dist.df <- dist.df[!is.na(dist.df$dist),]
    
    if(length(unique(M$value)) != length(clumps.out$Aggregates$ID)){
      warning("Not all clumps plotted: The clumps below the threshold are omitted as there is no distance information in them")
    }
    
    # Colours. Red is close to center, Blue is away from center
    cols <- myColorRamp(c("red", "blue"), dist.df$dist) 
    
    # Plot
    if(!is.null(pwidth) & !is.null(zstep)){
      rgl::plot3d(x = dist.df$Var1*pwidth, y = dist.df$Var2*pwidth, z = dist.df$Var3*zstep, col = cols, 
             aspect = c(img.dims[1]*pwidth,img.dims[2]*pwidth,img.dims[3]*zstep),
             xlab = expression(paste("x (",mu,m,")")), 
             ylab = expression(paste("y (",mu,m,")")), 
             zlab = expression(paste("z (",mu,m,")")))
    } else {
      rgl::plot3d(x = dist.df$Var1, y = dist.df$Var2, z = dist.df$Var3, col = cols, 
             xlab = "x", ylab = "y", zlab = "z",
             aspect = c(img.dims[1],img.dims[2],img.dims[3]))
      
    }
  }
  
  # Draw outline
  if(!is.null(pwidth) & !is.null(zstep)){
    if(!is.null(radius) & !is.null(center)){
      rgl::plot3d(x = cos(1:img.dims[1])*radius+center[1]*pwidth,
             y = sin(1:img.dims[2])*radius+center[2]*pwidth,
             z = rep(1,img.dims[1]), add = TRUE,
             aspect = c(img.dims[1]*pwidth,img.dims[2]*pwidth,img.dims[3]*zstep))
    }
  } else {
    if(!is.null(radius) & !is.null(center)){
      rgl::plot3d(x = cos(1:img.dims[1])*radius+center[1],
             y = sin(1:img.dims[2])*radius+center[2],
             z = rep(1,img.dims[1]), add = TRUE,
             aspect = c(img.dims[1],img.dims[2],img.dims[3]))
    }
  }

}

