This is a tutorial, showing examples of the various functions for analysis of confocal images.
==============================================================================================

The images has to be binary, and are assumed to have been thresholded
already The example image has four channels (named "xan","pan","ste" and
"mic").

First lets load the package and some packages for plotting

    library(RCon3D)
    library(ggplot2)
    library(reshape2)
    library(scatterplot3d)

Load the images. If the images have already been loaded we can use
findIMG to load in the images. The path should lead to folder with a
.tif for each image (with all z-stacks in one), or a folder with
subfolders in which the images is split in z-stacks and channels. An
internal function tiffToArray is partly borrowed from
github/rmnppt/iMage

    #myimg <- loadIMG("//a00143.science.domain/cmf483/Documents/PhD/Projects/Image Analysis/Data",c("xan","pan","ste","mic"),split=TRUE)
    myimg <- findIMG("//a00143.science.domain/cmf483/Documents/PhD/Projects/Image Analysis/Data")

Quantify pixels for each layer for the four channels. The naming
argument is optional but can be used to look through the names of the
images and add corresponding variables Here it looks for "24h" in the
image name, and makes a variable called Time. This is of course only
useful when there are several images with different metadata. (Eg.
Time=c("12h","24h"))

    myq <- Quantify(myimg,channels=c("xan","pan","mic","ste"),naming=list(Time=c("24h")))
    head(myq)

    ##                          Img Channel Count Layer Time
    ## 1 FourSpecies24h_xan_Array.R     xan  1583     1  24h
    ## 2 FourSpecies24h_xan_Array.R     xan  1985     2  24h
    ## 3 FourSpecies24h_xan_Array.R     xan  2225     3  24h
    ## 4 FourSpecies24h_xan_Array.R     xan  2542     4  24h
    ## 5 FourSpecies24h_xan_Array.R     xan  3012     5  24h
    ## 6 FourSpecies24h_xan_Array.R     xan  3508     6  24h

Plot quantification

    p <- ggplot(data=myq,aes(x=Layer,y=Count,colour=Channel,group=Channel)) +
      theme_classic() +
      geom_freqpoly(stat="identity",position=position_dodge(width = 0),size=1) +
      coord_flip()
    p

![](RCon3DTutorial_files/figure-markdown_strict/unnamed-chunk-4-1.png)

As the bottom of the specimen is in the high numbers of the layers, lets
reverse layers and plot again. Note that trim=TRUE. This is because we
think the layer with most fill is the actual bottom of the specimen, and
we therefore trim away all that is below that layer

    myq.std <- LayerStd(myq,layer.start = "Top",trim=TRUE)

    p <- ggplot(data=myq.std,aes(x=NewLayer,y=Count,colour=Channel,group=Channel)) +
      theme_classic() +
      geom_freqpoly(stat="identity",position=position_dodge(width = 0),size=1) +
      coord_flip()
    p

![](RCon3DTutorial_files/figure-markdown_strict/unnamed-chunk-5-1.png)

Lets calculate cross-correlation between channels "pan" and "xan". A
cross-correlation of 1 equals

First we find out how many microns we can scan. It has to be a multiple
of both zstep and pwidth

    pwidth <- 0.75
    zstep <- 0.25
    library(rootSolve)
    uniroot.all(function(x) x%%pwidth + x%%zstep,interval=c(0,100),n=1e6)

    ##   [1]  0.00  0.75  1.50  2.25  3.00  3.75  4.50  5.25  6.00  6.75  7.50
    ##  [12]  8.25  9.00  9.75 10.50 11.25 12.00 12.75 13.50 14.25 15.00 15.75
    ##  [23] 16.50 17.25 18.00 18.75 19.50 20.25 21.00 21.75 22.50 23.25 24.00
    ##  [34] 24.75 25.50 26.25 27.00 27.75 28.50 29.25 30.00 30.75 31.50 32.25
    ##  [45] 33.00 33.75 34.50 35.25 36.00 36.75 37.50 38.25 39.00 39.75 40.50
    ##  [56] 41.25 42.00 42.75 43.50 44.25 45.00 45.75 46.50 47.25 48.00 48.75
    ##  [67] 49.50 50.25 51.00 51.75 52.50 53.25 54.00 54.75 55.50 56.25 57.00
    ##  [78] 57.75 58.50 59.25 60.00 60.75 61.50 62.25 63.00 63.75 64.50 65.25
    ##  [89] 66.00 66.75 67.50 68.25 69.00 69.75 70.50 71.25 72.00 72.75 73.50
    ## [100] 74.25 75.00 75.75 76.50 77.25 78.00 78.75 79.50 80.25 81.00 81.75
    ## [111] 82.50 83.25 84.00 84.75 85.50 86.25 87.00 87.75 88.50 89.25 90.00
    ## [122] 90.75 91.50 92.25 93.00 93.75 94.50 95.25 96.00 96.75 97.50 98.25
    ## [133] 99.00 99.75

Ok. lets try 6 (microns) then

    #mycc <- CrossCorX(imgs=myimg,channels=c("xan","pan"),size=15,npixel=100,dstep=1,pwidth=0.75,zstep=0.25,R=5)

and plot the result

    #p <- ggplot(mycc,aes(x=Distance,y=CC,group=R)) +
    #  theme_classic() +
    #  geom_hline(yintercept=1) +
    #  geom_line() 
    #p

Lets define the Top from the top until 75 percent of the image is filled
with pixels. And the bottom as the bottom 20 layers

    myq.split <- LayerSplit(myq,side=512,pt=0.75,add.t=0,add.b=20,channel=NULL,trim=FALSE,layer.start="Top")

    p <- ggplot(data=myq.split,aes(x=rev(Layer),y=Count,colour=Channel,linetype=Split,group=interaction(Channel,Split))) +
      theme_classic() +
      geom_freqpoly(stat="identity",position=position_dodge(width = 0),size=1) +
      coord_flip() +
      scale_linetype_manual(values=c("dashed","solid","dotted"))
    p

![](RCon3DTutorial_files/figure-markdown_strict/unnamed-chunk-9-1.png)

Lets find aggregates of "mic"

    #my.agg <- Agg(myimg,"mic",pwidth=0.75,zstep=0.25)

Lets plot the 3D image of aggregates larger than 20000 pixels

    #M <- melt(my.agg[[2]][[1]])
    #M <- M[!is.na(M$value),]
    #tabl <- as.data.frame(table(M$value))
    #subtable <- tabl[tabl$Freq>20000,]
    #M <- M[M$value %in% subtable$Var1,]
    #scatterplot3d(M$Var1,M$Var2,M$Var3,color=M$value,angle=45,pch=19)
