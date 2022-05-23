########################################################################################################
# ---------------------------------------------------------------------------------------------------- #
#                                                                                                      #
#        Crestline Analysis: Automated Calcium Wave Data Analysis (R Script)                           #
#                                                                                                      #
# ---------------------------------------------------------------------------------------------------- #
########################################################################################################
# ---------------------------------------------------------------------------------------------------- #
# (c) 2022 Martin Zauser, Rik Brugman, Janos Loeffler, Milan, Zupunski, Guido Grossman, Juergen Pahle  #
#     Biological Information Processing Group, BIOMS / BioQuant, Heidelberg University, Germany        #
#     Centre for Organismal Studies, Heidelberg University, Germany                                    #
# Developer: Martin Zauser                                                                             #
# Corresponding author: Juergen Pahle                                                                  #
# Adress: Im Neuenheimer Feld 267, 69100 Heidelberg                                                    #
# Email:  juergen.pahle@bioquant.uni-heidelberg.de                                                     #
# ---------------------------------------------------------------------------------------------------- #
########################################################################################################
# ---------------------------------------------------------------------------------------------------- #
# tested with: R version 4.2.0 (2022-04-22)                                                            #
#              Platform: x86_64-w64-mingw32                                                            #
#              OS:       Windows 7 Enterprise Service Pack 1                                           #
#              R Studio: 2022.02.0                                                                     #
# ---------------------------------------------------------------------------------------------------- #
########################################################################################################
# ---------------------------------------------------------------------------------------------------- #
#  Note: if you want to clean up the workspace before executing the R script you can use               #
#  the following R command for deleting the workspace (all variables/objects !!!):                     #
#  rm(list=ls(all=TRUE))                                                                               #
# ---------------------------------------------------------------------------------------------------- #
########################################################################################################

########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# |                                            COLORMAP                                              | #
# +--------------------------------------------------------------------------------------------------+ #
# | create a color map using nine colors within the following range:                                 | #
# |          dark blue => blue => green => orange => yellow                                          | #
# | use this function by providing the number of color steps for example: colormap()(64)             | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
colormap <- function() colorRampPalette(c("#352A86", "#0362E0", "#1483D4", "#05A5C7", "#33B8A0", "#8CBE74", "#D2BA58", "#FDCA30", "#F8FA0D"))


########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | show image (reduce size by factor 0.8 due to long text)                                          | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# | img: image matrix                                                                                | #
# | title: title of the image                                                                        | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
imgshow <- function(img, title="", sub="", col=colormap()(256), titlesize=0.8, rotate=TRUE, timeSpaceDiagram=FALSE, stimulationStartMinutes=0,
                    x=NA, y=NA, pixelx=NA, pixely=NA, xlab="distance [mm]", ylab="time [min]", filenameOutput=NA, showif=TRUE, ...) {
  # do nothing at all if "showif" is FALSE
  if (!showif) {
    return()
  }
  # check matrix size
  if (length(dim(img)) == 0) {
    stop("Image must be a matrix with number of rows > 1 and number of columns > 1.")
  }
  # get image dimensions (important because dimensions change after rotation)
  dimx <- dim(img)[2]
  dimy <- dim(img)[1]
  # use identical axis labels if the image is no timeSpaceDiagram
  if (timeSpaceDiagram == FALSE) {
    ylab <- xlab
    pixely <- pixelx
  }
  # calculate image dimensions (must be done *before* image will be rotated)
  if (!is.na(pixelx)) {
    x <- dim(img)[2] * pixelx
  }
  if (!is.na(pixely)) {
    y <- dim(img)[1] * pixely
  }
  # -------------------------------------------
  # rotate image
  if (rotate) {
    img <- t(img[nrow(img):1,,drop=FALSE])
  }
  # draw image
  graphics::image(z=img, col=col, axes=FALSE, sub=sub, useRaster=TRUE, xlab=xlab, ylab=ylab, ...)
  # add title
  if (title != "") {
    graphics::title(title, cex.main=titlesize)
  }
  # -------------------------------------------
  # add axes
  if (!is.na(x)) {
    n <- 10
    axis(side=1, at=(0:round(x*n))/round(x*n), labels=as.character((0:round(x*n))/n))
  }
  if (!is.na(y)) {
    if (timeSpaceDiagram == FALSE) { # space axis
      n <- 10
      axis(side=2, at=(0:round(y*n))/round(y*n), labels=as.character((0:round(y*n))/n))
    } else { # time axis
      if (stimulationStartMinutes == 0) {
        step <- 2
        pos <- 1 - (0:(round(y/step))/(y/step))
        val <- as.character((0:(round(y/step)))*step)
        graphics::axis(side=2, at=pos[which(pos >= 0)], labels=val[which(pos >= 0)]) # remove negative values
      } else {
        step <- 3
        offset <- 3
        pos <- 1 - (0:(round(y/step))/(y/step)) - (stimulationStartMinutes - offset)/y
        val <- as.character((0:(round(y/step)))*step - offset)
        graphics::axis(side=2, at=pos, labels=val) # remove negative values
      }
    }
  }
  # -------------------------------------------
  # save image to disk
  if (!is.na(filenameOutput)) {
    grDevices::jpeg(filename=filenameOutput, width=1280, height=960)
    # draw image
    graphics::image(z=img, col=col, axes=FALSE, sub=sub, useRaster=TRUE, ...)
    # add title, reduce size by factor 0.8 due to long text
    if (title != "") {
      graphics::title(title, cex.main=titlesize)
    }
    # close device
    grDevices::dev.off()
  }
}
