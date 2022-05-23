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
# |                                       GENERAL FUNCTIONS                                          | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################

# (simple) ifelse replacement/modification
#   similar to ifelse but the condition is checked only once (not for every element)
#   compare the different behavior: ifelse(TRUE, 1:5) versus ifel(TRUE, 1:5)
ifel <- function(condition, yes, no="") {
  if (condition) {
    return(yes)
  } else {
    return(no)
  }
}

# limit function (methods: "cut" and "ring")
limit <- function(value, lowerLimit=0, upperLimit=1, method="cut")
{
  if (method == "cut") {
    return(pmax(pmin(value, upperLimit), lowerLimit))
  } else {
    if ((value >= lowerLimit) && (value <= upperLimit)) {
      retrun(value)
    }
    if (value < lowerLimit) {
      value <- upperLimit - (lowerLimit - value) + 1
    }
    if (value > upperLimit) {
      value <- lowerLimit + (value - upperLimit) - 1
    }
    # apply cut method additionally (in case the distance is too large for overrun)
    return(pmax(pmin(value, upperLimit), lowerLimit))
  }
}
limitLow    <- function(value, lowerLimit) pmax(value, lowerLimit)
limitHigh   <- function(value, upperLimit) pmin(value, upperLimit)
limitRange  <- function(start, end, lowerLimit, upperLimit) limitLow(start, lowerLimit):limitHigh(end, upperLimit)
limitRangeOffset <- function(start, end, offset) limitLow(start + offset, start):limitHigh(end + offset, end)

# calculate maximum difference within a vector
maxdiff <- function(a) return(max(a) - min(a))

########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# |                                        SLIDING WINDOWS                                           | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
# .................................................................................................................................................................
# calculate mean value for every position using a sliding window
# input parameters: a: input data
#                   windowSize: size of the window
#                   preserveSize: if TRUE size of the resulting vector is identical to the size of the input data ("a")
#                                    in this case the values at the borders depend on the parameter "keepWindowSizeConstant" (see below)
#                                 if FALSE the size of the resulting vector is smaller or equal to the size of the input data depending on the size of the window
#                                    in this case "keepWindowSizeConstant" is ignored and assumed to be TRUE
#                   keepWindowSizeConstant: if TRUE window at the borders is shifted inwards (meaning some values close to the border are identical)
#                                           if FALSE window is cut at the border (meaning window size shrinks down to size 1 for the values close to the border)
# .................................................................................................................................................................
########################################################################################################
slidingFunction <- function(a, windowSize=10, preserveSize=FALSE, keepWindowSizeConstant=TRUE, circular=FALSE, func=mean) {
  # ignore empty input data
  if (length(a) == 0) {
    return(a)
  }
  # ignore negative values and too large values for the window size
  if (windowSize < 1) {
    windowSize <- 1
  }
  if (windowSize > length(a)) {
    windowSize <- length(a)
  }
  # get data length
  len <- length(a)
  # calculate offset size for the sliding window and number of values
  if (preserveSize == FALSE) {
    # calculate offset size for the sliding window
    offset.left  <- 0
    offset.right <- windowSize - 1
    n <- len - windowSize + 1  # number of data elements
  } else {
    # calculate offset size for the sliding window
    if ((windowSize %% 2) == 0) {
      offset.left  <- windowSize / 2 - 1
      offset.right <- windowSize / 2
    } else {
      offset.left  <- (windowSize - 1) / 2
      offset.right <- (windowSize - 1) / 2
    }
    n <- len  # number of data elements
  }
  res <- rep(0, n)         # initialize result vector
  if (circular) {
    for(i in 1:n) {
      pos.left  <- i - offset.left
      pos.right <- i + offset.right
      # calculate mean or median
      if (pos.left < 1) {
        res[i] <- func(c(a[(len + pos.left):len], a[1:pos.right]))
      } else {
        if (pos.right > len) {
          res[i] <- func(c(a[pos.left:len], a[1:(pos.right - len)]))
        } else {
          res[i] <- func(a[pos.left:pos.right])
        }
      }
    }
  } else {
    for(i in 1:n) {
      pos.left  <- i - offset.left
      pos.right <- i + offset.right
      if (pos.left < 1) {
        if (keepWindowSizeConstant) {
          pos.right <- pos.right - pos.left + 1
        }
        pos.left  <- 1
      }
      if (pos.right > len) {
        if (keepWindowSizeConstant) {
          pos.left  <- len - windowSize + 1
        }
        pos.right <- len
      }
      # calculate mean or median
      res[i] <- func(a[pos.left:pos.right])
    }
  }
  return(res)
}

smean <- function(a, windowSize=10, preserveSize=FALSE, keepWindowSizeConstant=TRUE, circular=FALSE) {
  slidingFunction(func=mean, a=a, windowSize=windowSize, preserveSize=preserveSize, keepWindowSizeConstant=keepWindowSizeConstant, circular=circular)
}
smedian <- function(a, windowSize=10, preserveSize=FALSE, keepWindowSizeConstant=TRUE, circular=FALSE) {
  slidingFunction(func=median, a=a, windowSize=windowSize, preserveSize=preserveSize, keepWindowSizeConstant=keepWindowSizeConstant, circular=circular)
}
sminmax <- function(a, windowSize=10, preserveSize=FALSE, keepWindowSizeConstant=TRUE, circular=FALSE)
  slidingFunction(func=maxdiff, a=a, windowSize=windowSize, preserveSize=preserveSize, keepWindowSizeConstant=keepWindowSizeConstant, circular=circular)


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

########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | linear or quadratic regression                                                                   | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# parameters:
#   x: if x and y are provided: the X-coordinates of your data (e.g. time points)
#      if only x is provided: the function values (Y-coordinates) of your data
#           - in this case x is assumed to be running from 1 to (number of y)
#   y: function values of your data (Y-coordinates) if x is provided separately
#   vertex (two values):   x and y coordinates of the vertex; the vertex has the coordinates P(vertex[1]]|vertex[2]])
#   vertex (single value): index of the vertex;               the vertex has the coordinates P(x[vertex]|y[vertex])
#           if vertex is NA the center point of the data is assumed to be the vertex
#   fixVertex: if TRUE X- and Y-coordinate are assumed to be fix at P(x[vertex]|y[vertex])
#              if FALSE only the X-coordinate of the vertex is assumed to be fix at x[vertex]
#   resolution: (only in charge if plot is TRUE) number of data points for the plot
#   plot: if TRUE plot provided data and predicted values of the model and print model
#         characteristics to the console
#   return value: a list containing the model and the predicted values
# =--------------------------------------------------------------------------------------------------= #  
# |  examples:                                                                                       | #
# |     regression(c(1:50, 50:1))                                                                    | #
# |     regression(c(1:50, 50:1), fixVertex=FALSE)                                                   | #
# =--------------------------------------------------------------------------------------------------= #  
########################################################################################################
regression <- function(x, y=NA, vertex=NA, fixVertexXY=TRUE, fixVertexX=TRUE, resolution=NA, plot=TRUE,
                       type="quadratic") {
  # use default sequence if only y values are provided
  if (is.na(y[1])) {
    y <- x
    x <- 1:length(x)
  }
  # if vertex is a vector of size 2 => x and y value of the vertex are given
  if (length(vertex) == 2) {
    xs <- vertex[1]
    ys <- vertex[2]
  } else {
    # calculate vertex if not given or outside of the curve: midpoint of the data
    if (is.na(vertex) || (vertex < 1) || (vertex > length(x))) {
      vertex <- round((length(x) + 1) / 2)
    }
    xs <- x[vertex]
    ys <- y[vertex]
  }
  # fit model
  #   convert errors into warnings by using "control=nls.control(warnOnly=TRUE)"
  #   and suppress these warnings by using "suppressWarnings"
  if (type == "quadratic")
    # parabola equation: y = a * x^2 + b * x + c
    #                or: y = a * (x - xs)^2 + ys (this is the form we use)
    # xs is always assumed to be fix
    # ys as either assumed to be fix or a free parameter (fixVertex=TRUE/FALSE)
    suppressWarnings(model <- stats::nls(y ~ ifelse(x > xs, abs(a1) * (x - xs)^2 + ys, abs(a2) * (x - xs)^2 + ys),
                                         data=data.frame(x=x, y=y),
                                         start=ifel(fixVertexXY, list(a1=1, a2=1), ifel(fixVertexX, list(a1=1, a2=1, ys=ys), list(a1=1, a2=1, ys=ys, xs=xs))),
                                         control=stats::nls.control(warnOnly=TRUE)))
  if (type == "linear")
    # parabola equation: y = a * x^2 + b * x + c
    #                or: y = a * (x - xs)^2 + ys (this is the form we use)
    # xs is always assumed to be fix
    # ys as either assumed to be fix or a free parameter (fixVertex=TRUE/FALSE)
    suppressWarnings(model <- stats::nls(y ~ ifelse(x > xs, a1 * (x - xs) + ys, a2 * (x - xs) + ys),
                                         data=data.frame(x=x, y=y),
                                         start=ifel(fixVertexXY, list(a1=1, a2=1), ifel(fixVertexX, list(a1=1, a2=1, ys=ys), list(a1=1, a2=1, ys=ys, xs=xs))),
                                         control=stats::nls.control(warnOnly=TRUE)))
  # plot model and values
  if (plot == TRUE) {
    # print model
    print(model)
    # calculate x values for plotting the prediction:
    #  x.pred is just a copy of x if no resolution is given
    #  otherwise the x values are calculated as a sequence of size "resolution"
    x.pred <- ifel(is.na(resolution), x, seq(x[1], x[length(x)], length.out=resolution))
    # calculate predicted values
    y.pred <- stats::predict(model, data.frame(x=x.pred))
    # plot prediction (red line) and data points (black circles)
    plot(x.pred, y.pred, type="l", col=2, ylim=c(min(y, na.rm=TRUE),max(y, na.rm=TRUE)))
    graphics::points(x,y)
  }
  # dummy return value for model and predicted data for the given x values
  list(model, stats::predict(model, data.frame(x=x)))
}

########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | calculate "longest line" of a vector:                                                            | #
# |  find longest sequence of non-NAs within a vector                                                | #
# |  => all shorter non-NA lines are replaced by NA                                                  | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# | a: vector including NA values                                                                    | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
longestLine <- function(a) {
  longestLine <- rep(NA, length(a))
  startpos <- 0
  startmax <- 0
  endmax <- 0
  maxlength <- 0
  found <- FALSE
  for(i in 1:length(a)) {
    if (found) {
      if (is.na(a[i])) {
        if (i - startpos > maxlength) {
          startmax <- startpos
          endmax <- i - 1
          maxlength <- endmax - startmax + 1
        }
        found <- FALSE
        startpos <- 0
      }
    } else {
      if (!is.na(a[i])) {
        startpos <- i
        found <- TRUE
      }
    }
  }
  if (found) {
    if (length(a) + 1 - startpos > maxlength) {
      startmax <- startpos
      endmax <- length(a)
    }
  }
  if (startmax > 0) {
    longestLine[startmax:endmax] <- a[startmax:endmax]
  }
  return(longestLine)
}

