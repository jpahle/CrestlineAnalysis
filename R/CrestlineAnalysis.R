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
# |                                        LOADING FUNCTIONS                                         | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################

source("utils.R")

########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# |                                        LOADING LIBRARIES                                         | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
# load base packages
library(graphics)
library(grDevices)
library(utils)
library(stats)
library(datasets)
library(methods)
# load tiff package (package has to be installed before starting the script otherwise this command
#                    terminates with an error message)
library(tiff)
# optional libraries used only on demand
#library(cluster)
#library(ggplot2)
#library(plotly)
#library(R.matlab)


########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | linear correction of a vector                                                                    | #
# +--------------------------------------------------------------------------------------------------+ #
# | windowSize: window for calculating the first and the final value (NA => use windowSizeSeconds)   | #
# | periodSeconds: number of seconds for each time step                                              | #
# | windowSizeSeconds: window size (in seconds) for calculating the first and the final value        | #
# | onlyIfDownwards: TRUE: apply correction only if final value is smaller than starting value       | #
# | cutMinMAx: limit values if they exceed the original minimum/maximum after correction             | #
# | rescale: scale values to the range 0 to 1 after linear correction                                | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
linearCorrection <- function(x, windowSize=NA, windowSizeSeconds=60, periodSeconds=1, onlyIfDownwards=TRUE, cutMinMAx=FALSE, rescale=TRUE) {
  # calculate window size if not given
  if (is.na(windowSize)) {
    windowSize <- round(windowSizeSeconds / periodSeconds)
  }
  # get minimum and maximum
  minimum <- min(x, na.rm=TRUE)
  maximum <- max(x, na.rm=TRUE)
  # correct vector line linearly
  indexStart <- 1:windowSize
  indexEnd <- (length(x) - windowSize + 1):length(x)
  valueStart <- median(x[indexStart])
  valueEnd <- median(x[indexEnd])
  if (onlyIfDownwards && (valueStart <= valueEnd)) return(x)
  # apply linear correction
  for (i in 1:length(x)) {
    x[i] <- x[i] - (valueStart - valueEnd) / (length(x) - 1) * (length(x) - i)
  }
  # cut values which are below or above the original minimum/maximum
  if (cutMinMAx) {
    x[which(x < minimum)] <- minimum
    x[which(x > maximum)] <- maximum
  }
  # rescale vector
  if (rescale) {
    minimumNew <- min(x, na.rm=TRUE)
    maximumNew <- max(x, na.rm=TRUE)
    if (maximumNew - minimumNew > 0) { 
      x <- (x - minimumNew) / (maximumNew - minimumNew)
    }
    x[which(x < 0)] <- 0
    x[which(x > 1)] <- 1
  }  
  return (x)
}


########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | calculate separation point between two or more waves                                             | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# |   x: crestline vector                                                                            | #
# |   windowSizeSmedian: window size for smoothing                                                   | #
# |   rangemm: valid range for cutoff points in mm from the left and right border                    | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
cutoffPoint <- function(x, mmPerPixel, windowSizeSmedian=60, windowSizeSminmax=40, rangemm=c(0.2, 0.3), maxPoints=1,
                        minDistancePoints=200, minLimit=250, minLimitRelative=0.5, minSizeCutoffPoints=0.8) {
  # init variables
  res <- NULL
  if (maxPoints <= 0) {
    return (res)
  }
  # calculate range
  range <- c(round(rangemm[1] / mmPerPixel), length(x) - round(rangemm[2] / mmPerPixel))
  # smooth curve => calculate median
  xmedian <- smedian(x, windowSize=windowSizeSmedian, preserveSize=TRUE)
  # calculate minmax
  xminmax <- sminmax(xmedian, windowSize=windowSizeSminmax, preserveSize=TRUE)
  xminmaxmedian <- median(xminmax, na.rm=TRUE)
  # get maximum
  pos <- which.max(xminmax)[1]
  # check if maximum was found
  if (length(pos) > 0) {
    # get maximum
    m <- xminmax[pos]
    # check if maximum is within the range
    if ((pos >= range[1]) && (pos <= range[2])) {
      # check if maximum is large enough
      if ((m > minLimit) && (m >  (1 + minLimitRelative) * xminmaxmedian)) {
        res[length(res) + 1] <- pos
        while (TRUE) {
          if (length(res) >= maxPoints) {
            break
          }
          # erase maximum
          xminmax[limitLow(pos - minDistancePoints, 1):limitHigh(pos + minDistancePoints, length(xmm))] <- 0
          # search for additional cutoff points
          pos <- which.max(xminmax)[1]
          if (length(pos) > 0) {
            if ((pos >= range[1]) && (pos <= range[2])) {
              if (xminmax[pos] > minSizeCutoffPoints * m) {
                res[length(res) + 1] <- pos
              } else {
                break
              }
            } else {
              break
            }
          } else {
            break
          } 
        }
      }
    }
  }
  # return cutoff point(s) sorted by size
  return(sort(res))
}


########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | fit line data to a math function                                                                 | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# | line: line vector                                                                                | #
# | type: regression type of the fit ("linear" or "quadratic"), for "raw" line is only smoothed      | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
fitLine <- function(line, plot=FALSE, borderlevel=0.1, type="linear", mmPerPixel=NA, windowSizemm=0.05) {
  # fit line
  if (type == "raw") {
    # smooth line by applying a sliding median
    smedian(line, windowSize=round(windowSizemm / mmPerPixel), preserveSize=TRUE)
  } else {
    # calculate offset
    offsetLeft <- which(!is.na(line))[1] - 1
    offsetRight <-  which(!is.na(rev(line)))[1] - 1
    # calculate border
    border <- round((length(line) - offsetRight - offsetLeft) * borderlevel)
    # calculate vertex
    v <- round(median(which.min(line[(border + offsetLeft):(length(line) - border - offsetRight)]))) + border - 1 + offsetLeft
    # fit line
    regression(line, vertex=v, plot=plot, fixVertexX=TRUE, type=type)[[2]]
  }
}


########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | detect waves in the crestline and fit crestline data to a math function                          | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# | crestline: crestline vector                                                                      | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
fitCrestline <- function(crestline, ylim, mmPerPixel, maxNumberOfWaves=1, cutoffPointsmm=NULL, type="linear") {
  waveList <- list()
  # fit crestline
  if (maxNumberOfWaves <= 1) {
    cp <- NULL
  } else {
    # use given cutoff points(s)
    if (!is.null(cutoffPointsmm)) {
      cp <- round(cutoffPointsmm / mmPerPixel)
    } else {
      # calculate cutoff point(s) if not given
      cp <- cutoffPoint(crestline, mmPerPixel=mmPerPixel, maxPoints=(maxNumberOfWaves - 1))
    }
  }
  # fit lines
  for (i in 0:length(cp)) {
    section <- crestline
    if (i > 0) {
      section[1:(cp[i] - 1)] <- NA
    }
    if (i < length(cp)) {
      section[cp[i + 1]:length(section)] <- NA
    }
    line <- round(fitLine(section, mmPerPixel=mmPerPixel, windowSizemm=0.05, plot=FALSE, type=type))
    # cut borders
    line[which((line > ylim) | (line < 1))] <- NA
    waveList[[i + 1]] <- line
  }
  return (waveList)
}



########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | mark crestline in the crestline plot                                                             | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# | cre: crestline matrix                                                                            | #
# | line: crestline vector                                                                           | #
# | thickness: thickness of marked line in pixel                                                     | #
# | col: vector of intensities for above, center, and below the line                                 | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
markCrestline <- function(cre, line, thickness=5, col=c(1, 0, 1)) {
  ydim <- dim(cre)[1]
  for (i in 1:length(line)) {
    if (!is.na(line[i])) {
      if (!is.na(col[1])) cre[limitLow(line[i] - thickness - 1, 1):limitLow(line[i] - 1, 1), i] <- col[1]
      if (!is.na(col[2])) cre[line[i]:limitHigh(line[i] + thickness - 1, ydim), i] <- col[2]
      if (!is.na(col[3])) cre[limitHigh(line[i] + thickness, ydim):limitHigh(line[i] + 2 * thickness - 1, ydim), i] <- col[3]
    }
  }
  cre
}

########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | calculate crestline vector from crestline plot                                                   | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# |   crestlineImage: crestline plot, show: TRUE/FALSE => show profile plots and crestline           | #
# |   ignoreRows: ignore x rows from the beginning (i.e. ignore data before stimulation)             | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
calculateCrestline <- function(crestlineImage, show=FALSE, ignoreRows=0, ignoreColumns=0, removeOutliers=TRUE,
                               analyzeProfile=FALSE, level=0.5, windowSize=10) {
  # plot representative intensity profiles
  if (show) {
    xpos <- c(120, 580, 1000)
    color <- c("blue", "black", "green")
    for (i in 1:length(xpos)) {
      if (i == 1) {
        plot(crestlineImage[, xpos[i]], type="l", col="gray")
      } else {
        graphics::lines(crestlineImage[, xpos[i]], type="l", col="gray")
      }
      graphics::lines(smedian(crestlineImage[, xpos[i]], windowSize=10, preserveSize=TRUE),  type="l", col=color[i])
    }
  }
  # --------------------------------------------------------------------------------------------------- #
  xdim <- dim(crestlineImage)[2]
  ydim <- dim(crestlineImage)[1]
  maxline <- rep(NA, xdim)
  # --------------------------------------------------------------------------------------------------- #
  # algorithm I:
  # check for increase: find first local maximum following a continuous increase in the intensity profile
  if (analyzeProfile) {
    # adjust level
    level <- level * windowSize
    for (i in 1:xdim) {
      profile <- crestlineImage[(ignoreRows + 1):ydim,i]
      # calculate difference
      di <- profile[(windowSize + 1):(ydim - ignoreRows)] - profile[1:(ydim - ignoreRows - windowSize)]
      # find first local maximum following a continuous increase of at least 0.5 (parameter "level")
      sum <- 0
      for (j in 1:length(di)) {
        if (di[j] > 0) {
          sum <- sum + di[j]
        } else {
          if (sum > level * windowSize) {
            # maximum was found => store value and continue with next column
            jmax <- j + windowSize + 1
            if (jmax > ydim) {
              jmax <- ydim
            }
            maxline[i] <- ignoreRows + j - 1 + which.max(profile[j:jmax])[1]
            break
          }
          sum <- 0
        }
      }
    }
    if (show) {
      plot(maxline, type="l", col="black")
    }
    return(maxline)
  }
  # --------------------------------------------------------------------------------------------------- #
  # algorithm II:
  # calculate crestline as the maximum
  # determine start and end column
  # start column: 1 + number of columns to be ignored + remove columns for which all intensity values are zero
  # end column: remove columns for which all intensity values are zero
  colStart <- ignoreColumns + 1
  while ((colStart < xdim) && (sum(crestlineImage[,colStart]) == 0)) {
    colStart <- colStart + 1
  }
  colEnd <- xdim
  while ((colEnd > colStart) && (sum(crestlineImage[,colEnd]) == 0)) {
    colEnd <- colEnd - 1
  }
  # apply median to the maximum values if several maximum values exist
  for (i in colStart:colEnd) {
    maxline[i] <- median(which.max(crestlineImage[(ignoreRows + 1):ydim,i])) + ignoreRows
  }
  # plot crestline
  if (show) {
    plot(maxline, type="l", col="blue")
  }
  # remove outliers
  if (removeOutliers) {
    d <- sminmax(maxline)
    maxline2 <- maxline
    maxline2[which(d > 25)] <- NA
    maxline3 <- longestLine(maxline2)
  } else {
    maxline3 <- maxline
  }
  if (show) {
    graphics::lines(maxline3, type="l", col="green")
  }
  maxline4 <- smedian(maxline3, preserveSize=TRUE)
  if (show) {
    graphics::lines(maxline4, type="l", col="red")
  }
  return(maxline4)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate crestline vector from crestline plot    *** EXTENDED ***
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ccExt <- function(cre, show=FALSE, ignoreRows=0, level=0.5, windowSize=10, fillNA=FALSE, smooth=0) {
  xdim <- dim(cre)[2]
  ydim <- dim(cre)[1]
  maxline <- rep(NA, xdim)
  # adjust level
  level <- level * windowSize
  for (i in 1:xdim) {
    profile <- cre[(ignoreRows + 1):ydim,i]
    # calculate difference
    di <- profile[(windowSize + 1):(ydim - ignoreRows)] - profile[1:(ydim - ignoreRows - windowSize)]
    # find first local maximum following a continuous increase of at least 0.2 (parameter "level")
    found <- FALSE
    sum <- 0
    for (j in 1:length(di)) {
      if (di[j] > 0) {
        sum <- sum + di[j]
      } else {
        if (sum > level) {
          found <- TRUE
          break
        }
        sum <- 0
      }
    }
    if (found) {
      jmax <- j + windowSize + 1
      if (jmax > ydim) {
        jmax <- ydim
      }
      maxline[i] <- ignoreRows + j - 1 + which.max(profile[j:jmax])[1]
    }
  }
  # fill NAs between gaps
  if (fillNA) {
    maxline <- fillNA(maxline)
  }
  # smooth line
  if (smooth > 0) {
    maxline <- smedian(maxline, windowSize=smooth, preserveSize=TRUE)
  }
  if (show) {
    plot(maxline, type="l", col="black")
  }
  return(maxline)
}


########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | analyze crestline and calculate starting time and position and propagation speeds                | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# | line: crestline plot (matrix) or crestline (vector) - if plot is given calculate crestline first | #
# | correctLinearly: TRUE => apply linear correction to all columns of the crestline plot            | #
# | tipSizemm: tip size (in mm) = size of columns (in mm) to be ignored at the tip site              | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
crestlineAnalysis <- function(line, param=NULL, resolutionFactor=NA, stimulationStartMinutes=NA, periodSeconds=NA, mmPerPixel=NA, title="", 
                              consoleOutput=TRUE, imageOutput=FALSE, tipSizemm=0.1, maxNumberOfWaves=1, type="linear", cutoffPointsmm=NULL, correctLinearly=TRUE,
                              thickness=5, level=0.5, windowSize=10, useMinValue=FALSE, show=FALSE, analyzeProfile=FALSE, useFit=FALSE) {
  # check if mmPerPixelKymo or resolution factor is given
  if (((is.na(mmPerPixel)) && (is.na(resolutionFactor))) || ((is.null(param)) && (is.na(mmPerPixel)))) {
    stop("No resolution factor given. Please define resolution factor or mmPerPixel.")
  }
  # get parameters
  if (!is.null(param)) {
    filename <- param[[1]]
    periodSeconds <- param[[2]]
    stimulationStartMinutes <- param[[3]]
    # define spatial resolution
    if (!is.na(resolutionFactor)) {
      mmPerPixel <- param[[4]] * resolutionFactor
    }
    # use filename as title
    if (title == "") {
      title <- paste0("crestline plot","   -   ", filename)
    }
  }
  # init return value
  res <- list()
  # calculate crestline and return a list with data for each wave
  waveList <- list()
  if (class(line)[1] == "matrix") {
    # copy crestline into a separate variable named "cre" for crestline plot
    cre <- line
    # linear correction of all columns of the crestline plot
    if (correctLinearly) {
      for (i in 1:dim(cre)[2]) {
        cre[,i] <- linearCorrection(cre[,i], windowSizeSeconds=60, periodSeconds=periodSeconds)
      }
    }
    # calculate rows to be ignored before stimulation and columns to be ignored at tip site
    ignoreRows <- round(stimulationStartMinutes * 60 / periodSeconds)
    ignoreColumns <- round(tipSizemm / mmPerPixel)
    # calculate crestline
    line <- calculateCrestline(cre, ignoreRows=ignoreRows, ignoreColumns=ignoreColumns, show=show, removeOutliers=FALSE, analyzeProfile=analyzeProfile,
                               level=level, windowSize=windowSize)
    if (analyzeProfile) {
      waveList <- list(line)
    } else {
      # analyze distinct waves in the crestline plot
      waveList <- fitCrestline(line, ylim=dim(cre)[1], mmPerPixel=mmPerPixel, maxNumberOfWaves=maxNumberOfWaves, cutoffPointsmm=cutoffPointsmm, type=type)
    }
    # show line
    if (imageOutput) {
      cre <- markCrestline(cre, line, thickness=thickness, col=c(NA, 0, NA))
      for (i in 1:length(waveList)) {
        line <- waveList[[i]]
        cre <- markCrestline(cre, line, thickness=thickness, col=c(1, 0, 1))
      }
    }
  } else {
    imageOutput <- FALSE
  }
  # calculate speeds
  for (i in 1:length(waveList)) {
    line <- waveList[[i]]
    if (length(which(!is.na(line))) == 0) {
      starttime <- NA
    } else {
      starttime <- min(line, na.rm=TRUE)
      if (imageOutput) {
        cre[round(stimulationStartMinutes * 60) - 1,] <- 1
        cre[round(stimulationStartMinutes * 60)    ,] <- 1
        cre[round(stimulationStartMinutes * 60) + 1,] <- 1
        cre[(starttime - 1):(starttime + 1),] <- 1
      }
    }
    startpos <- median(which.min(line))
    if (useMinValue) {
      leftpos <- which.max(line[1:startpos])[1]
      rightpos <- which.max(line[startpos:length(line)])[1] + startpos - 1
    } else {
      leftpos <- firstIndex(line)
      rightpos <- lastIndex(line)
    }
    leftval <- line[leftpos]
    rightval <- line[rightpos]
    speedTip <- (startpos - leftpos) * mmPerPixel /  ((leftval - starttime) * periodSeconds)
    speedShoot <- (rightpos - startpos) * mmPerPixel /  ((rightval - starttime) * periodSeconds)
    startposmm <- startpos * mmPerPixel
    starttimeSeconds <- starttime * periodSeconds - (stimulationStartMinutes * 60)
    # store result
    
    currentResult <- c(starttimeSeconds / 60, startposmm * 1000, speedTip * 1000, speedShoot * 1000)
    res[[length(res) + 1]] <- currentResult
    # print result
    if (consoleOutput) {
      cat(paste0("start (time/pos): ", round(currentResult[1], 2), " min / ", round(currentResult[2]), " um", "   ",
                 "speed (tip/shoot): ", round(currentResult[3], 2), " um/s / " , round(currentResult[4], 2), " um/s", "\n"))
    }
    # mark result by lines
    if (imageOutput) {
      xdim <- dim(cre)[2]
      # define left and right column of the line
      left <- round((thickness - 1) / 2)
      right <- thickness - left - 1
      # add lines
      cre[, limitLow(startpos - left, 1):limitHigh(startpos + right, xdim)] <- 1
      cre[, limitLow(leftpos - left, 1):limitHigh(leftpos + right, xdim)] <- 1
      cre[, limitLow(rightpos - left, 1):limitHigh(rightpos + right, xdim)] <- 1
    }
  }
  # plot crestline plot and fitted lines
  if (imageOutput) {
    imgshow(cre, title=title, x=dim(cre)[2] * mmPerPixel, y=dim(cre)[1] * periodSeconds / 60, stimulationStartMinutes=stimulationStartMinutes, timeSpaceDiagram=TRUE)
  }
  # return result
  res
}


########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | plot kymograph or crestline plot (for plotting kymograph and crestline plot are treated equally) | #
# +--------------------------------------------------------------------------------------------------+ #
# |  img: time/space intensity matrix of the kymograph (rows: time, columns: space)                  | #
# |  kymograph parameters must be given EITHER by "param" and "resolutionFactor"                     | #
# |     OR by "stimulationStartMinutes", "periodSeconds" and "mmPerPixel"                            | #
# |  param: parameter set of the kymograph (includes period, stimulation start etc.)                 | #
# |  crestline: TRUE => img is a crestline matrix, FALSE => img is a kymograph matrix                | #
# |  show: TRUE => show kymograph, FALSE => do nothing                                               | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
plotKymograph <- function(img, param=NULL, resolutionFactor=NA, stimulationStartMinutes=NA,
                          periodSeconds=NA, mmPerPixel=NA, title="", show=TRUE, crestline=FALSE) {
  # show=FALSE => do nothing
  if (!show) {
    return
  }
  # check if mmPerPixelKymo or resolution factor is given
  if (((is.na(mmPerPixel)) && (is.na(resolutionFactor))) || ((is.null(param)) && (is.na(mmPerPixel)))) {
    stop("No resolution factor given. Please define resolution factor or mmPerPixel.")
  }
  # get parameters
  if (!is.null(param)) {
    filename <- param[[1]]
    periodSeconds <- param[[2]]
    stimulationStartMinutes <- param[[3]]
    # define spatial resolution
    if (!is.na(resolutionFactor)) {
      mmPerPixel <- param[[4]] * resolutionFactor
    }
    # use filename as title
    if (title == "") {
      title <- filename
    }
  }
  # add blanks and hyphen to title
  if (title != "") {
    title <- paste0("   -   ", title)
  }
  # define title
  if (crestline) {
    title <- paste("crestline plot", title)
  } else {
    title <- paste("kymograph", title)
  }
  # plot kymograph
  imgshow(img, title=title,
          x=dim(img)[2] * mmPerPixel, y=dim(img)[1] * periodSeconds / 60, stimulationStartMinutes=stimulationStartMinutes, timeSpaceDiagram=TRUE)
}



########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | create crestline plot                                                                            | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
plotCrestline <- function(...) {
  # create crestline plot
  plotKymograph(..., crestline=TRUE)
}

########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | show surface plot                                                                                | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# |   mmPerPixel: spatial resolution (for kymograph/crestline data: resolutionFactor is 0.5)         | #
# |               => divide by 2; for simulation data: analysis factor has default value 1)          | #
# |   digits=c(xdigits, ydigits): define digits for floating point rounding                          | #
# |   digits=c(xdigits, ydigits): define digits for floating point rounding                          | #
# |   ticksteps=c(xticksteps, yticksteps): define tick resolution                                    | #
# +--------------------------------------------------------------------------------------------------+ #
# | examples:                                                                                        | #
# | surfacePlot(chitin[[3]][[1]], "Kymograph for Dataset CH122", 1, 0.5*0.001/0.9091)                | #
# | surfacePlot(chitin[[3]][[2]], "Crestline Plot for Dataset CH122", 1, 0.5*0.001/0.9091)           | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
surfacePlot <- function(img, title="kymograph", periodSeconds=1, mmPerPixel=0.5*0.001/0.9091, plot=FALSE,
                        cmap=colormap(), stepsize=100, digits=c(2, 1), ticksteps=c(5, 10)) {
  # load plotly library (and the required library ggplot2)
  library(ggplot2)
  library(plotly)
  # get data dimensions
  ysize <- dim(img)[1]
  xsize <- dim(img)[2]
  # set NAs to the minimum value
  img <- replaceNAbyMin(img)
  # prepare image data
  # mirror matrix and scale data with factor (for example: divide all values by 1000 => zFactor = 0.001)
  dat <- apply(img, 2, rev)
  # show only each xth point
  if (stepsize > 1) {
    x <- round(seq(1, xsize, length=stepsize))
    y <- round(seq(1, ysize, length=stepsize))
    dat <- dat[y, x]
  }
  # calculate resolution
  if (stepsize > 1) {
    xres <- stepsize
    yres <- stepsize
  } else {
    xres <- xsize
    yres <- ysize
  }
  # calculate tickfactor = c(x, y, z)
  tickfactor <- c(xsize * mmPerPixel, ysize / 60 * periodSeconds, 1)
  # calculate ticks
  xtickvals <- seq(0, xres, xres / ticksteps[1])
  xticktext <- as.character(round(seq(0, 1, 1 / ticksteps[1]) * tickfactor[1], digits[1]))
  ytickvals <- seq(0, yres, yres / ticksteps[2])
  yticktext <- as.character(round(rev(seq(0, 1, 1 / ticksteps[2])) * tickfactor[2], digits[2]))   # apply reverse order to y axis => rev()
  # create axes and define plot parameters
  opacity <- 1
  fontTitle <- list(family="Arial,Raleway,sans-serif", size=14)
  fontAxis <- list(family="Arial,Raleway,sans-serif", size=10)
  fontColorbar <- list(family="Arial,Raleway,sans-serif", size=14)
  xaxis <- list(title="root axis [mm]", titlefont=fontAxis, showgrid=TRUE, tickmode="array", tickvals=xtickvals, ticktext=xticktext)
  yaxis <- list(title="time [min]", titlefont=fontAxis, showgrid=TRUE, tickmode="array", tickvals=ytickvals, ticktext=yticktext)
  zaxis <- list(title="intracellular Ca2+ [a.u.]", titlefont=fontAxis, showgrid=TRUE)
  colorbar <- list(title="[Ca2+]", titlefont=fontColorbar)
  cameraposition <- list(list(-0.35, 0.55, 0, -0.5), list(0, 0, -0.24), 2.8)
  scene <- list(xaxis=xaxis, yaxis=yaxis, zaxis=zaxis, cameraposition=cameraposition)
  # add frame rate and resolution to the title
  title <- paste0(title, "<br>", "(frame rate: ", round(1 / periodSeconds, 2), " Hz, ",
                  "resolution: data: ", xsize, "x", ysize, " view: ", xres, "x", yres, ")")
  # create surface plot (ignore warnings)
  p <- plotly::plot_ly(z=(dat * tickfactor[3]), type="surface", colors=cmap(256), color=0:1, opacity=opacity, colorbar=colorbar) %>%
    layout(title=list(text=title, font=fontTitle), scene=scene)
  # show plot if plot=TRUE
  if (plot) {
    print(p)
  }
  # return plotly htmlwidget
  return(p)
}


########################################################################################################
# +--------------------------------------------------------------------------------------------------+ #
# | create boxplots for publication                                                                  | #
# | (example: kymographs / crestline plots for "ATP" and "chitin")                                   | #
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+ #
# | kymographData: list of precalculated data (kymograph, crestline plot, parameters:                | #
# |                      period, start time of stimulation, spatial resolution of the kymograph)     | #
# | stimulus: name vector of the stimulus                                                            | #
# | ylimit: y-axis limits for starting time (sec), starting point (micrometer),                      | #
# |                           speed to the tip / shoot (micrometer/second)                           | #
# | type: fitting type of crestline analysis: "raw", "linear" or "quadratic"                         | #
# | allDatapoints: TRUE => plot all data points, FALSE => plot only outliers                         | #
# | pch, cex, lwd, col: if allDatapoints is TRUE: graphic parameters for data points                 | #
# +--------------------------------------------------------------------------------------------------+ #
########################################################################################################
showBoxplot <- function(kymographData, stimulus=c("ATP", "chitin"), cutoffPointsmm=NULL, consoleOutput=FALSE,  maxNumberWaves=2, type="quadratic",
                        ylimit=c(10, 1000, 30, 30), allDatapoints=TRUE, pch=1, cex=2, lwd=2, col="blue") {
  # define names
  caption <- c("Delay after stimulation", "Distance from tip", "Wave speed towards tip", "Wave speed towards shoot")
  ylab <- c("time [min]", "distance [µm]", "speed [µm/s]", "speed [µm/s]")
  # calculate waves
  waveData <- list()
  for (k in 1:length(kymographData)) {
    numberWaves <- 0
    waveData[[k]] <- list()
    for (i in 1:length(kymographData[[k]])) {
      periodSeconds <- kymographData[[k]][[i]][[4]][1]
      stimulationStartMinutes <- kymographData[[k]][[i]][[4]][2] / 60
      mmPerPixel <- kymographData[[k]][[i]][[4]][3]
      res <- crestlineAnalysis(kymographData[[k]][[i]][[2]], mmPerPixel=mmPerPixel, stimulationStartMinutes=stimulationStartMinutes, periodSeconds=periodSeconds,
                               type=type, analyzeProfile=FALSE, maxNumberOfWaves=maxNumberWaves, consoleOutput=consoleOutput, imageOutput=FALSE)
      # calculate maximum number of waves
      numberWaves <- max(numberWaves, length(res))
      if (length(waveData[[k]]) < numberWaves) {
        for (j in (length(waveData[[k]]) + 1):numberWaves) {
          waveData[[k]][[j]] <- matrix(NA, nrow=length(kymographData[[k]]), ncol=4)
          colnames(waveData[[k]][[j]]) <- caption
        }
      }
      for (j in 1:length(res)) {
        res[[j]][which(is.nan(res[[j]]))] <- NA
        res[[j]][which(res[[j]] == Inf)] <- NA
        waveData[[k]][[j]][i,] <- res[[j]]
      }
    }
  }
  # prepare data for boxplot
  listBoxplot <- list()
  for (i in 1:4) {
    count <- 0
    listBoxplot[[i]] <- list()
    for (k in 1:length(kymographData)) {
      for (j in 1:length(waveData[[k]])) {
        count <- count + 1
        listBoxplot[[i]][[count]] <- waveData[[k]][[j]][,i]
        if (length(waveData[[k]]) == 1) {
          names(listBoxplot[[i]]) <- stimulus
        } else {
          names(listBoxplot[[i]])[count] <- paste0(stimulus[k], " ", j)
        }
      }
    }
  }
  # create boxplot
  graphics::par(mfrow=c(1,4))
  for (i in 1:4) {
    # create subtitle
    subtitle <- "\n  "
    if (length(listBoxplot[[i]]) == 1) {
      subtitle <- paste0(stimulus, "\n  ")
    }
    for (k in 1:length(kymographData)) {
      for (j in 1:length(waveData[[k]])) {
        subtitle <- paste0(subtitle, ifel(j > 1, "      "), "(n = ", length(which(!is.na(waveData[[k]][[j]][,i]))),")  ")
      }
    }
    # show boxplot
    graphics::boxplot(listBoxplot[[i]], main=paste0(LETTERS[i], "\n", caption[i]), xlab=subtitle, ylab=ylab[i], ylim=c(0, ylimit[i]), outline=!allDatapoints)
    if(allDatapoints) {
      count <- 0
      for (k in 1:length(waveData)) {
        for (j in 1:length(waveData[[k]])) {
          count <- count + 1
          graphics::points(rep(count, length(waveData[[k]][[j]][,i])), waveData[[k]][[j]][,i], pch=pch, cex=cex, lwd=lwd, col=col)
        }
      }
    }
  }
  graphics::par(mfrow=c(1,1))
}
