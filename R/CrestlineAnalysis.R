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
                               type=type, usePeak=TRUE, maxNumberOfWaves=maxNumberWaves, consoleOutput=consoleOutput, imageOutput=FALSE)
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
