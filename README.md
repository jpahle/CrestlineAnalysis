# CrestlineAnalysis
Analysing kymographs using the crestline conversion.

## Overview

Development Stage: **Beta**

The CrestlinePackage provides functions for the automated evaluation of image series of calcium waves in plant roots.
The core of the package is the Crestline algorithm, which compensates for the weakening intensity of the wave by scaling the kymograph
columnwise, making the wave more visible and machine-readable.
The CrestlinePackage contains a collection of functions for generating a kymograph from the image data,
creating the crestline plot and determining characteristic parameters of the wave.

CrestlineAnalyis features:

-   TODO
-   TODO



The CrestlineAnalysis package comes with the Artistic License 2.0. By using the package you agree to
this license.

## Installation

Install the CrestlineAnalysis package directly from GitHub:

``` r
install.packages("remotes")
remotes::install_github("jpahle/CrestlineAnalysis")
```

The installation may take a few minutes.


## Usage

``` r
library(CoRC)
analyzeCrestline(ATPdata, stimulationStartMinutes=3, periodSeconds=1, mmPerPixel=0.001, resolutionFactor=0.5, type="quadartic", maxNumberOfWaves=2, cutoffPointsmm=NULL, analyzeProfile=FALSE, imageOutput=TRUE)
```
