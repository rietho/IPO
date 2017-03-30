# IPO

[![Build Status](https://travis-ci.org/rietho/IPO.svg?branch=master)](https://travis-ci.org/rietho/IPO)
[![Years in Biocs](https://bioconductor.org/shields/years-in-bioc/IPO.svg)](https://bioconductor.org/packages/IPO/)

IPO (‘Isotopologue Parameter Optimization’) is a tool for automated Optimization of XCMS Parameters. It is a fast and free of labeling steps, and applicable to data from different kinds of samples and data from different methods of liquid chromatography - high resolution mass spectrometry and data from different instruments.

Find out more about the usage in the [package vignette](https://bioconductor.org/packages/devel/bioc/vignettes/IPO/inst/doc/IPO.html). The publication is available from http://www.biomedcentral.com/1471-2105/16/118.

## Installation

Get the release version from Bioconductor:
```{r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("IPO")
````

Or the development version from github:

```R
# install.packages("devtools")
devtools::install_github("rietho/IPO") 
```
