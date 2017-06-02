# LESYMAP
Lesion to Symptom Mapping (R toolbox)  

*****  
## Description
[![Travis Build Status](https://travis-ci.org/dorianps/LESYMAP.png?branch=master)](https://travis-ci.org/dorianps/LESYMAP)
Version: 0.0.0.9001
License: Apache License 2.0
Author: Dorian Pustina

*****  
##  Requirements:  
* Linux or Mac (since Oct 2016 ANTsR [can run in Windows 10](https://github.com/stnava/ANTsR/wiki/Installing-ANTsR-in-Windows-10-(along-with-FSL,-Rstudio,-Freesurfer,-etc).))  
* [R v3.0 or above](http://www.r-project.org/) 
* The [ANTsR](http://stnava.github.io/ANTsR/) package

*****
## Install

The quickest way to install is:
```
install.packages('devtools') # if you don't have it yet.
devtools::install_github('dorianps/LESYMAP')
```
The above code will install all the dependencies, including ANTsR. It may take a while.

*****
## Test
```
library(LESYMAP)

# run an example analyses using data provided with lesymap
example(lesymap)

# All functions have appropriate documentation. Start by typing
?lesymap
```
  

