# LESYMAP
Lesion to Symptom Mapping (R toolbox)  

*****  
#### Package details  
Version:  0.0.0.9001  
Systems:  Linux, Mac or [Windows Linux Subsystem](https://github.com/stnava/ANTsR/wiki/Installing-ANTsR-in-Windows-10-(along-with-FSL,-Rstudio,-Freesurfer,-etc).)  
Language: R (version 3.0 or above)  
License:  Apache License 2.0  
Author:   Dorian Pustina  
  
[![Travis Build Status](https://travis-ci.org/dorianps/LESYMAP.png?branch=master)](https://travis-ci.org/dorianps/LESYMAP)  
  
*****
## Install

The quickest way to install is:
```
install.packages('devtools') # if you don't have it yet.
devtools::install_github('dorianps/LESYMAP', dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))
```
This will install all the dependencies, including ANTsR (may take one hour).  
  
*****
## Use
```
library(LESYMAP)

# run an example analyses using data provided with lesymap
example(lesymap)

# All functions have appropriate documentation. Start by typing
?lesymap
```
 
### Note
Package is under development, the behavior of some functions may change.
