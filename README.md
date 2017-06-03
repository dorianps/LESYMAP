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
If this doesn't work, check out the more detailed [Installation Wiki](https://github.com/dorianps/LESYMAP/wiki/Installation).
  
*****
## Use
```r
library(LESYMAP)

# All functions have appropriate documentation. Start by typing
?lesymap

# run an example analyses using data provided with lesymap
example(lesymap)
```
```
18:34:35 Running LESYMAP 0.0.0.9001   
18:34:35 Checking a few things...  
18:34:35 Loading behavioral data...131 scores found.  
18:34:35 Computing mask from average >= 10% ...  
18:35:23 Computing unique patches...  
18:36:01 Found 195102 patches in 326828 voxels - 1.7 times more voxels  
18:36:01 Using existing lesion matrix... 131x195102  
18:36:01 Running analysis: BMfast ...  
18:36:04 Correcting p-values: fdr ...  
18:36:04 Preparing images...  
18:36:05 Logging call details...  
18:36:05 Done! 1.5 mins   
Hit \<Return\> to see next plot: 
```

*****    
### Note
Package is under development, the behavior of some functions may change.
