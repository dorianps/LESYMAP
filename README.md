# LESYMAP
Lesion to Symptom Mapping (R toolbox)  

*****  
#### Package details  
[![Open Source Love](https://badges.frapsoft.com/os/v3/open-source.svg?v=103)](https://github.com/dorianps/LESYMAP)  
Version:  0.0.0.9210  
Systems:  Linux, Mac or [Windows Linux Subsystem](https://github.com/stnava/ANTsR/wiki/Installing-ANTsR-in-Windows-10-(along-with-FSL,-Rstudio,-Freesurfer,-etc).)  
Language: R (version 3.0 or above)  
License:  Apache License 2.0  
  
[![Travis](https://travis-ci.org/dorianps/LESYMAP.svg?branch=master)](https://travis-ci.org/dorianps/LESYMAP)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1452007.svg)](https://doi.org/10.5281/zenodo.1452007)  

  
*****
## What does it do  
Takes lesion maps and cognitive performance scores from patients with stroke, and maps the brain areas responsible for the cognitive deficit.  

  
*****
## Install

The quickest way to install is:
```
if (! 'devtools' %in% installed.packages()) install.packages('devtools')
devtools::install_github('dorianps/LESYMAP')
```
This will install all the dependencies, including ANTsR (may take one hour on some computers). Here is the [video](https://youtu.be/HSK2txFvbMU) of the installation process. If it doesn't work, check out the more detailed [installation instructions](https://github.com/dorianps/LESYMAP/wiki/Lesymap-Installation).
  
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
13:46:03 Running LESYMAP 0.0.0.9201 
13:46:03 Checking a few things...
13:46:03 Loading behavioral data...131 scores found.
13:46:03 Filenames as input, checking lesion values on 1st image...
13:46:04 Searching voxels lesioned >= 10% subjects...326828 found
13:46:36 Computing unique patches...
13:46:55 Found 195102 patches in 326828 voxels - 1.7 times more voxels
13:46:55 Using existing lesion matrix... 131x195102
13:46:55 Running analysis: BMfast ...
13:46:57 Correcting p-values: fdr ...
13:46:57 Preparing images...
13:46:57 Logging call details...
13:46:57 Done! 54.1 secs 
Hit <Return> to see next plot: 
```
Check out how this example looks on [the screen](https://youtu.be/0WQXEgip_zk).  
  
*****    
## Documentation
[What's new in v0.0.0.9201](https://github.com/dorianps/LESYMAP/wiki/What's-new-in-v0.0.0.9201)  
[What is lesion to symptom mapping](https://github.com/dorianps/LESYMAP/wiki/What-is-lesion-to-symptom-mapping)  
[Installing LESYMAP](https://github.com/dorianps/LESYMAP/wiki/Lesymap-Installation)  
[Registering lesions in template space](https://github.com/dorianps/LESYMAP/wiki/Registering-lesions-in-template-space)  
[Provided data](https://github.com/dorianps/LESYMAP/wiki/Data)  
[Understanding permutations](https://github.com/dorianps/LESYMAP/wiki/Understanding-permutations)  
[Fast lesion drawing in ITKsnap](https://www.youtube.com/watch?v=ZVmINdWk5R4)  
[All Videos](https://github.com/dorianps/LESYMAP/wiki/Videos)  
[SCCAN questions](https://github.com/dorianps/LESYMAP/wiki/SCCAN-questions)  
[Known limitations](https://github.com/dorianps/LESYMAP/wiki/Known-Limitations)  
[LESYMAP function reference](https://github.com/dorianps/LESYMAP/raw/master/LESYMAP.pdf)  

*****    
### Note
Package under development, the behavior of some functions may change.
