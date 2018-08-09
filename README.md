# LESYMAP
Lesion to Symptom Mapping (R toolbox)  

*****  
#### Package details  
[![Open Source Love](https://badges.frapsoft.com/os/v3/open-source.svg?v=103)](https://github.com/dorianps/LESYMAP)  
Version:  0.0.0.9008  
Systems:  Linux, Mac or [Windows Linux Subsystem](https://github.com/stnava/ANTsR/wiki/Installing-ANTsR-in-Windows-10-(along-with-FSL,-Rstudio,-Freesurfer,-etc).)  
Language: R (version 3.0 or above)  
License:  Apache License 2.0  
  
[![Travis](https://img.shields.io/travis/dorianps/LESYMAP.svg?branch=master)](https://travis-ci.org/dorianps/LESYMAP)  

  
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
Hit <Return> to see next plot: 
```
Check out how this example looks on [the screen](https://youtu.be/0WQXEgip_zk).  
  
*****    
## Documentation
[What is lesion to symptom mapping](https://github.com/dorianps/LESYMAP/wiki/What-is-lesion-to-symptom-mapping)  
[Installing LESYMAP](https://github.com/dorianps/LESYMAP/wiki/Lesymap-Installation)  
[Registering lesions in template space](https://github.com/dorianps/LESYMAP/wiki/Registering-lesions-in-template-space)  
[Provided data](https://github.com/dorianps/LESYMAP/wiki/Data)  
[Understanding permutations](https://github.com/dorianps/LESYMAP/wiki/Understanding-permutations)  
[Fast lesion drawing in ITKsnap](https://www.youtube.com/watch?v=ZVmINdWk5R4)  
[All Videos](https://github.com/dorianps/LESYMAP/wiki/Videos)  
***[New]*** [SCCAN questions](https://github.com/dorianps/LESYMAP/wiki/SCCAN-questions)  
[Known limitations](https://github.com/dorianps/LESYMAP/wiki/Known-Limitations)

*****    
### Note
Package under development, the behavior of some functions may change.
