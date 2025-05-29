<img width=400 src="https://github.com/predsci/OFT/assets/4073260/cc737381-833b-4b7b-b07a-b8e48c92cfdc" alt="OFT" />  
  
# OFT: Open Source Flux Transport  

[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  
  
## OVERVIEW  
Open-source Flux Transport (OFT) is a software suite for generating full-Sun magnetograms through acquiring & processing observational data, generating realistic convective flows, and running them through a surface flux transport (SFT) model.  

OFT includes the data acquisition/mapping code [MagMAP](https://github.com/predsci/magmap), the convective flow generation code [ConFlow](https://github.com/predsci/conflow), and the surface flux transport code [HipFT](https://github.com/predsci/hipft) as submodules.
  
--------------------------------  
  

## HOW TO BUILD OFT  

Check out the repository using git's `recursive` option:  

```
git clone https://github.com/predsci/oft --recursive
```
If forgotten, initialize the submodules with:  
```
git submodule update --init
```
When pulling updates, to ensure you get the submodule updates as well, use:  
```
git pull --recurse-submodules
```
  
Once you have the OFT repository including submodules, follow the installation instructions within each submodule folder.  

All executables/scripts need to be in your PATH for OFT to work.  
For BASH, we have provided a script that can be sourced to do this:
```
. ./load_oft_env.sh
```
The script can also be used as a reference for making an equivalent script for other shells.
  
--------------------------------  
  
## HOW TO RUN OFT  
  
First, make sure all OFT tools and scripts are in your PATH (see above).    

OFT has three main components:  
```
 - MagMAP:  Obtain, process, and map HMI magnetograms, and prepare them for data assimilation in HipFT  
 - ConFlow: Generate a sequence of convective flows for use in HipFT  
 - HipFT:   The flux transport model.  
```  
Depending on the use case, using MagMAP and/or ConFlow can be  optional (for example, when running HipFT without HMI data assimilation and/or without convective flows).  

See the run instructions in each submodule for further information.  

## OFTSWA for Space Weather Applications (OFTSWA)  

OFTSWA (`bin/oftswa.py`) is a flexible Python tool that links the various software components of OFT. It allows users to easily create ensembles of full-Sun magnetic field maps from start to finish, ready for ingestion into space weather models. A YAML-formatted text file is used to specify the requested details of the full OFT run including post-processing and re-binning of the maps. A default YAML file is included and used to allow user YAML files to be trimmed down for specific tasks.
  
A sample OFTSWA run including pre-computed MagMAP and Conflow datasets is provided in the zenodo package:  
  
[OFT Sample Run of Convective Flows and Data Assimilation](https://zenodo.org/records/14799033)
  
Documentation for OFTSWA is provided in [`doc/oftswa`](https://htmlpreview.github.io/?https://github.com/predsci/OFt/doc/oftswa/index.html).  
  
--------------------------------  

