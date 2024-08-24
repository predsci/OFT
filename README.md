<img width=400 src="https://github.com/predsci/OFT/assets/4073260/cc737381-833b-4b7b-b07a-b8e48c92cfdc" alt="OFT" />  
  
# OFT: Open Source Flux Transport  

[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  
  
## OVERVIEW  
Open-source Flux Transport (OFT) is a software suite for generating full-Sun magnetograms through acquiring & processing observational data, generating realistic convective flows, and running them through a surface flux transport (SFT) model.  

OFT includes the data acquisition/mapping code [MagMAP](https;//github.com/predsci/magmap), the convective flow generation code [ConFlow](https;//github.com/predsci/conflow), and the surface flux transport code [HipFT](https;//github.com/predsci/hipft) as submodules.
  
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
  
Once you have the OFT repository with its submodules, follow the instructions in each submodule to install them.  Then, add the OFT `bin` folder to your PATH.  For example, with bash:  
```
export PATH=<LOCATION_OF_OFT>/bin:$PATH
```
  
--------------------------------  
  

## HOW TO RUN OFT  
  
Coming soon...
  
--------------------------------  

