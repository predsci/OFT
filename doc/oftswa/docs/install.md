# OFTSWA Installation Instructions
# <font color="#0000FF">Linux Ubuntu-based (20.04/22.04) distributions</font>
### Install the required Ubuntu packages:  
  
```
sudo apt install bash make build-essential gfortran tar openmpi-bin libopenmpi-dev git libhdf5-dev libhdf5-openmpi-dev python3-pip hdf5-tools ffmpeg htop libblas-dev liblapack-dev qtbase5-dev xterm x11-apps
```
### Install Python Packages:  
```
pip install numpy matplotlib h5py scipy zenodo_get pandas astropy sunpy
```
### Create and Enter an OFT installation folder:
```
cd <PATH-TO-INSTALL-TO>
mkdir oft
mkdir oft/code
cd oft/code
```
### Clone the git repository:  
```
git clone --recursive https://github.com/predsci/OFT.git .
```
### Install OFTSWA:  
The main OFTSWA script called `oftswa.py` provided with this package is not yet part 
of the github repository.  
With the username and password given to you by PSI, you can download the script here: 
```
https://predsci.com/downloads/ccmc_oftswa_sbir_phase1/oftswa.py
```
After downloading the script, copy it into the OFT bin folder:
```
cp oftswa.py <PATH-TO-INSTALL-TO>/oft/code/bin/  
```
### Enter the HipFT directory and build HipFT:  
```
cd hipft
./build_examples/build_gcc_ubuntu20.04_cpu.sh
```
### Test the HipFT installation:  
```
cd testsuite
./run_test_suite.sh -mpicall='mpirun -bind-to socket -np' -hipftexe=${PWD}/../bin/hipft
```
### Return to the `oft` folder:  
```
cd ../../
```
### Enter the ConFlow directory and build ConFlow:  
```
cd conflow
./build_examples/build_gcc_ubuntu20.04_cpu.sh
```
### Return to the `oft` folder:  
```
cd ../
```
### Enter the MagMAP directory and install the MagMAP python package:  
```
cd magmap
pip install -e ${PWD}
```
### Return to the `oft` folder:  
```
cd ../
```
### Source the startup script with the command:  
```
. load_oft_env.sh
```
## Troubleshooting
If there are any issues in the above installation, or if you are installing on a different platform and/or on GPUs, see the documentation in the hipft, conflow, and magmap submodules for more information/instructions.
