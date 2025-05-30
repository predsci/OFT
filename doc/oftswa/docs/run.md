# OFTSWA Run Instructions  
  
OFTSWA is a high-level wrapper for the OFT model, specifically designed for space weather applications.  It can run all three major components of OFT, taking care of file names, folder structure, post processing, and map processing (remapping, smoothing, etc.)  
  
## Inputs

### Input YAML file

OFTSWA uses a text file written in a format called YAML to run the model.  A default YAML file is included in the OFT repository in `rsrc/default.yaml`.  This default YAML file contains all available input parameters, and is loaded during an OFTSWA run.  The user's custom YAML file is then loaded, and any parameters that differ from the default are updated, and the final combined YAML file is used for the run and saved in the output run folder.  The user's YAML file therefore only need specify parameters that differ from the default, allowing for streamlined input files.

For details on the format and available options in the input YAML file, see [here](yaml.md).

### Initial full-Sun magnetic map (optional)

When running flux transport with OFTSWA, a map file can be specified to be used as the initial full-Sun magnetic map in HipFT.  This map must match the grid and resolution of the requested HipFT run.  If no initial map is specified, a zero-value map will be used.

## Outputs

The outputs of an OFTSWA run is dependent on which model components were run and with what options set in the YAML file.  These include observational magnetic disk images and maps, convective flow fields, full-Sun magnetic maps (the raw maps from the flux transport model, processed versions of those maps, and/or processed versions on the grid style used in the WSA solar wind model), post processing quantities including run history diagnostics, images and movies of the maps and butterfly diagrams.  

## Folder Structure
  
OFTSWA uses a user specified main run directory to store and compute all OFTSWA runs.  
The location of the main run directory is set with the `-o` flag as shown below.  

For any OFTSWA run, a run folder is created within the main run directory.  The name of this run folder is set to the run name specified in the input YAML file.

For each component of OFTSWA being run or specified (MagMAP, ConFlow, and HipFT), a subfolder appears within the run folder.  For new MagMAP and/or ConFlow runs, `magmap` and/or `conflow` subfolders are created.  If previously run MagMAP and/or ConFlow runs are being used, soft links of the same name are created that point to the preexisting runs.  The HipFT component run is performed in a `hipft` subfolder.  

When a run is executed, the merged input YAML file is placed in the run folder.  If post processing is activated, the post processing directories are also placed in the run folder.  All standard output and standard error outputs from the various codes and tools are piped into the `oftswa.log` and `oftswa.err` files respectively.

An example folder structure with two runs, one that ran only ConFlow and MagMAP (`run1`), and one that ran HipFT using the original run (`run2`) is shown here:
```
oft_runs/

oft_runs/run1/
oft_runs/run1/conflow/
oft_runs/run1/magmap/
oft_runs/run1/magmap/magmap_data_disks/
oft_runs/run1/magmap/magmap_data_maps/
oft_runs/run1/run1.yaml
oft_runs/run1/oftswa.log
oft_runs/run1/oftswa.err
  
oft_runs/run2/
oft_runs/run2/conflow --> ../run1/conflow
oft_runs/run2/magmap --> ../run1/magmap
oft_runs/run2/hipft/
oft_runs/run2/hipft/input_map/
oft_runs/run2/hipft/output_maps/
oft_runs/run2/hipft/processed_maps/
oft_runs/run2/hipft/processed_maps_wsa/
oft_runs/run2/post_processing_raw/
oft_runs/run2/post_processing_processed/
oft_runs/run2/run2.yaml
oft_runs/run2/oftswa.log
oft_runs/run2/oftswa.err
```

## Loading OFT 
To load OFT into your terminal environment, enter the OFT code directory, and source the included loading script:
```
cd {OFT_DIR}  
. load_oft_env.sh
```
This assumes that all components of OFT have been properly installed (see the intallation guide) and that the `oftswa.py` has been copied into the `bin` folder.

## OFTSWA Run Command
Once the OFT code has been loaded, and the desired input YAML file is created, OFTSWA can then be run with the command:
```
oftswa.py -o <PATH-TO-MAIN-RUN-FOLDER> your_input_file.yaml
```

