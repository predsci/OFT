
# OFTSWA Input YAML File Details
  
OFTSWA uses a text file written in a format called YAML to run the model.  A default YAML file is included in the OFT repository in `rsrc/default.yaml`.  This default YAML file contains all available input parameters, and is loaded during an OFTSWA run.  The user's custom YAML file is then loaded, and any parameters that differ from the default are updated, and the final combined YAML file is used for the run and saved in the output run folder.  The user's YAML file therefore only need specify parameters that differ from the default, allowing for streamlined input files.
  
The YAML file is organized into the several sections.  The following is a detailed description of each available parameter, shown as they appear in the default YAML file.  
  
### General OFTSWA run information:

```
run:
  run_name: 'run1'
  preexisting_magmap_runname: ''
  preexisting_conflow_runname: ''
  date_start: '2022-01-01T00:00:00'
  date_end:   '2022-01-28T00:00:00'
  map_resolution:
    np: 1024
    nt: 512
  mpi_run_command: 'mpirun -bind-to none -np <NP>'
  mpi_num_procs: 1
  omp_num_threads: ''
```

 - `run_name`: This is the user's chosen run name.  It is independent of the name of the input YAML file, and will be used as the name of the run directory within the main OFTSWA output run folder.  It is also used to name the final merged YAML file stored in the run folder.

 - `preexisting_magmap_runname`: If not requesting a new MagMAP database of assimilation data, this specifies the run name of a previously existing MagMAP run.  This preexisting database is used with the OFTSWA run and/or is updated based on the parameter choices.

 - `preexisting_conflow_runname`: If not requesting a new ConFlow database of flow data, this specifies the run name of a previously existing ConFlow run.  This preexisting database is used with the OFTSWA run.

 - `date_start`: The UTC date/time for the start of the run.  This is specified in the format 'YYYY-MM-DDTHH:MM:SS'.  This date represents the start date for the HipFT run, and is also used as the initial requested date/time when creating a new MagMAP database.

 - `date_end`: The UTC date/time for the end of the run.  This is specified in the format 'YYYY-MM-DDTHH:MM:SS'.  This date represents the end date for the HipFT run, and is also used as the final requested date/time when creating a new MagMAP database.

 - `map_resolution`: The `np` and `nt` parameters in this entry determine or declare the desired resolution for the OFTSWA run.  Here, `np` is the number of longitudinal (phi) points, while `nt` is the number of latitudinal points in colatitude (theta).  For OFTSWA runs that are creating new data sets of ConFlow and/or MagMAP, this resolution will be used.  For runs using preexisting ConFlow/MagMAP data sets, the specified resolution must match that used to create the preexisting data sets.  Also, if one is using an initial magnetic map file for the HipFT run, the resolution between that file and that specified here must match.

 - `mpi_run_command`: The command to use when launching MPI codes (currently, only HipFT uses MPI).  The special tag `<NP>` is replaced by the number of MPI processes when the code is launched.  

 - `mpi_num_procs`: The number of MPI processes to use.  When running on GPUs, this should be set to the number of GPUs per compute node being used.  

 - `omp_num_threads`:  The number of OpenMP threads to use for codes and tools that are multi-threaded on the CPU.  This includes HipFT when running on CPUs, in which case it specifies the number of CPU threads per MPI process.

<img width=150 src="./conflow.png" alt="ConFlow" />

### ConFlow Convective Flow Options
```
conflow:
  run:  False
  overwrite:  False
```

 - `run`: Toggle to run ConFlow (`True`) or not (`False`). 
 - `overwrite`: Toggle to overwrite a previous existing ConFlow run within the run folder.  If set to `True`, the previous ConFlow run is deleted before the new one is run.  If set to `False`, if a previous run exists in the run folder, ConFlow will not run. 

<img width=300 src="./magmap.png" alt="MagMAP" />

### MagMAP Data Assimilation Data Options
 
```
magmap:
  run:  False
  overwrite:  False
  update:  False
  map_cadence_hr: 1
  map_source: 'hmi_los' 
```
 - `run`: Toggle to run MagMAP (`True`) or not (`False`). 
 - `overwrite`: Toggle to overwrite a previous existing MagMAP run within the run folder.  If set to `True`, the previous MagMAP run is deleted before the new one is run.  If set to `False`, if a previous run exists in the run folder, MagMAP will not run, unless `update` is set to `True` (see below).  
 - `update`: If set to `True`, an update command will be applied to a preexisting MagMAP run.  This will fill in any missing data in the date range that may not have been available during the original launch of the preexisting MagMAP run.
 - `map_cadence_hr`: Set the map cadence for data acquisition and mapping in units of hours.
 - `map_source`: Set the observational source of the MagMAP data to acquire and map.  Currently, only `hmi_los` is available.

<img width=300 src="./hipft.png" alt="HipFT" />

### HipFT Flux Transport Options  
```
hipft:
  run:  False
  overwrite:  False
  output_map_cadence_hr: 24
  initial_map_filename: ''
  initial_map_mult_fac: 1.0
  data_assimilation_mult_fac: 1.0
  random_flux_lifetime: 0.3
  realization_combination_mode: 'sweep1d'
  realization_parameters:
    diffusion_coefs: [300.0]
    flow_attenuation_values: [500.0]
    data_assimilation_mu_limits: [0.1]
    data_assimilation_mu_powers: [4.0]
    data_assimilation_lat_limits: [0.0]
    random_flux_amounts: [150.0]
    flow_meridional_p1: [22.0]
    flow_meridional_p3: [11.0]
    flow_meridional_p5: [-28.0]
    flow_differential_p0: [46.0]
    flow_differential_p2: [-262.0]
    flow_differential_p4: [-379.0]
```
 - `run`: Toggle to run HipFT (`True`) or not (`False`)
 - `overwrite`: Toggle to overwrite a previous existing HipFT run within the run folder.  If set to `True`, the previous HipFT run is deleted before the new one is run.  If set to `False`, if a previous run exists in the run folder, HipFT will not run. 
 - `output_map_cadence_hr`: Set the output map cadence for HipFT in units of hours. This can differ from the cadence of the MagMAP and ConFlow runs.
 - `initial_map_filename`: Set the full path to the chosen input full-Sun magnetic map.  The file must have the HipFT grid format and the same resolution as the current HipFT run.  If no map is specified, the run will start with a zero-value map.
 - `initial_map_mult_fac`:  Factor to multiply the initial map in HipFT.  
 - `data_assimilation_mult_fac`: Factor to multiply the MagMAP assimilated data before ingestion into HipFT.
 - `random_flux_lifetime`: Time span (in hours) between successive random flux generation.
#### Realization Parameter Specification:

A key feature of OFTSWA is specifying a range of realization parameters to run.  
Three modes of specifying the parameters are implemented:  
 1) `sweep1d`:  
For each list of parameters, span the range given while keeping all other parameters at default values.  For example, if one sets `par1=[a,b,c]` and `par2=[d,e,f]`, and assuming the default.yaml defines `par1=[a]` and `par2=[d]`, there would be 5 realizations: `(a,d),(b,d),(c,d),(a,e),(a,f)`.  The defaults can be changed by using a double list syntax, for example, `par1=[[g][a,b,c,g]]`.  This would set the default to `g` for all other parameters' 1D sweeps.  

 2) `crossall`:   
All possible combinations of the lists of parameters are combined and used.  Any parameter not specified uses the default for all realizations, and any parameter that is set to a single value, resets the default.  For example, if one sets `par1=[a,b,c]`, `par2=[d,e,f]`, and `par3=[g]`, there will be 9 realizations: `(a,d,g),(a,e,g),(a,f,g),(b,d,g),(b,e,g),(b,f,g),(c,d,g),(c,e,g),(c,f,g)`.  

 3) `manual`:  
This mode requires the user to manually specify the full list for each parameter one wants to vary.  Any parameter not specified uses the default for all realizations, and any parameter that is set to a single value, resets the default.  For example, to specify the realizations shown in the `sweep1d` example above, one would need to set: `par1=[a,b,c,a,a]` and `par2=[d,d,d,e,f]`.
  
 - `realization_combination_mode`: Set the combination mode described above.  Options are `sweep1d`, `crossall`, or `manual`.
 - `realization_parameters`: For each of the following, they can be specified as single lists (`[a,b,c]`) or double lists (`[[g][a,b,c]]`).  See above for the description of what each means.  For more details on what each parameter does in the model, see [OFT Paper I](https://arxiv.org/pdf/2501.06377).
     - `diffusion_coefs`: Diffusion coefficient in units of km^2/s.
     - `flow_attenuation_values`:  Field strength in Gauss to saturate the flow velocities.
     - `data_assimilation_mu_limits`: Data cutoff in mu=cos(theta) for assimilated data.
     - `data_assimilation_mu_powers`: Power for scaling weights for assimilated data.
     - `data_assimilation_lat_limits`: Data cutoff in degrees latitude for assimilated data.
     - `random_flux_amounts`: Amplitude of random flux emergence in units of 10^21 Mx/hr.
     - `flow_meridional_p1`: Factor for cos(theta) term in the analytic meridional flow profile in m/s.
     - `flow_meridional_p3`: Factor for cos^3(theta) term in the analytic meridional flow profile in m/s.
     - `flow_meridional_p5`: Factor for cos^5(theta) term in the analytic meridional flow profile in m/s.
     - `flow_differential_p0`: Factor for scalar term in the analytic differential rotation flow profile in m/s.
     - `flow_differential_p2`: Factor for cos^2(theta) term in the analytic differential rotation flow profile in m/s.
     - `flow_differential_p4`: Factor for cos^4(theta) term in the analytic differential rotation flow profile in m/s.

### Map Processing Options
```
output_map_processing:
  run: False
  overwrite: False
  map_resolution_delta_deg: 1.0
  smoothing_factor: 1.0
  map_multiplier: 1.0
  make_wsa_dataset: False
```  
 - `run`: Toggle to run map processing of the HipFT output (`True`) or not (`False`)
 - `overwrite`:  Toggle to overwrite a previous existing map processing.  If set to `True`, the previous map processing results are deleted before new ones are computed.  If set to `False`, if previous map processing results exist in the run folder, no map processing will be computed.
 - `map_resolution_delta_deg`:  Re-bin the maps to the chosen resolution using a flux-preserving method.  The resolution is specified by the delta-angle in degrees (e.g. setting to 1.0 will produce 361x181 maps).
 - `smoothing_factor`: If set, maps will be smoothed after any re-binning using the  surface diffusion advance in HipFT.  This factor determines the level of smoothing.
 - `map_multiplier`: If set, maps will be multiplied by this scalar.
 - `make_wsa_dataset`:  Toggle (`True`/`False`) to produce an additional set of maps with the same resolution as the processed maps, but on the grid format used in the WSA model (e.g. a 361x181 processed map becomes a 360x180 WSA input map).
  
### Post Processing Options
```
hipft_post_process:
  run:  True
  overwrite:  False
  history_plot_samples: 500
  history_plot_samples_markers: 20
```
 - `run`: Toggle to run post processing of the HipFT output and processed maps (`True`) or not (`False`)
 - `overwrite`: Toggle to overwrite a previous existing post processing within the run folder.  If set to `True`, the previous post processing results are deleted before new ones are computed.  If set to `False`, if previous post processing results exist in the run folder, no post processing will be computed.
 - `history_plot_samples`: Set the number of data points to use for the lines in the HipFT history plots.
 - `history_plot_samples_markers`: Set the number of data points to use for the markers in the HipFT history plots.



