#!/usr/bin/env python3
########################################################################
# ****** oftswa.py: Open-source Flux Transport model for 
#                   Space Weather Applications
#
#     Predictive Science Inc.
#     www.predsci.com
#     San Diego, California, USA 92121
########################################################################
# Copyright 2025 Predictive Science Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
########################################################################
#
#  Version 1.5.0 :  Added flux balancing options.
#  Version 1.4.4 :  Changed output to avoid reporting an error when no error was there.
#  Version 1.4.3 :  Fixed OMP_NUM_THREADS reporting.
#  Version 1.4.2 :  Now if smoothing is set to 0, no smoothing is done.
#  Version 1.4.1 :  Error checking and yaml output not overwritten (indexed).
#  Version 1.4.0 :  Added butterfly plot options.
#  Version 1.3.6 :  Initial version.
#
########################################################################
import signal
import sys
import argparse
import os
from pathlib import Path
from datetime import datetime
import pytz
import yaml
import shutil
import subprocess
import pandas as pd
import glob
import numpy as np
import h5py as h5
import re

version = '1.5.0'

def signal_handler(signal, frame):
        print('You pressed Ctrl+C! Stopping!')
        sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)


def argParsing():
    parser = argparse.ArgumentParser(description='OFTSWA: Run the OFT model for Space Weather Applications.')

    parser.add_argument('input_yaml_file',
                        help='Input yaml file.')

    parser.add_argument('-o',
                        help='Directory to store runs. (default: $PWD/oftswa_runs).',
                        dest='outdir',
                        required=False)
    
    parser.add_argument('-dry_run',
                      help='Dry run mode',
                      dest='dry_run',
                      required=False,
                      default=False,
                      action='store_true')    
    
    parser.add_argument('-overwrite',
                      dest='overwrite',
                      default=False,
                      action='store_true',
                      help=argparse.SUPPRESS,
                      required=False)  
    
    return parser.parse_args()


def check_oft_loaded(args):
  """
    Checks if the OFT environment is loaded by verifying that the necessary
    directories are present in the system's PATH. 
  """
  # Retrieve the system's PATH environment variable
  path_env = os.environ.get("PATH", "")
  path_dirs = path_env.split(os.pathsep)

  # Iterate through directories in PATH to check for an OFT installation  
  for directory in path_dirs:
    # Identify directories ending with "/bin" (indicating an executable location)  
    if directory.endswith("/bin"):
      # Extract the base OFT installation directory 
      oft_dir = directory[:-4]

      # Define required subdirectories for a valid OFT installation 
      required_subdirs = [
          oft_dir,
          os.path.join(oft_dir, "magmap", "bin"),
          os.path.join(oft_dir, "conflow", "bin"),
          os.path.join(oft_dir, "hipft", "bin"),
      ]

      # Verify that all required subdirectories are present in PATH
      if all(subdir in path_dirs for subdir in required_subdirs):
        print(f'    ==> Found OFT installation: {oft_dir}')
        # Store the OFT installation path in args and return
        args.oft_dir = oft_dir
        return

  # Exit with an error message if OFT is not properly loaded  
  sys.exit("ERROR: Please load the OFT environment by sourcing 'load_oft_env.sh'")


def represent_pair_of_lists(dumper, value):
  """
    Ensures that YAML output represents a pair of lists in the compact '[[], []]'
    format instead of the default block style.
  """
  return dumper.represent_sequence('tag:yaml.org,2002:seq', value, flow_style=True)


def read_run_yaml(args):
  """
    Reads and merges the default and user-provided YAML configuration files.
  """
  # Define the path to the default YAML file
  default_yaml = os.path.join(args.oft_dir, 'rsrc','default.yaml')
  print(f'    ==> Reading default yaml file from OFT rsrc directory.')

  # Load the default YAML configuration
  with open(default_yaml, 'r') as stream:
    default_data = yaml.safe_load(stream)
  
  # Load the user-provided YAML configuration
  print(f'    ==> Reading user input yaml file.')
  with open(args.input_yaml_file, 'r') as stream:
    input_data = yaml.safe_load(stream)

  # Check if 'sweep1d' mode is enabled in HipFT settings
  is_sweep1d = input_data.get('hipft', {}).get('realization_combination_mode') == 'sweep1d'

  # Merge the default and user-provided YAML configurations
  print(f'    ==> Merging default and user yaml files.')
  deep_update(default_data, input_data, is_sweep1d)

  # Replace empty string values with None to standardize missing data representation
  return replace_empty_strings_with_none(default_data)


def deep_update(dict1, dict2, is_sweep1d):
  """
    Recursively updates `dict1` (default YAML) with values from `dict2` (input YAML).
  """
  # Iterate through key-value pairs in the input YAML dictionary 
  for key, value in dict2.items():
    # If both values are dictionaries, merge recursively 
    if isinstance(value, dict) and key in dict1 and isinstance(dict1[key], dict):
      # Special handling for 'realization_parameters' 
      if key == 'realization_parameters':
        combine_update(dict1[key], value, is_sweep1d)
      else:
        deep_update(dict1[key], value, is_sweep1d)
    else:
      # Overwrite the value in the default YAML with the input YAML value
      dict1[key] = value


def combine_update(dict1, dict2, is_sweep1d):
  """
    Performs a special update for the 'realization_parameters' section of the YAML configuration.
  """
  # Iterate through key-value pairs in the realization_parameters dictionary 
  for key, value in dict2.items():
    # Check if the value is a list of numbers and 'sweep1d' mode is enabled
    if all(isinstance(x, (int, float)) for x in value) and is_sweep1d:
      # Combine into `[[], []]` format
      dict1[key] = [dict1[key],value]
    # If not in 'sweep1d' mode or value is already `[[], []]`, use input YAML data
    elif len(value) == 2 or not is_sweep1d:
      dict1[key] = value
    else:
      # Raise an error for invalid input
      check_err(1,f'Bad input format for {key}')


def replace_empty_strings_with_none(data):
  """
    Recursively replaces empty strings ("") with None in a nested dictionary or list.
  """
  # If the input is a dictionary, process each key-value pair recursively
  if isinstance(data, dict):
    return {k: replace_empty_strings_with_none(v) for k, v in data.items()}
  # If the input is a list, process each element recursively
  elif isinstance(data, list):
      return [replace_empty_strings_with_none(v) for v in data]
  # If the value is an empty string, replace it with None
  elif data == "":
    return None
  # Otherwise, return the value unchanged
  else:
    return data


def conflow_submodule(args, conflow_params, run_params):
  """
    Manages running conflow.
  """
  # Define the path for the 'conflow' directory
  conflow_dir = os.path.join(args.orun, 'conflow')
  # Check if conflow directory exists or if overwrite is requested
  if not os.path.isdir(conflow_dir) or conflow_params.get('overwrite') is True:
    # Handle overwrite logic
    if conflow_params.get('overwrite') is True:
      print(f'    ==> Overwriting existing conflow run folder...')
      # Remove symbolic link if it exists
      if os.path.islink(conflow_dir):
          os.unlink(conflow_dir)
      # Remove existing directory if it exists
      elif os.path.exists(conflow_dir):
          shutil.rmtree(conflow_dir)

    # Create the conflow directory and change to it
    os.makedirs(conflow_dir, exist_ok=True)
    print(f'    ==> Conflow output folder: {conflow_dir}')
    os.chdir(conflow_dir)

    # Copy conflow.dat from the OFT directory
    conflow_dat = os.path.join(args.oft_dir, "rsrc", "conflow.dat")
    shutil.copy(conflow_dat, conflow_dir)
 
    # Update conflow.dat with map resolution parameters
    map_resolution_params = run_params.get('map_resolution')
    if map_resolution_params:
      if map_resolution_params.get('np'):
        sed('n_long', map_resolution_params.get('np'), 'conflow.dat')
      if map_resolution_params.get('nt'):
        sed('n_lat', map_resolution_params.get('nt'), 'conflow.dat')
    
    # Print output to terminal
    print(f'    ==> Copied {conflow_dat} to {conflow_dir}')
    print(f'    ==> OMP_NUM_THREADS set to {os.getenv("OMP_NUM_THREADS")}')
    print(f'    ==> ACC_NUM_CORES set to {os.getenv("ACC_NUM_CORES")}')
    print(f'    ==> Running conflow...')

    # Run the conflow simulation if not in dry-run mode
    if not args.dry_run:
      command = 'conflow'
      ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, f'Failed: {command}')
      args.logout.flush()
      args.logerr.flush()
    
    # Return to the original directory
    os.chdir(args.orun)
  else:
    # Inform the user if a previous conflow run folder exists
    print(f'    ==> Previous conflow run folder found.')


def rdhdf_2d_dims(h5_filename):
  """
    Retrieves the NT and NP dimensions from a given 2D HDF5 file.
  """
  # Open the H5 file in read mode
  with h5.File(h5_filename, 'r') as h5file:
    # Access the 'Data' dataset from the file
    data = h5file['Data']
    # Get the number of dimensions of the dataset
    ndims = data.ndim
    # Check if the dimensions exist and extract NT and NP values
    NT = len(data.dims[0][0]) if ndims > 0 and len(data.dims[0].keys()) != 0 else 0
    NP = len(data.dims[1][0]) if ndims > 1 and len(data.dims[1].keys()) != 0 else 0
  # Return the NT and NP dimensions
  return NT, NP



def preexist_magmap_submodule(args, preexisting_magmap_runname, run_params):
  """
    Handles the setup and validation for a preexisting MagMAP run. 
  """
  first_obs_datetime = None

  # Get the absolute path of the preexisting MagMAP directory
  preexisting_magmap = os.path.join(args.outdir, preexisting_magmap_runname,'magmap')
  preexisting_magmap = os.path.abspath(preexisting_magmap)
  print(f'    ==> Preexisting magmap specified: {preexisting_magmap}')

  # Check if the specified MagMAP directory exists
  if os.path.isdir(preexisting_magmap):
    # Define the path where the symbolic link to the MagMAP directory will be created
    magmap_dir = os.path.join(args.orun, 'magmap')

    # If the symbolic link doesn't already exist, create one pointing to the preexisting MagMAP
    if not os.path.islink(magmap_dir):
      # Get the relative path from the current directory to the preexisting MagMAP directory
      relative_path = os.path.relpath(preexisting_magmap, args.orun)
      # Create the symbolic link
      os.symlink(relative_path, magmap_dir)
      print('    ==> Symbolic link to MagMAP created.')
    else:
      print('    ==> Symbolic link to MagMAP run already exists.')

    # Get the date range in the 'magmap_maps_all.csv
    csv_filepath = os.path.join(magmap_dir, 'magmap_data_maps/index_files', 'magmap_maps_all.csv')
    data = pd.read_csv(csv_filepath)
    first_obs_datetime = data['obs_datetime_utc'].iloc[0]
    last_obs_datetime = data['obs_datetime_utc'].iloc[-1]

    # Check the magmap start date to the input yaml start date
    period_start = datetime.strptime(run_params.get("date_start"), '%Y-%m-%dT%H:%M:%S')
    period_start_utc = pytz.utc.localize(period_start)
    magmap_start = datetime.strptime(first_obs_datetime, '%Y-%m-%dT%H:%M:%S')
    magmap_start_utc = pytz.utc.localize(magmap_start)
    if magmap_start_utc > period_start_utc:
      print('    ==> NOTE: Specified start date is before the first available MagMAP map:')
      print(f'        ==> Specified start date:  {period_start_utc}')
      print(f'        ==> MagMAP first map date: {magmap_start_utc}')

    # Check the magmap end date to the input yaml end date
    period_end = datetime.strptime(run_params.get("date_end"), '%Y-%m-%dT%H:%M:%S')
    period_end_utc = pytz.utc.localize(period_end)
    magmap_end = datetime.strptime(last_obs_datetime, '%Y-%m-%dT%H:%M:%S')
    magmap_end_utc = pytz.utc.localize(magmap_end)
    if magmap_end_utc < period_end_utc:
      print('    ==> NOTE: Specified end date is past the final available MagMAP map.')
      print(f'        ==> Specified end date:    {period_end_utc}')
      print(f'        ==> MagMAP final map date: {magmap_end_utc}')

    # Check the dimensions (NP, NT) of the map data from the first file in the preexisting MagMAP run
    file = os.path.join(magmap_dir,'magmap_data_maps',data['map_path'].iloc[0])
    NP, NT = rdhdf_2d_dims(file)
    # Get the map resolution parameters from the input run parameters
    map_resolution_params = run_params.get('map_resolution')
    # Compare the dimensions (NP, NT) of the preexisting MagMAP run with the input parameters
    if NP != map_resolution_params.get('np'):
      check_err(1, f'NP of {map_resolution_params.get("np")} differs from the preexisting magmap run: {NP}')
    if NT != map_resolution_params.get('nt'):
      check_err(1, f'NT of {map_resolution_params.get("nt")} differs from the preexisting magmap run: {NT}')
    
    # Return the first observation datetime from the preexisting MagMAP run
    return first_obs_datetime
  else:
    # If the MagMAP directory doesn't exist, raise an error
    check_err(1, '==> ERROR: Could not find the preexisting MagMAP run.')


def magmap_submodule(args, magmap_params, run_params):
  """
    Manages running magmap.
  """
  first_obs_datetime = None
  # Define the MagMAP directory path
  magmap_dir = os.path.join(args.orun, 'magmap')

  # Check if MagMAP directory exists or if overwrite or update flags are set
  if not os.path.isdir(magmap_dir) or magmap_params.get('overwrite') is True or magmap_params.get('update') is True:
    # Check for conflicts between overwrite and update flags
    if magmap_params.get('overwrite') is True and magmap_params.get('update') is True:
      check_err(1, 'Cannot update and overwite magmap at the same time!')
    # Handle overwriting the existing MagMAP directory
    if magmap_params.get('overwrite') is True:
      print(f'    ==> Overwriting existing magmap run folder')
      # Remove symbolic link if it exists
      if os.path.islink(magmap_dir):
          os.unlink(magmap_dir)
      # Remove existing directory if it exists
      elif os.path.exists(magmap_dir):
          shutil.rmtree(magmap_dir)
    # Handle updating an existing MagMAP directory
    if magmap_params.get('update') is True:
      print(f'    ==> Updating existing magmap run folder')

    # Create the MagMAP directory and change to it
    os.makedirs(magmap_dir, exist_ok=True)
    print(f'    ==> Created magmap output folder: {magmap_dir}')
    os.chdir(magmap_dir)

    # MagMAP data disks location
    magmap_data_disks = os.path.join(magmap_dir, 'magmap_data_disks')

    # This is a secret option for advance users...
    # If a raw data directory is provided, create a symbolic link to it
    if magmap_params.get('raw_data_dir'):
      print(f'    ==> Existing download directory provided for magmap_get_data.py :')
      existing_magmap_data_disks = os.path.abspath(magmap_params.get('raw_data_dir'))
      print(f'\t {existing_magmap_data_disks}')
      magmap_data_disks = os.path.join(magmap_dir, 'magmap_data_disks')
      # Create symbolic link to existing raw data directory
      if not os.path.islink(magmap_data_disks):
        relative_path = os.path.relpath(existing_magmap_data_disks, magmap_dir)
        os.symlink(relative_path, magmap_data_disks)
        print('    ==> Creating symbolic link to :')
      else:
        print('    ==> Symbolic link already exists :')
      print('    ==> Creating symbolic link to :')
      print(f'\t {magmap_data_disks}')
      # Skip fetching data if raw data directory is provided
      print(f'    ==> Skipping magmap_get_data.py') 
    # If no raw data directory is provided, fetch data using magmap_get_data.py
    else:
      print(f'    ==> Running magmap_get_data.py')
      # Construct the command to run the magmap_get_data.py
      command = 'magmap_get_data.py'
      print(f'    ==> Download directory for magmap_get_data.py is :')
      print(f'\t {magmap_data_disks}')
      # Output directory
      command += f' -odir {magmap_data_disks}'
      # Map cadence parameter
      if magmap_params.get('map_cadence_hr'):
        command += f' -cadence {magmap_params.get("map_cadence_hr")}'

      # Ensure the start and end dates are provided in the run parameters
      if not run_params.get('date_start'):
        check_err(1, 'No date_start under Run Parameters')
      if not run_params.get('date_end'):
        check_err(1, 'No date_end under Run Parameters')
      # Add start and end dates to the command
      command += f' {run_params.get("date_start")} {run_params.get("date_end")}'

      # If not a dry run, execute the command to fetch the data
      if not args.dry_run:
        ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, f'Failed: {command}')
        args.logout.flush()
        args.logerr.flush()

    # Run magmap_disk2map.py to convert the downloaded data into magnetogram maps
    # Construct the command to run the magmap_disk2map.py
    command = f'magmap_disk2map.py {magmap_data_disks}'
    print(f'    ==> Running magmap_disk2map.py')
    print(f'    ==> Magnetogram disk images data directory is :')
    # MagMAP data maps location
    magmap_data_maps = os.path.join(magmap_dir, 'magmap_data_maps')
    print(f'\t {magmap_data_maps}')
    # Output directory
    command += f' -odir {magmap_data_maps}'

    # Ensure the start and end dates are provided in the run parameters
    if not run_params.get('date_start'):
      check_err(1, 'No date_start under Run Parameters')
    if not run_params.get('date_end'):
      check_err(1, 'No date_end under Run Parameters')
    # Add start and end dates to the command
    command += f' -startdate {run_params.get("date_start")} -enddate {run_params.get("date_end")}'
    # Add the map resolution parameters from the input run parameters
    map_resolution_params = run_params.get('map_resolution')
    if map_resolution_params:
      print(f'    ==> Magnetogram disk images will have size:  {map_resolution_params.get("np")} x {map_resolution_params.get("nt")}')
      command += f' -npout {map_resolution_params.get("np")}'
      command += f' -ntout {map_resolution_params.get("nt")}'
    
    # Add the map_source
#   map_source = magmap_params.get('map_source'):
#   if (map_source == 'hmi_los' or map_source == 'hmi_vec'):
#      command += f' -obs {map_source}'
#   else:
#      print(f'ERROR!  MagMAP source ({map_source}) invalid!  Valid options: hmi_los, hmi_vec.'
    
    # Execute the command if not a dry run
    if not args.dry_run:
      ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, f'Failed: {command}')
      args.logout.flush()
      args.logerr.flush()

    # Return to the original directory
    os.chdir(args.orun)
  else:
    print(f'    ==> Previous magmap run folder exits!')

  # Get the date range in the 'magmap_maps_all.csv
  csv_filepath = os.path.join(magmap_dir, 'magmap_data_maps/index_files', 'magmap_maps_all.csv')
  data = pd.read_csv(csv_filepath)
  first_obs_datetime = data['obs_datetime_utc'].iloc[0]
  last_obs_datetime = data['obs_datetime_utc'].iloc[-1]

  # Check the magmap start date to the input yaml start date
  period_start = datetime.strptime(run_params.get("date_start"), '%Y-%m-%dT%H:%M:%S')
  period_start_utc = pytz.utc.localize(period_start)
  magmap_start = datetime.strptime(first_obs_datetime, '%Y-%m-%dT%H:%M:%S')
  magmap_start_utc = pytz.utc.localize(magmap_start)
  if magmap_start_utc > period_start_utc:
    print('    ==> NOTE: Specified start date is before the first available MagMAP map:')
    print(f'        ==> Specified start date:  {period_start_utc}')
    print(f'        ==> MagMAP first map date: {magmap_start_utc}')

  # Check the magmap end date to the input yaml end date
  period_end = datetime.strptime(run_params.get("date_end"), '%Y-%m-%dT%H:%M:%S')
  period_end_utc = pytz.utc.localize(period_end)
  magmap_end = datetime.strptime(last_obs_datetime, '%Y-%m-%dT%H:%M:%S')
  magmap_end_utc = pytz.utc.localize(magmap_end)
  if magmap_end_utc < period_end_utc:
    print('    ==> NOTE: Specified end date is past the final available MagMAP map.')
    print(f'        ==> Specified end date:    {period_end_utc}')
    print(f'        ==> MagMAP final map date: {magmap_end_utc}')

  # Return the first observation datetime from the preexisting MagMAP run
  return first_obs_datetime


def hipft_submodule(args, hipft_params, run_params):
  """
    Manages running HipFT. 
  """
  # Define the HipFT directory path
  hipft_dir = os.path.join(args.orun, 'hipft')

  # Check if HipFT directory exists or if overwrite is requested
  if not os.path.isdir(hipft_dir) or hipft_params.get('overwrite') is True:
    # Handle overwrite logic
    if hipft_params.get('overwrite') is True:
      print(f'    ==> Overwriting existing HipFT run folder')
      try:
        shutil.rmtree(hipft_dir)
      except FileNotFoundError:
          pass
      
    # Create the HipFT directory and change to it
    os.makedirs(hipft_dir, exist_ok=True)
    print(f'    ==> Created HipFT output folder')
    os.chdir(hipft_dir)

    # Copy hipft.in from the OFT directory
    hipft_in_new = os.path.join(hipft_dir, "hipft.in")
    hipft_in = os.path.join(args.oft_dir, "rsrc", "hipft.in")
    shutil.copy(hipft_in, hipft_in_new)  
    print(f'    ==> Copied HipFT template input file to run directory')

    # Update hipft.in with output map cadence
    if hipft_params.get('output_map_cadence_hr'):
      sed('output_map_time_cadence', hipft_params.get('output_map_cadence_hr'), 'hipft.in')
    
    # Update hipft.in with output flux balance
    if hipft_params.get('output_map_flux_balance') is not None:
      sed('output_map_flux_balance', hipft_params.get('output_map_flux_balance'), 'hipft.in')      
    
    # Check if initial_map_filename specified
    if hipft_params.get('initial_map_filename'):
      # Make the input_map folder
      input_map_dir = os.path.join(hipft_dir, 'input_map')
      os.makedirs(input_map_dir, exist_ok=True)
      # Copy initial_map_filename into the input_map folder
      input_map = Path(hipft_params.get('initial_map_filename')).resolve()
      command = f'cp {input_map} {input_map_dir}'
      ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, command)
      args.logout.flush()
      args.logerr.flush()
      # Check if file is in PT
      check_pt(input_map, run_params)

      # Create symbolic link from  the input_map_dir to the hipft run folder
      from_link = os.path.join(input_map_dir, input_map.name)
      to_link = os.path.join(input_map_dir, 'br_input_map.h5')
      if not os.path.islink(to_link):
        relative_path = os.path.relpath(from_link, input_map_dir)
        os.symlink(relative_path, to_link)
        print(f'    ==> Creating symbolic link to {input_map_dir}')
      else:
        print(f'    ==> Symbolic link already exists: {input_map_dir}')
      sed('initial_map_filename', "'input_map/br_input_map.h5'", 'hipft.in')
      print(f'    ==> Initial map found: {input_map}')
    else:
      print('    ==> WARNING: Initial map not specified, a zero-value map will be used')
      sed('initial_map_filename', '', 'hipft.in')

    # Update hipft.in with initial_map_mult_fac
    if hipft_params.get('initial_map_mult_fac'):
      sed('initial_map_mult_fac', hipft_params.get('initial_map_mult_fac'), 'hipft.in')

    # Update hipft.in with random_flux_lifetime
    if hipft_params.get('random_flux_lifetime'):
      sed('source_rfe_lifetime', hipft_params.get('random_flux_lifetime'), 'hipft.in')
    map_resolution_params = run_params.get('map_resolution')

    # Update hipft.in with data_assimilation_mult_fac
    if hipft_params.get('data_assimilation_mult_fac'):
      sed('assimilate_data_mult_fac', hipft_params.get('data_assimilation_mult_fac'), 'hipft.in')
     
    # Update hipft.in with data_assimilation_balance_flux
    if hipft_params.get('data_assimilation_balance_flux') is not None:
      sed('assimilate_data_balance_flux', hipft_params.get('data_assimilation_balance_flux'), 'hipft.in')      

    # Update hipft.in with np and nt
    if map_resolution_params:
      if map_resolution_params.get('np'):
        sed('res_np', map_resolution_params.get('np'), 'hipft.in')
      if map_resolution_params.get('nt'):
        sed('res_nt', map_resolution_params.get('nt'), 'hipft.in')

    # Ensure the start and end date are provided in the run parameters
    if not run_params.get('date_start'):
      check_err(1, 'No date_start under Run Parameters')
    if not run_params.get('date_end'):
      check_err(1, 'No date_end under Run Parameters')

    # Calculate the total physical time of run 
    print('    ==> Setting HipFT start and end times...')
    period_start = datetime.strptime(run_params.get("date_start"), '%Y-%m-%dT%H:%M:%S')
    min_datetime_utc = pytz.utc.localize(period_start)
    period_end = datetime.strptime(run_params.get("date_end"), '%Y-%m-%dT%H:%M:%S')
    max_datetime_utc = pytz.utc.localize(period_end)
    time_difference = (max_datetime_utc - min_datetime_utc).total_seconds() / 3600
    print(f'    ==> HipFT: Total physical time of run:  {time_difference} hours')    

    # Update hipft.in with time_start
    if args.first_obs_datetime is not None:
      first_obs_datetime = datetime.strptime(args.first_obs_datetime, '%Y-%m-%dT%H:%M:%S')
      first_obs_datetime_utc = pytz.utc.localize(first_obs_datetime)
      start_time = (min_datetime_utc - first_obs_datetime_utc).total_seconds() / 3600
      sed('time_start', start_time, 'hipft.in')
    else:
      start_time = 0.0 

    # Update hipft.in with time_end
    sed('time_end', start_time+time_difference, 'hipft.in')

    # Run hipft_input_lists.py to expand realizations into input file
    print(f'    ==> HipFT: Expanding realizations into input file...') 
    command = f'hipft_input_lists.py {args.output_yaml_file} {hipft_in_new}'
    ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
    check_err(ierr.returncode, f'Failed: {command}')
    args.logout.flush()
    args.logerr.flush()

    # Execute the HipFT if not a dry run
    if not args.dry_run:
      # Get the number of realizations from hipft.in
      result = os.popen(f'grep "  n_realizations = " "hipft.in"').read().strip()
      realizations = 0
      if result:
        try:
          realizations = int(result.split('=')[-1].strip())
        except ValueError:
          print(f"Integer not found in {result}.")
      else:
        print("n_realizations not found in hipft.in.")
      # Get number of mpi procs specified
      ranks = run_params.get("mpi_num_procs")

      # Check if there are more mpi procs then realizations
      if ranks > realizations:
        check_err(1, 'mpi_num_procs must be <= n_realizations')
      print(f'    ==> HipFT: Running code with {ranks} MPI process(es)...')

      # Execute HipFT
      command = f'{run_params.get("mpi_run_command", "").replace("<NP>", str(ranks))} hipft hipft.in 1>hipft.log 2>hipft.err'
      ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, f'Failed: {command}')
      args.logout.flush()
      args.logerr.flush()

      # Check if there was an error
      hipft_err = os.path.join(hipft_dir, "hipft.err")
      berr = check_if_file_is_empty(hipft_err)
      if not berr:
        print(f"Warning!  Error log of HipFT {hipft_err} not empty, something may have gone wrong.")


    # Change back to the original directory
    os.chdir(args.orun)
  else:
    print(f'    ==> Previous HipFT run folder exists')


def check_pt(file, run_params):
  """
    Checks if a file has the correct NP and NT dimensions and verifies that the file
    is in PT format by comparing its dimensions to those specified in the run parameters.
  """
  # Get the NP and NT dimensions from the file
  NP, NT = rdhdf_2d_dims(file)
  # Get the expected map resolution parameters from run_params
  map_resolution_params = run_params.get('map_resolution')
  # Check if the NP and NT values matches the expected value
  if NP != map_resolution_params.get('np'):
    check_err(1, f'NP of {map_resolution_params.get("np")} differs from the initial_map_filename: {NP}')
  if NT != map_resolution_params.get('nt'):
    check_err(1, f'NT of {map_resolution_params.get("nt")} differs from the initial_map_filename: {NT}')


def post_processing_submodule(args, run_params, hipft_post_process_params, step):
  """
    Manages running HipFT post-processing.
  """
  # If 'overwrite' is enabled, remove existing output folder
  if hipft_post_process_params.get('overwrite') is True:
    print(f'    ==> Overwriting existing  post_processing_{step} output folder')
    try:
      shutil.rmtree(f'post_processing_{step}')
    except FileNotFoundError:
      pass

  # Create the post-processing directory if it doesn't already exist and change to it
  post_processing_dir = os.path.join(args.orun, f'post_processing_{step}')
  os.makedirs(post_processing_dir, exist_ok=True)
  print(f'    ==> Created post_processing_{step} output folder: post_processing_{step}')
  os.chdir(post_processing_dir)

  # Construct the command to run the post-processing script
  command = 'hipft_post_process_all.py'
  # Set the output directory
  command += f' -output_dir {post_processing_dir}'
  # Set the run directory for HipFT data
  command += f' -rundir {os.path.join(args.orun, "hipft")}'
  # Samples for plotting
  command += f' -history_plot_samples {hipft_post_process_params.get("history_plot_samples")}'
  # Number plot markers
  command += f' -history_plot_samples_markers {hipft_post_process_params.get("history_plot_samples_markers")}'
  # Butterfly options
  command += f' -butterfly_make_options "-al {hipft_post_process_params.get("butterfly_averaging_width_hr")}"'

  # Determine the correct output path based on the step type ('raw' or 'processed')
  if step == 'raw':
    # Raw output path
    command += f' -outpath {os.path.join(args.orun, "hipft", "output_maps")}'
  elif step == 'processed':
    # Processed output path
    command += f' -outpath {os.path.join(args.orun, "hipft", "processed_maps")}'
    # History directory for processed data
    command += f' -hist_dir {os.path.join(args.orun, "hipft", "processed_maps")}'

  # Ensure the start date is provided in the run parameters and add it to the command
  if not run_params.get('date_start'):
    check_err(1, 'No date_start under Run Parameters')
  command += f' -utstart {run_params.get("date_start")}'

  # If 'overwrite' flag is set, include it in the command
  if hipft_post_process_params.get('overwrite') is True:
    command += f' -overwrite'
  
  # If dry run is not enabled, execute the post-processing command
  if not args.dry_run:
    print(f'    ==> Running post processing...')  
    print(f'    ==> Running post processing with command: {command}',file=args.logout)
    args.logout.flush()
    args.logerr.flush()    
    ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
    check_err(ierr.returncode, f'Failed: {command}')
    args.logout.flush()
    args.logerr.flush()

    # Copy the output map list files to the appropriate directories for both raw and processed data
    if step == 'raw':
      shutil.copy('hipft_output_map_list_tai.out',os.path.join(args.orun, "hipft", "output_maps"))
      shutil.copy('hipft_output_map_list_utc.out',os.path.join(args.orun, "hipft", "output_maps"))
    elif step == 'processed':
      shutil.copy('hipft_output_map_list_tai.out',os.path.join(args.orun, "hipft", "processed_maps"))
      shutil.copy('hipft_output_map_list_utc.out',os.path.join(args.orun, "hipft", "processed_maps"))
  
  # Change back to the original directory
  os.chdir(args.orun)


def map_process(args, output_map_processing_params, hipft_params):
  """
    This function handles the processing of output maps from HipFT. 
  """
  # Get map processing flags
  map_resolution_delta_deg = output_map_processing_params.get('map_resolution_delta_deg')

  # Initialize command string for map preparation
  command_map_prep=''

  # Adjust map resolution based on the specified delta
  if map_resolution_delta_deg:
    new_np = int(np.ceil(360.0/map_resolution_delta_deg) + 1.0)
    delta_angle_np = 360.0/(new_np-1.0)
    new_nt = int(np.ceil(180.0/map_resolution_delta_deg) + 1.0)
    delta_angle_nt = 180.0/(new_nt-1.0)
    print(f'    ==> Processed maps resolution: {new_np}x{new_nt}  dt_deg: {delta_angle_nt}  dp_deg: {delta_angle_np}')
    command_map_prep += f' -np {new_np}'
    command_map_prep += f' -nt {new_nt}'
  # Add smoothing factor
  smoothing_factor = output_map_processing_params.get('smoothing_factor')
  if smoothing_factor > 0.0:
    command_map_prep += f' -smoothfac {smoothing_factor}'
  else:
    command_map_prep += f' -nosmooth'
  # Add map multiplier
  if output_map_processing_params.get('map_multiplier'):
    command_map_prep += f' -mfac {output_map_processing_params.get("map_multiplier")}'
  #print(f'    ==> Map processing flags:  {command_map_prep}')
 
  # Define necessary directory paths
  hipft_dir               = os.path.join(args.orun, 'hipft')
  hipft_output_maps_dir   = os.path.join(hipft_dir, 'output_maps')
  processed_maps_dir      = os.path.join(hipft_dir, 'processed_maps')
  processed_maps_wsa_dir  = os.path.join(hipft_dir, 'processed_maps_wsa')
  processed_maps_temp_dir = os.path.join(processed_maps_dir, 'temp')
  processed_maps_wsa_temp_dir = os.path.join(processed_maps_wsa_dir, 'temp')

  # Skip execution if dry run is enabled
  if (args.dry_run):
    return

  # Check if previous run exists and stop if not overwriting
  if os.path.isdir(processed_maps_dir) and not output_map_processing_params.get('overwrite') is True:
     print(f'    ==> Previous processed_maps folder exits!')
  else:
    # Overwrite existing processed maps folder if specified
    if output_map_processing_params.get('overwrite') is True:
      print(f'    ==> Overwriting existing processed_maps folder')
      try:
        shutil.rmtree(processed_maps_dir)
      except FileNotFoundError:
          pass
      
    # Create the processed_maps folder and temporary folder for map processing
    os.makedirs(processed_maps_dir, exist_ok=True)
    print(f'    ==> Created processed_maps folder: processed_maps')
    os.makedirs(processed_maps_temp_dir, exist_ok=True)

    # List original (3D) output maps
    output_maps_pattern = os.path.join(hipft_output_maps_dir, 'hipft_brmap_idx*.h5')
    output_maps_list = sorted(glob.glob(output_maps_pattern))
    
    # Extract idx range from output maps
    idx_list = sorted({match.group(0) for f in output_maps_list if (match := re.search(r'idx(\d+)', f))})
    t0 = int(idx_list[0][3:])
    tf = int(idx_list[-1][3:])

    # Check the dimensionality of the first map file to see if 2D or 3D
    with h5.File(output_maps_list[0], 'r') as h5file:
      f = h5file['Data']
      ndims = np.ndim(f)
      if ndims == 3:
        is3d = True
      elif ndims == 2:
        is3d = False
      else:
        check_err(10,'Invalid number of dimensions ({ndims}) in {file_list[0]}') 
 
    # Change to temporary directory for processed maps
    print(f'    ==> Entering processed maps temporary directory: processed_maps/temp')
    os.chdir(processed_maps_temp_dir)

    # Iterate over output maps and extract or copy them
    print(f'    ==> Extracting/copying output maps:')
    for output_map in output_maps_list:
      if is3d:
        command = f'hipft_extract_realization.py {output_map}'
      else:
        temp_map=Path(output_map).stem+'.h5'
        command = f'cp {output_map} {temp_map}'
      idx_current = int(re.search(r'idx(\d+)', output_map).group(1))
      print(f'\r        ==> Extracting/copying map {idx_current}/{tf}', end = " ") 
      ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, f'Failed: {command}')
      args.logout.flush()
      args.logerr.flush()
    print('')
    
    # Get all maps in temporary folder
    maps_in_temp = os.path.join(processed_maps_temp_dir, 'hipft_brmap_idx*.h5')
    maps_in_temp_list = sorted(glob.glob(maps_in_temp))

    # Process the maps with psi_map_prep
    print(f'    ==> Processing output maps with psi_map_prep.py...')
    realization_list = sorted({match.group(0) for f in maps_in_temp_list if (match := re.search(r'r(\d+)', f))})
    rf = int(realization_list[-1][1:]) if realization_list else None
    # Process each map file based on the realization index
    if rf:
      for map in maps_in_temp_list:
        r_current = int(re.search(r'r(\d+)', map).group(1))
        idx_current = int(re.search(r'idx(\d+)', map).group(1))
        print(f'\r        ==> Processing map index: {idx_current}/{tf} realization: {r_current}/{rf}', end = " ") 
        fullCommand = f'psi_map_prep.py {map} {command_map_prep}'
        ierr = subprocess.run(["bash","-c",fullCommand],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, f'Failed : {fullCommand}')
        args.logout.flush()
        args.logerr.flush()
    else:
      for map in maps_in_temp_list:
        idx_current = int(re.search(r'idx(\d+)', map).group(1))
        print(f'\r        ==> Processing map index: {idx_current}/{tf}', end = " ") 
        fullCommand = f'psi_map_prep.py {map} {command_map_prep}'
        ierr = subprocess.run(["bash","-c",fullCommand],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, f'Failed : {fullCommand}')
        args.logout.flush()
        args.logerr.flush()
    print('')

    # If 3D files, pack them back into a cube
    if is3d:
      print('    ==> Packing processed map realizations into 3D map files:')
      for idx in idx_list:
        processed_maps_idx = os.path.join(processed_maps_temp_dir, f'processed_hipft_brmap_{idx}_r*.h5')
        processed_maps_idx_list = sorted(glob.glob(processed_maps_idx))
        command = 'hipft_pack_realizations.py'
        command += f' {",".join(processed_maps_idx_list)}'
        idx_current = int(re.search(r'idx(\d+)', idx).group(1))
        print(f'\r        ==> Packing map index {idx_current}/{tf}', end = " ")
        ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, f'Failed: {command}')
        args.logout.flush()
        args.logerr.flush()
      print('')        

    # Move final processed maps to the main directory
    for idx in idx_list:
      oname = os.path.join(processed_maps_temp_dir, f'processed_hipft_brmap_{idx}.h5')
      fname = os.path.join(processed_maps_dir, f'hipft_brmap_{idx}.h5')
      os.rename(oname, fname)

    # Clean up temporary directories
    print('    ==> Cleaning up')
    os.chdir(processed_maps_dir)
    shutil.rmtree(processed_maps_temp_dir)

    # Run history generation command
    command = 'hipft_get_histories_from_files.py'
    command += f' -folder "{processed_maps_dir}"'
    command += f' -t0 {t0}'
    command += f' -tf {tf}'
    if not is3d:
      command += f' -o hipft_history_sol_r000001.out'
    if hipft_params.get('output_map_cadence_hr'):
      command += f' -cadence { hipft_params.get("output_map_cadence_hr")}'
    command += f' -bfile hipft_brmap_idx'
    print('    ==> Generating history files from processed maps...')
    ierr = subprocess.run(["bash", "-c", command],stdout=args.logout,stderr=args.logerr)
    check_err(ierr.returncode, f'Failed: {command}')
    args.logout.flush()
    args.logerr.flush()

    # Create WSA grid data if selected:
    make_wsa_dataset = output_map_processing_params.get('make_wsa_dataset')
    if make_wsa_dataset:
      print('    ==> Making copy of processed maps interpolated to WSA grid:')

      # Create necessary directories for WSA data
      os.makedirs(processed_maps_wsa_temp_dir, exist_ok=True)
      os.makedirs(processed_maps_wsa_dir, exist_ok=True)
      print(f'    ==> Created processed_maps_wsa folder: processed_maps_wsa')
      # Change to temporary directory for WSA grid data
      print(f'    ==> Entering processed maps temporary directory: processed_maps_wsa/temp')
      os.chdir(processed_maps_wsa_temp_dir)

      # Adjust resolution for WSA grid
      new_np_wsa = new_np-1
      new_nt_wsa = new_nt-1

      # Generate WSA grid template
      print('    ==> Generating WSA grid template:')
      new_grid_command = f'oft_make_2d_inner_mesh_template.py -np {new_np_wsa} -nt {new_nt_wsa} -o {processed_maps_wsa_temp_dir}/oft_template_2d_inner_mesh_np{new_np_wsa}_nt{new_nt_wsa}.h5'
      ierr = subprocess.run(["bash","-c",new_grid_command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, 'Failed : '+new_grid_command)
      args.logout.flush()
      args.logerr.flush()      

      command_map_prep_wsa = f' -nosmooth -noflux -template {processed_maps_wsa_temp_dir}/oft_template_2d_inner_mesh_np{new_np_wsa}_nt{new_nt_wsa}.h5'

      # Make list of processed maps:
      processed_maps_pattern = os.path.join(processed_maps_dir, 'hipft_brmap_idx*.h5')
      processed_maps_list = sorted(glob.glob(processed_maps_pattern))

      # Iterate over files in hipft processed_maps
      print(f'    ==> Extracting/copying processed maps to temporary folder:')
      for processed_map in processed_maps_list:
        if is3d:
          command = f'hipft_extract_realization.py {processed_map}'
        else:
          temp_map=Path(processed_map).stem+'.h5'
          command = f'cp {processed_map} {temp_map}'
        idx_current = int(re.search(r'idx(\d+)', processed_map).group(1))
        print(f'\r        ==> Extracting/copying map {idx_current}/{tf}', end = " ")
        ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, f'Failed: {command}')
        args.logout.flush()
        args.logerr.flush()        
      print('')        

      #Get all maps in processed_maps_temp_dir
      maps_in_temp = os.path.join(processed_maps_wsa_temp_dir, 'hipft_brmap_idx*.h5')
      maps_in_temp_list = sorted(glob.glob(maps_in_temp))

      # Iterate over files in processed_maps_temp_dir and run psi_map_prep
      print(f'    ==> Processing processed maps with psi_map_prep.py ...')
      if rf:
        for map in maps_in_temp_list:
          r_current = int(re.search(r'r(\d+)', map).group(1))
          idx_current = int(re.search(r'idx(\d+)', map).group(1))
          print(f'\r        ==> Processing map index: {idx_current}/{tf} realization: {r_current}/{rf}', end = " ") 
          fullCommand = f'psi_map_prep.py {map} {command_map_prep_wsa}'
          ierr = subprocess.run(["bash","-c",fullCommand],stdout=args.logout,stderr=args.logerr)
          check_err(ierr.returncode, f'Failed : {fullCommand}')
          args.logout.flush()
          args.logerr.flush()
      else:
        for map in maps_in_temp_list:
          idx_current = int(re.search(r'idx(\d+)', map).group(1))
          print(f'\r        ==> Processing map index: {idx_current}/{tf}', end = " ") 
          fullCommand = f'psi_map_prep.py {map} {command_map_prep_wsa}'
          ierr = subprocess.run(["bash","-c",fullCommand],stdout=args.logout,stderr=args.logerr)
          check_err(ierr.returncode, f'Failed : {fullCommand}')
          args.logout.flush()
          args.logerr.flush()
      print('')

      # If 3d files pack back into a cube
      if is3d:
        print('    ==> Packing processed map realizations on WSA grid into 3D map files:')
        for idx in idx_list:
          processed_maps_idx = os.path.join(processed_maps_wsa_temp_dir, f'processed_hipft_brmap_{idx}_r*.h5')
          processed_maps_idx_list = sorted(glob.glob(processed_maps_idx))
          command = 'hipft_pack_realizations.py'
          command += f' {",".join(processed_maps_idx_list)}'
          idx_current = int(re.search(r'idx(\d+)', idx).group(1))
          print(f'\r        ==> Packing map index {idx_current}/{tf}',end = " ")
          ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
          check_err(ierr.returncode, f'Failed: {command}')
          args.logout.flush()
          args.logerr.flush()
      print('')

      # Move final processed maps to processed_maps_dir
      print(f'    ==> Moving processed maps on WSA grid from temp folder to processed_maps_wsa')
      for idx in idx_list:
        oname = os.path.join(processed_maps_wsa_temp_dir, f'processed_hipft_brmap_{idx}.h5')
        fname = os.path.join(processed_maps_wsa_dir, f'hipft_brmap_{idx}.h5')
        os.rename(oname, fname)

      # Clean up WSA temporary directory
      print('    ==> Cleaning up')
      os.chdir(processed_maps_wsa_dir)
      shutil.rmtree(processed_maps_wsa_temp_dir)


def check_err(ierr,message):
  """
    Checks the error code and prints an error message if the code is greater than 0.
    Exits the program if an error is encountered.
  """
  # If the error code is greater than 0, print the message and exit
  if ierr > 0:
    print(' ')
    print(message)
    print('Value of error code: '+str(ierr))
    sys.exit(1)


def sed(match, value, file):
  """
    Performs a search and replace in the given file using the 'sed' command.
    It replaces the entire line containing the 'match' string with a new line 
    containing 'match = value'.
  """
  
  # Set logicals to the correct Fortran value.
  if (value is True):
    value = ".true."
  if (value is False):
    value = ".false."
  
  os.system(f'sed -i "s|.*{match}.*|  {match} = {value}|" "{file}"')


def generate_unique_filename(base_name, extension):
  """
    Generate a unique filename by appending a counter if the file already exists.
    The counter starts at 0 and increments until a non-existing filename is found.
  """
  counter = 0
  # Loop until a unique filename is found
  while True:
    # Create the filename with or without a counter based on the current value
    filename = f"{base_name}{f'_{counter}' if counter > 0 else ''}{extension}"
    # Check if the file already exists
    if not os.path.exists(filename):
      # If the counter is greater than 0, print a warning that the base file existed
      if counter > 0:
        print(f'    ==> WARNING:  File exists:     {base_name}{extension}')
        print(f'    ==>           Writing logs to: {filename}')
      return filename
    # Increment the counter to check for the next possible filename
    counter += 1


def collect_output_submodule(args):
  """
    Collects the meta data files and copies them to the run_info directory.
  """
  # Define the path for the 'run_info' directory and create it
  run_info_dir = os.path.join(args.orun, 'run_info')
  os.makedirs(run_info_dir, exist_ok=True)

  # Define the path for the 'hipft' files 
  hipft_dir = os.path.join(args.orun, 'hipft')
  hipft_out = os.path.join(hipft_dir, 'hipft_timing.out')

  # Copy the HipFT meta data files to the 'run_info' directory
  shutil.copy(hipft_out,hipft_out.replace(hipft_dir, run_info_dir))


def get_hipft_add_dates_to_map_output_list(args,run_params, step):
  """
    Runs the `hipft_add_dates_to_map_output_list.py` script to add date information to 
    the HipFT output map list.
  """
  # Define the working directory and change to it
  output_map_dir = os.path.join(args.orun, "hipft", step)
  os.chdir(output_map_dir)

  # Construct the command to run the HipFT script
  command =  'hipft_add_dates_to_map_output_list.py'
  # Add the start date parameter
  command += f' {run_params.get("date_start")}' 
  # Input file
  command += f' -maplist {os.path.join(args.orun, "hipft", "hipft_output_map_list.out")}'
  # Output file name
  command += f' -o {"hipft_output_map_list"}'

  # Execute the command using subprocess and check for errors
  ierr = subprocess.run(["bash","-c",command],stdout=args.logout,stderr=args.logerr)
  check_err(ierr.returncode, f'Failed: {command}')


def check_if_file_is_empty(file_path):
    """
    Checks if the specified file is empty.
    """
    try:
        return os.path.exists(file_path) and os.stat(file_path).st_size == 0
    except Exception as e:
        print(f"Error checking file: {e}")
        return False
    

def main():
  """
    This is the main master call function. 
  """
  print('')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                       ╔═╗╔═╗╔╦╗╔═╗╦ ╦╔═╗')
  print('                       ║ ║╠╣  ║ ╚═╗║║║╠═╣')
  print('                       ╚═╝╚   ╩ ╚═╝╚╩╝╩ ╩')
  print(f'                          Version {version}')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')

  # Get input parameters
  args = argParsing()
  if args.dry_run:
    print('')
    print('                   Dry run mode activated!!                       ')
    print('')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('')

  print(f'==> Checking for OFT installation.')
  
  # Checks if the OFT environment is loaded
  check_oft_loaded(args)

  print(f'==> Reading input file: {args.input_yaml_file}')
  #print(f'==> Reading input file: {args.input_yaml_file}',file=args.logput)

  # Reads and merges the default and user-provided YAML files
  yaml_data = read_run_yaml(args)

  #Set date if "NOW" sepcified for end date
  current_date = yaml_data.get('run', {}).get('date_end').strip()
  if current_date.upper() == "now":
    current_date = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
    yaml_data['run']['date_end'] = current_date
    print(f"==> Date end set to current date: {current_date}")
  else:
    try:
        datetime.strptime(current_date, '%Y-%m-%dT%H:%M:%S')
    except ValueError:
        check_err(1,f"Invalid date format for date_end: {current_date}")

  # Get the runame and check if it is specified
  run_params = yaml_data.get('run', {})
  runname = run_params.get('run_name')
  if not runname:
    sys.exit("ERROR: A runname must be provided in the input yaml file")
  print(f'==> Run name: {runname}')

  # Get the path to the output directory
  args.outdir = args.outdir or os.path.join(os.getcwd(), 'oftswa_runs')
  args.outdir = os.path.abspath(args.outdir)

  # Check if the run exists and create or overwrite it
  args.orun = os.path.join(args.outdir, runname)
  if args.overwrite:
    print(f'==> Removing preexisting run folder (overwrite requested)')
    try:
      shutil.rmtree(args.orun)
    except FileNotFoundError:
      pass
  if (os.path.isdir(args.orun)):
    print(f'==> WARNING: A run called {runname} already exists.')
  os.makedirs(args.orun, exist_ok=True)
  print(f'==> Output directory: {args.orun}')

  # Output the run yaml file
  args.output_yaml_file = generate_unique_filename(os.path.join(args.outdir, runname, f'{runname}'),'.yaml')
  
  
  with open(args.output_yaml_file, 'w') as stream:
      yaml.add_representer(list, represent_pair_of_lists)
      yaml.dump(yaml_data, stream)
  print(f'==> Run parameters written to: {args.output_yaml_file}')
  
  # specify the logout and logerr files
  log_base_name = os.path.join(args.orun,'oftswa')
  log_filename = generate_unique_filename(log_base_name, ".log")
  err_filename = generate_unique_filename(log_base_name, ".err")
  args.logout = open(log_filename, 'a')
  args.logerr = open(err_filename, 'a')

  # Get the OMP_NUM_THREADS from the YAML and environment
  original_omp_num_threads = os.getenv('OMP_NUM_THREADS')
  omp_num_threads = run_params.get('omp_num_threads')

  if omp_num_threads:
    # Set OMP_NUM_THREADS from YAML
    print(f'==> OMP_NUM_THREADS set by input yaml file: {omp_num_threads}')
    os.environ['OMP_NUM_THREADS'] = str(int(omp_num_threads))
  else:
    # Check and output the enironment OMP_NUM_THREADS
    if original_omp_num_threads:
      print(f'==> OMP_NUM_THREADS set by environment:  {original_omp_num_threads}')
    else:
      print(f'==> OMP_NUM_THREADS not set, using default threads')

  # Get the number of mpi_num_procs from the YAML
  tmp = run_params.get('mpi_num_procs')
  print(f'==> MPI processes/ranks:                 {tmp}')
  
  print('')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                            CONFLOW')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')
  
  # Check if an preexisting conflow was specified and get the conflow parameters
  preexisting_conflow_runname = run_params.get('preexisting_conflow_runname')
  conflow_params = yaml_data.get('conflow', {})

  if preexisting_conflow_runname:
    print('==> Preexisitng ConFlow run selected')
    # Get path to preexisting conflow run folder
    preexisting_conflow = os.path.join(args.outdir, preexisting_conflow_runname,'conflow')
    preexisting_conflow = os.path.abspath(preexisting_conflow)
    print(f'    ==> Preexisting ConFlow specified: {preexisting_conflow}')
    if os.path.isdir(preexisting_conflow):
      # Get the preexisting conflow.dat
      preexisting_conflow_dat = os.path.join(preexisting_conflow, 'conflow.dat')
      # Check if the preexisting conflow run has the correct NT and NP
      read_result = os.popen(f'grep "  n_long = " "{preexisting_conflow_dat}"').read().strip()
      result = int(re.search(r"n_long = (\d+)", read_result).group(1)) 
      map_resolution_params = run_params.get('map_resolution')
      if result != map_resolution_params.get('np'):
        check_err(1, f'NP of {map_resolution_params.get("np")} differs from the preexisting conflow run: {result}')
      read_result = os.popen(f'grep "  n_lat = " "{preexisting_conflow_dat}"').read().strip()
      result = int(re.search(r"n_lat = (\d+)", read_result).group(1)) 
      if result != map_resolution_params.get('nt'):
        check_err(1, f'NT of {map_resolution_params.get("nt")} differs from the preexisting conflow run: {result}')
      
      # Get path to current output conflow folder
      conflow_dir = os.path.join(args.orun, 'conflow')
      # Create symbolic link if it does not exist
      if not os.path.islink(conflow_dir):
        relative_path = os.path.relpath(preexisting_conflow, args.orun)
        os.symlink(relative_path, conflow_dir)
        print(f'    ==> Symbolic link to ConFlow run created.')
      else:
        print(f'    ==> Symbolic link to ConFlow run already exists.')
    else:
      check_err(1, '==> ERROR: Could not find the preexisting ConFlow run.')
  elif conflow_params.get('run') is True:
    # Run conflow if preexisting conflow not specific and run set to True
    print('==> Running ConFlow...')
    conflow_submodule(args, conflow_params, run_params)
  else:
    print('==> WARNING: ConFlow disabled.')
  
  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                            MAGMAP')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  
  # Check if an preexisting MagMAP was specified and get the MagMAP parameters
  preexisting_magmap_runname = run_params.get('preexisting_magmap_runname')
  magmap_params = yaml_data.get('magmap', {})
  args.first_obs_datetime = None
  if preexisting_magmap_runname:
    # Check the preexisting MagMAP folder and create a symbolic link
    print('==> Preexisitng MagMAP run selected')
    args.first_obs_datetime = preexist_magmap_submodule(args, preexisting_magmap_runname, run_params)
  elif magmap_params.get('run') is True:
    # Run MagMAP if preexisting MagMAP not specific and run set to True
    print('==> Starting MagMAP run...')
    args.first_obs_datetime = magmap_submodule(args, magmap_params, run_params)
  else:
    print('==> WARNING: MagMAP disabled')

  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                            HIPFT')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  

  # Get the HipFT parameters
  hipft_params = yaml_data.get('hipft', {})
  if hipft_params.get('run') is True:
    # Run HipFT if run set to True
    print('==> Starting HipFT ...')
    hipft_submodule(args, hipft_params, run_params)
  else:
    print('==> HipFT is disabled.')
    
  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                         POST PROCESSING')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  
  # Get the HipFT post-processing parameters
  hipft_post_process_params = yaml_data.get('hipft_post_process', {})
  if hipft_post_process_params.get('run') is True:
    # Run HipFT post-processing if run set to True
    if os.path.isdir(os.path.join(args.orun, "hipft", "output_maps")):
      # If output_maps exists run post-processing
      print('==> Starting hipft post-processing...')
      post_processing_submodule(args, run_params, hipft_post_process_params, 'raw')
    else:
      # If output_maps does not exists give warning
      print('==> HipFT post-processing could not run. Path did not exist:')
      print(f'\t {os.path.join(args.orun, "hipft", "output_maps")}')
  else:
    # If run set to False only add date information to the HipFT output map list
    print('==> Adding dates to map output list...')
    get_hipft_add_dates_to_map_output_list(args, run_params, "output_maps")

  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                         MAP PROCESSING')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  
  # Get the HipFT map-processing parameters
  output_map_processing_params = yaml_data.get('output_map_processing', {})
  if output_map_processing_params.get('run') is True:
    # Run HipFT map-processing if run set to True
    print('==> Starting HipFT map processing...')
    map_process(args, output_map_processing_params, hipft_params)
  else:
    print('==> Map processing disabled.')

  # Run HipFT post-processing if run set to True
  if (hipft_post_process_params.get('run') is True):
  
    print('')  
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('                  PROCESSED MAPS POST PROCESSING')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('')  
    
    # Get the HipFT post-processing parameters
    hipft_post_process_params = yaml_data.get('hipft_post_process', {})
    if os.path.isdir(os.path.join(args.orun, "hipft", "processed_maps")):
      # If processed_maps exists run post-processing
      print('==> Starting HipFT post-processing...')
      post_processing_submodule(args, run_params, hipft_post_process_params, 'processed')
    else:
      # If processed_maps does not exists give warning
      print('==> ERROR! HipFT post process on processed maps failed. Path did not exist:')
      print(f'\t {os.path.join(args.orun, "hipft", "processed_maps")}')
  else:
    # If run set to False only add date information to the HipFT output map list
    print('==> Adding dates to map output list...')
    get_hipft_add_dates_to_map_output_list(args, run_params, "output_maps")

#  print('')  
#  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
#  print('                         COLLECTING OUTPUT')
#  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') 
  
#  # TODO:  Copy all log files to run_info folder.
#  if not args.dry_run:
#    collect_output_submodule(args)

  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')

  # Set the environment OMP_NUM_THREADS back to it's original value
  if original_omp_num_threads and omp_num_threads:
    print(f'==> Resetting environment OMP_NUM_THREADS back to: {original_omp_num_threads}')
    os.environ['OMP_NUM_THREADS'] = original_omp_num_threads

  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                      OFTSWA RUN COMPLETE!')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  
  args.logout.close()
  args.logerr.close()

if __name__ == '__main__':
  main()

# TODO NEW:   Insert timers for each major section - Dsilay the total time in days, hours, minute,s seconds for each section and a timing summary at the end???

# Conflow:
# Run the code (check OMP_NUM_THREADS, ACC_NUM_CORES)
# Make sure it ran correctly 


# Meta data collection:   grab hipt.out, etc etc etc etc and put into "run_info" folder.
#Error checking throughout....
