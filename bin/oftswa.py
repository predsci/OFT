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
#  Version 1.0.0 :  Initial "working" version.  Still in development.
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
import psi_io


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


# Check if oft loaded
def check_oft_loaded(args):
  # Get the path variable
  path_env = os.environ.get("PATH", "")
  path_dirs = path_env.split(os.pathsep)

  # Loop through path directories to see if magmap, conflow, and hipft loaded
  for directory in path_dirs:
    if directory.endswith("/bin"):
      oft_dir = directory[:-4]
      required_subdirs = [
          oft_dir,
          os.path.join(oft_dir, "magmap", "bin"),
          os.path.join(oft_dir, "conflow", "bin"),
          os.path.join(oft_dir, "hipft", "bin"),
      ]
      if all(subdir in path_dirs for subdir in required_subdirs):
        print(f'    ==> Found OFT installation: {oft_dir}')
        args.oft_dir = oft_dir
        break
  else:
    sys.exit("ERROR: Please load the OFT environment by sourcing 'load_oft_env.sh'")


# Make the output yaml have [[],[]] instead of other format
def represent_pair_of_lists(dumper, value):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', value, flow_style=True)


# Read Yaml
def read_run_yaml(args):

  # Read the default yaml
  default_yaml = os.path.join(args.oft_dir, 'rsrc','default.yaml')
  print(f'    ==> Reading default yaml file from OFT rsrc directory.')
  with open(default_yaml, 'r') as stream:
    default_data = yaml.safe_load(stream)
  
  # Read the input yaml
  print(f'    ==> Reading user input yaml file.')
  with open(args.input_yaml_file, 'r') as stream:
    input_data = yaml.safe_load(stream)

  # Check if sweep1d is on
  isSweep1d = input_data.get('hipft', {}).get('realization_combination_mode') == 'sweep1d'

  # Combine the yamls
  print(f'    ==> Merging default and user yaml files.')
  deep_update(default_data, input_data, isSweep1d)

  # Replace empty string values with None
  default_data = replace_empty_strings_with_none(default_data)

  return default_data


# Recursively update the default yaml with the input yaml
def deep_update(dict1, dict2, isSweep1d):
  # iterate through dictionary items
  for key, value in dict2.items():
    # if dictionary recursively call
    if isinstance(value, dict) and key in dict1 and isinstance(dict1[key], dict):
      # if a realization_parameters combine specially
      if key == 'realization_parameters':
        combine_update(dict1[key], value, isSweep1d)
      else:
        deep_update(dict1[key], value, isSweep1d)
    # update the default value with the input yaml
    else:
      dict1[key] = value


# Special update of the default yaml for realization_parameters
def combine_update(dict1, dict2, isSweep1d):
  for key, value in dict2.items():
    # if Sweep1d and a list of numbers, combine default and input yaml to construct [[],[]]
    if all(isinstance(x, (int, float)) for x in value) and isSweep1d:
      dict1[key] = [dict1[key],value]
    # if not Sweep1d or [[],[]] use the input yaml data
    elif len(value) == 2 or not isSweep1d:
      dict1[key] = value
    else:
      check_err(1,f'Bad input format for {key}')


# Replace empty string values recursively with None
def replace_empty_strings_with_none(data):
    if isinstance(data, dict):
        return {k: replace_empty_strings_with_none(v) for k, v in data.items()}
    elif data == "":
        return None
    else:
        return data


# Run conflow
def conflow_submodule(args, conflow_params, run_params):
  conflow_dir = os.path.join(args.orun, 'conflow')
  if not os.path.isdir(conflow_dir) or conflow_params.get('overwrite') is True:
    if conflow_params.get('overwrite') is True:
      print(f'    ==> Overwriting existing conflow run folder...')
      if os.path.islink(conflow_dir):
          os.unlink(conflow_dir)
      elif os.path.exists(conflow_dir):
          shutil.rmtree(conflow_dir)
    os.makedirs(conflow_dir, exist_ok=True)
    print(f'    ==> Conflow output folder: {conflow_dir}')
    os.chdir(conflow_dir)
    conflow_dat = os.path.join(args.oft_dir, "rsrc", "conflow.dat")
    Command = f'cp {conflow_dat} {conflow_dir}'
    ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
    check_err(ierr.returncode, Command)
    args.logout.flush()
    args.logerr.flush()    
    map_resolution_params = run_params.get('map_resolution')
    if map_resolution_params:
      if map_resolution_params.get('np'):
        sed('n_long', map_resolution_params.get('np'), 'conflow.dat')
      if map_resolution_params.get('nt'):
        sed('n_lat', map_resolution_params.get('nt'), 'conflow.dat')
    print(f'    ==> Copied {conflow_dat} to {conflow_dir}')
    print(f'    ==> OMP_NUM_THREADS set to {os.getenv("OMP_NUM_THREADS")}')
    print(f'    ==> ACC_NUM_CORES set to {os.getenv("ACC_NUM_CORES")}')
    print(f'    ==> Running conflow...')
    if not args.dry_run:
      Command = 'conflow'
      ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, 'Failed : '+Command)
      args.logout.flush()
      args.logerr.flush()
    os.chdir(args.orun)
  else:
     print(f'    ==> Previous conflow run folder found.')


# Get the NT, NP dimension of a h5 file
# This shoudl really not need to read data, just the dimension sizes
def rdhdf_2d_dims(h5_filename):
    with h5.File(h5_filename, 'r') as h5file:
        data = h5file['Data']
        ndims = data.ndim
        x = data.dims[0][0] if ndims > 0 and len(data.dims[0].keys()) != 0 else np.array([])
        y = data.dims[1][0] if ndims > 1 and len(data.dims[1].keys()) != 0 else np.array([])
        NT = len(x)
        NP = len(y)
    return NT, NP

# If preexisting magmap specified
def preexist_magmap_submodule(args, preexisting_magmap_runname, run_params):
    first_obs_datetime = None

    # Get the relative preexisting_magmap directory
    preexisting_magmap = os.path.join(args.outdir, preexisting_magmap_runname,'magmap')
    preexisting_magmap = Path(preexisting_magmap).resolve()
    print(f'    ==> Preexisting magmap specified: {preexisting_magmap}')

    # Create/ Check the symbolic link to the preexisting_magmap directory
    if os.path.isdir(preexisting_magmap):
      magmap_dir = os.path.join(args.orun, 'magmap')
      if not os.path.islink(magmap_dir):
        os.symlink(preexisting_magmap, magmap_dir)
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
        check_err(1, f'date_start not in the existing magmap directory')

      # Check the magmap end date to the input yaml end date
      period_end = datetime.strptime(run_params.get("date_end"), '%Y-%m-%dT%H:%M:%S')
      period_end_utc = pytz.utc.localize(period_end)
      magmap_end = datetime.strptime(last_obs_datetime, '%Y-%m-%dT%H:%M:%S')
      magmap_end_utc = pytz.utc.localize(magmap_end)
      if magmap_end_utc < period_end_utc:
        check_err(1, f'date_end not in the existing magmap directory')

      #Check if the magmap data maps are in pt format and the dimensions specified in the input yaml
      file = os.path.join(magmap_dir,'magmap_data_maps',data['map_path'].iloc[0])
      NP, NT = rdhdf_2d_dims(file)
      map_resolution_params = run_params.get('map_resolution')
      if NP != map_resolution_params.get('np'):
        check_err(1, f'NP of {map_resolution_params.get("np")} differs from the preexisting magmap run: {NP}')
      if NT != map_resolution_params.get('nt'):
        check_err(1, f'NT of {map_resolution_params.get("nt")} differs from the preexisting magmap run: {NT}')

      return first_obs_datetime
    else:
      check_err(1, '==> ERROR: Could not find the preexisting MagMAP run.')


# Run magmap
def magmap_submodule(args, magmap_params, run_params):
  magmap_dir = os.path.join(args.orun, 'magmap')
  if not os.path.isdir(magmap_dir) or magmap_params.get('overwrite') is True or magmap_params.get('update') is True:
    if magmap_params.get('overwrite') is True and magmap_params.get('update') is True:
      check_err(1, 'Cannot update and overwite magmap at the same time!')
    if magmap_params.get('overwrite') is True:
      print(f'    ==> Overwriting existing magmap run folder')
      if os.path.islink(magmap_dir):
          os.unlink(magmap_dir)
      elif os.path.exists(magmap_dir):
          shutil.rmtree(magmap_dir)
    if magmap_params.get('update') is True:
      print(f'    ==> Updating existing magmap run folder')
    os.makedirs(magmap_dir, exist_ok=True)
    print(f'    ==> Created magmap output folder: {magmap_dir}')
    os.chdir(magmap_dir)
    # Run magmap_get_data.py
    magmap_data_disks = os.path.join(magmap_dir, 'magmap_data_disks')

    # This is a secret option for advance users...
    if magmap_params.get('raw_data_dir'):
      print(f'    ==> Existing download directory provided for magmap_get_data.py :')
      existing_magmap_data_disks = Path(magmap_params.get('raw_data_dir')).resolve()
      print(f'\t {existing_magmap_data_disks}')
      magmap_data_disks = os.path.join(magmap_dir, 'magmap_data_disks')
      if not os.path.islink(magmap_data_disks):
        os.symlink(existing_magmap_data_disks, magmap_data_disks)
        print('    ==> Creating symbolic link to :')
      else:
        print('    ==> Symbolic link already exists :')
      print('    ==> Creating symbolic link to :')
      print(f'\t {magmap_data_disks}')
      print(f'    ==> Skipping magmap_get_data.py') 
    else:
      print(f'    ==> Running magmap_get_data.py') 
      Command = 'magmap_get_data.py'
      print(f'    ==> Download directory for magmap_get_data.py is :')
      print(f'\t {magmap_data_disks}')
      Command += f' -odir {magmap_data_disks}'
      if magmap_params.get('map_cadence_hr'):
        Command += f' -cadence {magmap_params.get("map_cadence_hr")}'
      if not run_params.get('date_start'):
        check_err(1, 'No date_start under Run Parameters')
      if not run_params.get('date_end'):
        check_err(1, 'No date_end under Run Parameters')
      Command += f' {run_params.get("date_start")} {run_params.get("date_end")}'
      if not args.dry_run:
        ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, 'Failed : '+Command)
        args.logout.flush()
        args.logerr.flush()
    Command = f'magmap_disk2map.py {magmap_data_disks}'
    print(f'    ==> Running magmap_disk2map.py')
    print(f'    ==> Magnetogram disk images data directory is :')
    magmap_data_maps = os.path.join(magmap_dir, 'magmap_data_maps')
    print(f'\t {magmap_data_maps}')
    Command += f' -odir {magmap_data_maps}'
    Command += f' -startdate {run_params.get("date_start")} -enddate {run_params.get("date_end")}'
    map_resolution_params = run_params.get('map_resolution')
    if map_resolution_params:
      print(f'    ==> Magnetogram disk images will have size:  {map_resolution_params.get("np")} x {map_resolution_params.get("nt")}')
      Command += f' -npout {map_resolution_params.get("np")}'
      Command += f' -ntout {map_resolution_params.get("nt")}'
    if not args.dry_run:
      ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, 'Failed : '+Command)
      args.logout.flush()
      args.logerr.flush()
      os.chdir(args.orun)
  else:
    print(f'    ==> Previous magmap run folder exits!')


# Run hipft
def hipft_submodule(args, hipft_params, run_params):
  hipft_dir = os.path.join(args.orun, 'hipft')
  if not os.path.isdir(hipft_dir) or hipft_params.get('overwrite') is True:
    if hipft_params.get('overwrite') is True:
      print(f'    ==> Overwriting existing HipFT run folder')
      try:
        shutil.rmtree(hipft_dir)
      except FileNotFoundError:
          pass
    os.makedirs(hipft_dir, exist_ok=True)
    print(f'    ==> Created HipFT output folder')
    os.chdir(hipft_dir)
    hipft_in_new = os.path.join(hipft_dir, "hipft.in")
    hipft_in = os.path.join(args.oft_dir, "rsrc", "hipft.in")
    Command = f'cp {hipft_in} {hipft_in_new}'
    ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
    check_err(ierr.returncode, Command)
    args.logout.flush()
    args.logerr.flush()    
    print(f'    ==> Copied HipFT template input file to run directory')
    if hipft_params.get('output_map_cadence_hr'):
      sed('output_map_time_cadence', hipft_params.get('output_map_cadence_hr'), 'hipft.in')
    if hipft_params.get('initial_map_filename'):
      input_map_dir = os.path.join(hipft_dir, 'input_map')
      os.makedirs(input_map_dir, exist_ok=True)
      input_map = Path(hipft_params.get('initial_map_filename')).resolve()
      Command = f'cp {input_map} {input_map_dir}'
      ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, Command)
      args.logout.flush()
      args.logerr.flush()
      check_pt(input_map, run_params)
      from_link = os.path.join(input_map_dir, input_map.name)
      to_link = os.path.join(input_map_dir, 'br_input_map.h5')
      if not os.path.islink(to_link):
        os.symlink(from_link, to_link)
        print(f'    ==> Creating symbolic link to {input_map_dir}')
      else:
        print(f'    ==> Symbolic link already exists: {input_map_dir}')
      print(f'    ==> Initial map found: {input_map}')
    else:
        check_err(1, 'initial_map_filename not provided')
#
    if hipft_params.get('initial_map_mult_fac'):
      sed('initial_map_mult_fac', hipft_params.get('initial_map_mult_fac'), 'hipft.in')
#
    if hipft_params.get('random_flux_lifetime'):
      sed('source_rfe_lifetime', hipft_params.get('random_flux_lifetime'), 'hipft.in')
    map_resolution_params = run_params.get('map_resolution')
#
    if hipft_params.get('data_assimilation_mult_fac'):
      sed('assimilate_data_mult_fac', hipft_params.get('data_assimilation_mult_fac'), 'hipft.in')

    if map_resolution_params:
      if map_resolution_params.get('np'):
        sed('res_np', map_resolution_params.get('np'), 'hipft.in')
      if map_resolution_params.get('nt'):
        sed('res_nt', map_resolution_params.get('nt'), 'hipft.in')
#
    if not run_params.get('date_start'):
      check_err(1, 'No date_start under Run Parameters')
    if not run_params.get('date_end'):
      check_err(1, 'No date_end under Run Parameters')
#
    print('    ==> Setting HipFT start and end times...')
    period_start = datetime.strptime(run_params.get("date_start"), '%Y-%m-%dT%H:%M:%S')
    min_datetime_utc = pytz.utc.localize(period_start)
    period_end = datetime.strptime(run_params.get("date_end"), '%Y-%m-%dT%H:%M:%S')
    max_datetime_utc = pytz.utc.localize(period_end)
    time_difference = (max_datetime_utc - min_datetime_utc).total_seconds() / 3600
    print(f'    ==> HipFT: Total physical time of run:  {time_difference} hours')    
#
    if args.first_obs_datetime is not None:
      first_obs_datetime = datetime.strptime(args.first_obs_datetime, '%Y-%m-%dT%H:%M:%S')
      first_obs_datetime_utc = pytz.utc.localize(first_obs_datetime)
      start_time = (min_datetime_utc - first_obs_datetime_utc).total_seconds() / 3600
      sed('time_start', start_time, 'hipft.in')
    else:
      start_time = 0.0 
#
    sed('time_end', start_time+time_difference, 'hipft.in')
#
    # Run hipft_input_lists.py
    print(f'    ==> HipFT: Expanding realizations into input file...') 
    Command = f'hipft_input_lists.py {args.output_yaml_file} {hipft_in_new}'
    ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
    check_err(ierr.returncode, 'Failed : '+Command)
    args.logout.flush()
    args.logerr.flush()
    if not args.dry_run:
      result = os.popen(f'grep "  n_realizations = " "hipft.in"').read().strip()
      realizations = 0
      if result:
        try:
          realizations = int(result.split('=')[-1].strip())
        except ValueError:
          print(f"Integer not found in {result}.")
      else:
        print("n_realizations not found in hipft.in.")
      ranks = run_params.get("mpi_num_procs")
      if ranks > realizations:
        check_err(1, 'mpi_num_procs must be <= n_realizations')
      print(f'    ==> HipFT: Running code with {ranks} MPI process(es)...')
      Command = f'{run_params.get("mpi_run_command", "").replace("<NP>", str(ranks))} hipft hipft.in 1>hipft.log 2>hipft.err'
      ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, 'Failed : '+Command)
      args.logout.flush()
      args.logerr.flush()
    os.chdir(args.orun)
  else:
    print(f'    ==> Previous HipFT run folder exists')


# Check if a file has the right np and nt dimensions and is in PT format
def check_pt(file, run_params):
  NP, NT = rdhdf_2d_dims(file)
  map_resolution_params = run_params.get('map_resolution')
  if NP != map_resolution_params.get('np'):
    check_err(1, f'NP of {map_resolution_params.get("np")} differs from the initial_map_filename: {NP}')
  if NT != map_resolution_params.get('nt'):
    check_err(1, f'NT of {map_resolution_params.get("nt")} differs from the initial_map_filename: {NT}')


# Run post_processing
def post_processing_submodule(args, run_params, hipft_post_process_params, step):
  if hipft_post_process_params.get('overwrite') is True:
    print(f'    ==> Overwriting existing  post_processing_{step} output folder')
    try:
      shutil.rmtree(f'post_processing_{step}')
    except FileNotFoundError:
      pass

  #Get map post processing directory and maje
  post_processing_dir = os.path.join(args.orun, f'post_processing_{step}')
  os.makedirs(post_processing_dir, exist_ok=True)
  print(f'    ==> Created post_processing_{step} output folder: post_processing_{step}')

  # Change to directry
  os.chdir(post_processing_dir)

  # Construct command for hipft_post_process_all.py
  Command = 'hipft_post_process_all.py'
  Command += f' -output_dir {post_processing_dir}'
  Command += f' -rundir {os.path.join(args.orun, "hipft")}'
  Command += f' -history_plot_samples {hipft_post_process_params.get("history_plot_samples")}'
  Command += f' -history_plot_samples_markers {hipft_post_process_params.get("history_plot_samples_markers")}'

  # If raw output from hipft make correct outpath
  if step == 'raw':
    Command += f' -outpath {os.path.join(args.orun, "hipft", "output_maps")}'
  # If map_process output make correct paths
  elif step == 'processed':
    Command += f' -outpath {os.path.join(args.orun, "hipft", "processed_maps")}'
    Command += f' -hist_dir {os.path.join(args.orun, "hipft", "processed_maps")}'

  # Set utstart
  if not run_params.get('date_start'):
    check_err(1, 'No date_start under Run Parameters')
  Command += f' -utstart {run_params.get("date_start")}'

  # Call overwrite if specified
  if hipft_post_process_params.get('overwrite') is True:
    Command += f' -overwrite'
  
  if not args.dry_run:
    print(f'    ==> Running post processing...')  
    print(f'    ==> Running post processing with command: {Command}',file=args.logout)
    ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
    check_err(ierr.returncode, 'Failed : '+Command)
    args.logout.flush()
    args.logerr.flush()    

    # Copy TAI and UTC output map list .out text files into the output_maps and processed_maps directories
    if step == 'raw':
      shutil.copy('hipft_output_map_list_tai.out',os.path.join(args.orun, "hipft", "output_maps"))
      shutil.copy('hipft_output_map_list_utc.out',os.path.join(args.orun, "hipft", "output_maps"))
    elif step == 'processed':
      shutil.copy('hipft_output_map_list_tai.out',os.path.join(args.orun, "hipft", "processed_maps"))
      shutil.copy('hipft_output_map_list_utc.out',os.path.join(args.orun, "hipft", "processed_maps"))
  
  os.chdir(args.orun)


def map_process(args, output_map_processing_params, hipft_params):

  #Get map processing flags
  map_resolution_delta_deg = output_map_processing_params.get('map_resolution_delta_deg')

  Command_map_prep=''

  if map_resolution_delta_deg:
    new_np = int(np.ceil(360.0/map_resolution_delta_deg) + 1.0)
    delta_angle_np = 360.0/(new_np-1.0)
    new_nt = int(np.ceil(180.0/map_resolution_delta_deg) + 1.0)
    delta_angle_nt = 180.0/(new_nt-1.0)
    print(f'    ==> Processed maps resolution: {new_np}x{new_nt}  dt_deg: {delta_angle_nt}  dp_deg: {delta_angle_np}')
    Command_map_prep += f' -np {new_np}'
    Command_map_prep += f' -nt {new_nt}'
  if output_map_processing_params.get('smoothing_factor'):
    Command_map_prep += f' -smoothfac {output_map_processing_params.get("smoothing_factor")}'
  if output_map_processing_params.get('map_multiplier'):
    Command_map_prep += f' -mfac {output_map_processing_params.get("map_multiplier")}'
  #print(f'    ==> Map processing flags:  {Command_map_prep}')
 
  # Needed directory paths:
  hipft_dir               = os.path.join(args.orun, 'hipft')
  hipft_output_maps_dir   = os.path.join(hipft_dir, 'output_maps')
  processed_maps_dir      = os.path.join(hipft_dir, 'processed_maps')
  processed_maps_wsa_dir  = os.path.join(hipft_dir, 'processed_maps_wsa')
  processed_maps_temp_dir = os.path.join(processed_maps_dir, 'temp')
  processed_maps_wsa_temp_dir = os.path.join(processed_maps_wsa_dir, 'temp')

  if (args.dry_run):
    return

  # Check if previous run exists and stop if not overwriting
  if os.path.isdir(processed_maps_dir) and not output_map_processing_params.get('overwrite') is True:
     print(f'    ==> Previous processed_maps folder exits!')
  # Run map processing
  else:
    # Check if overwriting and overwrite
    if output_map_processing_params.get('overwrite') is True:
      print(f'    ==> Overwriting existing processed_maps folder')
      try:
        shutil.rmtree(processed_maps_dir)
      except FileNotFoundError:
          pass
    # Create processed_maps folder
    os.makedirs(processed_maps_dir, exist_ok=True)
    print(f'    ==> Created processed_maps folder: processed_maps')
    # Make temporary folder for map processing
    os.makedirs(processed_maps_temp_dir, exist_ok=True)
    # Make list of original (3D) output maps:
    output_maps_pattern = os.path.join(hipft_output_maps_dir, 'hipft_brmap_idx*.h5')
    output_maps_list = sorted(glob.glob(output_maps_pattern))
    # Get idx range:
    idx_list = sorted({match.group(0) for f in output_maps_list if (match := re.search(r'idx(\d+)', f))})
    t0 = int(idx_list[0][3:])
    tf = int(idx_list[-1][3:])
    #check if 3d h5 file:
    with h5.File(output_maps_list[0], 'r') as h5file:
      f = h5file['Data']
      ndims = np.ndim(f)
      if ndims == 3:
        is3d = True
      elif ndims == 2:
        is3d = False
      else:
        check_err(10,'Invalid number of dimensions ({ndims}) in {file_list[0]}') 
 
    # Change to processed_maps_temp_dir
    print(f'    ==> Entering processed maps temporary directory: processed_maps/temp')
    os.chdir(processed_maps_temp_dir)

    # Iterate over files in hipft output_maps
    print(f'    ==> Extracting/copying output maps:')
    for output_map in output_maps_list:
      if is3d:
        Command = f'hipft_extract_realization.py {output_map}'
      else:
        temp_map=Path(output_map).stem+'.h5'
        Command = f'cp {output_map} {temp_map}'
      print(f'\r        ==> Extracting/copying map {Path(output_map).stem}', end = " ") 
      ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, 'Failed : '+Command)
      args.logout.flush()
      args.logerr.flush()
    print('')
    
    #Get all maps in processed_maps_temp_dir
    maps_in_temp = os.path.join(processed_maps_temp_dir, 'hipft_brmap_idx*.h5')
    maps_in_temp_list = sorted(glob.glob(maps_in_temp))

    # Iterate over files in processed_maps_temp_dir and run psi_map_prep
    print(f'    ==> Processing output maps with psi_map_prep.py...')
    for map in maps_in_temp_list:
#TEST!    
      print(f'\r        ==> Processing map {Path(map).stem} ...', end = " ") 
      FullCommand = f'psi_map_prep.py {map} {Command_map_prep}'
      ierr = subprocess.run(["bash","-c",FullCommand],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, 'Failed : '+FullCommand)
      args.logout.flush()
      args.logerr.flush()
    print('')

    # If 3d files pack back into a cube
    if is3d:
      print('    ==> Packing processed map realizations into 3D map files:')
      for idx in idx_list:
        processed_maps_idx = os.path.join(processed_maps_temp_dir, f'processed_hipft_brmap_{idx}_r*.h5')
        processed_maps_idx_list = sorted(glob.glob(processed_maps_idx))
        Command = 'hipft_pack_realizations.py'
        Command += f' {",".join(processed_maps_idx_list)}'
        print(f'\r        ==> Packing map index {idx} ...', end = " ")
        ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, 'Failed : '+Command)
        args.logout.flush()
        args.logerr.flush()
      print('')        

    # Move final processed maps to processed_maps_dir
    for idx in idx_list:
      oname = os.path.join(processed_maps_temp_dir, f'processed_hipft_brmap_{idx}.h5')
      fname = os.path.join(processed_maps_dir, f'hipft_brmap_{idx}.h5')
      os.rename(oname, fname)

    # Change directory and remove temp directory
    print('    ==> Cleaning up')
    os.chdir(processed_maps_dir)
    shutil.rmtree(processed_maps_temp_dir)

    # Run hipft_get_histories_from_files.py
    Command = 'hipft_get_histories_from_files.py'
    Command += f' -folder "{processed_maps_dir}"'
    Command += f' -t0 {t0}'
    Command += f' -tf {tf}'
    if hipft_params.get('output_map_cadence_hr'):
      Command += f' -cadence { hipft_params.get("output_map_cadence_hr")}'
    Command += f' -bfile hipft_brmap_idx'
    print('    ==> Generating history files from processed maps...')
    ierr = subprocess.run(["bash", "-c", Command],stdout=args.logout,stderr=args.logerr)
    check_err(ierr.returncode, 'Failed : '+Command)
    args.logout.flush()
    args.logerr.flush()

    # Create WSA grid data if selected:
    make_wsa_dataset = output_map_processing_params.get('make_wsa_dataset')
    if make_wsa_dataset:
    
      print('    ==> Making copy of processed maps interpolated to WSA grid:')

      os.makedirs(processed_maps_wsa_temp_dir, exist_ok=True)
      os.makedirs(processed_maps_wsa_dir, exist_ok=True)
      print(f'    ==> Created processed_maps_wsa folder: processed_maps_wsa')
     # Change to processed_maps_temp_dir
      print(f'    ==> Entering processed maps temporary directory: processed_maps_wsa/temp')
      os.chdir(processed_maps_wsa_temp_dir)

      new_np_wsa = new_np-1
      new_nt_wsa = new_nt-1

      print('    ==> Generating WSA grid template:')

      new_grid_command = f'oft_make_2d_inner_mesh_template.py -np {new_np_wsa} -nt {new_nt_wsa} -o {processed_maps_wsa_temp_dir}/oft_template_2d_inner_mesh_np{new_np_wsa}_nt{new_nt_wsa}.h5'
      ierr = subprocess.run(["bash","-c",new_grid_command],stdout=args.logout,stderr=args.logerr)
      check_err(ierr.returncode, 'Failed : '+new_grid_command)
      args.logout.flush()
      args.logerr.flush()      

      Command_map_prep_wsa = f' -nosmooth -noflux -template {processed_maps_wsa_temp_dir}/oft_template_2d_inner_mesh_np{new_np_wsa}_nt{new_nt_wsa}.h5'

      # Make list of processed maps:
      processed_maps_pattern = os.path.join(processed_maps_dir, 'hipft_brmap_idx*.h5')
      processed_maps_list = sorted(glob.glob(processed_maps_pattern))

      # Iterate over files in hipft processed_maps
      print(f'    ==> Extracting/copying processed maps to temporary folder:')
      for processed_map in processed_maps_list:
        if is3d:
          Command = f'hipft_extract_realization.py {processed_map}'
        else:
          temp_map=Path(processed_map).stem+'.h5'
          Command = f'cp {processed_map} {temp_map}'
        print(f'\r        ==> Extracting/copying map {Path(processed_map).stem}', end = " ")
        ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, 'Failed : '+Command)
        args.logout.flush()
        args.logerr.flush()        
      print('')        

      #Get all maps in processed_maps_temp_dir
      maps_in_temp = os.path.join(processed_maps_wsa_temp_dir, 'hipft_brmap_idx*.h5')
      maps_in_temp_list = sorted(glob.glob(maps_in_temp))

      # Iterate over files in processed_maps_temp_dir and run psi_map_prep
      print(f'    ==> Processing processed maps with psi_map_prep.py ...')
      for map in maps_in_temp_list:
        print(f'\r        ==> Processing map {Path(map).stem} ...', end = " ")
        FullCommand = f'psi_map_prep.py {map} {Command_map_prep_wsa}'
        ierr = subprocess.run(["bash","-c",FullCommand],stdout=args.logout,stderr=args.logerr)
        check_err(ierr.returncode, 'Failed : '+FullCommand)
        args.logout.flush()
        args.logerr.flush()
      print('')

      # If 3d files pack back into a cube
      if is3d:
        print('    ==> Packing processed map realizations on WSA grid into 3D map files:')
        for idx in idx_list:
          processed_maps_idx = os.path.join(processed_maps_wsa_temp_dir, f'processed_hipft_brmap_{idx}_r*.h5')
          processed_maps_idx_list = sorted(glob.glob(processed_maps_idx))
          Command = 'hipft_pack_realizations.py'
          Command += f' {",".join(processed_maps_idx_list)}'
          print(f'\r        ==> Packing map index {idx} ...',end = " ")
          ierr = subprocess.run(["bash","-c",Command],stdout=args.logout,stderr=args.logerr)
          check_err(ierr.returncode, 'Failed : '+Command)
          args.logout.flush()
          args.logerr.flush()
      print('')

      # Move final processed maps to processed_maps_dir
      print(f'    ==> Moving processed maps on WSA grid from temp folder to processed_maps_wsa')
      for idx in idx_list:
        oname = os.path.join(processed_maps_wsa_temp_dir, f'processed_hipft_brmap_{idx}.h5')
        fname = os.path.join(processed_maps_wsa_dir, f'hipft_brmap_{idx}.h5')
        os.rename(oname, fname)

      print('    ==> Cleaning up')
      os.chdir(processed_maps_wsa_dir)
      shutil.rmtree(processed_maps_wsa_temp_dir)


def check_err(ierr,message):
  if ierr > 0:
    print(' ')
    print(message)
    print('Value of error code: '+str(ierr))
    sys.exit(1)


def sed(match, value, file):
  os.system(f'sed -i "s/.*{match}.*/  {match} = {value}/" "{file}"')


def generate_unique_filename(base_name, extension):
    """
    Generate a unique filename based on the base name and extension.
    Appends a counter if the file already exists.
    """
    counter = 0
    while True:
        filename = f"{base_name}{f'_{counter}' if counter > 0 else ''}{extension}"
        if not os.path.exists(filename):
            if counter > 0:
              print(f'    ==> WARNING:  File exists:     {base_name}{extension}')
              print(f'    ==>           Writing logs to: {filename}')
            return filename
        counter += 1


def collect_output_submodule(args):
  run_info_dir = os.path.join(args.orun, 'run_info')
  hipft_dir = os.path.join(args.orun, 'hipft')
  os.makedirs(run_info_dir, exist_ok=True)

  hipft_out = os.path.join(hipft_dir, 'hipft_timing.out')
  shutil.copy(hipft_out,hipft_out.replace(hipft_dir,run_info_dir))


# Main Script
def main():

  version = '1.0.0'

  print('')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                       ╔═╗╔═╗╔╦╗╔═╗╦ ╦╔═╗')
  print('                       ║ ║╠╣  ║ ╚═╗║║║╠═╣')
  print('                       ╚═╝╚   ╩ ╚═╝╚╩╝╩ ╩')
  print(f'                          Version {version}')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')

  args = argParsing()
  if args.dry_run:
    print('')
    print('                   Dry run mode activated!!                       ')
    print('')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('')

  print(f'==> Checking for OFT installation.')
  
  check_oft_loaded(args)

  print(f'==> Reading input file: {args.input_yaml_file}')
  #print(f'==> Reading input file: {args.input_yaml_file}',file=args.logput)

  yaml_data = read_run_yaml(args)
  
  # Get the runame
  run_params = yaml_data.get('run', {})
  runname = run_params.get('run_name')
  if not runname:
    sys.exit("ERROR: A runname must be provided in the input yaml file")
  print(f'==> Run name: {runname}')

  args.outdir = args.outdir or os.path.join(os.getcwd(), 'oftswa_runs')
  args.outdir = Path(args.outdir).resolve()

  # Check if run exists and create or overwrite it
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
  args.output_yaml_file = os.path.join(args.outdir, runname, f'{runname}.yaml')
  with open(args.output_yaml_file, 'w') as stream:
      yaml.add_representer(list, represent_pair_of_lists)
      yaml.dump(yaml_data, stream)
  print(f'==> Run parameters written to: {args.output_yaml_file}')
  
  log_base_name = os.path.join(args.orun,'oftswa')
  log_filename = generate_unique_filename(log_base_name, ".log")
  err_filename = generate_unique_filename(log_base_name, ".err")
  args.logout = open(log_filename, 'a')
  args.logerr = open(err_filename, 'a')

  # Grab OMP_NUM_THREADS from shell, 
  # if not empty, save as "old_omp_num_threads"
  # After run is all done, reset the ENV to the old value.

  omp_num_threads = run_params.get('omp_num_threads')
  if omp_num_threads:
    print(f'==> Setting OMP_NUM_THREADS to:  {omp_num_threads}')
    os.environ['OMP_NUM_THREADS'] = str(int(omp_num_threads))
  else:
    print(f'HERE')  #  Write out the user's current OMP_NUM_THREADS if set.
                    # Overwise, say "OMP not set, using default threads"
  tmp = run_params.get('mpi_num_procs')
  print(f'==> MPI processes/ranks:         {tmp}')
  
  print('')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                            CONFLOW')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')
  
  preexisting_conflow_runname = run_params.get('preexisting_conflow_runname')
  conflow_params = yaml_data.get('conflow', {})
  if preexisting_conflow_runname:
    print('==> Preexisitng ConFlow run selected')
    preexisting_conflow = os.path.join(args.outdir, preexisting_conflow_runname,'conflow')
    preexisting_conflow = Path(preexisting_conflow).resolve()
    print(f'    ==> Preexisting ConFlow specified: {preexisting_conflow}')
    if os.path.isdir(preexisting_conflow):
      preexisting_conflow_dat = os.path.join(preexisting_conflow, 'conflow.dat')
      read_result = os.popen(f'grep "  n_long = " "{preexisting_conflow_dat}"').read().strip()
      result = int(re.search(r"n_long = (\d+)", read_result).group(1)) 
      map_resolution_params = run_params.get('map_resolution')
      if result != map_resolution_params.get('np'):
        check_err(1, f'NP of {map_resolution_params.get("np")} differs from the preexisting conflow run: {result}')
      read_result = os.popen(f'grep "  n_lat = " "{preexisting_conflow_dat}"').read().strip()
      result = int(re.search(r"n_lat = (\d+)", read_result).group(1)) 
      if result != map_resolution_params.get('nt'):
        check_err(1, f'NT of {map_resolution_params.get("nt")} differs from the preexisting conflow run: {result}')
      conflow_dir = os.path.join(args.orun, 'conflow')  
      if not os.path.islink(conflow_dir):
        os.symlink(preexisting_conflow, conflow_dir)
        print(f'    ==> Symbolic link to ConFlow run created.')
      else:
        print(f'    ==> Symbolic link to ConFlow run already exists.')
    else:
      check_err(1, '==> ERROR: Could not find the preexisting ConFlow run.')
  elif conflow_params.get('run') is True:
    print('==> Running ConFlow...')
    conflow_submodule(args, conflow_params, run_params)
  else:
    print('==> WARNING: ConFlow disabled.')
  
  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                            MAGMAP')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  
  preexisting_magmap_runname = run_params.get('preexisting_magmap_runname')
  magmap_params = yaml_data.get('magmap', {})
  args.first_obs_datetime = None
  if preexisting_magmap_runname:
    print('==> Preexisitng MagMAP run selected')
    args.first_obs_datetime = preexist_magmap_submodule(args, preexisting_magmap_runname, run_params)
  elif magmap_params.get('run') is True:
    print('==> Starting MagMAP run...')
    magmap_submodule(args, magmap_params, run_params)
  else:
    print('==> WARNING: MagMAP disabled')

  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                            HIPFT')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  hipft_params = yaml_data.get('hipft', {})
  if hipft_params.get('run') is True:
      print('==> Starting HipFT ...')
      hipft_submodule(args, hipft_params, run_params)
  else:
    print('==> HipFT is disabled.')
    
  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                         POST PROCESSING')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  
  hipft_post_process_params = yaml_data.get('hipft_post_process', {})
  if hipft_post_process_params.get('run') is True:
    if os.path.isdir(os.path.join(args.orun, "hipft", "output_maps")):
      print('==> Starting hipft post-processing...')
      post_processing_submodule(args, run_params, hipft_post_process_params, 'raw')
    else:
      print('==> hipft post process could not run. Path did not exist:')
      print(f'\t {os.path.join(args.orun, "hipft", "output_maps")}')
  else:
    print('==> Post processing is disabled.')

  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                         MAP PROCESSING')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  
  output_map_processing_params = yaml_data.get('output_map_processing', {})
  if output_map_processing_params.get('run') is True:
    print('==> Starting HipFT map processing...')
    map_process(args, output_map_processing_params, hipft_params)
  else:
    print('==> Map processing disabled.')

  if (hipft_post_process_params.get('run') is True):
  
    print('')  
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('                  PROCESSED MAPS POST PROCESSING')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('')  
    
    hipft_post_process_params = yaml_data.get('hipft_post_process', {})
    if os.path.isdir(os.path.join(args.orun, "hipft", "processed_maps")):
      print('==> Starting HipFT post-processing...')
      post_processing_submodule(args, run_params, hipft_post_process_params, 'processed')
    else:  
      print('==> ERROR! HipFT post process on processed maps failed. Path did not exist:')
      print(f'\t {os.path.join(args.orun, "hipft", "processed_maps")}')
  else:
    print('==> Post processing is disabled')

  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                         COLLECTING OUTPUT')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') 
  
  # TODO:  Copy all log files to run_info folder.
  collect_output_submodule(args)

  print('')  
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('                      OFTSWA RUN COMPLETE!')
  print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  print('')  
  
  args.logout.close()
  args.logerr.close()

if __name__ == '__main__':
  main()

# For progres changing text lines (see "\r"), change (for example):
# "Processing map hipft_brmap_idx000005_r000001"
# to
# "Processing map index: 5/10 realization: 1/5"


# Run hipft_add_dates_to_map_output_list.py as part of hipft run step in case post processing is disabled.  Put output in all map folders and all post processing folders



# NEW:   Insert timers for each major section - Dsilay the total time in days, hours, minute,s seconds for each section and a timing summary at the end???




#  Directory structure:    
#  oftswa_runs/<runname>/conflow
#  oftswa_runs/<runname>/magmap
#  oftswa_runs/<runname>/hipft
#  oftswa_runs/<runname>/post_processing_raw
#  oftswa_runs/<runname>/post_processing_processed
#  oftswa_runs/<runname>/<runname>.yaml
#  Read DEFAULT YAML file from "rsrc" folder
#  Read USER YAML file and use any fields provided
#  Write "used" YAML to run folder.

# Conflow:
# In rsrc, we have a conflow.dat
# Copy that - and run the code (check OMP_NUM_THREADS, ACC_NUM_CORES)
# Make sure it ran correctly 

# MagMAP:
# Run magmap disk script based on uttimes and cadence
# Run magmap map script based on uttimes and cadence, etc.

# HipFT
# We need a hipft.in in rsrc
# Run realization script bnased on yaml selections (cross-mult, seq, manual)
# Use yaml dicutionary to sed things into hipft.in
# Run HipFT
# Run post processing 
# Run map processing  "/processed_maps"
# (Run make_histories on procesed maps, and run post proicessing)
# Meta data collection:   grab hipt.out, etc etc etc etc and put into "run_info" folder.
#Error checking throughout....
