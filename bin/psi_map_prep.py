#!/usr/bin/env python3
########################################################################
# ****** psi_remap_mm.py: PSI magnetogram remapping
#
#     Predictive Science Inc.
#     www.predsci.com
#     San Diego, California, USA 92121
########################################################################
# Copyright 2024 Predictive Science Inc.
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
########################################################################
#
#  Version 1.0.0
#
########################################################################

import signal
import numpy as np
import os
import subprocess
import argparse
import sys
from pathlib import Path
#
import psi_io

def signal_handler(signal, frame):
        print('You pressed Ctrl+C! Stopping!')
        sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)

def argParsing():

    parser = argparse.ArgumentParser(description='PSI Map Prep:  Process a hdf5 phi-theta magnetogram.  This script requires HIPFT (github.com/predsci/hipft) to be installed and "hipft" in the user path.  The psi_io.py must be in the current directory as well.')

    parser.add_argument(help='FULL PATH for input magnetogram file (h5).',
                        dest='inputfilename')

    parser.add_argument('-o',
                        help='Output file name.  Default is br_final.h5.',
                        dest='outputfilename',
                        required=False)

    parser.add_argument('-nt',
                        help='Uniform theta resolution.',
                        dest='nt',
                        default='181',
                        required=False)

    parser.add_argument('-np',
                        help='Uniform phi resolution.',
                        dest='np',
                        default='361',
                        required=False)

    parser.add_argument('-template',
                        help='2D h5 file with desired grid scales.',
                        dest='template',
                        required=False)

    parser.add_argument('-mfac',
                        help='Factor to multiply the Br map by (e.g. 1.4 for HMI to make it more like MDI)',
                        dest='multfac',
                        default=1.0,
                        type=float,
                        required=False)

    parser.add_argument('-noflux',
                        help='Do not flux balance the map.',
                        dest='flux',
                        action='store_false',
                        default=True,
                        required=False)

    parser.add_argument('-nosmooth',
                        help='Do not smooth the map.',
                        dest='smooth',
                        action='store_false',
                        default=True,
                        required=False)

    parser.add_argument('-smoothfac',
                        help='Smoothing factor.',
                        dest='smoothfac',
                        default='0.5',
                        type=float,
                        required=False)

    parser.add_argument('-noremap',
                        help='Do not remap the map to a different resolution.',
                        dest='remap',
                        action='store_false',
                        default=True,
                        required=False)

    parser.add_argument('-hipftexe',
                        help='Full path to hipft executable (otherwise assume hipft is in the path)',
                        dest='hipftexe',
                        type=str,
                        required=False)

    return parser.parse_args()


## Get input agruments:
args = argParsing()

ext='.h5'

# Set location of hipft if not specified.
# Should check if hipft exists in specified local, or in standard path, or in user path.
# If none, error out.
if args.hipftexe is None:
  args.hipftexe=sys.path[0]+'/../hipft/bin/hipft'

# Get base name of input file:
fileroot = Path(args.inputfilename).stem
fname=fileroot+ext

# Get full path of input file:
inputfilename = str(Path(args.inputfilename).resolve())

# Set default output name if not specified:
if args.outputfilename is None:
  args.outputfilename = fileroot+'_processed'+ext

print('')
print('############################################################')
print('#####  PSI MAP PREP')
print('#####  Input map name: '+args.inputfilename)
print('############################################################')
print('')

print('=> Making temporary directory \'map_processing\'')
os.makedirs('map_processing', exist_ok=True)
os.chdir('map_processing')
print('=> Making soft link br'+ext+' to input file.')
os.system('ln -sf '+inputfilename+' br'+ext)
curr_name='br'


if (args.multfac != 1.0):
    print('=> Reading in file to apply multiplicative factor ('+str(args.multfac)+')')
    x,y,brmap = psi_io.rdhdf_2d(inputfilename)
    brmap = args.multfac*brmap
    curr_name = curr_name + '_scaled'
    fname = curr_name + ext
    psi_io.wrhdf_2d(fname,x,y,brmap)
    print('=> Wrote file: '+fname)


if args.remap:
    print('')
    print('=> Remapping map to specified grid.')
    if args.template is not None:
      strtmp = '-template '+args.template
    else:
      strtmp = '-nt '+args.nt+' -np '+args.np
    command = sys.path[0]+'/psi_remap_mm.py '+strtmp+' '+curr_name+ext+' '+curr_name+'_remapped'+ext
    subprocess.run(["bash","-c",command])
    curr_name = curr_name+'_remapped'
    print('=> Wrote file: '+curr_name+ext)


if args.smooth or args.flux:
    print('=> Starting to prepare HipFT run...')
    if args.smooth and not args.flux:
        nm='smoothed'
    if args.flux and not args.smooth:
        nm='flxbal'
    if args.flux and args.smooth:
        nm='smoothed_flxbal'
    hipft_input_string = '&hipft_input_parameters\n'                        +\
                         '  initial_map_filename=\''+curr_name+ext+'\'\n'   +\
                         '  output_map_root_filename=\''+curr_name+'_'+nm+'\'\n'
    curr_name = curr_name+'_'+nm
    if args.smooth:
        print('    => Adding smoothing options to HipFT input file.')
        hipft_input_string = hipft_input_string                                 +\
                             '  diffusion_coef_factor='+str(args.smoothfac)+'\n'+\
                             '  time_end=1.0\n'                                 +\
                             '  advance_diffusion=.true.\n'                     +\
                             '  diffusion_coef_grid=.true.\n'
    if args.flux:
        print('    => Adding flux balancing to HipFT input file.')
        hipft_input_string = hipft_input_string +\
                             '  output_map_flux_balance=.true.\n'
    if args.smooth or args.flux:
        print('    => Writing HipFT input file to disk.')
        hipft_input_string = hipft_input_string + '/\n'
        f = open("hipft.in", "w")
        f.write(hipft_input_string)
        f.close()
#     - Run Hipft
#       ERROR CHECK - make sure hipft is in the path!
        command='mpiexec -np 1 '+args.hipftexe+' 1>hipft.log 2>hipft.err'
        print('    => Running HipFT with command: '+command)
        subprocess.run(["bash","-c",command])
        os.system('mv '+curr_name+'_final'+ext+' '+curr_name+ext)


# - Get result back (final step).
print('=> Copying final map.')
os.chdir('../')
os.system('cp map_processing/'+curr_name+ext+' '+args.outputfilename)

print('')
print('############################################################')
print('#####  Final map ready!')
print('#####  Final map name: '+args.outputfilename)
print('############################################################')
print('')
