#!/usr/bin/env python3
########################################################################
# ****** psi_map_prep.py: PSI magnetogram processing
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
#  Version 1.1.0
#
########################################################################

import numpy as np
import scipy.interpolate as sp
import argparse
from pathlib import Path
#
import downsample_grid
import psi_io as ps

# This can down-sample or up-sample, but only works with main-main mesh maps.
# Whether the input maps are pt or tp, the output is currently pt.

def argParsing():
    parser = argparse.ArgumentParser(\
             description='PSI_REMAP_MM: Remap a full-Sun magnetogram to a new resolution.  This only works with maps on a main-main psi mesh.')
    parser.add_argument(help='Input map filename.',
                        dest='input_file',
                        type=str)
    parser.add_argument('-nt',
                        help='Uniform theta resolution.',
                        dest='nt',
                        default='181',
                        type=int,
                        required=False)
    parser.add_argument('-np',
                        help='Uniform phi resolution.',
                        dest='np',
                        default='361',
                        type=int,
                        required=False)
    parser.add_argument('-template',
                        help='2D h5 file with desired grid scales.',
                        dest='template_file',
                        type=str,
                        required=False)
    parser.add_argument('-resample_template',
                        help='Resample the template to the resolution specified by nt and np.',
                        dest='resample_template',
                        default=False,
                        action='store_true',
                        required=False)
    parser.add_argument(help='Output map filename.',
                        dest='output_file',
                        type=str,
                        )
    return parser.parse_args()

args = argParsing()



# Get full path of input file:
args.input_file = str(Path(args.input_file).resolve())

# Read input map and make into tp:
xvec,yvec,input_data = ps.rdhdf_2d(args.input_file)
if (np.max(xvec) > 3.5):
  tvec_centers_input = yvec
  pvec_centers_input = xvec
  input_data = np.transpose(input_data)
  tp = False
else:
  tvec_centers_input = xvec
  pvec_centers_input = yvec
  tp = True

# Get area metric:
# Modify boundary points of theta for correct near-pole calcualtions.
tvec_tmp = tvec_centers_input.copy()
tvec_tmp[0]  = (tvec_centers_input[0]  + tvec_centers_input[1] )*0.25
tvec_tmp[-1] = tvec_centers_input[-1] - np.abs((tvec_centers_input[-1] - tvec_centers_input[-2])*0.25)
da_input_map = input_data.copy()
for i in range(len(pvec_centers_input)):
  da_input_map[i,:] = np.sin(tvec_tmp[:])

# Make scales for cell edges:
# As per Jamie's instructions, set boundary edges to original center positions.
# This assumes both maps are on a main-main mesh!

tvec_edges_input       = np.zeros(len(tvec_centers_input) + 1)
tvec_edges_input[1:-1] = (tvec_centers_input[0:-1] + tvec_centers_input[1:]) * 0.5
tvec_edges_input[0]    = tvec_centers_input[0]
tvec_edges_input[-1]   = tvec_centers_input[-1]

pvec_edges_input       = np.zeros(len(pvec_centers_input) + 1)
pvec_edges_input[1:-1] = (pvec_centers_input[0:-1] + pvec_centers_input[1:]) * 0.5
pvec_edges_input[0]    = pvec_centers_input[0]
pvec_edges_input[-1]   = pvec_centers_input[-1]

# Read template file if selected and make into tp, otherwise make uniform grid.
if (args.template_file is not None):
  xvec,yvec,_ = ps.rdhdf_2d(args.template_file)
  if (np.max(xvec) > 3.5):
    tvec_centers_output = yvec
    pvec_centers_output = xvec
  else:
    tvec_centers_output = xvec
    pvec_centers_output = yvec
    
# Now rescale the template to the requested nt and np if desired:
  if (args.resample_template):
    tvec_interpolator = sp.interp1d(np.arange(len(tvec_centers_output)), tvec_centers_output, fill_value="extrapolate")
    tvec_centers_output = tvec_interpolator(np.arange(args.nt)/(args.nt - 1)*(len(tvec_centers_output) - 1))
    pvec_interpolator = sp.interp1d(np.arange(len(pvec_centers_output)), pvec_centers_output, fill_value="extrapolate")
    pvec_centers_output = pvec_interpolator(np.arange(args.np)/(args.np - 1)*(len(pvec_centers_output) - 1))
else:
# Make uniform grid.
  tvec_centers_output = np.linspace(0,  np.pi,args.nt)
  pvec_centers_output = np.linspace(0,2*np.pi,args.np)

# Make output cell edges scales:
tvec_edges_output       = np.zeros(len(tvec_centers_output) + 1)
tvec_edges_output[1:-1] = (tvec_centers_output[0:-1] + tvec_centers_output[1:]) * 0.5
tvec_edges_output[0]    = tvec_centers_output[0]
tvec_edges_output[-1]   = tvec_centers_output[-1]

pvec_edges_output       = np.zeros(len(pvec_centers_output) + 1)
pvec_edges_output[1:-1] = (pvec_centers_output[0:-1] + pvec_centers_output[1:]) * 0.5
pvec_edges_output[0]    = pvec_centers_output[0]
pvec_edges_output[-1]   = pvec_centers_output[-1]

# Remesh the map:
output_data = downsample_grid.downsamp_reg_grid(input_data, tvec_edges_input,  pvec_edges_input,  \
                                                            pvec_edges_output, tvec_edges_output, \
                                                            da=da_input_map)
# Write out result:
if (tp):
  ps.wrhdf_2d(args.output_file,tvec_centers_output,pvec_centers_output,output_data)
else:
  ps.wrhdf_2d(args.output_file,pvec_centers_output,tvec_centers_output,np.transpose(output_data))
  
  
