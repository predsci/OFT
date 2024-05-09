#!/usr/bin/env python3
import numpy as np
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
else:
  tvec_centers_input = xvec
  pvec_centers_input = yvec

# Get area metric:
# Modify boundary points of theta for correct near-pole calcualtions.
tvec_tmp = tvec_centers_input.copy()
tvec_tmp[0]  = (tvec_centers_input[0]  + tvec_centers_input[1] )*0.25
tvec_tmp[-1] = tvec_centers_input[-1] - np.abs((tvec_centers_input[-1] - tvec_centers_input[-2])*0.25)
da_input_map = input_data
for i in range(len(pvec_centers_input)):
  da_input_map[i,:] = np.sin(tvec_tmp[:])

print(da_input_map.shape)
print(np.min(da_input_map))
print(np.max(da_input_map))

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

#print(np.all(np.diff(tvec_edges_input) > 0))
#print(np.all(np.diff(pvec_edges_input) > 0))
#print(np.all(np.diff(tvec_edges_output) > 0))
#print(np.all(np.diff(pvec_edges_output) > 0))

print(np.min(tvec_edges_input))
print(np.max(tvec_edges_input))
print(np.min(pvec_edges_input))
print(np.max(pvec_edges_input))

print(np.min(tvec_centers_input))
print(np.max(tvec_centers_input))
print(np.min(pvec_centers_input))
print(np.max(pvec_centers_input))

# Remesh the map:
output_data = downsample_grid.downsamp_reg_grid(input_data, tvec_edges_input,  pvec_edges_input,  \
                                                            pvec_edges_output, tvec_edges_output, \
                                                            da=da_input_map)
# Write out result:
ps.wrhdf_2d(args.output_file,pvec_centers_output,tvec_centers_output,np.transpose(output_data))
