#!/usr/bin/env python3
########################################################################
# ****** oft_make_2d_inner_mesh_template.py:
#        Make 2D longitude-colatitude inner mesh template.
#
#     Predictive Science Inc.
#     www.predsci.com
#     San Diego, California, USA 92121
########################################################################
# Copyright 2025 Predictive Science Inc.
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
import numpy as np
import argparse
import psi_io

def argParsing():
    parser = argparse.ArgumentParser(\
             description='Make 2D longitude-colatitude inner mesh template.  This is a grid that does not include points at phi=0 and phi=2*pi, nor theta=0 and theta=pi.  The points are cell centers, whose edges fill the domain.  This is the grid that the WSA model uses.')
    parser.add_argument('-np',
                        help='phi (longitude) resolution (default 360).',
                        dest='NP',
                        default='360',
                        type=int,
                        required=False)
    parser.add_argument('-nt',
                        help='theta (colatitude) resolution (default 180).',
                        dest='NT',
                        default='180',
                        type=int,
                        required=False)
    parser.add_argument("-o",
                        help='Output hdf5 file name (default is oft_template_2d_inner_mesh_np<NP>_nt<NT>.h5)',
                        dest='output_file',
                        required=False,
                        type=str)
    return parser.parse_args()

def main():

    # Get input agruments:
    args = argParsing()

    # Get desired resolution:
    NT = int(args.NT)
    NP = int(args.NP)

    # Set output name:
    if args.output_file is not None:
      oname = args.output_file
    else:
      oname = f'oft_template_2d_inner_mesh_np{NP}_nt{NT}.h5'

    dt = 180.0/NT
    dp = 360.0/NP

    tvec = np.linspace(dt/2.0,180-dt/2.0,num=NT)
    pvec = np.linspace(dp/2.0,360-dp/2.0,num=NP)

#    print('NT:         '+str(NT))
#    print('NP:         '+str(NP))
#    print('DT:         '+str(dt))
#    print('DP:         '+str(dp))
#    print('DT0:        '+str(tvec[1]-tvec[0]))
#    print('DP0:        '+str(pvec[1]-pvec[0]))
#    print('DTN:        '+str(tvec[NT-1]-tvec[NT-2]))
#    print('DPN:        '+str(pvec[NP-1]-pvec[NP-2]))    
#    print('min(theta): '+str(np.min(tvec)))
#    print('max(theta): '+str(np.max(tvec)))
#    print('min(phi):   '+str(np.min(pvec)))
#    print('max(phi):   '+str(np.max(pvec)))

    tvec_rad = tvec*(np.pi/180.0)
    pvec_rad = pvec*(np.pi/180.0)

#    print('min(theta) (rad): '+str(np.min(tvec_rad)))
#    print('max(theta) (rad): '+str(np.max(tvec_rad)))
#    print('min(phi)   (rad): '+str(np.min(pvec_rad)))
#    print('max(phi)   (rad): '+str(np.max(pvec_rad)))
 
    data = np.zeros((NT,NP))

    psi_io.wrhdf_2d(oname,pvec_rad,tvec_rad,data)

if __name__ == '__main__':
    main()


