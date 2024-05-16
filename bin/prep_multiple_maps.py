#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
from pathlib import Path


# INPUT:  - Map directory (use Path from pathlib to get full path)
#         - Output directory (default is new local folder called "output_swig")

def argParsing():
    parser = argparse.ArgumentParser(description='Run multiple maps through PSI Map Prep')

    parser.add_argument(help='Input directory of magnetogram files (h5).',
                        dest='input_directory')

    parser.add_argument('-psi_map_prep',
                        help='Path to psi_map_prep.py.',
                        dest='psi_map_prep',
                        required=False,
                        type=str)

    parser.add_argument('-outdir',
                        help='Directory where output will go.',
                        dest='outdir',
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

def run(args):

  # Get full path of input directory:
  input_directory = Path(args.input_directory).resolve()

  # Make rundir and go there
  if args.outdir is None:
      args.outdir = str(Path('.').resolve())+'/processed_maps'



  # Get all files in input directory
  h5_files = list(input_directory.glob('*.h5'))

  if len(h5_files) < 1:
    print(' ')
    print('No h5 files found.')
    sys.exit(1)

  if args.psi_map_prep is None:
    pmpdir = sys.path[0]
    args.psi_map_prep = pmpdir+'/psi_map_prep.py'

  args.psi_map_prep = str(Path(args.psi_map_prep).resolve())

  print(Path(args.psi_map_prep))

  os.makedirs(args.outdir, exist_ok=True)
  os.chdir(args.outdir)

  empty_idx=99999
  for h5_file in h5_files:
    h5_file=str(h5_file)
    oidx=' -o processed_'+h5_file.split('/')[-1]

    print('=> Running map : ' +h5_file.split('/')[-1])

    print(args.psi_map_prep)

    Command=args.psi_map_prep+' '+h5_file+oidx+ \
      ' -nt '+str(args.nt)+' -np '+str(args.np)+' -template '+str(args.template)+\
      ' -mfac '+str(args.multfac)+' -smoothfac '+str(args.smoothfac)

    if args.hipftexe is not None:
        Command=Command+' -hipftexe'

    if (args.flux):
      Command=Command+' -noflux'
    if (args.smooth):
      Command=Command+' -nosmooth'
    if (args.remap):
      Command=Command+' -noremap'

    ierr = subprocess.run(["bash","-c",Command])
    check_error_code_NON_CRASH(ierr.returncode,'Failed : '+Command)

def check_error_code_NON_CRASH(ierr,message):
  if ierr > 0:
    print(' ')
    print(message)
    print('Error code of fail : '+str(ierr))


def check_error_code(ierr,message):
  if ierr > 0:
    print(' ')
    print(message)
    print('Error code of fail : '+str(ierr))
    sys.exit(1)


def main():
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()
