#!/bin/bash

# Script to load OFT environment.
# This script sets all the paths to the OFT and its submodule executables.

# To use this script, source it in your BASH terminal with:
#   . load_oft_env.sh

# You should then have all executables and scripts in your path,
# allowing you to run OFT.
# Note - this assumes you have followed the instructions in each submodule
#        to build/compile each submodule.

cR="\033[1;31m"
cG="\033[32m"
cX="\033[0m"
cM="\033[35m"
echo="echo -e"

${echo} "${cG}==> OFT Environment Setup${cX}"

(return 0 2>/dev/null) && sourced=1 || sourced=0

if [ $sourced -eq "0" ]
then
  ${echo} "${cR}==> ERROR! It seems this script was executed instead of being sourced.${cX}"
  ${echo} "${cR}    Please source this script with: '. load_oft_env.sh'${cX}"
  exit 1
fi

oft_dir="$( dirname -- "$( readlink -f -- "${BASH_SOURCE[0]}"; )"; )"

echo "==> Checking that submodule components are installed..."

if [ ! -e ${oft_dir}/magmap/bin ]
then
  ${echo} "${cR}==> ERROR! MAGMAP submodule seems to be missing!${cX}"
  ${echo} "${cR}    Please ensure the submodule was checked out${cX}"
  ${echo} "${cR}    with: 'git clone --recursive' ${cX}"
  ${echo} "${cR}    or activate submodules with: 'git submodule update --init'.${cX}"
  return
fi

if [ ! -e ${oft_dir}/conflow/bin ]
then
  ${echo} "${cR}==> ERROR! CONFLOW submodule seems to be missing!${cX}"
  ${echo} "${cR}    Please ensure the submodule was checked out${cX}"
  ${echo} "${cR}    with: 'git clone --recursive' ${cX}"
  ${echo} "${cR}    or activate submodules with: 'git submodule update --init'.${cX}"
  return
fi

if [ ! -e ${oft_dir}/hipft/bin ]
then
  ${echo} "${cR}==> ERROR! HIPFT submodule seems to be missing!${cX}"
  ${echo} "${cR}    Please ensure the submodule was checked out${cX}"
  ${echo} "${cR}    with: 'git clone --recursive' ${cX}"
  ${echo} "${cR}    or activate submodules with: 'git submodule update --init'.${cX}"
  return
fi

if [ ! -e ${oft_dir}/conflow/bin/conflow ]
then
  ${echo} "${cM}==> WARNING! CONFLOW does not seem to be built.${cX}"
  ${echo} "${cM}    If you plan to use ConFlow, please ensure the submodule${cX}"
  ${echo} "${cM}    was checked out and build the code.${cX}"
fi

if [ ! -e ${oft_dir}/hipft/bin/hipft ]
then
  ${echo} "${cM}==> WARNING! HIPFT does not seem to be built!${cX}"
  ${echo} "${cM}    If you plan to use HipFT, please ensure the submodule${cX}"
  ${echo} "${cM}    was checked out and build the code.${cX}"
fi

echo "==> Appending OFT and its submodules located at [${oft_dir}] to PATH..."

export PATH=${oft_dir}:${oft_dir}/bin:${oft_dir}/magmap/bin:${oft_dir}/conflow/bin:${oft_dir}/hipft/bin:$PATH

${echo} "${cG}==> Done!${cX}"

