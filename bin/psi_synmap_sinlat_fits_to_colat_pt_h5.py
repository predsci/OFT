#!/usr/bin/env python3
import argparse
import sys
import os
import numpy as np
from astropy.io import fits

import psi_io  # pip install psi-io

def main():
    parser = argparse.ArgumentParser(
        description="Convert a Synoptic map in a sinlatitude-longitude FITS file into a PSI grid colatitude-longitude map in H5 format ready for use with psi_map_prep.py."
    )
    parser.add_argument("fitsfile", help="Input FITS file")
    parser.add_argument("-hdffile", help="Output H5 file name",required=False)
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Open and find the image data (supports primary HDU or first extension)
    # ------------------------------------------------------------------
    try:
        verify_mode = "exception"
        with fits.open(args.fitsfile, mode="readonly", checksum=verify_mode) as hdul:
            # Primary HDU has data?
            if hdul[0].data is not None and hdul[0].data.ndim >= 2:
                img_hdu = hdul[0]
            # Otherwise look for the first image extension (common for SDO/AIA etc.)
            elif len(hdul) > 1 and hdul[1].data is not None and hdul[1].data.ndim >= 2:
                img_hdu = hdul[1]
            else:
                raise ValueError("No 2-D image data found in the FITS file.")

            data = img_hdu.data.squeeze()          # remove degenerate dimensions
            if data.ndim != 2:
                raise ValueError(f"Found {data.ndim}-D data; only 2-D images are supported.")

            ny, nx = data.shape                        # Fortran order in the original code: (nx,ny)

    except Exception as e:
        print("\n### ERROR: Could not read image from FITS file")
        print(f"    File : {args.fitsfile}")
        print(f"    Reason: {e}")
        sys.exit(1)

    # ------------------------------------------------------------------
    # BITPIX <= 32  →  32-bit float in the original Fortran code
    # ------------------------------------------------------------------
    bitpix = abs(img_hdu.header.get("BITPIX", 32))  # Default to 32 if missing
    dtype = np.float32 if bitpix <= 32 else np.float64
    data = data.astype(dtype)

    # Replace NaNs with zero (same behaviour as the original FITSIO nullval handling)
    if np.any(np.isnan(data)):
        print("### WARNING: Image contains NaNs → replaced with 0.0")
        data = np.nan_to_num(data, nan=0.0)

    # ------------------------------------------------------------------
    # 2. Create scales and flip data to match 0->pi colat direction.
    #    Assume Sinlat-long on a "inner half-half" mesh
    # ------------------------------------------------------------------

    dx = (np.pi*2)/nx
    x = np.linspace(dx/2,2*np.pi-dx/2,nx,dtype=dtype)

    dy = 2/ny
    y = np.linspace(-1+dy/2,1-dy/2,ny,dtype=dtype)

    y = np.arccos(y)
    y = np.flipud(y)

    data = np.flipud(data)
    
    # ------------------------------------------------------------------
    # 3. Add edge points to data and scales to make 
    #    the map on the PSI "main-main" grid.
    # ------------------------------------------------------------------
    
    # Take average around pole:
    pole_n = np.sum(data[0,:])/nx
    pole_s = np.sum(data[ny-1,:])/nx
   
    #Extend data and scales:
    
    data = np.pad(data, ((1,1),(0,0)), mode='constant', constant_values=(pole_n,pole_s))
    y = np.pad(y, (1,1), mode='constant', constant_values=(0,np.pi))
    ny = ny+2
    
    data = np.pad(data, ((0,0),(1,1)), mode='constant', constant_values=-9999)
    x = np.pad(x, (1,1), mode='constant', constant_values=(0,2*np.pi))
    nx = nx+2

    data[:,0]    = 0.5*(data[:,1]+data[:,nx-2])
    data[:,nx-1] = data[:,0]

    # ------------------------------------------------------------------
    # 4. Write out map with non-uniform colat-long scales 
    # ------------------------------------------------------------------
    if not args.hdffile:
      root, _ = os.path.splitext(args.fitsfile)
      hdffile=root + '.h5'
    else:
      hdffile=args.hdffile
    try:
        psi_io.wrhdf_2d(hdffile, x, y, data)
        print(f"Successfully wrote map to {hdffile}")
    except Exception as e:
        print("\n### ERROR: Failed to write HDF file")
        print(f"    File : {hdffile}")
        print(f"    Reason: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()