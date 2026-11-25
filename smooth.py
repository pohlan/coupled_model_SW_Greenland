from pathlib import Path
import os
import geoutils as gu
import numpy as np
from scipy.ndimage import gaussian_filter

# define data paths relative to this script's location
script_dir     = Path(__file__).parent
project_dir    = script_dir.parent
data_dir       = project_dir / "Greenland_data"
bedmachine_dir = data_dir / "bedmachine"
bedmachine_dir = bedmachine_dir.relative_to(project_dir)

def load_and_smooth(name, sig=9):
    min_x = -2.5e5
    min_y = -2.675e6
    max_x = 6e4
    max_y = -2.450e6

    r = gu.Raster("NETCDF:"+str(bedmachine_dir/"BedMachineGreenland-v5.nc")+":"+name)
    r.crop((min_x,min_y,max_x,max_y), inplace=True)

    dat_new = gaussian_filter(r.data.data, sigma=sig)
    r_new = gu.Raster.from_array(data=dat_new, transform=r.transform, crs=r.crs)
    r_new.save(bedmachine_dir/f"BedMachineGreenland-v5_{name}_smooth_sig{sig}.nc")

load_and_smooth("bed")
load_and_smooth("thickness")
