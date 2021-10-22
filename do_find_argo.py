# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 17:38:16 2021

@author: gavin
"""

from pathlib import Path
from netCDF4 import Dataset
import numpy.ma as ma
import numpy as np
from cftime import num2date
import pandas as pd

dataDir = Path(r'C:\Users\gavin\Data - not synced\temp\argo_download')
outputDir = Path(r'C:\Users\gavin\OneDrive - Aqualyd Limited\Documents\Aqualyd\Projects\2021-05 SIO ORH survey analysis\Data')

profiles = dataDir.glob('**/*_prof.nc')

first = True
for p in profiles:
    print(p.name)
    rootgrp = Dataset(p, "r")
    if first:
        lats = rootgrp['LATITUDE'][:]
        lons = rootgrp['LONGITUDE'][:]
        units = rootgrp['JULD'].units        
        t = rootgrp['JULD'][:]
        time = num2date(t, units = units, calendar='julian')
        cycle_number = rootgrp['CYCLE_NUMBER'][:]
        filepath = [str(p)] * lats.size
        first = False
    else:
        lats = ma.append(lats, rootgrp['LATITUDE'][:])
        lons = ma.append(lons, rootgrp['LONGITUDE'][:])
        t = rootgrp['JULD'][:]
        time = ma.append(time, num2date(t, units=units, calendar='julian'))
        cycles = rootgrp['CYCLE_NUMBER'][:]
        cycle_number = ma.append(cycle_number, cycles)
        filepath = filepath + [str(p)] * len(cycles)
        
    rootgrp.close()
    
# Now find profiles that are with the months of June to August
start_month = 6
end_month = 8
start_year = 2018
end_year = 2021

months = np.array([t.month for t in time.data])
years = np.array([t.year for t in time.data])

m = (months >= start_month) & (months <= end_month) & (years >= start_year) & (years <= end_year)

# and write to a csv for loading into a GIS
d = pd.DataFrame({'timestamp': time.data[m], 'month': months[m], 'year': years[m], 
                  'lat': lats[m], 'lon': lons[m], 'cycle_number': cycle_number[m], 
                  'filepath': np.array(filepath)[m]})
d.to_csv(outputDir.joinpath('argo_positions.csv'), index=False)
