# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 12:45:29 2021

@author: gavin
"""

# This script copies data from the per-survey directories into one directory
# Only surveys that meet the quality and time period criteria are copied.
# These copied files are what the pre-processing and echo-integration is done
# on.

from pathlib import Path
import pandas as pd
import shutil
import datetime as dt
from pathvalidate import sanitize_filename
import echopype as ep
from haversine import haversine, Unit
from netCDF4 import Dataset
from netCDF4 import chartostring
import gsw
from scipy.stats import hmean

newDrive = 'E:'

projectDir = Path(r'C:\Users\gavin\OneDrive - Aqualyd Limited\Documents\Aqualyd\Projects\2021-05 SIO ORH survey analysis')

originalDataDir = Path(r'E:\Aqualyd\SIO_ORH\Data')
surveyDataDir = originalDataDir.joinpath('Survey_data')

metaDataFile = projectDir.joinpath(r'Data\Metadata.xlsx')

# Get the transect/survey metadata
md = pd.read_excel(metaDataFile, sheet_name='Transects', parse_dates=['start_time', 'end_time']).fillna('')

# Get the first transect in each survey
t1 = md.loc[md['transect'] == 1]

# work out whether to accept the survey
start_month = t1['start_time'].dt.month
t1 = t1.assign(accepted=(t1['noise'] != 'high') & (t1['motion'] != 'high') & (start_month >= 6) & (start_month <= 8) & (t1['pulse_type'] == 'CW'))

# sort on area, then feature name, then start time, then shapshot
t1 = t1.sort_values(by=['area', 'feature', 'start_time', 'snapshot'])

print(f"There are {t1.shape[0]} potential surveys, {t1['accepted'].sum()} of which are accepted.")

# a dataset with just the accepted surveys
t2 = t1[t1.accepted]

# save this selected set for use by other code
accepted_transects = md[md['survey_id'].isin(t2.survey_id.values)]
accepted_transects.to_csv(metaDataFile.parent.joinpath('Metadata_selected_transects_generated.csv'), index=False)


# for each survey that passes the filter copy the files to a per-year directory and make
# up a LSSS calibration file that contains the appropriate sound speed and absorption, 
# derived from the Argo data.

# Get the Argo CTD cast times and locations
ctds = pd.read_csv(projectDir.joinpath('Data').joinpath('argo_positions.csv'), parse_dates=[0])
# and reformat into a list of lat,lon pairs
ctd_positions = np.array(list(zip(ctds.lat.values, ctds.lon.values)))

calibration_start_time = []
calibration_end_time = []
calibration_sound_speed = []
calibration_alpha = []
calibration_argo_distance = []
calibration_argo_time_period = []


for survey in t2.itertuples():
    filesToCopy = list(surveyDataDir.joinpath(survey.data_directory).glob('*.raw'))
     
    # find the closest Argos CTD profile for the current dataset.
    # to do that we need a position, so take the first such data from the first
    # file in the directory and use the echopype library to parse the file.
    firstFile = filesToCopy[0]
    
    print('\tGetting position from first file in directory.')
    file = open(firstFile, "rb")
    headers = file.read(8)
    file.close()

    first_datagram = headers[4:].decode()

    if first_datagram == 'CON0': # EK/ES60 and ES70 have a CON0 datagram first
        sonar_model = 'EK60'
    elif first_datagram == 'XML0': # ES80/EK80 file have a XML0 datagram first
        sonar_model = 'EK80'
    else:
        print(f'Unknown first datagram type: {first_datagram}')
        sonar_model = 'unknown'
            
    ed = ep.open_raw(firstFile, sonar_model=sonar_model)
        
    # Pull out the lat/lon for the current file and time of first ping
    lats = ed.platform.latitude.to_series()
    longs = ed.platform.longitude.to_series()
    first_ping_time = ed.beam.backscatter_r.ping_time[0].data
    del ed
    
    lat = lats.values.mean()
    lon = longs.values.mean()

    # Get the ping time of the last ping in the last file too
    ed = ep.open_raw(filesToCopy[-1], sonar_model=sonar_model)    
    last_ping_time = ed.beam.backscatter_r.ping_time[-1].data
    del ed
    
    # Work out the distance and times between the first transect in the directory
    # and all ctd locations
    dist = lambda ctd_pos: haversine(ctd_pos, (lat, lon), unit=Unit.KILOMETERS)
    dists = np.array([dist(i) for i in ctd_positions]) # [km]
    times = ctds.timestamp - pd.to_datetime(transects.iloc[0].start_time)
    time_diffs = times.apply(lambda x: abs(x.total_seconds())/3600) # [hours]
    
    dist_weight = 0.1 # weight applied to the distance in km
    time_weight = 1 # weight applied to the time period in hours
    # we'll minimise this sum:
    spatial_temporal_distance = dists*dist_weight + time_diffs*time_weight
    i_min = spatial_temporal_distance.idxmin()
    print(f'\tClosest CTD was {dists[i_min]:.0f} km and {time_diffs[i_min]:.0f} hours away')
    
    # now get the CTD data for the selected cast and calculate water properties.
    ctd_to_use = ctds.iloc[i_min]
    
    rootgrp = Dataset(ctd_to_use.filepath, "r")
    platform_num = chartostring(rootgrp['PLATFORM_NUMBER'][0].data)
    cycle_number = rootgrp['CYCLE_NUMBER'][:]
    i = np.where(cycle_number == ctd_to_use.cycle_number)[0][0] # assume there is always and only one value to find
    pres = rootgrp['PRES'][i]
    temp = rootgrp['TEMP'][i]
    psal = rootgrp['PSAL'][i]
    rootgrp.close()
    
    # Assuming that most of the data is about 1000 m deep, calculate mean sound 
    # speed and absorption between the surface and that depth.
    max_depth = 1000.0 # [m]
    sa = gsw.SA_from_SP(psal, pres, lon, lat)
    ct = gsw.CT_from_t(sa, temp, pres)

    c = gsw.sound_speed(sa, ct, pres)
    rho = gsw.rho(sa, ct, pres)
    alpha = sw_absorption(38.0, psal, temp, pres, 'doonan', 8.0)
    
    depth = pres.data # assuming mBar equals metres...
    depth_mask = depth <= max_depth
    
    mean_c = hmean(c[depth_mask])
    mean_abs = np.mean(alpha[depth_mask])
    mean_t = np.mean(temp)
    mean_s = np.mean(psal)
    
    calibration_start_time.append(str(first_ping_time)[:19]+'Z')
    calibration_end_time.append(str(last_ping_time)[:19]+'Z')
    calibration_sound_speed.append(mean_c)
    calibration_alpha.append(mean_abs)
    calibration_argo_distance.append(dists[i_min])
    calibration_argo_time_period.append(time_diffs[i_min])
    
    # Make up and create, if necessary, the survey directory
    dirName = surveyDataDir.joinpath('AllSurveys')
    dirName.mkdir(parents=True, exist_ok=True)
    
    # Copy the raw files to the new survey directory
    print(f'Copying files from {survey.data_directory} to {dirName.name}')
    # for f in filesToCopy:
    #     to_path = dirName.joinpath(f.name)
    #     shutil.copy(f, to_path)

# Get the calibration gains to use and use them along with the per survey
# sound speed and absorption estimates

# and make up the calibration.xml file for LSSS
xml = '<?xml version="1.0" encoding="UTF-8"?>\n<calibration>\n'
for i_cal in range(len(calibration_start_time)):
    xml = xml + f'   <calibration begin="{calibration_start_time[i_cal]}" end="{calibration_end_time[i_cal]}">\n'
    xml = xml + f'      <case channel="1" g="26.5" SA="1024:0.0, 2048:0.0" abs="{calibration_alpha[i_cal]*0.001:.5f}" sound="{calibration_sound_speed[i_cal]:.1f}"/>\n'
    xml = xml + '   </calibration>\n'
xml = xml + '</calibration>\n'

dirNameNMEA = dirName.parent / (dirName.name + '-NMEA')
with open(dirNameNMEA.joinpath('calibration.xml'), 'w') as cal_file:
    cal_file.write(xml)
    
    
    

