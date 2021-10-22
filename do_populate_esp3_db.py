# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 17:33:37 2021

@author: gavin
"""

# Code to take the start/end times for identified transects and populate an
# esp3 log_book. database with that data. This is done on a per-directory
# basis.

from pathlib import Path
import pandas as pd
import numpy as np
import shutil
import sqlite3
import datetime as dt
import echopype as ep
from haversine import haversine, Unit
from netCDF4 import Dataset
import gsw
from scipy.stats import hmean

from sw_absorption import sw_absorption

baseDir = Path(r'E:\Aqualyd\SIO_ORH\Data')
surveyDataDir = baseDir.joinpath('Survey_data')

logbookFilename = 'echo_logbook.db'
logbookTemplate = baseDir.joinpath('template_' + logbookFilename)

surveyFilename = 'survey_options.xml'
surveyTemplate = baseDir.joinpath('template_' + surveyFilename)

projectDir = Path(r'C:\Users\gavin\OneDrive - Aqualyd Limited\Documents\Aqualyd\Projects\2021-05 SIO ORH survey analysis')
metaDataFile = projectDir.joinpath(r'Data\Metadata.xlsx')

# Get the survey directories
dataDirs = surveyDataDir.glob('*')
# -NEMA directories are for LSSS only
dataDirs = [p for p in dataDirs if not p.name.endswith('-NMEA')]

# Get the transect metadata
md = pd.read_excel(metaDataFile, sheet_name='Transects').fillna('')
# massage the timestamp format to suit what esp3 (matlab) wants
md.start_time = md.start_time.str.replace('T', ' ')
md.end_time = md.end_time.str.replace('T', ' ')

# Get the Argo CTD cast times and locations
ctds = pd.read_csv(projectDir.joinpath('Data').joinpath('argo_positions.csv'), parse_dates=[0])
# and reformat into a list of lat,lon pairs
ctd_positions = np.array(list(zip(ctds.lat.values, ctds.lon.values)))

for d in dataDirs:
    logbook = d.joinpath(logbookFilename)
    # copy the template database and survey file to the directory
    shutil.copy(logbookTemplate, logbook)
    shutil.copy(surveyTemplate, d.joinpath(surveyFilename))
    
    # find the transect start/end times for this directory
    transects = md.loc[md['data_directory'] == d.name]
    transect_num = transects.transect
    # sort these by feature, then snapshot, then start_time
    transects = transects.sort_values(by=['feature', 'snapshot', 'start_time'])
    
    # find the closest Argos CTD profile for the current dataset.
    # to do that we need a position, so take the first such data from the first
    # file in the directory and use the echopype library to parse the file.
    firstFile = list(d.glob('*.raw'))[0]
    
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
            
    print(f'  File type is {sonar_model} - first datagram type is {first_datagram}')
    ed = ep.open_raw(firstFile, sonar_model=sonar_model)
        
    # Pull out the lat/lon for the current file
    lats = ed.platform.latitude.to_series()
    longs = ed.platform.longitude.to_series()
    
    del ed
    lat = lats.values.mean()
    lon = longs.values.mean()
    
    # Work out the distance and times between the first transect in the directory
    # and all ctd locations
    d = lambda ctd_pos: haversine(ctd_pos, (lat, lon), unit=Unit.KILOMETERS)
    dists = np.array([d(i) for i in ctd_positions]) # [km]
    times = ctds.timestamp - pd.to_datetime(transects.start_time[0])
    time_diffs = times.apply(lambda x: abs(x.total_seconds())/3600) # [hours]
    
    # 1 km of distance has equal weight as 1 hour of time
    spatial_temporal_distance = dists * time_diffs
    i_min = spatial_temporal_distance.idxmin()
    print(f'Closest CTD was {dists[i_min]:.0f} km and {time_diffs[i_min]:.0f} hours away')
    
    # now get the CTD data for the selected cast and give that to esp3.
    ctd_to_use = ctds.iloc[i_min]
    
    rootgrp = Dataset(ctd_to_use.filepath, "r")
    cycle_number = rootgrp['CYCLE_NUMBER'][:]
    i = np.where(cycle_number == ctd_to_use.cycle_number)[0][0] # assume there is always and only one value to find
    pres = rootgrp['PRES'][i]
    temp = rootgrp['TEMP'][i]
    psal = rootgrp['PSAL'][i]
    rootgrp.close()
    
    # Assuming that most of the data is about 1000 m deep, calculate mean sound 
    # speed and absorption between the surface and that depth for use by esp3.
    max_depth = 1000.0 # [m]
    sa = gsw.SA_from_SP(psal, pres, lon, lat)
    ct = gsw.CT_from_t(sa, temp, pres)

    c = gsw.sound_speed(sa, ct, pres)
    rho = gsw.rho(sa, ct, pres)
    alpha = sw_absorption(38.0, psal, temp, pres, 'fandg', 8.0)
    
    depth_mask = pres.data <= max_depth
    
    mean_c = hmean(c[depth_mask])
    mean_abs = np.mean(alpha[depth_mask])
    
    # populate the database
    con = sqlite3.connect(logbook)
    
    # we redo the transect numbering sequentially based on when the snapshot changes
    current_snapshot = -1
    
    for i, t in transects.iterrows():
        feature = t.feature.replace("'", "''")

        if t.snapshot != current_snapshot:
            current_snapshot = t.snapshot
            new_transect = 1
        else:
            new_transect += 1

        if t.start_filename == t.end_filename:
            # transect is entirely within one file
            cmd = ("INSERT INTO logbook (Filename, Snapshot, Stratum, Type, Transect, StartTime, EndTime, Comment) "
                   f"VALUES('{t.start_filename}', {t.snapshot}, '{feature}', 'Acoustic', {new_transect}, '{t.start_time}', '{t.end_time}', '{t.comment}');")
            print(cmd)
            con.execute(cmd)
            con.commit()
        else: # transect covers more than one file
            # get sorted list of files in this survey and work from them...
            rawFiles = sorted(d.glob('*.raw'))
            lastFile = False
            for j, rf in enumerate(rawFiles):
                if rf.name == t.start_filename: 
                    filename = t.start_filename
                    start_time = t.start_time
                    # end time of current file is the start time of the next file, less a little
                    tt = dt.datetime.strptime(rawFiles[j+1].name[-21:], 'D%Y%m%d-T%H%M%S.raw') - dt.timedelta(seconds=1)
                    end_time = tt.strftime('%Y-%m-%d %H:%M:%S.000')
                    withinTransect = True
                elif withinTransect:
                    if rf.name == t.end_filename: # file that contains the end of the transect
                        filename = t.end_filename
                        start_time = end_time # start_time is the end time of previous file plus a second or so
                        end_time = t.end_time
                        lastFile = True
                    else: # a file between the start and end file
                        start_time = end_time # start_time is the end time of the previous file
                        # end time of current file is the start time of the next file, less a little
                        tt = dt.datetime.strptime(rawFiles[j+1].name[-21:], 'D%Y%m%d-T%H%M%S.raw') - dt.timedelta(seconds=1)
                        end_time = tt.strftime('%Y-%m-%d %H:%M:%S.000')
                        filename = rf.name
                if withinTransect:
                    cmd = ("INSERT INTO logbook (Filename, Snapshot, Stratum, Type, Transect, StartTime, EndTime, Comment) "
                          f"VALUES('{filename}', {t.snapshot}, '{feature}', 'Acoustic', {new_transect}, '{start_time}', '{end_time}', '{t.comment}');")
                    print(cmd)
                    con.execute(cmd)
                    con.commit()
                if lastFile:
                    withinTransect = False
    con.close()
    
    # create script file to process the directory
    
