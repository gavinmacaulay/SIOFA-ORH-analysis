# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 12:45:29 2021

@author: gavin
"""

# This script copies raw files from their original directories into a separate
# location, into directories named after the feature and date of the survey.
# This is done to provide a simplier directory hierarchy for subsequent 
# processing.

from pathlib import Path
import pandas as pd
import shutil
import datetime as dt
from pathvalidate import sanitize_filename

newDrive = 'E:'

projectDir = Path(r'C:\Users\gavin\OneDrive - Aqualyd Limited\Documents\Aqualyd\Projects\2021-05 SIO ORH survey analysis')

originalDataDir = Path(r'E:\Aqualyd\SIO_ORH\Data')
surveyDataDir = originalDataDir.joinpath('Survey_data')

metaDataFile = projectDir.joinpath(r'Data\Metadata.xlsx')

# read in the metadata that tells us what time periods apply to surveys
md = pd.read_excel(metaDataFile, sheet_name='Activities')

# for each row in the metadata, create a suitable directory then copy raw files
# into it, based on the start and end times for that row in the metadata
for row in md.itertuples():
    if not pd.isna(row.Directory):
        name = sanitize_filename(f'{row.Start_date}_{row.Feature}')
        
        # Make up and create, if necessary, the survey directory
        dirName = surveyDataDir.joinpath(name)
        dirName.mkdir(parents=True, exist_ok=True)

        # Get the survey start/end times as datatime variables        
        start_t = dt.datetime.strptime(row.Start_timestamp[0:19], '%Y-%m-%dT%H:%M:%S')
        end_t = dt.datetime.strptime(row.End_timestamp[0:19], '%Y-%m-%dT%H:%M:%S')
        
        # Get a list of all .raw files in the directory for the current row
        rawDir = newDrive + row.Directory[2:]
        rawFiles = list(Path(rawDir).glob('*.raw'))
        fileStartTime = list([])

        # Get the start time for each .raw file        
        for f in rawFiles:
            t = dt.datetime.strptime(f.stem[-17:], 'D%Y%m%d-T%H%M%S')
            fileStartTime.append(t)
        fileEndTime = fileStartTime[1:]
        fileEndTime.append(fileEndTime[-1] + dt.timedelta(hours=10)) # 10 is a guess
        # and put these into a DataFrame
        df = pd.DataFrame(data={'file': rawFiles, 'time_start': fileStartTime, 'time_end': fileEndTime}) 
        
        # Awkward code to select first file that includes timestamp start_t and the 
        # last file that includes end_t, but also cope with start_t being before
        # the first file and end_t being after the last file.
        firstFileIn = (start_t >= df.time_start) & (start_t <= df.time_end)
        if ((~firstFileIn).all()) & (start_t <= df.time_start.iloc[0]): # before start of first file
            firstFileIn.iloc[0] = True
        
        lastFileIn = (end_t >= df.time_start) & (end_t <= df.time_end)
        if ((~lastFileIn).all()) & (end_t >= df.time_end.iloc[-1]): # after end of last file
            lastFileIn.iloc[-1] = True

        firstFile = df.loc[firstFileIn]
        lastFile = df.loc[lastFileIn]
            
        filesToCopy = df.iloc[firstFile.index[0]:lastFile.index[0]+1]
        # and then copy the selected files to the new survey directory
        print(f'Copying {len(filesToCopy)} files from {f.parent.name} to {dirName.name}')
        for f in filesToCopy.itertuples():
            to_path = dirName.joinpath(f.file.name)
            shutil.copy(f.file, to_path)
        
    
