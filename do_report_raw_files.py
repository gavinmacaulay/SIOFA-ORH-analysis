# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 14:15:16 2021

@author: gavin
"""

# Find all .raw files in a directory hirerachy and list them and some metadata
# abou them.

from pathlib import Path
import pandas as pd
import datetime as dt
import re
import echopype as ep
import numpy as np


p_re = re.compile(r'D\d{8}-T\d{6}')

baseDir = Path(r'E:\Aqualyd\SIO_ORH\Data')
#baseDir = Path(r'E:\SIO_ORH\Data')
docsDir = Path(r'C:\Users\gavin\OneDrive - Aqualyd Limited\Documents\Aqualyd\Projects\2021-05 SIO ORH survey analysis\Docs')
dataDir = baseDir.joinpath('WW ES80 2018-2021')

dataDir = baseDir.joinpath('2005-17 selected')

# Choose specific directories for testing and re-running
#dataDir = Path(r'D:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 78\AOS\5m Calibrations 17.8.18')
#dataDir = Path(r'D:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 81 survey')
#dataDir = dataDir.joinpath(r'Trip 77/All trip')
dataDir = Path(r'C:\Users\gavin\Data - not synced\temp\Files with errors\NMEA checksum failed')



# The name of the file that the positions get saved to.
posFileName = 'pos.csv'
posAllFileName = 'pos_all.csv'
posPathFileName = 'pos_paths.csv'

rawFiles = dataDir.rglob('*.raw')
rawFiles = sorted(rawFiles)

print(f'There are {len(rawFiles)} .raw files.')

parentDirs = []
total_rf_dirs = 0
first_file = []
last_file = []
duration = []
n = []

pos_dfs = []

p = []
for f in rawFiles:
    p.append(f.parent)
parentDirs = np.unique(p)

# Directories to skip cause they are large or cause unresolved errors
ignoreDirs = [dataDir.joinpath(r'Trip 78/AOS/Calibration Trawl 18.08.18/Acoustic'),
              dataDir.joinpath(r'Trip 83/Calibration')]
#ignoreDirs = []

# New list of directories without the ignored ones
parentDirsTrimmed = [p for p in parentDirs if p not in ignoreDirs]

for p_i in range(len(parentDirsTrimmed)):
    
    p = parentDirsTrimmed[p_i]
    rf = sorted(p.glob('*.raw'))
    
    m = p_re.search(rf[0].name)
    t = dt.datetime.strptime(m.group(), 'D%Y%m%d-T%H%M%S')
    first_file.append(t)
    
    m = p_re.search(rf[-1].name)
    t = dt.datetime.strptime(m.group(), 'D%Y%m%d-T%H%M%S')
    last_file.append(t)
    
    duration.append((last_file[-1] - first_file[-1]) / dt.timedelta(hours=1))
    n.append(len(rf))
    
    # Load the raw file and pull out the lat/lon and store that. But need to know 
    # the file type for the echopype library, so work that out from the first
    # byte in the file.
    lats = pd.Series(None, dtype=np.float64)
    longs = pd.Series(None, dtype=np.float64)

    # A directory of files
    for rawfile in rf:
        print(f'Reading {rawfile}')

        file = open(rawfile, "rb")
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
        ed = ep.open_raw(rawfile, sonar_model=sonar_model)
        
        # Pull out the lat/lon for the current file
        lats = lats.append(ed.platform.latitude.to_series())
        longs = longs.append(ed.platform.longitude.to_series())
        
        del ed
    # end of one directory of files
    # save the positions to a csv file
    pos_df = lats.to_frame(name='lat').join(longs.to_frame(name='lon'))
    pos_df.index.name ='timestamp'
    pos_df['directory'] = str(p)
    
    # save GPS positions to a per directory file
    pos_df.to_csv(p.joinpath(posFileName), index=True)

##############################################################################
# now concatenate all the pos.csv files into one file to load into a GIS. Best
# suited for files that are small enough to load into the computer's memory.

fileSet1 = list(dataDir.rglob(posFileName))
# split into smaller sections (too large for qgis in the end)
fileSet2 = Path(r'D:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 83\All trip\pos.csv')
fileSet3 = Path(r'D:\Aqualyd\SIO_ORH\Data\WW ES80 2018-2021\Trip 81b\pos.csv')
if fileSet2 in fileSet1:
    fileSet1.remove(fileSet2)
if fileSet3 in fileSet1:
    fileSet1.remove(fileSet3)
    
posFileSet = [(fileSet1, 'pos1.csv')]
posFileSet.append(([fileSet2], 'pos2.csv'))
posFileSet.append(([fileSet3], 'pos3.csv'))

for s in posFileSet:
    if s[0][0].exists():
        with open(dataDir.joinpath(s[1]), 'w') as outfile:
            for fname in s[0]:
                with open(fname) as infile:
                    print(f'Adding {fname} to {s[1]}')
                    outfile.write(infile.read())

# do this again, but subsample the points (one point for every minute)
posFiles = dataDir.rglob(posFileName)
header = True
with open(dataDir.joinpath('points_subsampled.csv'), 'w') as outfile:
    for fname in posFiles:
        print(f'Adding {fname}')
        pos_df = pd.read_csv(fname)
        if len(pos_df) > 1000:
            pos_df = pos_df.iloc[::60,:]
        pos_df.to_csv(fname.parent.joinpath('pos_subsampled.csv'), index=False)
        pos_df.to_csv(outfile, mode='a', header=header, index=False)
        if header: # only first to_csv provides a header
            header = False

# problem here is that some files have too many points for qgis to have in 
# a single wkt linstring (seems to be a limit of 4 million charaters mer wkt)
# posFiles = dataDir.rglob(posFileName)
# with open(dataDir.joinpath(posPathFileName), 'w') as outfile:
#     outfile.write('directory;start_time;end_time;wkt\n')
#     for fname in posFiles:
#         print(f'Adding {fname}')
#         pos_df = pd.read_csv(fname)
        
#         outfile.write(str(fname.parent) + ';')
#         outfile.write(pos_df.iloc[0,:].timestamp + ';')
#         outfile.write(pos_df.iloc[-2,:].timestamp + ';') # sometimes last result is nan, so use -2
#         outfile.write('LINESTRING(')
#         with open(fname) as infile:
#             for index, row in pos_df.iterrows():
#                 outfile.write(f'{pos_df.lon[index]:.4f} {pos_df.lat[index]:.4f},')
#         outfile.write(')\n')
        

# Now concatenate the dataframes from each directory into one dataframe and then save
# that to csv for import into a GIS.
#pos_all_df = pd.concat(pos_dfs)
#pos_all_df.to_csv(baseDir.joinpath('pos.csv'), index=True)


#df = pd.DataFrame(list(zip(duration, n, first_file, last_file, parentDirsTrimmed)), columns=['duration', 'num', 'first', 'last', 'dir'])
#df.sort_values(by=['first'], inplace=True)
##df.to_csv(baseDir.joinpath('summary.csv'), index=False)
#print(df)


###########################################################################
# do a table for the ToR 1 report

if False:
    metadataFile = docsDir.joinpath('Metadata.xlsx')
    md = pd.read_excel(metadataFile, sheet_name='Activities')
    md.sort_values(inplace=True, by=['ORH area', 'Timestamp at start'], axis='index')
    
    duration = []
    date_start = []
    for index, row in md.iterrows():
        if isinstance(row['Timestamp at end'], str) and isinstance(row['Timestamp at start'], str):
            endt = dt.datetime.strptime(row['Timestamp at end'], '%Y-%m-%dT%H:%M:%S')
            startt = dt.datetime.strptime(row['Timestamp at start'], '%Y-%m-%dT%H:%M:%S')
            duration.append((endt - startt) / dt.timedelta(hours=1))
            date_start.append(startt.date())
        else:
            duration.append(dt.timedelta(seconds=0))
            date_start.append(None)
            
    md['duration'] = duration
    md['Date start'] = date_start
    
    md.to_csv(docsDir.joinpath('table.csv'), index=False, 
              columns=['ORH area', 'Date start', 'duration', 'No. transects', 'Explicit survey'])


