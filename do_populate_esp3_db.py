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
import shutil
import sqlite3
import datetime as dt

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
ctds = pd.read_csv(projectDir.joinpath('Data').joinpath('argo_positions.csv'))

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
    spatial_temporal_distance = 
    
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
    
