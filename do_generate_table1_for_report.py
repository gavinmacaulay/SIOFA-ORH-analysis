# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 14:36:09 2021

@author: gavin
"""

# Produce a Word table for the ToR 2 report that lists all assessed surveys

from pathlib import Path
import pandas as pd
from docx import Document


projectDir = Path(r'C:\Users\gavin\OneDrive - Aqualyd Limited\Documents\Aqualyd\Projects\2021-05 SIO ORH survey analysis')
metaDataFile = projectDir.joinpath(r'Data\Metadata.xlsx')

# Get the transect/survey metadata
md = pd.read_excel(metaDataFile, sheet_name='Transects', parse_dates=['start_time', 'end_time']).fillna('')

# Pick out all rows that have a non-blank transect_type column. This is 
# always transect 1 of a survey.
t1 = md.loc[md['survey_type'] != '']

# work out whether to accept the survey
start_month = t1['start_time'].dt.month
t1 = t1.assign(accepted=(t1['noise'] != 'high') & (t1['motion'] != 'high') & (start_month >= 6) & (start_month <= 8) & (t1['pulse_type'] == 'CW'))

# add a column with just the date and with with just the year
t1 = t1.assign(start_date = t1['start_time'].dt.date, year = t1['start_time'].dt.year)

# sort on area, then feature name, then start time, then shapshot
t1 = t1.sort_values(by=['area', 'feature', 'start_time', 'snapshot'])

# remove columns and change the order of the remainder in preparation for generating the Word table
t1_word = t1[['area', 'feature', 'start_date', 'snapshot', 'survey_type', 'noise', 'motion', 'pulse_type', 'accepted']]

print(f"There are {len(t1['accepted'])} potential surveys, {t1['accepted'].sum()} of which are accepted.")

# create a summary dataset also showing the number of surveys within each area by year
t1_summary = t1.groupby(by=['area', 'year', 'accepted'], sort=True).count()

# and output to a Word table
tableFile = projectDir.joinpath('docs').joinpath('ToR 2 report').joinpath('ToR 2 table 1.docx')
    
document = Document()
table = document.add_table(rows=t1_word.shape[0]+1, cols=t1_word.shape[1])
for i, column in enumerate(t1_word):
    prev_cell_str = ''
    table.cell(0, i).text = column.replace('_', ' ').capitalize()
    for row in range(t1_word.shape[0]):
        cell_str = str(t1_word[column].iloc[row])
        if (column == 'area'): # drop area name when it is in multiple consecutive rows
            if (cell_str == prev_cell_str):
                cell_str = ''
            else:
                prev_cell_str = cell_str

        if (column == 'feature'): # drop feature name when it is in multiple consecutive rows
            if (cell_str == prev_cell_str):
                cell_str = ''
            else:
                prev_cell_str = cell_str
            
        table.cell(row+1, i).text = cell_str
document.save(tableFile)




