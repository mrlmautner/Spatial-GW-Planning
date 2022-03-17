# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:17:39 2019

@author: MM
"""
import tools.model_functions as mf
import tools.satools as sa
from pathlib import Path
import numpy as np

'''
The model run will create an object with the following results:

    self.mthlyleak - array with the total quantity of leaks in the system per month in m3
    self.wells - dictionary of well objects input into the MODFLOW model which includes positive flow from wastewater treatment plants, leaks, and recharge basins and negative flow from pumping wells
    self.landuse - dictionary that contains a raster and list of the percentage of land use type (NATURAL, URBAN, or WATER) per model cell for each model phase
'''

# MODFLOW File location
exefile = r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe'

# Load alternatives
altpath = Path.cwd() / 'input' / 'decisions' / 'alternatives.csv'
alternatives = {}
altkeys = ['names','c1','c2','c3','c4','p1','p2','p3','labels']
with open(altpath) as a:
    lines = a.read().splitlines() 
    for i, line in enumerate(lines):
        templist = line.split(',')
        alternatives[altkeys[i]] = templist

sarun = 1
safolder = '20220316'
#model, sa_loc, error = mf.single_run(model_name=alternatives['names'][0], exefile=exefile, safolder=safolder, sarun=sarun, error_calc=True, verbose=True)
model_error, objectives, objectives_subregion, model_dict = mf.run_alternatives(alternatives, exefile, safolder, sarun=sarun, verbose=True)
#objectives, objectives_subregion, heads = mf.process_objectives(model_name=model.name, safolder=safolder, model=model, sa_loc=sa_loc, sarun=sarun, verbose=True)