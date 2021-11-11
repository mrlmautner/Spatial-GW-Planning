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

    self.wwtps - list of randomly selected wastewater treatment plants where reuse has been implemented
    self.basins - list of the row and column where each infiltration basin has been implemented
    self.mthlyleak - array with the total quantity of leaks in the system per month in m3
    self.wells - dictionary of well objects input into the MODFLOW model which includes positive flow from wastewater treatment plants, leaks, and recharge basins and negative flow from pumping wells
    self.landuse - dictionary that contains a raster and list of the percentage of land use type (NATURAL, URBAN, or WATER) per model cell for each model phase
'''

# MODFLOW File location
exefile = r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe'

# Parameters
model_name = 'Test_v1'
hydrographloc = 'Test_v1'

# Load alternatives
altpath = Path.cwd() / 'input' / 'decisions' / 'alternatives.csv'
alternatives = {}
altkeys = ['names','wwtps','basins','leakrepair','labels']
with open(altpath) as a:
    lines = a.read().splitlines() 
    for i, line in enumerate(lines):
        templist = line.split(',')
        alternatives[altkeys[i]] = templist

nsamples = 1000
soswrlim = 1000000000
sarun = 2
params = sa.gen_param_vals(nsamples)
safolder = '20211110'
objectives = mf.SA_mode(alternatives=alternatives, params=params, exefile=exefile, safolder=safolder, sarun=sarun, soswrlim=soswrlim, verbose=True)
    
### Assign supply source quantities
#cutz = np.loadtxt(Path.cwd() / 'input' / 'decisions' / 'new_cutz.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Imports from Cutzamala reservoir system
#lerm = np.loadtxt(Path.cwd() / 'input' / 'decisions' / 'new_lerm.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Imports from Lerma groundwater system
#pai = np.loadtxt(Path.cwd() / 'input' / 'decisions' / 'new_pai.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Imports from PAI groundwater system external to the model
#int_sw = np.loadtxt(Path.cwd() / 'input' / 'decisions' / 'new_int_sw.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Surface water sources within the basin
#int_ww = np.loadtxt(Path.cwd() / 'input' / 'decisions' / 'new_int_ww.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Wastewater reuse within the basin
#new_other = cutz + lerm + pai + int_sw + int_ww # Total of all other water supplies except local groundwater (m3/s)
#new_other = new_other.sum(axis=0) # Total of all other supplies (m3/s)
