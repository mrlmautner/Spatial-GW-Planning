# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:17:39 2019

@author: MM
"""
import tools.model_functions as mf
import tools.measureobjectives as mo
from pathlib import Path
import numpy as np

'''

'''
scale = ['MUN_CLUSTER-CLIP', 'CLUSTER', 'AGEB']

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

run = 0
folder = '20220323'
#model, sa_loc, error = mf.single_run(model_name=alternatives['names'][0], exefile=exefile, folder=folder, run=run, error_calc=True, verbose=True)
model_error, objectives, objectives_subregion, model_dict = mf.run_alternatives(alternatives, exefile, folder, run=run, verbose=True, scale=scale[run])
#objectives, objectives_subregion, heads = mf.process_objectives(model_name=model_dict['Historical'].name, safolder=safolder, model=model_dict['Historical'], sa_loc='N/A', sarun=sarun, verbose=True)
