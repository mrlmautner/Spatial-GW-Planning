# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:05:02 2022

@author: MM
"""

import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import tools.plot_results as plot

folder = '20220316'
scale_names = ['MUN_CLUSTER-CLIP', 'CLUSTER', 'AGEB']

obj_short = ['A','W','F']
obj_long = ['GW Availability','Water Quality','Urban Flooding']
obj_label = ['$\Delta$ GW Availability\n[L/capita]','$\Delta$ Water Quality\n[% cells]','$\Delta$ Urban Flooding\n[% cells]']

# Load social indicators
indpath = Path.cwd() / 'input' / 'socialindicators'
ind_array = np.loadtxt(indpath.joinpath('INDICATORS-WTD-BY-POP_GRID_2010.csv'), skiprows=1, delimiter=',', usecols=(0,1,2,3,4,5,6,7,8))
ind_short = ['ID', 'edu', 'wat', 'hlth', 'pov', 'percentile_edu', 'percentile_wat', 'percentile_hlth', 'percentile_pov']
ind_long = ['Subregion ID', 'Education Level', 'Houshold Water Access', 'Health Insurance', 'Above Poverty Line', 'Education (percentile)', 'Water (percentile)', 'Health (percentile)', 'Poverty (percentile)']
ind_label = ['Grade', '% Households', '% Population', '% Population']

# Load alternatives
altpath = Path.cwd() / 'input' / 'decisions' / 'alternatives.csv'
alternatives = {}
altkeys = ['names','c1','c2','c3','c4','p1','p2','p3','labels']
with open(altpath) as a:
    lines = a.read().splitlines() 
    for i, line in enumerate(lines):
        templist = line.split(',')
        alternatives[altkeys[i]] = templist

alt_short = alternatives['names']
alt_long = alternatives['labels']
#
#plot.plt_IndvObj_ind(folder, obj_short, obj_long, obj_label, alt_short, alt_long, ind_short, ind_long, ind_label)
#plot.plt_IndvObj_scale(folder, obj_short, obj_long, obj_label, alt_short, alt_long, ind_short, ind_long, ind_label)

plot.plot_nondom_heatmap(folder, alt_short, scale_names, obj_short, ind_short, ind_long)
