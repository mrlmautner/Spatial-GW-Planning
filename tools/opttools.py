# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:49:58 2019

@author: MM
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors as mcolors
import seaborn as sns

sns.set_theme()

def nondom_sort(solutions):
    objs = solutions.copy()
    num_solutions = len(objs)

    # use a boolean index to keep track of nondominated solns
    keep = np.zeros(num_solutions, dtype = bool)

    for i in range(num_solutions):
        for j in range(num_solutions):
            a = objs[i]
            b = objs[j]
            if np.all(a <= b) & np.any(a < b):
                keep[i] = True
                keep[j] = False

    return keep

def parallel_axis(data, figname, color_indicator, label_indicator, i_min, i_max):
    
    i_colormap = cm.get_cmap('viridis', 1000)
    color_indicator_range = ((color_indicator - i_min) / (i_max - i_min))*1000
    color_indicator_range[color_indicator_range > 1000] = 1000
    color_indicator_range[color_indicator_range < 1] = 1
    
    # first normalize each objective to [0,1] range
    data = (data - np.nanmin(data,axis=0)) / (np.nanmax(data,axis=0) - np.nanmin(data,axis=0))
    
    plt.figure()
    for i, objs in enumerate(data):
#        # filtering
#        if objs[2] < 0.5:
#            plt.plot(range(4), objs, color='steelblue', zorder=2)
#        else:
#            plt.plot(range(4), objs, color='0.8', zorder=1)
        plt.plot(range(3), objs, color=i_colormap(int(color_indicator_range[i]-1)), zorder=int(color_indicator_range[i]-1), alpha=0.2)
    
    cmap = i_colormap
    norm = mcolors.Normalize(vmin=i_min, vmax=i_max)
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), orientation='vertical', label=label_indicator)
    
    plt.gca().set_xticks(range(3))
    plt.gca().set_xticklabels(['Energy','Quality','Flooding'])

    plt.savefig(figname, dpi=600)
    plt.close()

alternatives = ['hist','wwtp','basin','leak']
indicator = 19
for alt in range(4):
    alt_text = alternatives[alt]
    for f in ['ageb','cluster','cluster-ed','mun','mun-clip']:
        
        raw_data = np.loadtxt('G:/My Drive/Davis/Research/Papers/Subregional_Objectives_Management/Data/objectives/objective_'+f+'_joh-params.csv',delimiter=',',skiprows=1)
        
        raw_data = raw_data[~np.isnan(raw_data[:,indicator]), :]
        objectives = np.append(np.append(raw_data[:,13+alt],raw_data[:,5+alt],axis=0),raw_data[:,9+alt],axis=0).reshape((raw_data.shape[0],3), order='F')
        
        save_file = 'G:/My Drive/Davis/Research/Papers/Subregional_Objectives_Management/Data/indicator vs objective/Figures/'+alt_text+'-'+ f +'-education.png'
        i_label = 'Average highest educational grade\nof adult population'
        # [Mun, E-HIST, E-WWTP, E-Basin, E-Leak, NE-HIST, NE-WWTP, NE-Basin, NE-Leak, W-HIST, W-WWTP, W-Basin, W-Leak, F-HIST, F-WWTP, F-Basin, F-Leak, pobtotal, VPH_AGUADV_% (18), GRAPROES_C (19), PDER_SS_% (20), ARR_LBIEN_% (21), ARR_LBIENMIN_%]
#        parallel_axis(data=objectives, color_indicator=raw_data[:,indicator], figname=save_file, label_indicator=i_label, i_min=np.nanmin(raw_data[:,indicator],axis=0), i_max=np.nanmax(raw_data[:,indicator],axis=0))
#        parallel_axis(data=objectives, color_indicator=raw_data[:,indicator], figname=save_file, label_indicator=i_label, i_min=0, i_max=np.nanmax(raw_data[:,indicator],axis=0))
        parallel_axis(data=objectives, color_indicator=raw_data[:,indicator], figname=save_file, label_indicator=i_label, i_min=6, i_max=14)
