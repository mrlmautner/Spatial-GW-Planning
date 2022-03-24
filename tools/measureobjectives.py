# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:59:47 2018

@author: MM
"""
import numpy as np
import pandas as pd
import math
import calendar
from pathlib import Path

## OLD
def measureEnergy(heads,supply_dict,dem,bottom,subregion_list):
    E = 0
    E_subregion = np.zeros(subregion_list.shape[0])
    
    efficiency = 0.7
    MJc = 9.81 # MegaJoules to lift 1 MegaLiter 1 meter
    kWhc = 3.6 # MegaJoules per kWh
    
#    pd.DatetimeIndex(freq='M',start='01/31/1986',end='12/31/2013')
    coords = np.zeros(2)
    
    for i, p in supply_dict.items():
        # Start in year 2
        if i > 23:
            d = calendar.monthrange(1984+math.floor(i/12),(i%12)+1)[1] # number of days in stress period
            h = heads.get_data(kstpkper=(8,i),mflay=1) # heads binary file for last time step in stress period
            
            # loop through all well values
            for n in p:
                # Only measure energy use for pumping wells (negative flow)
                if n[3] < 0:
                    coords[0] = n[1]
                    coords[1] = n[2]
                    wellHead = max(h[int(coords[0]),int(coords[1])], bottom[1][int(coords[0]),int(coords[1])])
                    surface = dem[int(coords[0]),int(coords[1])]
                    depthtowater = max(surface-wellHead,1)
                    pumping = n[3] * (-1) * 0.001 # m3/d --> ML/d
                    E += (efficiency * depthtowater * pumping * MJc / kWhc)*d # kWh per month (stress period) of pumping
                    
                    if n[4] > 0:
                        current_subregion = int(np.where(subregion_list == n[4])[0])
                    E_subregion[current_subregion] += (efficiency * depthtowater * pumping * MJc / kWhc)*d # kWh per month (stress period) of pumping
    
    return E, E_subregion

def measureWaterQuality(heads,dem,active,bottom,subregion_array,subregion_list):    
    cells_clay = np.sum(active[0]) # Number of cells in Clay layer
    cells_clay_subregion = np.zeros(subregion_list.shape[0]) # Number of cells in Clay layer in each subregion
    for m, current_subregion in enumerate(subregion_list):
        cells_clay_subregion[m] = np.sum(np.multiply(active[0], (subregion_array == current_subregion).astype(float)))
    
#    wqual = np.zeros(360) # total difference between bottom of clay layer and head in layer 2
    h = list(np.zeros(360)) # head map for each stress period
    cells_below = np.zeros(360) # Number of cells with head below bottom of clay layer
    cells_below_subregion = np.zeros((360, subregion_list.shape[0]))
    
    # loop through all stress periods
    for t in range(348,360):
        # check during the last time step of each stress period
        h2 = np.multiply(active[0],heads.get_data(kstpkper=(8,t),mflay=1)) # heads binary file layer 2
        
        b0 = np.multiply(active[0],bottom[0])
        b1 = np.multiply(active[0],bottom[1])
        
        h[t] = np.multiply(np.less(h2,b1),b1) + np.multiply(np.multiply(np.less(h2,b0), np.greater(h2,b1)),h2) + np.multiply(np.greater(h2,b0),b0) # Matrix containing bottom of layer 2 if layer 2 is dry, plus the head in layer 2 if less than the bottom of layer 1, plus the bottom of layer 1 if the head is greater than layer 1
        
        h_lessthan_b = np.multiply(np.less(h[t],b0), active[0])
        cells_below[t] = np.sum(h_lessthan_b)
                
        for m, current_subregion in enumerate(subregion_list):
            cells_below_subregion[t,m] = np.sum(np.multiply(h_lessthan_b, (subregion_array == current_subregion).astype(float)))

#        wq_temp = b0 - h[t] # Subtract the head matrix from the bottom of the active clay layer
        
#        wqual[t] = np.sum(wq_temp) # Sum the total depth to groundwater below the bottom of layer 1 in layer 2 over all cells

    annual_wqclaycells = np.sum(cells_below)/(cells_clay*12)
    annual_wqclaycells_subregion = np.sum(cells_below_subregion, axis=0)/(cells_clay_subregion*12)
            
    return annual_wqclaycells, annual_wqclaycells_subregion, h

## NEW
def measureAvailability(supply_dict,subregion_array,subregion_list):
    # Load population information
    modelGridIndicators = np.loadtxt(Path.cwd() / 'input' / 'socialindicators' / 'INDICATORS-WTD-BY-POP_GRID_2010.csv', delimiter=',', skiprows=1, usecols=[0,1,2,7])
    modelGridIndicators[:,1] = np.array([int(yi) for yi in modelGridIndicators[:,1] - 1]) # zero-index model grid column
    modelGridIndicators[:,2] = np.array([int(xi) for xi in modelGridIndicators[:,2] - 1])
    
    modelGridIndicators = np.c_[ modelGridIndicators, np.zeros((modelGridIndicators.shape[0],1)) ] # add a column that holds population within 10 grid cells array
    for row in modelGridIndicators:
        row[0] = subregion_array[int(row[1]),int(row[2])] # reassign first column to subregion
        row[4] = np.sum(modelGridIndicators[np.where((abs(modelGridIndicators[:,1]-row[1])<=10) & (abs(modelGridIndicators[:,2]-row[2])<=10)),3])
    
    modelGridIndicators = np.c_[ modelGridIndicators, np.zeros((modelGridIndicators.shape[0],1)) ] # add a column that holds cummulative monthly pumping within 10 grid cells array
    for i, p in supply_dict.items():
        # Start in year 2
        if i > 23:
            d = calendar.monthrange(1984+math.floor(i/12),(i%12)+1)[1] # number of days in stress period
            
            # loop through all well values
            for n in p:
                # Only measure energy use for pumping wells (negative flow)
                if n[3] < 0:
                    pumping = n[3] * (-1) * d # m3/month
                    modelGridIndicators[np.where((abs(modelGridIndicators[:,1]-n[1])<=10) & (abs(modelGridIndicators[:,2]-n[2])<=10)),5] += pumping
                    
    modelGridIndicators = np.c_[ modelGridIndicators, np.zeros((modelGridIndicators.shape[0],1)) ] # add a column that holds cummulative pumping divided by population within 10 grid cells array
    modelGridIndicators[:,6] = modelGridIndicators[:,5] / modelGridIndicators[:,4]
    
    A = np.nanmean(modelGridIndicators[:,6])
    
    A_subregion = np.zeros(subregion_list.shape[0]) # initialize subregion population list
    for s, subregion in enumerate(subregion_list):
         A_subregion_temp = modelGridIndicators[np.where(modelGridIndicators[:,0] == subregion), 6]
         A_subregion[s] = np.nanmean(A_subregion_temp)
    
    return A, A_subregion

def measurAvgDtoGW(heads,dem,active,bottom,subregion_array,subregion_list):
    dem = np.multiply(active[1],dem) # DEM only for active cells
    
    cells_active = np.sum(active[1]) # Number of cells in active area
    cells_active_subregion = np.zeros(subregion_list.shape[0]) # Number of active cells in each subregion
    for m, current_subregion in enumerate(subregion_list):
        cells_active_subregion[m] = np.sum(np.multiply(active[1],(subregion_array == current_subregion).astype(float)))
    
    d2gw = [np.zeros(12), np.zeros(12)] # depth to groundwater map for each stress period
    d2gw_subregion = [np.zeros((12, subregion_list.shape[0])), np.zeros((12, subregion_list.shape[0]))]
    
    third_last_per = [[24,36],[348,360]]
    
    for a in range(2):
        # loop through third or last year of stress periods
        for t in range(third_last_per[a][0],third_last_per[a][1]):
            # check during the last time step of each stress period
            h2 = np.multiply(active[1],heads.get_data(kstpkper=(8,t),mflay=1)) # heads binary file layer 2 for active area
            h2[h2==-0] = 0
            h2 = np.maximum(h2, np.multiply(active[1],bottom[1])) # set minimum head to bottom of layer 2
            
            d2gw_temp = dem - h2 # calculate depth to groundwater
            d2gw_temp = np.multiply(active[1],np.maximum(d2gw_temp, d2gw_temp*0)) # set negative values to 0
                   
            for m, current_subregion in enumerate(subregion_list):
                d2gw_subregion[a][t-third_last_per[a][0],m] = np.sum(np.multiply(d2gw_temp, (subregion_array == current_subregion).astype(float))) / cells_active_subregion[m]
    
            d2gw[a][t-third_last_per[a][0]] = np.sum(d2gw_temp) / cells_active # Sum the total depth to groundwater over all cells averaged over total active cells
    
    delta_average_d2gw = np.mean(d2gw[1]) - np.mean(d2gw[0])
    delta_average_d2gw_subregion = np.mean(d2gw_subregion[1], axis=0) - np.mean(d2gw_subregion[0], axis=0)
            
    return delta_average_d2gw, delta_average_d2gw_subregion, d2gw, d2gw_subregion

def measureMound(heads,dem,active,LU,subregion_array,subregion_list,PhasePer):
#    mound = 0 # cumulative head above DEM during model period
    urban_cells = 0
    mound_cells = 0
    
    urban_cells_subregion = np.zeros(subregion_list.shape[0])
    mound_cells_subregion = np.zeros(subregion_list.shape[0])
    
    hmatrix = {}
    hmatrix['h1'] = np.zeros((12,dem.shape[0],dem.shape[1]))
    hmatrix['h2'] = np.zeros((12,dem.shape[0],dem.shape[1]))
    i = 0
    
    for t in range(348,360):
        h1 = heads.get_data(kstpkper=(8,t),mflay=0) # heads binary file for last time step of stress period
        h1[h1 < -1] = -1
        hmatrix['h1'][i] = np.round(h1.astype(np.float32),3)
        h2 = heads.get_data(kstpkper=(8,t),mflay=1) # heads binary file for last time step of stress period
        h2[h2 < -1] = -1
        hmatrix['h2'][i] = np.round(h2.astype(np.float32),3)
        
        if t < PhasePer[0]:
            LUset = '1990'
        elif t < PhasePer[1]:
            LUset = '2000'
        else:
            LUset = '2010'
        
        lu_temp = np.multiply(LU[LUset]['ARRAY']['URBAN'],active[1])
        
        lu_active = np.sum(lu_temp) # Total urban area of active cells
        urban_cells += lu_active
                
        h = np.multiply(active[0],h1) + np.multiply((active[1]-active[0]),h2) # Choose topmost active layer for each cell
        above_dem = np.greater(h,np.multiply(dem,active[1]))
        
        lu_above = np.sum(np.multiply(lu_temp,above_dem)) # Total area of cells with mounding above dem
        mound_cells += lu_above
        
        for m, current_subregion in enumerate(subregion_list):
            lu_temp_subregion = np.multiply(lu_temp, (subregion_array == current_subregion).astype(float))
            lu_active_subregion = np.sum(lu_temp_subregion) # Total urban area of active cells
            urban_cells_subregion[m] += lu_active_subregion
            
            above_dem_subregion = np.multiply(above_dem, (subregion_array == current_subregion).astype(float))
            lu_above_subregion = np.sum(np.multiply(lu_temp,above_dem_subregion)) # Total area of cells with mounding above dem
            mound_cells_subregion[m] += lu_above_subregion
        
        i += 1
        
        mound_per = mound_cells/urban_cells
        mound_per_subregion = mound_cells_subregion/urban_cells_subregion
        
    return mound_per, mound_per_subregion, hmatrix

def get_oldobjectives(heads, wellinfo, landuse, dem, active, bottom, subregion_array):
    
    subregion_list = np.unique(np.unique(subregion_array)) #List of subregions
    subregion_list = subregion_list[subregion_list>0]
    subregions = subregion_list.shape[0]
    
    try:
        energy, energy_subregion = measureEnergy(heads, wellinfo, dem, bottom, subregion_list)
    except:
        energy, energy_subregion = np.nan, np.ones(subregions)*np.nan
#    try:
#        sub, sub_cells = measureSubidence(heads, dem, active, bottom, subregion_array, subregion_list)
#    except:
#        sub, sub_cells = np.nan, np.nan
    try:
        wq_cells, wq_cells_subregion, hwq = measureWaterQuality(heads, dem, active, bottom, subregion_array, subregion_list)
    except:
        wq_cells, wq_cells_subregion, hwq = np.nan, np.ones(subregions)*np.nan, np.nan
    try:
        mound_per, mound_per_subregion, h = measureMound(heads, dem, active, landuse, subregion_array, subregion_list, [132,252])
    except:
        mound_per, mound_per_subregion, h = np.nan, np.ones(subregions)*np.nan, np.nan

    return energy, energy_subregion, wq_cells, wq_cells_subregion, mound_per, mound_per_subregion

def get_newobjectives(heads, wellinfo, landuse, dem, active, bottom, subregion_array):
    
    subregion_list = np.unique(np.unique(subregion_array)) #List of subregions
    subregion_list = subregion_list[subregion_list>0]
    subregions = subregion_list.shape[0]
    
    try:
        a, a_subregion = measureAvailability(wellinfo,subregion_array,subregion_list)
    except:
        a, a_subregion = np.nan, np.ones(subregions)*np.nan
    try:
        wq_cells, wq_cells_subregion, hwq = measureWaterQuality(heads, dem, active, bottom, subregion_array, subregion_list)
    except:
        wq_cells, wq_cells_subregion, hwq = np.nan, np.ones(subregions)*np.nan, np.nan
    try:
        mound_per, mound_per_subregion, h = measureMound(heads, dem, active, landuse, subregion_array, subregion_list, [132,252])
    except:
        mound_per, mound_per_subregion, h = np.nan, np.ones(subregions)*np.nan, np.nan

    return a, a_subregion, wq_cells, wq_cells_subregion, mound_per, mound_per_subregion

def calculate_SOSWR(heads, stats):
    '''
    Calculates the sum-of-squared, weighted residual given weights (stats) and array (heads) with three columns: simulated head, observed head, observation ID
    '''
    soswr = 0
    maxerror = 0
    for i in range(len(stats)):
        currenterror = (heads[i+1][1] - heads[i+1][0])
        maxerror = max([maxerror, np.abs(currenterror)])
        soswr += ((1/stats[i]) * currenterror)**2
    
    return soswr, maxerror