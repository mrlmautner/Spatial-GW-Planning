# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:27:21 2019

@author: MM
"""

import ValleMexico_setup as vm
import tools.plot_results as pltvm
import tools.measureobjectives as mo
import flopy.utils.binaryfile as bf
import numpy as np
import time
import pandas as pd
from pathlib import Path
import shutil
#import pickle

def single_run(model_name, exefile, safolder, sarun=0, error_calc=False, verbose=False):
    # Create a folder for the sarun parameter set
    sa_loc = Path.cwd() / 'modflow' / safolder / str(sarun)
    sa_loc.mkdir(parents=True, exist_ok=True)
    
    # Run the simulation model under historical conditions
    hist_time = time.time()
    model = vm.model(name=str(Path(sa_loc) / model_name), sarun=sarun, exe_file=exefile)
    model.run_simulation_model(incl_obs=True, verbose=verbose)
    #model.run_simulation_model(alt_pumping=np.array([[1,1],[2,1],[3,1],[4,1]]),incl_obs=True, verbose=verbose)
    print(model_name + ' for model ' + '{:05d}'.format(sarun) + ' completed in: ' + str(time.time() - hist_time) + ' seconds', flush=True)
    
    if error_calc:
        # Load head observation information for historical model run
        stats = np.loadtxt(Path.cwd() / 'modflow' / 'OBS_JH_stats.csv',skiprows=1)
        df = pd.read_fwf(Path.cwd().joinpath('modflow').joinpath(safolder).joinpath(str(sarun)).joinpath(model_name+'.hob.out'),widths=[22,19,22])
        heads_obs = [df.columns.values.tolist()] + df.values.tolist()
        
        soswr, maxerror = mo.calculate_SOSWR(heads_obs, stats)
                
        # Move the head observation file into the SA experiment folder (safolder)
        sa_hob_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('hob').joinpath(safolder)
        try:
            shutil.move(str(sa_loc.joinpath(model_name + '.hob.out')), str(sa_hob_loc.joinpath('{:05d}'.format(sarun) + '.hob_out')))
        except:
            sa_hob_loc.mkdir(parents=True, exist_ok=True)
            shutil.move(str(sa_loc.joinpath(model_name + '.hob.out')), str(sa_hob_loc.joinpath('{:05d}'.format(sarun) + '.hob_out')))
    
        # Save the sum-of-squared, weighted residual error
        error = np.array([soswr, maxerror])
        hist_err_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('err').joinpath(safolder)
        try:
            np.savetxt(hist_err_loc.joinpath('{:05d}'.format(sarun) + '.csv'), error, delimiter=',')
        except:
            hist_err_loc.mkdir(parents=True, exist_ok=True)
            np.savetxt(hist_err_loc.joinpath('{:05d}'.format(sarun) + '.csv'), error, delimiter=',')
        
        return model, sa_loc, error
    
    else:
        return model, sa_loc, np.nan
    
def process_objectives(model_name, safolder, model, sa_loc, sarun=0, verbose=False, delfolder=False):        
    
    # Load subregion raster
    subregion_list = np.unique(np.unique(model.subregions)) # List of subregions
#    subregion_list = subregion_list[subregion_list>0]
    subregions = subregion_list.shape[0]
    
    # Calculate objectives
    try:
        heads = pltvm.get_heads([model_name], sarun)
        a, a_subregion, mound, mound_subregion = mo.get_2objectives(heads[model_name], model.wells, model.landuse, model.dem, model.actv, model.botm, model.subregions)
    except:
        a, a_subregion, mound, mound_subregion = np.nan, np.ones(subregions)*np.nan, np.nan, np.ones(subregions)*np.nan
    
    mound = mound*100
    mound_subregion = [i * 100 for i in mound_subregion]
        
    objectives = [a, mound]
    hist_obj_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('obj').joinpath(safolder)
    try:
        np.savetxt(hist_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')
    except:
        hist_obj_loc.mkdir(parents=True, exist_ok=True)
        np.savetxt(hist_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')

    objectives_subregion = np.array([a_subregion, mound_subregion])
    hist_subregion_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('subregion').joinpath(safolder).joinpath(str(sarun))

    for m, current_subregion in enumerate(subregion_list):
        try:
            np.savetxt(hist_subregion_loc.joinpath('{:05.0f}'.format(current_subregion) + '.csv'), objectives_subregion[:,m], delimiter=',')
        except:
            hist_subregion_loc.mkdir(parents=True, exist_ok=True)
            np.savetxt(hist_subregion_loc.joinpath('{:05.0f}'.format(current_subregion) + '.csv'), objectives_subregion[:,m], delimiter=',')
    
    # Delete the model directory and files for this SA parameter set
    if delfolder:
        shutil.rmtree(str(sa_loc), ignore_errors=True)
    
    return objectives, objectives_subregion, heads

def process_3objectives(model_name, safolder, model, sa_loc, sarun=0, verbose=False, delfolder=False):        
    
    # Load subregion raster
    subregion_list = np.unique(np.unique(model.subregions)) # List of subregions
#    subregion_list = subregion_list[subregion_list>0]
    subregions = subregion_list.shape[0]
    
    # Calculate objectives
    try:
        heads = pltvm.get_heads([model_name], sarun)
        energy, energy_subregion, wq, wq_subregion, mound, mound_subregion = mo.get_3objectives(heads[model_name], model.wells, model.landuse, model.dem, model.actv, model.botm, model.subregions)
    except:
        energy, energy_subregion, wq, wq_subregion, mound, mound_subregion = np.nan, np.ones(subregions)*np.nan, np.nan, np.ones(subregions)*np.nan
    
    mound = mound*100
    wq = wq*100
    mound_subregion = [i * 100 for i in mound_subregion]
    wq_subregion = [i * 100 for i in wq_subregion]
        
    objectives = [energy, wq, mound]
    hist_obj_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('obj').joinpath(safolder)
    try:
        np.savetxt(hist_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')
    except:
        hist_obj_loc.mkdir(parents=True, exist_ok=True)
        np.savetxt(hist_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')

    objectives_subregion = np.array([energy_subregion, wq_subregion, mound_subregion])
    hist_subregion_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('subregion').joinpath(safolder).joinpath(str(sarun))

    for m, current_subregion in enumerate(subregion_list):
        try:
            np.savetxt(hist_subregion_loc.joinpath('{:05.0f}'.format(current_subregion) + '.csv'), objectives_subregion[:,m], delimiter=',')
        except:
            hist_subregion_loc.mkdir(parents=True, exist_ok=True)
            np.savetxt(hist_subregion_loc.joinpath('{:05.0f}'.format(current_subregion) + '.csv'), objectives_subregion[:,m], delimiter=',')
    
    # Delete the model directory and files for this SA parameter set
    if delfolder:
        shutil.rmtree(str(sa_loc), ignore_errors=True)
    
    return objectives, objectives_subregion, heads

# Function to run SA: sarun is the index of the set of parameter values to use given a previously generated set of parameter values in SALib, weighted residuals allowed to proceed with an evaluation of the managed aquifer recharge alternatives
def run_alternatives(alternatives, exefile, safolder, sarun=0, verbose=False):
    
    num_alts = len(alternatives['names'])
    
    # Create objective arrays    
    a, mound = (np.ones(num_alts)*np.nan for i in range(2))
    a_subregion, mound_subregion = ([[] for _ in range(num_alts)] for i in range(2))
    
    # Create dictionary to store each model for all alternatives
    model_dict = {}

    # Execute the MODFLOW model for each alternative and collect results
    for i, name in enumerate(alternatives['names'][1:]):
        if verbose: print(name, 'Alternative', flush=True)
        error_calc = False
        if name == 'historical': error_calc = True
        
        alt_model, sa_loc, error = single_run(model_name=name, exefile=exefile, safolder=safolder, sarun=sarun, error_calc=error_calc, verbose=verbose)
        
        if name == 'historical': model_error = error
        model_dict[name] = alt_model
        
        objectives, objectives_subregion, heads = process_objectives(model_name='test', safolder=safolder, model=alt_model, sa_loc=sa_loc, sarun=sarun, verbose=verbose)
        
        # Load subregion raster
        subregion_list = np.unique(np.unique(alt_model.subregions)) # List of subregions
    #    subregion_list = subregion_list[subregion_list>0]
        subregions = subregion_list.shape[0]
            
        # Calculate objective performance
        if verbose: print('Calculating alternative performance under objectives', flush=True)
        
        try:
            heads = pltvm.get_heads([name], sarun)
            a[i], a_subregion[i], mound[i], mound_subregion[i] = mo.get_2objectives(heads[name], alt_model.wells, alt_model.landuse, alt_model.dem, alt_model.actv, alt_model.botm, alt_model.subregions)
        except:
            a[i], a_subregion[i], mound[i], mound_subregion[i] = np.nan, np.ones(subregions)*np.nan, np.nan, np.ones(subregions)*np.nan
        
    mound = mound*100
    mound_subregion = [i * 100 for i in mound_subregion]
        
    objectives = [a, mound]
    sa_obj_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('obj').joinpath(safolder)
    try:
        np.savetxt(sa_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')
    except:
        sa_obj_loc.mkdir(parents=True, exist_ok=True)
        np.savetxt(sa_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')

    objectives_subregion = np.array([a_subregion, mound_subregion])
    sa_subregion_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('subregion').joinpath(safolder).joinpath(str(sarun))
    

    for m, current_subregion in enumerate(subregion_list):
        try:
            np.savetxt(sa_subregion_loc.joinpath('{:05.0f}'.format(current_subregion) + '.csv'), objectives_subregion[:,:,m], delimiter=',')
        except:
            sa_subregion_loc.mkdir(parents=True, exist_ok=True)
            np.savetxt(sa_subregion_loc.joinpath('{:05.0f}'.format(current_subregion) + '.csv'), objectives_subregion[:,:,m], delimiter=',')

#    # Save last year of heads
#    sa_heads_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('heads').joinpath(safolder)
#    try:
#        f = open(sa_heads_loc.joinpath('{:05d}'.format(sarun) + '.pkl'),"wb")
#        pickle.dump(h,f)
#        f.close()
#    except:
#        sa_heads_loc.mkdir(parents=True, exist_ok=True)
#        f = open(sa_heads_loc.joinpath('{:05d}'.format(sarun) + '.pkl'),"wb")
#        pickle.dump(h,f)
#        f.close()

    return model_error, objectives, objectives_subregion, heads, model_dict