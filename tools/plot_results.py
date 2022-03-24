# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:02:08 2018

@author: MM
"""
import flopy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import flopy.utils.binaryfile as bf
from pathlib import Path
import pandas as pd
import matplotlib as mplt

sns.set(style="white", palette="muted", color_codes=True)
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.titlesize'] = 24
plt.rcParams.update({'font.size': 20})

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

STRT_YEAR = 1984
END_YEAR = 2014

ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns

indpath = Path.cwd() / 'input' / 'socialindicators'
plot_loc = Path.cwd() / 'output' / 'plots'

def get_heads(alt_list, folder=-1):
    '''
    Generate a dictionary of the piezometric head for each model in list alt_list
    '''
    S_heads = {}
    if folder >= 0: 
        headpath = Path.cwd().joinpath('modflow').joinpath(str(folder))
    else:
        headpath = Path.cwd().joinpath('modflow')
    for name in alt_list:
        S_heads[name] = bf.HeadFile(str(headpath.joinpath(name+'.hds')))

    return S_heads

def alt_selection(solutions, alt_short):
    objs = solutions.copy()
    objs = objs[:,~np.isnan(objs).any(axis=0)]

    # use a boolean index to keep track of nondominated solns
    keep = np.zeros(len(alt_short), dtype = bool)

    for i in range(len(alt_short)):
        for j in range(len(alt_short)):
            a = objs[i]
            b = objs[j]
            if np.all(a <= b) & np.any(a < b):
                keep[i] = True
                keep[j] = False

    return keep

def nondom_alternatives(folder, run, scale, alt_short, obj_short, ids):
    '''
    Non-dominated alternatives by municipality
    '''
    # FINISH: save to output
    
    obj_dict = {}
    for o, obj in enumerate(obj_short):
        obj_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('subregion').joinpath(folder).joinpath(str(run))
        obj_array = np.loadtxt(obj_loc.joinpath(obj+'_subregions.csv'), delimiter=',').T
        obj_dict[obj] = obj_array
    
    nondom_alternative = {}
    nondom_alt_array = np.ones((ids.shape[0],len(alt_short)))*np.nan
    for u, units in enumerate(ids):
        solutions = np.ones((len(alt_short),len(obj_short)))*np.nan
        for j in range(len(alt_short)):
            solutions[j,0] = obj_dict[obj_short[0]][u, j]*-1
            solutions[j,1] = obj_dict[obj_short[1]][u, j]
            solutions[j,2] = obj_dict[obj_short[2]][u, j]
        
        nondom_alternative[int(ids[u])] = alt_selection(solutions, alt_short)
        nondom_alt_array[u,:] = alt_selection(solutions, alt_short)
    
    return nondom_alternative, nondom_alt_array

def plot_nondom_heatmap(folder, alt_short, scale_names, obj_short, ind_short, ind_long):
    '''
    
    '''
    for s, scale in enumerate(scale_names):
        plot_loc_scale = plot_loc / 'Alternative Heatmap' / folder / scale
        plot_loc_scale.mkdir(parents=True, exist_ok=True)
        
        ind_array = np.genfromtxt(indpath.joinpath('soc-ind_'+scale+'.csv'), delimiter=',', skip_header=1)
        nondom_alternative, nondom_alt_array = nondom_alternatives(folder, s, scale, alt_short, obj_short, ind_array[:,0])
        
        nondom_heatmap = nondom_alt_array
        for r, row in enumerate(nondom_heatmap):
            if row.sum() == 1: nondom_heatmap[r,:] += row
        
        nondom_heatmap = pd.DataFrame(nondom_heatmap, columns=alt_short)
        for i, ind in enumerate(ind_short[1:5]):
            nondom_heatmap_temp = nondom_heatmap.copy()
            nondom_heatmap_temp[ind_long[i+5]] = np.round(ind_array[:,i+5]*100, decimals=1)
            nondom_heatmap_temp = nondom_heatmap_temp.set_index(ind_long[i+5]).sort_index().transpose()
            
            fig, ax = plt.subplots(figsize=(25,5))
            sns.heatmap(nondom_heatmap_temp, cmap="PuBuGn", ax=ax, cbar=False)
            
            plt.tight_layout()
            plt.savefig(plot_loc_scale.joinpath(ind+'.png'), dpi=300)
            plt.savefig(plot_loc_scale.joinpath(ind+'.svg'))
            plt.close()
    
def plt_IndvObj_ind(folder, obj_short, obj_long, obj_label, alt_short, alt_long, ind_short, ind_long, ind_label, scale_names=['MUN_CLUSTER-CLIP', 'CLUSTER', 'AGEB']):
    '''
    Given folder, scale_names, and objectives, calculate the change in objective values between each alternative and the historical. Then plot the social indicators (x-axis) vs management objectives (y-axis) as a scatter plot. Each column in the subplot is a different social indicator. There is a folder for each scale and a figure for each objective.
    '''
    for s, scale in enumerate(scale_names):
        
        plot_loc_scale = plot_loc / 'Indicator vs Objectives' / folder / scale
        plot_loc_scale.mkdir(parents=True, exist_ok=True)
        
        df_ind = pd.read_csv(indpath.joinpath('soc-ind_'+scale+'.csv'), names=ind_short, skiprows=1)
    
        for o, obj in enumerate(obj_short):
            
            obj_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('subregion').joinpath(folder).joinpath(str(s))
            obj_array = np.loadtxt(obj_loc.joinpath(obj+'_subregions.csv'), delimiter=',').T
            
            # Create dataframe with objective values
            df_obj = pd.DataFrame(data=obj_array, index=np.arange(0,len(obj_array)), columns=alt_short)
                                
            fig, axes = plt.subplots(nrows=1, ncols=4, sharey=True, figsize=(20, 6))
            fig.suptitle('Indicators vs ' + obj_long[o] + ' Objective')
#            fig.subplots_adjust(left=0.12, top=0.7)
            
            for i, ind in enumerate(ind_short[1:5]):

                df_deltaHist = pd.DataFrame(df_obj[alt_short[1:]].sub(df_obj[alt_short[0]], axis='index'))
                df_deltaHist[ind] = df_ind[ind]
                df_deltah_melt = pd.melt(df_deltaHist, id_vars=[ind], var_name='Alternative', value_name=obj)
                
                sns.scatterplot(data=df_deltah_melt, x=ind, y=obj, hue='Alternative', linewidth=0, ax=axes[i], legend=False)
                
#                axes[a,i].set_ylim(bottom=yllim[o], top=yulim[o])
                axes[i].set_ylabel('')
                axes[i].set_title(ind_long[i+1])
                axes[i].set_xlabel(ind_label[i])
                    
                if i == 0:
                    axes[i].set_ylabel(obj_label[o])
                    
            plt.legend(bbox_to_anchor=(1.25, 1), borderaxespad=0)
            plt.savefig(plot_loc_scale.joinpath(obj+'.png'), dpi=600)
            plt.savefig(plot_loc_scale.joinpath(obj+'.svg'))
            plt.close()

def plt_IndvObj_scale(folder, obj_short, obj_long, obj_label, alt_short, alt_long, ind_short, ind_long, ind_label, scale_names=['MUN_CLUSTER-CLIP', 'CLUSTER', 'AGEB'], scale_label=['Municipality\n(n = 42)\nMedian unit size: 135', 'Cluster\n(n = 200)\nMedian unit size: 32', 'Census Block\n(n = 3307)\nMedian unit size: 1']):
    '''
    Given folder, scale_names, and objectives, calculate the change in objective values between each alternative and the historical. Then plot the social indicators (x-axis) vs management objectives (y-axis) as a scatter plot. Each column in the subplot is a different social indicator. There is a folder for each scale and a figure for each objective.
    '''
    for i, ind in enumerate(ind_short[1:5]):
        
        plot_loc_ind = plot_loc / 'Indicator vs Objectives' / folder / ind
        plot_loc_ind.mkdir(parents=True, exist_ok=True)
    
        for o, obj in enumerate(obj_short):
            
            fig, axes = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(20, 6))
            fig.suptitle('Indicators vs ' + obj_long[o] + ' Objective')
            fig.subplots_adjust(left=0.12, top=0.7)
            
            for s, scale in enumerate(scale_names):
                
                df_ind = pd.read_csv(indpath.joinpath('soc-ind_'+scale+'.csv'), names=ind_short, skiprows=1)

                obj_loc = Path.cwd().joinpath('output').joinpath('spatial').joinpath('subregion').joinpath(folder).joinpath(str(s))
                obj_array = np.loadtxt(obj_loc.joinpath(obj+'_subregions.csv'), delimiter=',').T
                
                # Create dataframe with objective values
                df_obj = pd.DataFrame(data=obj_array, index=np.arange(0,len(obj_array)), columns=alt_short)
                                    
                df_deltaHist = pd.DataFrame(df_obj[alt_short[1:]].sub(df_obj[alt_short[0]], axis='index'))
                df_deltaHist[ind] = df_ind[ind]
                df_deltah_melt = pd.melt(df_deltaHist, id_vars=[ind], var_name='Alternative', value_name=obj)
                
                sns.scatterplot(data=df_deltah_melt, x=ind, y=obj, hue='Alternative', linewidth=0, ax=axes[s], legend=False)
                
#                axes[a,i].set_ylim(bottom=yllim[o], top=yulim[o])
                axes[s].set_ylabel('')
                axes[s].set_title(scale_label[s])
                axes[s].set_xlabel(ind_label[i])
                    
                if i == 0:
                    axes[s].set_ylabel(obj_label[o])
            
            plt.legend(bbox_to_anchor=(1.25, 1), borderaxespad=0)
            plt.savefig(plot_loc_ind.joinpath(obj+'.png'), dpi=600)
            plt.savefig(plot_loc_ind.joinpath(obj+'.svg'))
            plt.close()
            
def plt_alt_objectives(alt_names, num_alt, objectives):
    '''
    objectives is a list with an array of length number of alternatives for each objective
    '''
    print('Plotting alternative performance under objectives...')
    plt.rcParams['legend.fontsize'] = 20
    plt.rcParams['axes.titlesize'] = 22
    plt.rcParams['axes.labelsize'] = 22
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['figure.titlesize'] = 24
    plt.rcParams.update({'font.size': 20})
    
    c = ['k','goldenrod','blue','darkgreen']
    barWidth = 0.1
    r = np.arange(num_alt)*0.1 # bar position
    y_label = ['Pumping Energy (kWh)','Depth to Groundwater in Clay (m)','Percent of Urban Cells Flooded']
    obj_title = ['Energy Use','Subsidence Avoidance','Urban Flooding']

    normalized_o = np.zeros((num_alt, len(objectives)))

    fig, axes = plt.subplots(nrows=1, ncols=3,figsize=(12,7.2))

    for o, obj in enumerate(objectives):
        normalized_o[:,o] = obj / (obj.max(axis=0) - obj.min(axis=0))
        plt.subplot(1,len(objectives),o+1)
        plt.bar(r, obj, width=barWidth, edgecolor='white', color=c)
        plt.xticks(r, alt_names, rotation=35, ha='right')
        plt.ylabel(y_label[o])
        plt.title(obj_title[o], fontweight='bold', y=1.04)

    # Flip subsidence measure to be minimizing
    plt.gcf().subplots_adjust(wspace=0.45,left=0.09,right=.97,bottom=0.15,top=.9)
    plt.savefig(Path.cwd() / 'output' / 'plots' / 'Objectives.eps', dpi=600)
    plt.savefig(Path.cwd() / 'output' / 'plots' / 'Objectives.png', dpi=600)
    plt.close()

def parallel_axis(nondom_results, obj_labels, opt_run):
    # Plots a normalized parallel axis
    plt.figure()
    for ppoint in nondom_results:
        ppoint = (ppoint - nondom_results.min(axis=0)) / (nondom_results.max(axis=0) - nondom_results.min(axis=0))
        plt.plot(range(len(obj_labels)), ppoint, 'steelblue')

    plt.gca().set_xticks(range(len(obj_labels)))
    plt.gca().set_xticklabels(obj_labels)
    plt.savefig(Path.cwd().joinpath('parallelaxis_' + opt_run + '.png'))
    plt.close()

def get_budgets(alt_list, mapTitles, s_heads):
    # Budget
    df_Bdget = {}
    df_CumSum = {}
    
    for name in alt_list:
        mf_list = flopy.utils.MfListBudget(Path.cwd().joinpath('modflow').joinpath(name+".list"))
        incremental, cumulative = mf_list.get_budget()
    
        df_Bdget[name], df_CumSum[name] = mf_list.get_dataframes(start_datetime="04-30-1984")
        
        mthly_Bdget = df_Bdget[name].drop(['CONSTANT_HEAD_IN', 'TOTAL_IN', 'CONSTANT_HEAD_OUT', 'RECHARGE_OUT', 'TOTAL_OUT', 'IN-OUT', 'PERCENT_DISCREPANCY'], axis=1)
        
        mthly_Bdget['STORAGE_OUT'] = mthly_Bdget['STORAGE_OUT'].apply(lambda x: x*-1)
        mthly_Bdget['WELLS_OUT'] = mthly_Bdget['WELLS_OUT'].apply(lambda x: x*-1)
#        mthly_Bdget['DRAINS_OUT'] = mthly_Bdget['DRAINS_OUT'].apply(lambda x: x*-1)
        mthly_Bdget = mthly_Bdget.multiply(30/1000000)
        cols = mthly_Bdget.columns.tolist()
        # reorder columns
        cols = [cols[1]] + [cols[2]] + [cols[0]] + [cols[4]] + [cols[3]] 
        # "commit" the reordering
        mthly_Bdget = mthly_Bdget[cols]
        
        im = axes[i].imshow(hist_compare,vmin=0,vmax=20)
        CS = axes[i].contour(ACTIVE_LYR1, colors='k', linewidths=2)
        axes[i].xaxis.set_visible(False)
        axes[i].yaxis.set_visible(False)
        axes[i].set_title(mapTitle[i].format(i+1))
        
        fig.subplots_adjust(right=0.8)
        fig.colorbar(im, cax=cbar_ax, label='Change in Groundwater Head (m)')
        plt.savefig(Path.cwd() / 'model_output' / 'plots' / 'Hist_change_all.svg')
        plt.savefig(Path.cwd() / 'model_output' / 'plots' / 'Hist_change_all.png', dpi=600)
        plt.show()

    for name in alt_list:
        df_CumSum[name]['IN'] = df_CumSum[name]['RECHARGE_IN'].divide(1000000) + df_CumSum[name]['WELLS_IN'].divide(1000000)
        df_CumSum[name]['OUT'] = df_CumSum[name]['WELLS_OUT'].divide(1000000) #+ df_CumSum[name]['DRAINS_OUT'].divide(1000000)
        df_CumSum[name]['INOUTCUMSUM'] = df_CumSum[name]['IN'] - df_CumSum[name]['OUT']
    
    return df_Bdget, mthly_Bdget, df_CumSum
    
def plt_cum_sum(filename, alt_list, mapTitles, df_CumSum, start='01-31-1985', end='12-31-2013'):
    # Plotting defaults
    l = [4,2,2,2]
    c = ['k','goldenrod','blue','darkgreen']
    mark = ['-','-','-','--']
    
    for i,name in enumerate(alt_list):
        df_CumSum[name].INOUTCUMSUM[start:end].plot(linewidth=l[i],color=c[i],style=mark[i])
        
    plt.ylabel(r'Volume (million m$^3$)')
    #plt.xlabel(r'Year')
    #plt.ylim(-8,10)
    #plt.title('Cumulative In - Out')
    plt.legend(mapTitles,loc='lower left')
    plt.gcf().subplots_adjust(left=0.15,right=.95,bottom=0.1,top=.95)
    
    plt.savefig(filename+'.png', dpi=600)
    plt.savefig(filename+'.eps', dpi=600)
    
    plt.show()

def plt_obs_cluster(df, km_cluster, filename):
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    xs = df['LON']
    ys = df['LAT']
    zs = df['Z']
    
#    if time_marker:
#        markers = ['o','^','s']
#        for i in range(3):
#            ax.scatter(xs[df['t_group'].values==i], ys[df['t_group'].values==i], zs[df['t_group'].values==i], s=50, alpha=0.6, edgecolors=None, c=km_cluster[df['t_group'].values==i], cmap='Set1', marker=markers[i])
#    else:
#        ax.scatter(xs, ys, zs, s=50, alpha=0.6, c=km_cluster, cmap='Set1', edgecolors=None)
    
    ax.scatter(xs, ys, zs, s=50, alpha=0.6, c=km_cluster, cmap='Set1', edgecolors=None)
    
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Elevation (masl)')
    
    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    ax.set_zlim(zs.min(), zs.max())
    
    fig.tight_layout()

    fig.savefig(filename)
    plt.close()
