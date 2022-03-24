# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:23:15 2018

@author: MM
"""

import flopy
import numpy as np
import time
#import pickle
import calendar
from pathlib import Path

class model():

    # Initializer / Instance attributes
    def __init__(self, name, exe_file=str(Path('C:') / 'WRDAPP' / 'MF2005.1_12' / 'bin' / 'mf2005.exe'), modlims=[455000, 2107000, 539000, 2175000], cellsize=500, strt_yr=1984, end_yr=2014, ACTIVE=[Path.cwd() / 'input' / 'ACTIVE_VM_LYR1.asc', Path.cwd() / 'input' / 'ACTIVE_VM_LYR2.asc'], THICKNESS=[Path.cwd() / 'input' / 'THICK1_VM.asc', Path.cwd() / 'input' / 'THICK2_VM.asc'], GEO=[Path.cwd() / 'input' / 'GEO_VM_LYR1.asc', Path.cwd() / 'input' / 'GEO_VM_LYR2.asc'], DEM=Path.cwd() / 'input' / 'DEM_VM.asc', IH=Path.cwd() / 'input' / 'IH_1984_LT2750.asc', SUBR=Path.cwd() / 'input' / 'CLUSTER_VM.asc', MUN=Path.cwd() / 'input' / 'MUN_VM.asc', PAR=Path.cwd() / 'modflow' / 'params.pval', run=0):
        self.name = name # Assign name
        self.xll = modlims[0] # X coordinate of the lower left corner
        self.yll = modlims[1] # Y coordinate of the lower left corner
        self.xur = modlims[2] # X coordinate of the upper right corner
        self.yur = modlims[3] # Y coordinate of the upper right corner
        self.cellsize = cellsize # Grid size
        self.ncol = int((self.xur - self.xll) / self.cellsize) # Number of rows
        self.nrow = int((self.yur - self.yll) / self.cellsize) # Number of columns
        self.strt_yr = strt_yr
        self.end_yr = end_yr
        self.actv = [np.loadtxt(i,skiprows=6) for i in ACTIVE]  # Extent of model layers 1 through n
        self.thck = [np.loadtxt(i,skiprows=6) for i in THICKNESS] # Thickness of model layers 1 through n
        self.geo = [np.loadtxt(i,skiprows=6) for i in GEO] # Geologic formations in layers 1 through n
        self.dem = np.loadtxt(DEM,skiprows=6) # Digital elevation model of the basin (model top)
        self.ih = np.loadtxt(IH,skiprows=6) # Initial hydraulic head in layer 1 and layer 2
        self.subregions = np.loadtxt(SUBR,skiprows=6) # Geographic extent of each subregion
        self.mun = np.loadtxt(MUN,skiprows=6) # Geographic extent of each municipality
        self.nlay = 2 # This model only accepts 2 layers
        self.exe = exe_file
        self.unknown_pumping = []
        self.r_multiplier = []
        self.total_mthly_pumping = []
        
        # Create adjustable parameter dictionary
        # If PAR is a dictionary, extract parameter names and values
        if isinstance(PAR, dict):
            self.params = {}
            for i, name in enumerate(PAR['names']):
                p = PAR['values'][run][i]
                if int(PAR['transform'][i]) == 1:
                    p = np.exp(p)
                try:
                    self.params[name[:name.rfind('_')]].append(p)
                except:
                    try:
                        self.params[name[:name.rfind('_')]] = [p]
                    except:
                        pass
        # Otherwise, PAR is a filename: read file and create parameters
        else:
            self.params = {}
            with open(PAR) as f:
                pval = f.read().splitlines()
            
            for i in pval:
                try:
                    self.params[i[:i.rfind('_')]].append(float(i[i.rfind('   '):]))
                except:
                    try:
                        self.params[i[:i.rfind('_')]] = [float(i[i.rfind('   '):])]
                    except:
                        pass
                    
    def initializeFM(self):
        # modelname to set the file root 
        mf = flopy.modflow.Modflow(str(Path('modflow') / self.name), exe_name=str(self.exe))
        
        # Model domain and grid definition
        botm = np.zeros((self.nlay, self.nrow, self.ncol))
        sumthck = np.zeros((self.nrow, self.ncol))
        for b, thickness in enumerate(self.thck):
            sumthck = np.add(sumthck,thickness)
            botm[b,:,:] = self.dem - sumthck # Current layer bottom elevation
        self.botm = botm
        
        # Time discretization
        nper = (self.end_yr - self.strt_yr)*12 # Number of stress periods
        nstp = []
        for y in range(self.strt_yr,self.end_yr):
            for m in range(1,13):
                nstp.append(calendar.monthrange(y,m)[1])
        nstp = np.array(nstp)
        steady = np.zeros((nper),dtype=bool)
        
        dis = flopy.modflow.ModflowDis(mf, nlay=self.nlay, nrow=self.nrow, ncol=self.ncol, nper=nper, delr=self.cellsize, delc=self.cellsize, top=self.dem, botm=botm, perlen=nstp, nstp=9, tsmult=1.3, steady=steady, start_datetime='01/01/1984')
            
        # Model Boundaries & initial conditions
        # Active areas
        ibound = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.int32)
        for a, active in enumerate(self.actv):
            ibound[a,:,:] = active # Current layer active area
        
        # Variables for the BAS package
        strt = np.array([self.ih]*2)
        
        bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, ifrefm=True, ichflg=True, stoper=3)
        
        # Layer properties
        # Add LPF package to the MODFLOW model
        # Create a dictionary of arrays with geologic characteristics: HK = Hydraulic conductivity, VANI = Vertical anisotropy (H:V) of hydraulic conductivity, SS = Specific storage, SY = Specific yield
        geoarrays = {}
        
        # Loop through the layers and formations for each layer to apply the geologic parameters to each array
        for p in ['HK', 'VANI', 'SS', 'SY']:   
            geoarrays[p] = np.zeros((self.nlay,self.nrow,self.ncol))
            
            for l in range(self.nlay):
                for f, fval in enumerate(self.params[p]):
                    geoarrays[p][l,:,:] += (self.geo[l] == f+1) * fval
            
        layvka = [1]*self.nlay # Indicates that VANI represents the ratio of H:V hydraulic conductivity
        
        lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, laytyp=[1,1], layvka=layvka, hk=geoarrays['HK'], vka=geoarrays['VANI'], ss=geoarrays['SS'], sy=geoarrays['SY'], laywet=[1,1]) # Wetting, convertible both layers
#        lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, laytyp=[0,0], layvka=layvka, hk=geoarrays['HK'], vka=geoarrays['VANI'], ss=geoarrays['SS'], sy=geoarrays['SY']) # No wetting, confined both layers
        
        return mf, dis, bas, lpf
    
    def addNewWells(self, New_WEL, LYR, WEL_Dict=0, INFO_Dict=0, WEL_mult=1, start=0, end=0, dateType='per', coordType='xy', pumpwell=False, changepumping=False):
        '''
        New_WEL is an np array of the following format: X (or C), Y (or R), Start Year, End Year, Flow (m3/d)
        WEL_mult is a scalar multiplier to be applied to all wells in the data set New_WEL
        WEL_Dict is a dictionary that contains dictionary for each stress period, each dictionary contains an entry for each well with the layer, row, column, and pumping rate
        coordType is a marker that should be either 'xy' or 'rc' depending on the coordinate definition in the New_WEL array
        dateType is a marker that should be either 'yr' or 'per' depending on whether a value of year or stress period is passed in the Start and End time columns
        pumpwell
        changepumping
        '''
        
        # Initialize dictionary    
        if WEL_Dict == 0:
            WEL_Dict = {}
        if INFO_Dict == 0:
            INFO_Dict = {}
        
        # Assign start period and end period if defined in input
        if start > 0:
            New_WEL[:,2] = np.ones((New_WEL.shape[0]))*start
            
        if end > 0:
            New_WEL[:,3] = np.ones((New_WEL.shape[0]))*end
        
        # Convert X and Y to Column and Row
        if coordType == 'xy':
            cconvert = lambda x: int(np.floor((x - self.xll) / self.cellsize))
            New_WEL[:,0] = np.array([cconvert(xi) for xi in New_WEL[:,0]])
            rconvert = lambda y: int(np.floor((self.yur - y) / self.cellsize))
            New_WEL[:,1] = np.array([rconvert(yi) for yi in New_WEL[:,1]])
        
        if coordType == 'rc':
            New_WEL[:,0] = np.array([int(xi) for xi in New_WEL[:,0] - 1])
            New_WEL[:,1] = np.array([int(yi) for yi in New_WEL[:,1] - 1])
        
        # Convert data in year format to stress period format (months)
        if dateType == 'yr':
            New_WEL[:,2] = (New_WEL[:,2] - self.strt_yr) * 12 + 1
            New_WEL[:,3] = (New_WEL[:,3] - self.strt_yr) * 12 + 1
        
        # Loop through all wells in the dataset to fill dictionary
        for w in range(0,New_WEL.shape[0]):
            r = New_WEL[w,1]
            c = New_WEL[w,0]
            wellsubregion = self.subregions[int(r),int(c)]
                    
            # Reduce the pumping amount by a percentage by cluster
            if changepumping:
                wellcluster = New_WEL[w,5]
                P = float(self.altpump[np.where(self.altpump[:,0]==wellcluster)[0],1]) # the ratio of new pumping to old pumping
            else:
                P = 1
                    
            # Assign flow rate for each well to all stress periods indicated by start and end years
            for per in range(int(New_WEL[w,2] - 1),int(New_WEL[w,3] - 1)):
                
                try:
                    WEL_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult*P])
                    INFO_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult*P,wellsubregion]) # layer, row, column, volume (m3/d), subregion
                except:
                    WEL_Dict[per] = [[LYR,r,c,New_WEL[w,4]*WEL_mult*P]]
                    INFO_Dict[per]= [[LYR,r,c,New_WEL[w,4]*WEL_mult*P,wellsubregion]]
                    
        return WEL_Dict,INFO_Dict
    
    def addRecharge(self,LU_arrays, PRECIP, start=0, end=0, RCH_Dict=0, RCH_mult=[1,1,1], dateType='per'):
        '''
        Outputs a dictionary of recharge arrays based on land use multiplier, land use cover, and precipitation input
        LU_arrays: dictionary with 3 eantries, one for each land use type which contains gridded percent amounts for each land use type
        PRECIP: dictionary with 361 entries, one for each stress period which contains gridded precipitation
        RCH_Dict: existing dictionary holding recharge data or 0 if the dictionary must be initialized dateType: the date format for the start and end variables
        '''
        
        # Initialize dictionary: if there is no exisiting dictionary, create dictionary with no entries
        if RCH_Dict == 0:
            RCH_Dict = {}
        
        # If the recharge is for the first time step, apply only to the first time step
        if start == 0:
            for l, landuse in enumerate(['URBAN','NATURAL','WATER']):
                # If there is not already an entry for the selected stress period, create a new array
                try:
                    RCH_Dict[0] += PRECIP[0] * LU_arrays[landuse] * RCH_mult[l]
                except:
                    RCH_Dict[0] = PRECIP[0] * LU_arrays[landuse] * RCH_mult[l]
        
        # Convert data in year format to stress period format (months)
        if dateType == 'yr':
            start = (start - self.strt_yr) * 12 + 1
            end = (end - self.strt_yr) * 12 + 1
        
        # Loop through all stress periods between S_YR and E_YR
        else:
            for per in range(int(start - 1), int(end - 1)):
                
                # Apply recharge amounts for each land use type                
                for l, landuse in enumerate(['URBAN','NATURAL','WATER']):                    
                    # If there is not already an entry for the selected stress period, create a new array
                    try:
                        RCH_Dict[per] += PRECIP[per]*LU_arrays[landuse]*RCH_mult[l]
                    except:
                        RCH_Dict[per] = PRECIP[per]*LU_arrays[landuse]*RCH_mult[l]
    
        return RCH_Dict
    
    def outputControl(self,mf):
        ''' 
        Generate Output Control and Solver packages
        Add OC package to the MODFLOW model
        '''
        spd = {}
        data2record = ['save head', 'save drawdown', 'save budget', 'print budget']
        for y in range(0,30):
            for m in range(1,13):
                spd[y * 12 + m - 1, 8] = data2record.copy()
#                spd[y * 12 + m - 1, calendar.monthrange(self.strt_yr + y, m)[1] - 1] = data2record.copy() # If time steps in month is equal to number of days
#                for d in range(0,calendar.monthrange(self.strt_yr + y, m)[1]):
#                    spd[y * 12 + m - 1, d] = data2record.copy()
        spd[26,8] = ['save head', 'save drawdown', 'save budget', 'print budget', 'ddreference']
        oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

        # Add PCG package to the MODFLOW model
        pcg = flopy.modflow.ModflowPcg(mf,mxiter=20, iter1=20)
        
        return oc, pcg
    
    def run_simulation_model(self, alt_pumping=np.array([[1,1],[2,1],[3,1],[4,1]]), alt_pumping_reduction=[0,0,0], incl_obs=False, seed=1, verbose=False):
        '''
        
        '''
        
        np.random.seed(seed)
        timestart = time.time()
        if verbose: print('Processing data...')
        
        # Phase starting stress period
        PHASE_PER = [0, 132, 252, 360]
        phases = len(PHASE_PER) - 1
        S_per = PHASE_PER[0:len(PHASE_PER) - 1]
        E_per = PHASE_PER[1:len(PHASE_PER)]
        # Phase land use dataset year
        LU_PAR = ['1990', '2000', '2010']
        
        # Model internal variables
        sec2day = 60*60*24 # Second to day conversion
        
        self.altpump = alt_pumping
        
        # Water supply data
        hist_water_use = np.loadtxt(Path.cwd() / 'input' / 'decisions' / 'twu.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Initial (original) all other supplies before alternatives, matrix of size municipalities by phases (m3/s)
        total_water_use = hist_water_use*self.params['TWU']*sec2day # Multiply by total water use parameters (m3/d)
        self.twateruse = total_water_use.sum(axis=0) # Total water use for each model phase (m3/d)
        
        i_other = np.loadtxt(Path.cwd() / 'input' / 'decisions' / 'initial_supply.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Initial (original) all other supplies before alternatives (m3/s)
        i_other = i_other.sum(axis=0)*sec2day # Initial other water supply for each model phase (m3/d)
        
        # Calculate historical quantity of leaks in each municipality
        LEAK_MUN = np.loadtxt(Path('input') / 'leak' / 'LEAK_TOT_MUN.csv',delimiter=',',skiprows=1) # Total recharge percent per municipality: equal to percent of total water use (1997 values) x percent leak (~35%) x recharge percent (15%)
        leaks = np.zeros((LEAK_MUN.shape[0],phases+1))
        leaks[:,0] = LEAK_MUN[:,0]
        
        # Initialize the modflow model with the boundary conditions input above
        mf, dis, bas, lpf = self.initializeFM()
        
        if verbose: print('Basic, Discretization, and Layer packages generated in', str(time.time() - timestart), 'seconds')
        
        '''
        Land Use Type
        Fill a land use dictionary with the ARRAYs that represent the % of each land use cover in each cell and the LISTs that contain all the cells and percentages of each land use type
        '''
        LU = {}
        for i, LUset in enumerate(LU_PAR):
            LU[LUset] = {'ARRAY':{},'LIST':{}}
            
            for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
                filename = Path('input').joinpath('landuse').joinpath('LU-' + LUset + '-' + LUtype + '.asc')
                perarea =  np.loadtxt(filename,skiprows=6)
                LU[LUset]['ARRAY'][LUtype] = perarea
                
                LU[LUset]['LIST'][LUtype] = np.zeros((perarea.shape[0]*perarea.shape[1],5))
                
                l = 0
                for row in range(0,perarea.shape[0]):
                    for col in range(0,perarea.shape[1]):
                        if perarea[row,col] > 0.001:
                            LU[LUset]['LIST'][LUtype][l,2] = perarea[row,col]
                        LU[LUset]['LIST'][LUtype][l,0] = col
                        LU[LUset]['LIST'][LUtype][l,1] = row
                        LU[LUset]['LIST'][LUtype][l,3] = 1 - self.geo[0][row,col] # 0 if clay layer, 1 if no clay layer
                        LU[LUset]['LIST'][LUtype][l,4] = self.mun[row,col]
                        l += 1
                LU[LUset]['LIST'][LUtype] = LU[LUset]['LIST'][LUtype][LU[LUset]['LIST'][LUtype][:,2]>0,:]
    
#        # Save land use database for use in mounding objective
#        winfofile = r'model_files\optimization_data\objectives\LU_' + self.name + '.pickle'
#        with open(winfofile, 'wb') as handle:
#            pickle.dump(LU, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        '''
        Recharge
        Create recharge dictionary for MODFLOW RCH package based on land use multipliers and interpolated precipitation rasters
        '''
        newtime = time.time()
            
        RCH_DICT = {}
        Precip_Dict = {}
            
        for year in range(int(self.strt_yr),int(self.end_yr)):
            for month in range(1,13):
                per = (year - self.strt_yr) * 12 + month - 1
            
                filename = Path('input').joinpath('recharge').joinpath('claymult').joinpath('PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc')
                Precip_Dict[per] = np.loadtxt(filename,skiprows=6)
        
        for i, LUset in enumerate(LU_PAR):
            RCH_DICT = self.addRecharge(LU_arrays=LU[LUset]['ARRAY'], PRECIP=Precip_Dict, start=S_per[i]+1, end=E_per[i]+1, RCH_Dict=RCH_DICT, RCH_mult=self.params['RCH'])
        
        # Create MODFLOW RCH package
        rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_DICT)
        
        if verbose: print('RCH_Dict generated in', str(time.time() - newtime), 'seconds')
        
        '''
        Well objects: supply wells, distribution leaks
        '''
#        newtime = time.time()
        
        WEL_DICT = {}
        WEL_INFO = {}
        
        # Add supply wells
        # Import CONAGUA and SACM pumping datasets
        CAEM_array = np.loadtxt(Path('input') / 'wells' / 'PUMP_C.csv', delimiter=',', skiprows=1, usecols=[1,2,7,8,11,15]) # pumping in m3 per day
        WEL_DICT, WEL_INFO = self.addNewWells(CAEM_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, pumpwell=True, changepumping=True)
            
        SACM_array = np.loadtxt(Path('input') / 'wells' / 'PUMP_S.csv',delimiter=',', skiprows=1, usecols=[1,2,7,8,11,15]) # pumping in m3 per day
        WEL_DICT, WEL_INFO = self.addNewWells(SACM_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, pumpwell=True, changepumping=True)
        
        REPDA_array = np.loadtxt(Path('input') / 'wells' / 'PUMP_RC_Q.csv', delimiter=',', skiprows=1, usecols=[1,2,4,5,11,15]) # pumping in m3 per day
        self.r_multiplier = np.loadtxt(Path('input') / 'wells' / 'R_MULT.csv', delimiter=',', skiprows=1, usecols=[0]) # REPDA multiplier to ensure that total REPDA pumping is equal to the unknown pumping calibrated according to JoH

        # Include only municipalities with urban land cover
        mun = np.unique(self.mun)[1:].copy()
        
        # Loop through model phases
        for i, l in enumerate(LU_PAR):
            # Calculate unaccounted for water supply by subtraction to determine pumping in REPDA dataset
            total_mthly_pumping = max(0, self.twateruse[i] - i_other[i]) # Monthly pumping is equal to the total water use minus other supply (m3/d)
            self.total_mthly_pumping.append(total_mthly_pumping)
            
            # Generate monthly pumping datasets for REPDA data in single pumping value format
            # Loop through monthly periods in each model phase
            for p in range(S_per[i]+1,E_per[i]+1):
                        
                # Urban wells
                WEL_DICT, WEL_INFO = self.addNewWells(New_WEL=REPDA_array[REPDA_array[:,5]==1,:5], LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=self.params['Q'][i]*self.r_multiplier[p-1], start=p, end=(p + 1), pumpwell=True)
                
                # Peri-urban wells
                WEL_DICT, WEL_INFO = self.addNewWells(New_WEL=REPDA_array[REPDA_array[:,5]==0,:5], LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=self.r_multiplier[p-1], start=p, end=(p + 1), pumpwell=True)
        
        # Loop through model phases
        for i, l in enumerate(LU_PAR):
            leaks[:,i+1] = self.params['LK'][i]*LEAK_MUN[:,1]*(total_water_use[:,i]-alt_pumping_reduction[i]) # Total leak per municipality by model phase (m3/d)☻
            
            LEAK_array = np.zeros((LU[l]['LIST']['URBAN'].shape[0] * (PHASE_PER[i + 1] - PHASE_PER[i]),5))
            j = 0
            
            # Loop through monthly periods in each model phase
            for p in range(S_per[i]+1,E_per[i]+1):
                '''
                Leak Repair
                
                Create a well dictionary for all the leak cells. The leak cells will be treated
                as an injection well at each cell in which there is urban land cover. MODFLOW
                distributes injection wells evenly across the area of the cell. The leak
                percentage is based on the leak percentage determined by municipality. Then
                the leak amount determined for each cell is multiplied by the percent of
                urban land cover in that cell. Finally, leaks in cells that are located in the
                lacustrine zone are reduced by 90% assuming that the low hydraulic
                conductivity does not allow for high levels of infiltration and the sewer
                provides a preferential flow path out of the basin
                '''
            
                for n, m in enumerate(mun):
                    # Create an array for all urban model cells in this municipality
                    tempLeak = LU[l]['LIST']['URBAN'][(LU[l]['LIST']['URBAN'][:,4]==m),:4]
                    if len(tempLeak) > 0:
                        u_cells = tempLeak[:, 2].sum() # Number of urban model cells in municipality m
                        # Use total pumping for each stress period to determine leak quantities
                        LperCell = float(leaks[np.where(leaks==m)[0],i+1]) * self.params['IN'][0] / u_cells # Reference the leaks in the municipality (m3/d), multiply by infiltration rate, divide by the number of urban cells
                        tempLeak[:,2] *= LperCell
                        
                        # apply 90% returns to sewer under clay layer (Geologic formation 1)
                        tempLeak[tempLeak[:,3]==1, 2] *= 0.1
                        
                        # Get rows of all cells of urban land use type from list
                        LEAK_array[j:(j + tempLeak.shape[0]), 0] = tempLeak[:, 0]
                        
                        # Get columns of all cells of urban land use type from list
                        LEAK_array[j:(j + tempLeak.shape[0]), 1] = tempLeak[:, 1]
                        
                        # Set the period to the current stress period for all urban cells
                        LEAK_array[j:(j + tempLeak.shape[0]), 2] = p
                        
                        # Set the end of the period to the next stress period
                        LEAK_array[j:(j + tempLeak.shape[0]), 3] = p + 1
                        
                        # Set the multiplier to the percentage of urban land use type stored in list
                        LEAK_array[j:(j + tempLeak.shape[0]), 4] = tempLeak[:,2]
                        
                        # Set the new index to the previous index plus the number of cells added
                        j += tempLeak.shape[0]
          
            WEL_DICT, WEL_INFO = self.addNewWells(LEAK_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=self.params['LK'][i], coordType='rc')

        total_mthly_leak = np.zeros(PHASE_PER[3])
        for i in range(total_mthly_leak.shape[0]):
            total_mthly_leak[i] = np.sum(list(zip(*WEL_INFO[i]))[3])
        total_mthly_leak = (total_mthly_leak - total_mthly_pumping)
        
        wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT, options=['NOPRINT'])
            
        if verbose: print('WEL_Dict generated in', str(time.time() - newtime), 'seconds')
#        
#        ## Create pickle file of Well Data to be available for post processing of well energy use objective
#        winfofile = r'model_files\optimization_data\objectives\WEL_INFO_' + self.name + '.pickle'
#        with open(winfofile, 'wb') as handle:
#            pickle.dump(WEL_INFO, handle, protocol=pickle.HIGHEST_PROTOCOL)
#        print('WEL_Dict saved in',str(time.time()-newtime),'seconds')
        
        # Generate output control and solver MODFLOW packages 
        oc, pcg = self.outputControl(mf)
        
        if incl_obs:
            mf.add_existing_package(str(Path.cwd().joinpath('modflow').joinpath('OBS_JH.ob_hob')),ptype='HOB', copy_to_model_ws=False)
            mf.add_output(str(Path.cwd().joinpath('modflow').joinpath(self.name + '.hob.out')),unit=2001)
            
#        hob = flopy.modflow.ModflowHob.load(r'modflow\OBS.ob_hob', mf)
#        winfofile = r'modflow\OBS.pickle'
#        import pickle
#        with open(winfofile, 'wb') as handle:
#            pickle.dump(hob.obs_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
#    
        # Run Model and post processing
        ## Write the MODFLOW model input files
        if verbose:
            print('Data processed in', str(time.time() - timestart), 'seconds')
        
            newtime = time.time()
            print('Writing input file...')
            
        mf.write_input()
        
        if verbose:
            print('Input file written in', str(time.time() - newtime), 'seconds')
            
            newtime = time.time()
            print('Running MODFLOW model...')
            
        # Run the MODFLOW model
        success, buff = mf.run_model(silent=not verbose)
        
        if verbose:
            print('MODFLOW model completed run in', str(time.time() - newtime), 'seconds')
        
            print('Wrapper completed run in', str(time.time() - timestart), 'seconds')
        
        self.mthlyleak = total_mthly_leak
        self.wells = WEL_INFO
        self.landuse = LU
