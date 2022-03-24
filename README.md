# Spatial-GW-Planning

Run | Name | N | Avg grid cells | Median grid cells
----- | ------ | ------ | ---------------- | -------------------
0 | MUN_CLUSTER-CLIP | 41 | 188 | 135
1 | CLUSTER | 200 | 38.5 | 32
2 | AGEB | 3307 | 2.3 | 1

OLD
This repository takes raster and csv input data for a groundwater model of the Valley of Mexico and runs a simulation model to compare or optimize recharge alternatives under multiple objectives. The portfolio of intervention options includes repairs to the distribution network to reduce leaks, increased wastewater treatment and infiltration, and increased infiltration of imported water using recharge basins. These policies can be optimized using the NSGAII MOEA according to four objectives: minimize energy use for pumping, maximize hydraulic head in high subsidence regions, minimize groundwater mounding in urban areas, and minimize the cost of the interventions.

Dependencies: ```flopy``` ```osgeo``` ```numpy``` ```pickle``` ```mpi4py``` ```SALib``` ```pandas```

To run new simulations you must have MODFLOW installed in ```C:\WRDAPP\MF2005.1_12\bin\mf2005.exe``` which can be downloaded at https://water.usgs.gov/ogw/modflow/mf2005.html#downloads. However, you can still plot saved files when ```run_alternatives``` and ```run_optimization``` are set to ```False```.

## Data Inputs
These scripts process and reformat the data inputs listed below for use in the simulation model.
```Data_main.py```
```recharge```
![Data](/images/data_processing.png)

The data for the model must be in the following format:

Data Set Required | Format
-------------------- | --------------------
Digital elevation model of model area | raster
Active model cells by model layer | raster
Layer thickness by model layer | raster
Initial head for model area | raster
Percent of grid cell for each land use type (natural, urban, water) | raster
Interpolated monthly precipitation | raster
Pumping data by month or year | csv
Area covered by each municipality | raster
Leaks as percentage of total usage by municipality | csv
List of all potential recharge basin sites | csv
List and characteristics of existing wastewater treatment plants | csv

## Simulation Options
```main.py``` script provides a number of modes for running single simulations with various recharge alternative and parameter options.

The user must input the data above and integer choices for interventions to run the simulation model. The choices for interventions are defined as: the number of wastewater treatment plants to rehabilitate for wastewater injection into the aquifer from the existing plants, the number of infiltration basins that will recharge the aquifer using imported water from the potential sites, and the percent of fixed leaks when compared to historical leaks (0 indicates the same level as historical leaks and 100 indicates all leaks are fixed).

Each model run will create an object with the following results:

 ```self.wwtps``` is a list of randomly selected wastewater treatment plants where reuse has been implemented. ```self.basins``` is a list of the row and column where each infiltration basin has been implemented. ```self.mthlyleak``` is an array with the total quantity of leaks in the system per month in m3 cost is the total number of interventions time their weights defined in the model set-up. ```self.wells``` is a dictionary of well objects input into the MODFLOW model which includes positive flow from wastewater treatment plants, leaks, and recharge basins and negative flow from pumping wells. ```self.landuse``` is a dictionary that contains a raster and list of the percentage of land use type (NATURAL, URBAN, or WATER) per model cell for each model phase.

## Sensitivity Analysis Tools
Parallel runs of the simulation model may be necessary to develop large samples of simulation output for global sensitivity analysis. The following tools are used to run the simulation model efficiently and in an organized fashion for use on a computing cluster. It is assumed that the user will be running a global sensitivity analysis with each of the four default management alternatives for each parameter set out of a sample of the user's choosing.

The singular function ```SA_mode``` in ```model_functions.py``` takes a set of previously defined management alternatives using the options described above and runs a simulation for each alternative using a provided parameter set.

```SA_main.py``` uses ```SA_mode``` in parallel to run simulations across a set number of processors for subsets of a global parameter set sample. The sample can be predefined or generated by the script.

The main outputs are:
Path | Contents
-------------------- | --------------------
model_files\output\sa\hob | head observations at observation locations for the historical management alternative
model_files\output\sa\err | sum of squared weighted residual for the observations and maximum error across all observations
model_files\output\sa\obj | three management objectives calculated for each recharge alternative simulation
model_files\output\sa\mun | three management objectives calculated for each recharge alternative simulation at a subregional scale (municipal)

These data can then be processed into plottable datasets using ```process_SA.py```, which includes sensitivity analysis of error and objective values using ```SALib```.

```SA_main.py``` and ```process_SA.py``` have accompanying bash shell scripts with the same name to implement on a cluster environment.

## Plotting Tools
```plot_results.py```
This script contains functions for plotting results from single simulation runs. The current default saves these figures to ```model_files\output\plots```.

```sa_analysis_main.py``` uses ```sa_analysis_functions.py``` to analyze and plot processed sensitivity analysis data. The current default saves these figures to ```images\sa```.
