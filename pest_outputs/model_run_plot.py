# Author: Cecile Coulon
#
# ---------------------------------- Readme -----------------------------------
# 1) Check requirements (see below), install missing softwares
# 2) Edit the script preamble (just below)
# 3) Run the script below
#
#------------------------------ 1) Requirements -------------------------------
# -- MODFLOW-2005 groundwater model
# http://www.usgs.gov/software/modflow-2005-usgs-three-dimensional-finite-difference-ground-water-model
#
#----------------------------- 2) Script preamble -----------------------------
import os

# Enter working directory
working_dir = 'C:/Users/CECOU50/Documents/Cecile Coulon/4_Model Grande Entree/Model cecile/GitHub/pest_outputs/'
os.chdir(working_dir)

# Enter path to MODFLOW-2005
modflow_path = os.path.expanduser('C:/WRDAPP/MF2005.1_12/bin/mf2005')


# Choose what parameter set to work with (uncomment one of the following)
#mode = 'prior' # Prior to calibration (prior parameter set)
mode = 'posterior' # Post calibration (best parameter set)

#--------------------------------- 3) Script ----------------------------------

#%%------------------------- load paths and libraries -------------------------

# Path to model input files
path_bin_in = 'preproc_IDM_20m_qgis_v6.bin' # Path to the bin file generated with QGridder
path_mnw2_data = 'muni_wells_mnw2_data.csv' # Path to the csv containing information for MNW2
path_krig_factors = 'pp_30max.fac' # Path to the file containing the kriging factors

# Path to GIS folder
gis_path = os.path.join(working_dir, '../gis/')

# Load python libraries
import pickle, copy, pyemu, flopy
import pandas as pd
import numpy as np

# Load functions
import functions_model_plot as fplot

#%%----------------------- Load GIS data from bin file ------------------------

# Retrieve data from the bin file generated with QGridder
with open(path_bin_in, 'rb') as handle:
     objects = pickle.load(handle)
     (nrow, ncol, delr, delc, ibound, sea, dem, geology, thk_valley, cx, cy, 
      muni_wells_row_col, muni_pumping, domestic_wells_row_col, ind_old_wells_row_col,
      obs_wells_row_col, obs_head, obs_zeta, tdem_row_col, ert_row_col) = objects

# Create row_col dictionaries for h and zeta observations
Dict_z = {key: value for key, value in obs_zeta.items() if value == 1} # Create dictionary of zeta observation wells
obs_z_row_col = {key: value for key, value in obs_wells_row_col.items() if key in Dict_z.keys()} # Narrow obs_wells_row_col dictionary to obs_z_row_col (containing only zeta observations)
obs_z_row_col.update(muni_wells_row_col) # Add municipal wells to the list to create fictional zeta observations

Dict_h = {key: value for key, value in obs_head.items() if value == 1} # Create dictionary of head observation wells 
obs_h_row_col = {key: value for key, value in obs_wells_row_col.items() if key in Dict_h.keys()}  # Narrow obs_wells_row_col dictionary to obs_h_row_col (containing only head observations)
obs_h_row_col.update(muni_wells_row_col) # Add municipal wells to the list they all have h observations

#%%----------------------------- Model parameters -----------------------------

# Determine which file to read depending on the parameter set to use

if mode == 'prior': # Use prior parameter set (before calibration)
    
    add_to_name = '_priorcalib'
    folder = 'model_run' + add_to_name + '/'
    path_hk_pp = folder + 'hk1pp' + add_to_name + '.dat' # Path to the file containing the hk values at pilot points
    param_df =  pyemu.Pst('parameter_estimation/calreg_v16_pp_hp.pst').parameter_data.parval1.to_frame() # Path to the file containing the other parameter values
    param_df.columns = ['value']
    
if mode == 'posterior': # Use best parameter set, post calibration
    
    add_to_name = '_bpa'
    folder = 'model_run' + add_to_name + '/'
    path_hk_pp = folder + 'hk1pp' + add_to_name + '.dat' # Path to the file containing the hk values at pilot points
    param_df = pd.read_csv('parameter_estimation/calreg_v16_pp_hp.par', skiprows=1, header=None, usecols=[0,1], 
                           names=['parnme','value'], index_col=0, delim_whitespace=True) # Path to the file containing the other parameter values

# Read parameters

hk_dunes = param_df.loc['hk_dunes', 'value']
hk_valley = param_df.loc['hk_valley', 'value']
hk_seabed = param_df.loc['hk_seabed', 'value']
rech_mmyr = param_df.loc['rech_mmyr', 'value']
alpha_T = param_df.loc['alpha_t', 'value']
B = {well: param_df.loc['b' + well, 'value'] for well in muni_wells_row_col.keys()}
sy = param_df.loc['sy', 'value']
ss = param_df.loc['ss', 'value']
ne = param_df.loc['ne', 'value']
slvl = param_df.loc['slvl', 'value']

# Create hk grid from hk values determined at pilot points ('hk1pp.dat' pp_file) & kriging factors ('pp.fac' factors_file)

hk_EDC_array = pyemu.utils.fac2real(pp_file = path_hk_pp, factors_file = path_krig_factors, out_file=None, fill_value=np.nan)

#%%================================ 1) No SWI =================================

# Discretization package
nlay = 1 # Number of layers (along height)
top= 100 # Top elevation of layer
botm = -300 # Bottom elevation

nper = 1 # Number of stress periods
perlen = 2*24*3600 # Length of stress periods (2 days, converted to seconds)
nstp = 1 # Number of time steps in each stress period (scalar or array)
steady = True # Steady state simulation
tsmult = 1 # Time step multiplier
itmuni = 1 # Time units (1=seconds)
lenuni = 2 # Length units (2=meters)

# Basic package
h0 = 0 # Initial heads (m)

# Layer-Property Flow package
laytyp = 1 # Convertible layer (T = f(head) throughout simulation)
thk_dunes = 10 # m Thickness of sand dunes

# General-Head Boundary package
seabed_thk = 150 # Seabed thickness (m) = aquifer thickness/2
rho_f = 1000 # Freshwater density (kg/m3)
rho_s = 1025 # Saltwater density (kg/m3)

# Well package WEL
domestic_use_m3d = 80 # Total amount of pumped water from individual wells (Q3/day)

# Revised Multi-Node Well package MNW2
nnodes = 0 # Number of cells associated with the well
losstype = 'general' # Model for well loss
ppflag = 0 # Partial penetration flag

#%%--------------------------- Parameter processing ---------------------------

# Basic package
hstart = h0*np.ones((nlay, nrow, ncol)) # Grid of initial heads

# Layer-Property Flow package
hk_array = hk_EDC_array # Initialize
hk_valleyH = ( hk_EDC_array * (dem-botm-thk_valley) + hk_valley * thk_valley ) / (dem-botm)
hk_dunesH = ( hk_EDC_array * (dem-botm-thk_dunes) + hk_dunes * thk_dunes ) / (dem-botm)
hk_array[geology == 3] = hk_valleyH[geology == 3] # Paleoglacial valley
hk_array[geology == 4] = hk_dunesH[geology == 4] # Sand dunes

# Recharge package
recharge_ms = rech_mmyr/(1000*365.25*24*3600) # Convert recharge values from mm/year to m/s
rech = np.zeros((nrow,ncol)) # Recharge grid
rech[ sea == 0 ] = recharge_ms # For onshore cells, assign recharge

# General-Head Boundary package
epsilon = (rho_s - rho_f)/rho_f # Density ratio
corr_factor = 1 - ( alpha_T/(dem.mean()-botm) ) ** (1/4) # Density ratio correction factor (Lu & Werner 2013)
epsilon_corr = epsilon * corr_factor # Corrected density ratio

ghbc_EDC = (hk_seabed/seabed_thk) * delr * delc # Boundary conductance for EDC geology (m2/s)
ghb_cond = np.ones((nrow,ncol)) * ghbc_EDC # Default geology = EDC
ghb_stage = (epsilon+1) * slvl + (-epsilon) * dem # Stage (FW head at the ocean bottom) (m)
ghb_stage[sea != 1] = np.nan # GHB stage = nan for onshore cells

idx = np.logical_and(sea == 1, ibound == 1) # Boolean (indexes of active & ocean cells)
ghb_sp_data = [] # GHB grid
for row_col in np.argwhere(idx): # For indices where idx is True
    row, col = row_col
    ghb_sp_data.append( [0, row , col, ghb_stage[row,col], ghb_cond[row,col]] )
    # Each ghb cell defined by: [layer (here only 1), row, column, stage, conductance]

# Well package WEL
domestic_use_m3s = domestic_use_m3d/(24*3600) # Total domestic use: conversion from m3/day to m3/s
domestic_use_m3s_pc = domestic_use_m3s/len(domestic_wells_row_col) # Use per domestic well = total use/number of wells, in m3/s

wel_sp_data = {} # Build WEl dictionary = {nper: [layer, row, column, pumping rate]}
lay = 0 # Layer in which wells are pumping
for i in range(nper): # Compute for each stress period (here only 1)
    wells_data = []
    for well_id in domestic_wells_row_col.keys(): # Create one well_data array per well
        well_data = [lay] # Layers in which well are pumping
        well_data.append(domestic_wells_row_col[well_id][0][0]) # Append row to list
        well_data.append(domestic_wells_row_col[well_id][0][1]) # Append column to list
        well_data.append(-domestic_use_m3s_pc) # Append Q to list (m3/s, <0 for withdrawal)
        wells_data.append(well_data) # Append each individual well_data array to wells_data array
    wel_sp_data[i] = wells_data # Add wells_data array to WEL list   

# Revised Multi-Node Well package MNW2
node_data_df = pd.read_csv(path_mnw2_data) # Table with MNW2 data by node
node_data_df['nnodes'] = nnodes
node_data_df['ppflag'] = ppflag
node_data_df['losstype'] = losstype
node_data_df['B'] = node_data_df['wellid'].map(B)
node_data = node_data_df.to_records() # Convert dataframe to a rec array for compatibility with flopy

# Stress period information for MNW2
per = [0] # List of stress periods
active_mnw2_wells = len(node_data_df) # Number of active MNW2 wells
wel_sp_data_mnw2_df = pd.DataFrame(list(zip(per*active_mnw2_wells, node_data_df.wellid, -node_data_df.Qm3_s)), 
               columns =['per', 'wellid', 'qdes']) # qdes = actual volumetric pumping rate at the well (m3/s, <0 for withdrawal)
pers = wel_sp_data_mnw2_df.groupby('per') # Group wel_sp_data_mnw2_df by periods
wel_sp_data_mnw2 = {i: pers.get_group(i).to_records() for i in range(nper)} # Convert df to dictionary

# Multi-Node Well Information Package MNWI
unit_mnwi, qndflag, qbhflag = 35, 0, 0 # Unit number for the output files, flags for writing additional flow information in the output files
mnwi_list = [list(x) for x in zip(node_data_df.wellid, [unit_mnwi]*active_mnw2_wells, \
             [qndflag]*active_mnw2_wells, [qbhflag]*active_mnw2_wells)]

# Output Control package
oc_list = ['save head', 'save budget']
spd = {(i,0) : oc_list for i in range(nper)} # Save head & budget for (stress period i, timestep 1)

#%%------------------------------- Build model --------------------------------

name_no_swi = folder + 'idm_steady_no_swi_allwells' + add_to_name # Model name

# Create modflow model object
mf = flopy.modflow.Modflow(modelname=name_no_swi, exe_name = modflow_path)

# Discretization package
discret = flopy.modflow.ModflowDis(mf, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, 
                                   delr=delr, delc=delc, top=top, botm=botm, 
                                   perlen=perlen, nstp=nstp, tsmult=tsmult, 
                                   steady=steady, itmuni=itmuni, lenuni=lenuni)
# Basic package
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=hstart)
# Recharge package
rch = flopy.modflow.ModflowRch(mf, rech=rech)
# Layer-Property Flow package ## adapted to pilot points
lpf = flopy.modflow.ModflowLpf(mf, laytyp=laytyp, hk=hk_array, sy=sy, ss=ss)
# General-Head Boundary package
ghb = flopy.modflow.ModflowGhb(mf, stress_period_data = ghb_sp_data)
# Well package
wel = flopy.modflow.ModflowWel(mf, stress_period_data = wel_sp_data)
# Revised Multi-Node Well package
mnw2 = flopy.modflow.ModflowMnw2(mf, mnwmax = active_mnw2_wells, node_data = node_data,
                                 stress_period_data = wel_sp_data_mnw2,
                                 itmp = [active_mnw2_wells]*nper)
# Multi-Node Well Information Package
mnwi = flopy.modflow.ModflowMnwi(mf, byndflag=36, mnwobs=active_mnw2_wells, # bynd = output file for MNW2 information
                                 wellid_unit_qndflag_qhbflag_concflag = mnwi_list)
# Output Control package
oc = flopy.modflow.ModflowOc(mf, stress_period_data = spd)
# Preconditioned Conjugate-Gradient package
pcg = flopy.modflow.ModflowPcg(mf)

# Write input files
mf.write_input()

#%%-------------------------------- Run model ---------------------------------

#mf.run_model(silent=False)

#%%------------------------------- Read outputs -------------------------------

hds_noswi = flopy.utils.binaryfile.HeadFile(name_no_swi + '.hds') # Read modflow binary output file containing head values
h_noswi = hds_noswi.get_alldata()[0,0,:,:] # Get all the data from the binary file
h_noswi[ h_noswi < -888 ] = 0 # Assign nan to no data value

#%%================================== 2) SWI ==================================

# Discretization package
nper = 1 # Number of stress periods
perlen_yr = 1000 # Length of stress periods (years)
stplen_obj_day = 200 # Desired length of time steps (days)
steady = False # Type of simulation (transient)
tsmult = 1 # Time step multiplier

perlen = int(perlen_yr*365.25*24*3600) # Convert perlen from years to seconds
stplen_obj = stplen_obj_day*24*3600 # Convert stplen_obj from days to seconds
nstp = round(perlen/stplen_obj) # Number of time steps = length of stress period / length of time steps
stplen = round( (perlen/nstp)/(24*3600) , 2) # Real length of time steps (days)

# Saltwater Intrusion package
nsolver = 2 # (2 = use PCG solver)
toeslope = 0.16
tipslope = toeslope
nu = [0, epsilon_corr] # [0, density ratio]

zstart = -40 * h_noswi # Ghyben-Herzberg of h_no_SWI
zstart[ zstart < botm ] = botm # Avoid SWI to go below aquifer bottom
zstart = [zstart] # List of initial FW-SW interface elevations
isource = np.zeros((nrow, ncol), np.int) # Sources/sinks have same fluid density as active zone at the top of the aquifer 
isource[ sea == 1 ] = -2 # Ocean bottom: infiltrating water = SW and exfiltrating water = same type as water at the top of the aquifer

# Output Control package
output_nstp = 50 # Output frequency for zeta 
spd = {}
for istp in range(0,nstp,output_nstp): # Save head for istp ranging from 0 to nstp, with step = output_nstp
    spd[ (0, istp) ] = ['save head'] # Save head for: (stress period 1, timestep i) (kper i, kstp i)

#%%------------------------------- Build model --------------------------------

name_swi = folder + 'idm_transient_swi_allwells' + add_to_name # Model name
ipakcb = None # Flag used to determine if cell-by-cell budget data should be saved

# New model
mf = flopy.modflow.Modflow(modelname=name_swi, exe_name = modflow_path) 

# Discretization package
discret = flopy.modflow.ModflowDis(mf, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, 
                                   delr=delr, delc=delc, top=top, botm=botm, 
                                   perlen=perlen, nstp=nstp, tsmult=tsmult, 
                                   steady=steady, itmuni=itmuni, lenuni=lenuni)
# Basic package
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=h_noswi) 
# Recharge package
rch = flopy.modflow.ModflowRch(mf, rech = rech, ipakcb=ipakcb)
# Layer-Property Flow package ## adapted to pilot points
lpf = flopy.modflow.ModflowLpf(mf, laytyp=laytyp, hk=hk_array, sy=sy, ss=ss)
# General-Head Boundary package
ghb = flopy.modflow.ModflowGhb(mf, stress_period_data = ghb_sp_data, ipakcb=ipakcb)
# Well package
wel = flopy.modflow.ModflowWel(mf, stress_period_data = wel_sp_data, ipakcb=ipakcb)
# Revised Multi-Node Well package
mnw2 = flopy.modflow.ModflowMnw2(mf, mnwmax = active_mnw2_wells, node_data = node_data, 
                                 stress_period_data = wel_sp_data_mnw2, ipakcb = ipakcb,
                                 itmp = [active_mnw2_wells]*nper)
# Multi-Node Well Information Package
mnwi = flopy.modflow.ModflowMnwi(mf, byndflag=36, mnwobs=active_mnw2_wells,
                                 wellid_unit_qndflag_qhbflag_concflag = mnwi_list)
# Output control package
oc = flopy.modflow.ModflowOc(mf, stress_period_data = spd)
# Preconditioned Conjugate-Gradient package
pcg = flopy.modflow.ModflowPcg(mf)
# Saltwater-Intrusion package
swi = flopy.modflow.ModflowSwi2(mf, nsrf=1, istrat=1, ipakcb=ipakcb, 
                                toeslope=toeslope, tipslope=tipslope, nu=nu, 
                                zeta=zstart, alpha=0.1, beta=0.1, ssz=ne, 
                                isource=isource, nsolver=nsolver, iswizt=55)
# Write input files
mf.write_input()

#%%-------------------------------- Run model ---------------------------------

#mf.run_model(silent=False)

#%%------------------------------- Read outputs -------------------------------

#---- Freshwater heads

hds_swi = flopy.utils.binaryfile.HeadFile(name_swi + '.hds') # Read the .hds modflow binary output file, heads are saved here
times = hds_swi.get_times() # List of simulation times for which outputs are saved (in seconds)
hswi = hds_swi.get_alldata()[:,:,:,:] # Either get ALL head data from EVERY simulation time

#---- MNW2 data

mnw2_heads = pd.read_table(name_swi + '.bynd', delimiter='\s+', usecols = [0,7,8], 
                           names = ['name', 'hwell', 'hcell'], header = 0, index_col=0) # Read the .bynd modflow output file as a dataframe
mnw2_heads.index = mnw2_heads.index.str.lower()
mnw2_last = mnw2_heads.tail(9) # Retrieve information on last timestep

#---- Interface data

zta = flopy.utils.binaryfile.CellBudgetFile(name_swi + '.zta') # Read the .zta modflow binary output file, zeta values are saved here
kstpkper = zta.get_kstpkper() #  Get a list of unique stress periods and timesteps in the file (nper,nstp) = dictionary keys of spd (defined in Output Control package)
zeta = []
for kk in kstpkper: # Get ALL zeta data from EVERY simulation time (EVERY (nper,nstp) couple)
    zeta.append(zta.get_data(kstpkper=kk, text='ZETASRF  1')[0])
zeta = np.array(zeta)

#%%--------------------------- Post-process outputs ---------------------------

#---- Post-process head and zeta outputs

hswiP = copy.deepcopy(hswi) # Freshwater heads
zetaP = copy.deepcopy(zeta) # Interface elevations

# Assign Nan value to inactive model cells
for i in range(len(spd)):
    hswiP[i,0,:,:][ ibound == 0 ] = np.nan
    zetaP[i,0,:,:][ ibound == 0 ] = np.nan

#---- Compute thickness of the freshwater lens (m)
thk_FW_lens = (hswiP-zetaP)

#---- Compute total freshwater volume (m3)
FW_storage = []
for i in range(len(hswi)):
    fw = np.nansum(thk_FW_lens[i,0,:,:])*delr*delc # FW storage in the FW lens = sum(hswi-zeta)*delr*delc
    FW_storage.append(fw)

#%% Plots

# Locate the model in a real world coordinate reference system
xmin, xmax, ymin, ymax = cx[0,0] - delr/2, cx[0,-1] + delr/2, cy[-1,0] - delc/2, cy[0,0] + delc/2
mf.modelgrid.set_coord_info(xoff=xmin, yoff=ymin, angrot=0, epsg=2946)
mf.sr = flopy.utils.reference.SpatialReference(delr = delr*np.ones(ncol, dtype=float), delc = delc*np.ones(nrow, dtype=float), xll=xmin, yll=ymin, rotation=0, epsg=2946)

#%%------------------------------- Plot Fig. 6A -------------------------------
''' Map: Hydraulic conductivity field '''

fplot.plot_fig6a(mf, path_hk_pp, mode=mode)

#%%------------------------------- Plot Fig. 6B -------------------------------
''' Map: Freshwater-seawater interface elevation '''

fplot.plot_fig6b(mf, zetaP[-1,0,:,:], gis_path, mode=mode)
