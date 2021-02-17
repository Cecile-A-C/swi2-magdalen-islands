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
import os, pyemu

# Enter working directory
working_dir = 'C:/Users/CECOU50/Documents/Cecile Coulon/4_Model Grande Entree/Model cecile/GitHub/pest_outputs/'
os.chdir(working_dir)

# Enter path to MODFLOW-2005
modflow_path = os.path.expanduser('C:/WRDAPP/MF2005.1_12/bin/mf2005')

#--------------------------------- 3) Script ----------------------------------

#%%------------------------- load paths and libraries -------------------------

# PEST control file used for parameter estimation (containing the prior parameter set)

pst_priorcalib = pyemu.Pst('parameter_estimation/' + 'calreg_v16_pp_hp.pst') # Name of the PEST control file used for parameter estimation
namemodel_priorcalib = 'model_run_priorcalib/' + 'idm_transient_swi_allwells_priorcalib' # Name of the model run with the prior parameter set

# PEST control file used for the linear analysis (containing the best parameter set)

pst_bpa = pyemu.Pst('linear_analysis_bpa/' + 'calreg_v16_pp_hp_la_bpa.pst') # Name of the PEST control file used for the linear analysis
jco = os.path.join('linear_analysis_bpa/', 'calreg_v16_pp_hp_la_bpa.jco') # Load the .jco file
namemodel_bpa = 'model_run_bpa/' + 'idm_transient_swi_allwells_bpa' # Name of the model run with the best parameter set

# Path to GIS folder
gis_path = os.path.join(working_dir, '../gis/')

# Load python libraries
import copy
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# Load functions
import functions_linear_analysis as fpst

#%%------------------------------- Plot Fig. 5 --------------------------------

''' Scatter plots of simulated to observed data '''

fpst.plot_fig5a(pst_priorcalib) # Fig. 5A
fpst.plot_fig5b(pst_priorcalib) # Fig. 5B

#%%----------------------- Linear uncertainty analysis ------------------------

#---- Load the jacobian matrix to a dataframe

jco_mat = pyemu.Jco.from_binary(jco) # Load jacobian matrix
jco_df = jco_mat.to_dataframe() # Write to dataframe

#---- List of model forecasts

forecasts = jco_df[jco_df.index.str.startswith(('zmuni', 'fw_storage'))].index.to_list() # All observations starting with 'zmuni' or 'fw_storage'

#---- Create linear analysis object

la = pyemu.Schur(jco=jco, pst=pst_bpa, parcov=None, obscov=None, forecasts=forecasts, sigma_range=4.0)

print('# observations: ' + str(la.jco.shape[0]) + ', # adjustable parameters: ' + str(la.jco.shape[1]) + ', # forecasts: ' + str(len(la.forecast_names)))

#%%---------------- Prior vs posterior forecast distributions -----------------

# Retrieve prior & posterior forecast uncertainty (as variances)
for_sum = la.get_forecast_summary() # Dataframe containing prior & posterior forecasts variances (and percent reduction in variance)

# Compute the standard deviation from the variance

for_sum['prior_stdev'] = np.sqrt(for_sum.prior_var)
for_sum['post_stdev'] = np.sqrt(for_sum.post_var)

# Retrieve model forecasts from the simulation output files

fcsts_prior = fpst.get_forecasts_from_sim(namemodel_priorcalib)
fcsts_post = fpst.get_forecasts_from_sim(namemodel_bpa)

# Create a dataframe containing information on prior & posterior forecast distributions
# Dataframe with index = forecast name and columns = mean, stdev

prior_for = pd.DataFrame({'Mean': fcsts_prior.value, 'Stdev': for_sum['prior_stdev']})
post_for = pd.DataFrame({'Mean': fcsts_post.value, 'Stdev': for_sum['post_stdev']})

#%%------------------------------- Plot Fig. 7 --------------------------------

''' PLot prior & posterior probability distributions of the model forecasts,
defined by the mean +- 2 x standard deviation '''

fpst.plot_fig7(prior_for, post_for, 'fw_storage') # Fig 7A
fpst.plot_fig7(prior_for, post_for, 'zmuni1') # Fig 7B

#%%--------- Contribution of parameter groups to forecast uncertainty ---------

#---- Define parameter groups

hk = [par_name for par_name in la.pst.par_names if par_name.startswith('hk') \
              and par_name in pst_bpa.adj_par_names] #                                    1) Hydraulic conductivities
recharge = [par_name for par_name in la.pst.par_names if par_name.startswith('rech')] #   2) Recharge
alpha = [par_name for par_name in la.pst.par_names if par_name.startswith('alpha')] #     3) Transerve dispersivity (as a correction factor)
par_dict = {'Hydraulic conductivities': hk, 'Recharge': recharge, 'alpha_T': alpha}

# Get the posterior forecast uncertainties (as variance) as a result of parameter groups becoming perfectly known

par_contrib_gp_var = la.get_par_contribution(parlist_dict=par_dict) # Dataframe containing posterior forecast variances as a result of parameter groups becoming perfectly known

# Compute the standard deviation from the variance
par_contrib_gp_stdev = np.sqrt(par_contrib_gp_var)

#%%------------------------------- Plot Fig. 8 --------------------------------

''' Bar plot: Percent decrease in posterior forecast uncertainty (as standard deviation) 
when one parameter group is considered fully known '''

fpst.plot_fig8(post_for, par_contrib_gp_stdev.T, 'fw_storage') # Fig 8A
fpst.plot_fig8(post_for, par_contrib_gp_stdev.T, 'zmuni1') # Fig 8B

# Note: this function can plot all forecasts ('fw_storage', 'zmuni1', ... 'zmuni9')

#%%----------------- Forecast worth of existing observations ------------------

#---- Create observation groups (lists)

# Head observation groups

hshallow_gp = [obs for obs in la.get_obs_group_dict()['hshallow'] if obs in la.pst.nnz_obs_names] # From shallow wells
hdeep_gp = [obs for obs in la.get_obs_group_dict()['hdeep'] if obs in la.pst.nnz_obs_names] #       From deep wells
hmuni_gp = [obs for obs in la.get_obs_group_dict()['hmuni'] if obs in la.pst.nnz_obs_names] #       From municipal pumping wells

# Interface observation groups

zwells_gp = [obs for obs in la.get_obs_group_dict()['zwells'] if obs in la.pst.nnz_obs_names] # From deep wells
ztdem_gp = la.get_obs_group_dict()['ztdem'] #                                                   From the TDEM survey
zert_gp = la.get_obs_group_dict()['zert'] #                                                     From the ERT survey


#---- Create dictionaries with the observation groups

# Divide observations in 6 groups

obs_dict_6gps = {'hshallow': hshallow_gp, #  1) Head observations from shallow wells
                 'hdeep': hdeep_gp, #        2) Head observations from deep wells
                 'hmuni': hmuni_gp, #        3) Head observations from municipal pumping wells
                 'zwells': zwells_gp, #      4) Interface observations from deep wells
                 'ztdem': ztdem_gp, #        5) Interface observations from the TDEM survey
                 'zert': zert_gp} #          6) Interface observations from the ERT survey

# Divide observations in 3 groups

obs_dict_3gps = {'h': hshallow_gp + hdeep_gp + hmuni_gp, #  1) All head observations
                 'zwells': zwells_gp, #                     2) Interface observations from deep wells
                 'zgeo': ztdem_gp + zert_gp} #              3) Interface observations from all geophysical surveys

# Divide observations in 2 groups

obs_dict_2gps = {'h': hshallow_gp + hdeep_gp + hmuni_gp, #  1) All head observations
                 'zeta': zwells_gp + ztdem_gp + zert_gp} #  2) All interface observations


#---- Compute the decrease in prior forecast standard deviation when observations / observation groups are added individually to the calibration dataset

# Compute the decrease in prior forecast variance
# Values are stored in dataframes

added_obs_all_var = la.get_added_obs_importance(obslist_dict=la.pst.nnz_obs_names, #          1) When observations are added individually (no groups)
                                                base_obslist=None, reset_zero_weight=False)

added_obs_6gps_var = la.get_added_obs_importance(obslist_dict=obs_dict_6gps, #                2) When observation groups are added individually (6 groups)
                                                 base_obslist=None, reset_zero_weight=False)

added_obs_3gps_var = la.get_added_obs_importance(obslist_dict=obs_dict_3gps, #                3) When observation groups are added individually (3 groups)
                                                 base_obslist=None, reset_zero_weight=False)

added_obs_2gps_var = la.get_added_obs_importance(obslist_dict=obs_dict_2gps, #                4) When observation groups are added individually (2 groups)
                                                 base_obslist=None, reset_zero_weight=False)


# Compute the standard deviation from the variance
# Dataframes now contain standard deviations

added_obs_all_std = np.sqrt(added_obs_all_var)
added_obs_6gps_std = np.sqrt(added_obs_6gps_var)
added_obs_3gps_std = np.sqrt(added_obs_3gps_var)
added_obs_2gps_std = np.sqrt(added_obs_2gps_var)


# Add the posterior standard deviation to dataframes, resulting from the inclusion of ALL observations in the calibration

def add_post_std_df(df):
    ''' Function that adds a column to an existing dataframe (df), 
    containing the forecast's posterior standard deviation '''
    
    df = df.T
    df['post'] = for_sum['post_stdev']
    df = df.T
    
    return df

added_obs_all_std = add_post_std_df(added_obs_all_std)
added_obs_6gps_std = add_post_std_df(added_obs_6gps_std)
added_obs_3gps_std = add_post_std_df(added_obs_3gps_std)
added_obs_2gps_std = add_post_std_df(added_obs_2gps_std)


# Combine dataframes to obtain Fig. 9

added_obs_9gps_std = copy.deepcopy(added_obs_6gps_std) # Copy the dataframe with the 6 observations groups (added_obs_6gps_std)
added_obs_9gps_std = added_obs_9gps_std.append(added_obs_3gps_std.loc[['h', 'zgeo'], :]) # Append the results from 2 other groups (All head observations + All geophysical interface observations)
added_obs_9gps_std = added_obs_9gps_std.append(added_obs_2gps_std.loc[['zeta'], :]) # Append the results from 1 other group (All interface observations)

#%%------------------------------- Plot Fig. 9 --------------------------------

''' Bar plot: Percent decrease in prior forecast uncertainty (standard deviation) 
when one or several observation groups is added to the initially empty calibration dataset '''

fpst.plot_fig9(post_for, added_obs_9gps_std.T, 'fw_storage') # Fig. 9A
fpst.plot_fig9(post_for, added_obs_9gps_std.T, 'zmuni1') # Fig. 9B

# Note: this function can plot all forecasts ('fw_storage', 'zmuni1', ... 'zmuni9')

#%%------------------------------ Plot Fig. A.1 -------------------------------

''' Map: Percent decrease in prior forecast uncertainty (standard deviation) 
when an individual observation is added to the initially empty calibration dataset '''

#---- Processing

# Normalize absolute decrease (in prior forecast uncertainty) to percent decrease
added_obs_all_std_pct = fpst.normalize_df(added_obs_all_std, 'pct_decrease')

# Add a 'name_model' column to the added_obs dataframe, to be able to use the merge function in fpst.plot_figA1
added_obs_all_std_pct['name_model'] = added_obs_all_std_pct.index

#---- Plot

fpst.plot_figA1(added_obs_all_std_pct, 'zmuni1', gis_path, hshallow_gp, hdeep_gp, hmuni_gp, zwells_gp) # Fig. A.1

# Note: this function can plot all forecasts ('fw_storage', 'zmuni1', ... 'zmuni9')