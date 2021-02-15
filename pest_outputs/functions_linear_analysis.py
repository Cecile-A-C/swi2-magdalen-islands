# Author: Cecile Coulon
#
# ---------------------------------- Readme -----------------------------------
# 1) Edit the script preamble (just below)
# 2) No other modifications needed to the script below
#
#----------------------------- 1) Script preamble -----------------------------
# Absolute path to the bin file generated with QGridder
path_bin_in = 'C:/Users/CECOU50/Documents/Cecile Coulon/4_Model Grande Entree/Model cecile/GitHub/pest_outputs/preproc_IDM_20m_qgis_v6.bin' # Path to bin file generated with GIS

#--------------------------------- 2) Script ----------------------------------

# Load python libraries
import pickle, fiona, copy, flopy, pyemu
import pandas as pd
import numpy as np
import geopandas as gpd
from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Load functions
from functions_model_run import get_row_col_id, sim_to_df

#---- Load relevant data

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

# Create df name, row, col dataframes for municipal wells & observations
muni_df = get_row_col_id(muni_wells_row_col)
hwells_df = get_row_col_id(obs_h_row_col)
zwells_df = get_row_col_id(obs_z_row_col)
tdem_df = get_row_col_id(tdem_row_col)
ert_df = get_row_col_id(ert_row_col)

#---- Define functions

def cm2inch(*tupl):
    ''' From a tuple containing centimeters, return a tuple containing inches '''
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def plot_fig5a(pst):
    
    # Export parameters
    fontsze = 11
    figsze = cm2inch(9.5, 9.5)
    plt.rcParams['font.size'] = fontsze
    plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Arial']})
    
    # Load residuals & observations dataframes
    res = pst.res # Residuals
    obs = pst.observation_data # Observations
    
    # Group observations
    grouper = obs.groupby(obs.obgnme).groups
    
    # Drop the observation groups that are not to be plotted
    [grouper.pop(key, None) for key in ['fw_storage', 'hdomest']]
    [grouper.pop(key, None) for key in ['zwells', 'ztdem', 'zert']]
    
    # For each remaining observation group, define legend and marker type, size, color, edgecolor, linewidth
    grouper_lgd = {'hdeep': r'$h_{\rm f\ deep\ wells}$', 'hshallow': r'$h_{\rm f\ shallow\ wells}$','hmuni': r'$h_{\rm f\ pumping\ wells}$'}    
    grouper_mkr = {'hdeep': '^', 'hshallow': 'o','hmuni': 's'}
    grouper_sz = {'hdeep': 120, 'hshallow': 70,'hmuni': 50}
    grouper_clr = {'hdeep': 'red', 'hshallow': 'lime', 'hmuni': 'black'}
    grouper_edgeclr = {'hdeep': 'black', 'hshallow': 'black','hmuni': None}
    grouper_linewdth = {'hdeep': 0.5, 'hshallow': 0.5, 'hmuni': None}
    
    # Start plotting
    fig, ax = plt.subplots(figsize=figsze)
    
    # For each observation group
    for groupname, obsnames in grouper.items():
        
        obs_g = obs.loc[obsnames, :]
        obs_g.loc[:, 'sim'] = res.loc[obsnames, 'modelled']
        obs_g.loc[:,'res'] = obs_g.sim - obs_g.obsval # Define residual = sim - obs
        obs_g = obs_g.loc[obs_g.weight > 0, :] # Exclude observations with weight = 0
        
        ax.scatter([obs_g.obsval], [obs_g.sim], 
                   label=grouper_lgd[groupname], 
                   marker=grouper_mkr[groupname], s=grouper_sz[groupname], 
                   color=grouper_clr[groupname], edgecolors=grouper_edgeclr[groupname],
                   linewidths=grouper_linewdth[groupname])        
    
    # Modify limits of the scatter plot
    Max = 3
    Min = 0
    ax.plot([Min,Max], [Min,Max], 'k--', lw=1.0) # Plot the 1:1 diagonal line
    ax.set_xlim(Min,Max)
    ax.set_ylim(Min,Max)
        
    # Other figure parameters
    ax.axis('square')
    ax.grid()
    ax.set_axisbelow(True)
    
    # Order the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [2,1,0]
    ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper left')
    
    ax.set_xlabel('Observed (m)', labelpad=0.1)
    ax.set_ylabel('Simulated (m)', labelpad=0.1)
    ax.set_title('Freshwater heads ' + r'$h_{\rm f}$', loc='left', fontsize=fontsze)
    
    plt.tight_layout()
    plt.show()
    plt.savefig('Fig5A.png', dpi = 600)


def plot_fig5b(pst):
    
    # Export parameters
    fontsze = 11
    figsze = cm2inch(9.5, 9.5)
    plt.rcParams['font.size'] = fontsze
    plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Arial']})
    
    # Load residuals & observations dataframes
    res = pst.res # Residuals
    obs = pst.observation_data # Observations
    
    # Group observations
    grouper = obs.groupby(obs.obgnme).groups
    
    # Drop the observation groups that are not to be plotted
    [grouper.pop(key, None) for key in ['fw_storage', 'hdomest']]
    [grouper.pop(key, None) for key in ['hdeep', 'hshallow', 'hmuni']]
    
    # For each remaining observation group, define legend and marker type, size, color, edgecolor, linewidth
    grouper_lgd = {'zwells': r'$\zeta_{\rm\ deep\ wells}$', 'ztdem': r'$\zeta_{\rm\ TDEM}$','zert': r'$\zeta_{\rm\ ERT}$'}
    grouper_mkr = {'zwells': '^', 'ztdem': 'p','zert': 'D'}
    grouper_sz = {'zwells': 120, 'ztdem': 70,'zert': 8}
    grouper_clr = {'zwells': 'cyan', 'ztdem': 'cornflowerblue','zert': 'navy'}
    grouper_edgeclr = {'zwells': 'black', 'ztdem': 'black','zert': None}
    grouper_linewdth = {'zwells': 0.5, 'ztdem': 0.5, 'zert': None}
    
    # Start plotting
    fig, ax = plt.subplots(figsize=figsze)
    
    # For each observation group
    for groupname, obsnames in grouper.items():
        
        obs_g = obs.loc[obsnames, :]
        obs_g.loc[:, 'sim'] = res.loc[obsnames, 'modelled']
        obs_g.loc[:,'res'] = obs_g.sim - obs_g.obsval # Define residual = sim - obs
        obs_g = obs_g.loc[obs_g.weight > 0, :] # Exclude observations with weight = 0
        
        ax.scatter([obs_g.obsval], [obs_g.sim], 
                   label=grouper_lgd[groupname], 
                   marker=grouper_mkr[groupname], s=grouper_sz[groupname], 
                   color=grouper_clr[groupname], edgecolors=grouper_edgeclr[groupname],
                   linewidths=grouper_linewdth[groupname])
    
    # Modify limits of the scatter plot
    Max = 0
    Min = -90
    ax.plot([Min,Max], [Min,Max], 'k--', lw=1.0) # Plot the 1:1 diagonal line
    ax.set_xlim(Min,Max)
    ax.set_ylim(Min,Max)

    # Other figure parameters
    ax.axis('square')
    ax.grid()
    ax.set_axisbelow(True)
    
    # Order the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [2,1,0]
    ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper left')
    
    # Other figure parameters
    ax.set_xlabel('Observed (m)', labelpad=0.1)
    ax.set_ylabel('Simulated (m)', labelpad=0.1)
    ax.set_title('Interface elevations ' + r'$\zeta$', loc='left', fontsize=fontsze)
    
    plt.tight_layout()
    plt.show()
    plt.savefig('Fig5B.png', dpi = 600)


def get_forecasts_from_sim(name_model):
    
    ''' Read output files of a MODFLOW-SWI2 simulation and return a dataframe containing all model forecasts
    Input: name_model = model name (string)
    Output: df_fcsts = dataframe with index = forecast name and column = forecast value
    '''
    
    #----------------------------- Read outputs ------------------------------
    hds_swi = flopy.utils.binaryfile.HeadFile(name_model + '.hds') # Read the .hds modflow binary output file, heads are saved here
    hswi = hds_swi.get_alldata()[:,:,:,:] # Either get ALL head data from EVERY simulation time
    zta = flopy.utils.binaryfile.CellBudgetFile(name_model + '.zta') # Read the .zta modflow binary output file, zeta values are saved here
    kstpkper = zta.get_kstpkper() #  Get a list of unique stress periods and timesteps in the file (nper,nstp) = dictionary keys of spd (defined in Output Control package)
    zeta = []
    for kk in kstpkper: # Get ALL zeta data from EVERY simulation time (EVERY (nper,nstp) couple)
        zeta.append(zta.get_data(kstpkper=kk, text='ZETASRF  1')[0])
    zeta = np.array(zeta)
    
    #------------------------- Post-process outputs --------------------------
    hswiP = copy.deepcopy(hswi) # Post-process hswi
    for i in range(len(hswiP)):
        hswiP[i,0,:,:][ ibound == 0 ] = np.nan # Assign hswi = 0 to non active cells
    zetaP = copy.deepcopy(zeta) # Post-process zeta
    for i in range(len(zetaP)):
        zetaP[i,0,:,:][ ibound == 0 ] = np.nan # Assign zeta = 0 to non active cells
    
    #------------- Correct simulated h & z at muni pumping wells --------------
    
    #-------------------------- Compute FW storage ---------------------------
    thk_FW_lens = (hswiP-zetaP) # Thickness of FW lens (m)
    FW_storage = [] # Freshwater storage in lens (m3)
    for i in range(len(hswi)):
        fw = np.nansum(thk_FW_lens[i,0,:,:])*delr*delc # FW storage in the FW lens = ∑(hswi-zeta)*delr*delc
        FW_storage.append(fw)
    
    #-------------------- Get h, zeta at municipal wells ---------------------
    # Create dataframes containing the data simulated at observation points
    hmuni = sim_to_df(muni_wells_row_col, hswiP, -1)
    zmuni = sim_to_df(muni_wells_row_col, zetaP, -1)
    # Name of forecasts
    hmuni.name = 'h' + hmuni.name
    zmuni.name = 'z' + zmuni.name
    
    #------- Return dataframe with all model forecasts from simulation -------
    df_fcsts = pd.concat([hmuni, zmuni])
    df_fcsts.index = df_fcsts.name
    df_fcsts.drop(columns=['row', 'col', 'name'], inplace=True)
    df_fcsts.loc['fw_storage'] = FW_storage[-1]
    df_fcsts.sort_index(inplace=True)
    
    return df_fcsts


def plot_fig7(df_prior, df_post, fcstnme):
    
    # Export parameters
    fontsze = 11
    figsze = cm2inch(9.5,9.5)
    plt.rcParams['font.size'] = fontsze
    plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Arial']})
    
    # Prepare dataframes
    df_prior_subset = df_prior[df_prior.index.str.startswith(fcstnme)]
    df_post_subset = df_post[df_post.index.str.startswith(fcstnme)]
        
    # Set up labels
    if 'zmuni' in fcstnme:
        well_no = fcstnme.split('muni')[1] # Get well number
        clean_nme = r'$\zeta$'
        unit = 'masl'
        Fmt = '%0.0f'
        xlabel = 'Interface elevation ' + r'$\zeta_{\rm muni} $' + well_no + ' (' + unit + ')'
    
    elif fcstnme == 'fw_storage':
        df_prior_subset = df_prior_subset / (10**9) # From m3 to km3
        df_post_subset = df_post_subset / (10**9) #   From m3 to km3
        clean_nme = r'$V_{\rm f}$'
        unit = 'km$^3$'
        Fmt = '%0.1f'
        xlabel = 'Freshwater volume ' + r'$V_{\rm f}$' + ' (' + unit + ')'
    
    # Start plot
    plt.figure(figsize=figsze)
    ax = plt.subplot(1,1,1)
    
    # Plot prior distribution (Gaussian)
    mean = df_prior_subset.Mean[0]
    stdev = df_prior_subset.Stdev[0]
    
    x1, y1 = pyemu.plot_utils.gaussian_distribution(mean, stdev)
    label = 'Prior ' + clean_nme + ': ' + Fmt % (mean) + ' \u00B1 ' + Fmt % (2 * stdev) + ' ' + unit
    ax.plot(x1, y1, linestyle='dashed', color = 'grey', label = label)
    
    # Plot posterior distribution (Gaussian)
    mean = df_post_subset.Mean[0]
    stdev = df_post_subset.Stdev[0]
    
    x2, y2 = pyemu.plot_utils.gaussian_distribution(mean, stdev)
    label = 'Posterior ' + clean_nme + ': ' + Fmt % (mean) + ' \u00B1 ' + Fmt % (2 * stdev) + ' ' + unit
    ax.fill_between(x2, 0, y2, alpha=0.5, color = 'grey', label = label)
    
    
    # Figure limits
    Max = max(y1.max(), y2.max())
    ax.set_ylim(0, Max + Max * 0.3)
    
    # Other figure parameters
    ax.legend(loc = 'upper left')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Probability density')
    
    plt.tight_layout()
    plt.show()
    plt.savefig('Fig7_' + fcstnme + '.png', dpi = 600)


def plot_fig8(post_for, df_abs, fcstnme):
    
    # Export parameters
    fontsze = 11
    figsize = cm2inch(9.5, 8)
    plt.rcParams['font.size'] = fontsze
    plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Arial']})
    
    # Prepare dataframe
    df_abs_subset = df_abs[df_abs.index.str.startswith(fcstnme)]
    
    # Replace column names with the following list (for the legend)
    cols = [r'$\sigma_{\rm\ post}$', #          Posterior forecast stdev
            r'If $K$ field' + '\nknown', #      hcond parameter group
            r'If $R$ known', #                  recharge parameter group
            r'If $\alpha_{\rm T}$ known'] #     alpha parameter group
    df_abs_subset.columns = cols
    
    # Sort dataframe in descending order
    df_abs_subset = df_abs_subset.T
    df_abs_subset = df_abs_subset.sort_values(by=fcstnme, ascending=False)
    df_abs_subset = df_abs_subset.T
    
    # Set up a different color for each parameter group
    colors=[]
    for col in df_abs_subset.columns:
        if col == cols[0]: #    posterior forecast stdev
            c = 'k'
        if col == cols[1]: #    hcond parameter group
            c = 'tab:red'
        elif col == cols[2]: #  recharge parameter group
            c = 'tab:blue'
        elif col == cols[3]: #  alpha parameter group
            c = 'tab:green'
        colors.append(c)
    
    # Set up labels
    if 'muni' in fcstnme:
        well_no = fcstnme.split('muni')[1] # Get well number
        unit_stdev = 'm'
        unit_fcst = 'masl'
        Fmt = '%0.0f'
        fcst = post_for.loc[fcstnme, 'Mean']
        title = 'Interface elevation ' + r'$\zeta_{\rm muni}$ ' + well_no + \
            ' (%0.1f '  % (fcst) + unit_fcst + ')'
    
    elif fcstnme == 'fw_storage':
        df_abs_subset = df_abs_subset / (10**9)
        unit_stdev = 'km$^3$'
        unit_fcst = unit_stdev
        Fmt = '%0.1f'
        fcst = post_for.loc[fcstnme, 'Mean']  / (10**9)
        title = 'Freshwater volume ' + r'$V_{\rm f}$ (%0.1f '  % (fcst) + unit_fcst + ')'
    
    # Bar plot
    ax = df_abs_subset.plot(kind = 'bar', figsize = figsize, rot = 45, grid = True, 
                         color = colors, legend = None, width = 15, edgecolor = 'white', linewidth = 5)
    
    # Annotate the first bar differently (base posterior forecast uncertainty)
    for p in ax.patches[0:1]:
        sigma_base = p.get_height()
        ann = df_abs_subset.columns[0]
        ax.annotate(ann, (p.get_x() * 1.005, p.get_height() * 0.980), va='bottom', fontsize = fontsze - 1)
    
    # Annotate the other bars with 1) legend and 2) percent reduction in forecast uncertainty
    for i, p in enumerate(ax.patches[1:]):
        sigma = p.get_height()
        label = 100 * (sigma_base - sigma) / sigma_base # Compute percent reduction in forecast uncertainty
        ann = df_abs_subset.columns[i+1] + '\n' + '- ' + Fmt % (label) + ' %%'
        ax.annotate(ann, (p.get_x() * 1.005, p.get_height() * 0.980),  va='bottom', fontsize = fontsze - 1)
    
    # Figure limits
    ax.set_ylim(0, sigma_base + sigma_base * 0.15)
    
    # Other figure parameters
    ax.set_axisbelow(True)
    x_axis = ax.axes.get_xaxis()
    x_axis.set_visible(False)
    ax.set_ylabel('Standard deviation ' + r'$\sigma$' + ' (' + unit_stdev + ')')
    plt.title(title, loc='left', fontsize = fontsze)
    
    plt.tight_layout()
    plt.show()
    plt.savefig('Fig8_' + fcstnme + '.png', dpi = 600)


def plot_fig9(post_for, df_stdev, fcstnme):
    
    # Export parameters
    fontsze = 11
    plt.rcParams['font.size'] = fontsze
    plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Arial']})
    figsze = cm2inch(19, 6)
    
    # Prepare dataframe
    df_stdev_subset = df_stdev[df_stdev.index.str.startswith(fcstnme)]
    
    # Replace column names with the following list (for the legend)
    cols = ['No data ' + r'($\sigma_{\rm prior}$)', #    Prior forecast stdev (no observations)
            r'+ $h_{\rm shallow}$', #                    h shallow observation group
            r'+ $h_{\rm deep}$', #                       h deep observation group
            r'+ $h_{\rm muni}$', #                       h muni observation group
            r'+ $\zeta_{\rm deep}$', #                   zeta deep well observation group  
            r'+ $\zeta_{\rm TDEM}$', #                   zeta TDEM observation group
            r'+ $\zeta_{\rm ERT}$', #                    zeta ERT observation group
            'All data\n' + r'($\sigma_{\rm post}$)', #   Posterior forecast stdev (all observations)
            '+ All ' + r'$h$', #                         All h observations
            '+ All\n' + r'$\zeta_{\rm geophy}$', #       All geophysical zeta observations
            '+ All ' + r'$\zeta$'] #                     All zeta observations
    df_stdev_subset.columns = cols
    
    # Sort dataframe in descending order
    df_stdev_subset = df_stdev_subset.T
    df_stdev_subset = df_stdev_subset.sort_values(by=fcstnme, ascending=False)
    df_stdev_subset = df_stdev_subset.T
    
    # Set up a different color for each observation group
    colors=[]
    for col in df_stdev_subset.columns:
        if col == cols[0]: # prior forecast stdev
            c = 'tab:orange'
        if col == cols[1]: # hshallow
            c = 'chartreuse'
        elif col == cols[2]: # hdeep
            c = 'red'
        elif col == cols[3]: # hmuni
            c = 'grey'
        elif col == cols[4]: # zdeep
            c = 'aquamarine'
        elif col == cols[5]: # ztdem
            c = 'cornflowerblue'
        elif col == cols[6]: # zert
            c = 'tab:blue'
        elif col == cols[7]: # posterior forecast stdev
            c = 'k'
        elif col == cols[8]: # hobs
            c = 'tab:green'
        elif col == cols[9]: # zgeo
            c = 'tab:purple'
        elif col == cols[10]: # zobs
            c = 'navy'
        colors.append(c)
    
    # Set up labels
    if 'muni' in fcstnme:
        well_no = fcstnme.split('muni')[1]
        unit_stdev = 'm'
        unit_fcst = 'm/sea level'
        fcst = post_for.loc[fcstnme, 'Mean']
        title = 'Interface elevation ' + r'$\zeta_{\rm muni} $ ' + well_no + \
            ' (%0.1f '  % (fcst) + unit_fcst + ')'
    
    elif fcstnme == 'fw_storage':
        df_stdev_subset = df_stdev_subset / (10**9)
        unit_stdev = 'km$^3$'
        unit_fcst = unit_stdev
        fcst = post_for.loc[fcstnme, 'Mean']  / (10**9)
        title = 'Freshwater volume ' + r'$V_{\rm f}$ (%0.1f '  % (fcst) + unit_fcst + ')'
    
    
    # Bar plot
    ax = df_stdev_subset.plot(kind = 'bar', figsize = figsze, rot = 45, grid = True, 
                         color = colors, legend = None, 
                         width = 15, edgecolor = 'white', linewidth = 5) # Cycle through default colors
    
     # Annotate the first bar differently (base prior forecast uncertainty)
    for p in ax.patches[0:1]:
        sigma_base = p.get_height()
        ann = df_stdev_subset.columns[0] #stdev_symbol
        ax.annotate(ann, (p.get_x() * 1.005, p.get_height() * 0.960), va='bottom', fontsize = fontsze - 1)
    
    # Annotate the other bars with 1) legend and 2) percent reduction in prior forecast uncertainty
    for i, p in enumerate(ax.patches[1:]):
        sigma = p.get_height()
        label = 100 * (sigma_base - sigma) / sigma_base
        ann = df_stdev_subset.columns[i+1] + '\n' + '-%0.1f %%' % (label)
        ax.annotate(ann, (p.get_x() * 1.005, p.get_height() * 0.900), va='bottom', fontsize = fontsze - 1)
    
    # Figure limits
    ax.set_ylim(0, sigma_base + sigma_base * 0.1)
    ax.set_axisbelow(True)
    
    # Other figure parameters
    x_axis = ax.axes.get_xaxis()
    x_axis.set_visible(False)
    ax.set_ylabel('Standard deviation ' + r'$\sigma$' + ' (' + unit_stdev + ')')
    plt.title(title, loc='left', fontsize = fontsze)
    
    plt.tight_layout()
    plt.show()
    plt.savefig('Fig9_' + fcstnme + '.png', dpi = 600)


def normalize_df(df, Type):
    ''' Normalize all the standard deviations in a dataframe (df) to the maximum standard deviation 
    Type = 'pct_decrease' or 'pct_increase' (string): whether to compute percent decrease or increase '''
    
    df_base = df.loc['base', :].copy() # Base row = prior forecast variance resulting from all
    
    # Compute percent decrease
    if Type == 'pct_decrease':
        df_percent = 100.0 * (df.apply(lambda x: (df_base - x) / df_base, axis=1))
    
    # Compute percent increase
    elif Type == 'pct_increase':
        df_percent = 100.0 * (df.apply(lambda x: (x - df_base) / df_base, axis=1))
    
    # Drop the base row
    df_percent.drop('base', inplace=True)
    
    return df_percent


def records(filename, usecols, **kwargs):
    
    ''' Read only a specific attribute column of a shapefile '''
    
    with fiona.open(filename, **kwargs) as source:
        for feature in source:
            f = {k: feature[k] for k in ['id', 'geometry']}
            f['properties'] = {k: feature['properties'][k] for k in usecols}
            yield f


def get_coordinates(df_row_col):
    
    ''' Get the x and y coordinates for the cells specified in (name, row, col) dataframes
    Input: df_row_col = dataframe with columns = name, row, col, for certain locations
    Output: coord_dict = dictionary of tuples with {name: (x-coordinate, y-coordinate)} '''
    
    xcoord = [cx[row, col] for row,col in zip(df_row_col.row, df_row_col.col)]
    ycoord = [cy[row, col] for row,col in zip(df_row_col.row, df_row_col.col)]
    
    coord_dict = dict( zip(df_row_col.name, zip(xcoord, ycoord)) )
    
    return coord_dict


def merge_df_with_gpd(pst, df, gis_path):
    
    # Read wells_obs.shp
    hwells = gpd.GeoDataFrame.from_features(records(gis_path + 'wells_obs' + '.shp', ['obs_name _', 'obs_name_m']))
    hwells.rename(columns={'obs_name_m': 'name_model', 'obs_name _': 'name_ori'}, inplace=True)
    zwells = copy.deepcopy(hwells)
    
    # Create dataframe for all h observations
    hwells.name_model = 'h' + hwells.name_model.astype(str)
    
    # Extract hshallow, hdeep and hmuni groups from hwells
    list_hshallow = pst.observation_data.loc[pst.observation_data.obgnme == 'hshallow']['obsnme'].tolist()
    list_hdeep = pst.observation_data.loc[pst.observation_data.obgnme == 'hdeep']['obsnme'].tolist()
    hshallow = hwells[hwells.name_model.apply(lambda x : x in list_hshallow)]
    hdeep = hwells[hwells.name_model.apply(lambda x : x in list_hdeep)]
    hmuni = hwells.loc[hwells.name_model.str.startswith('hmuni')]
    
    # Extract and create dataframe for zeta well observations
    zwells.name_model = 'z' + zwells.name_model.astype(str)
    list_zwells = pst.observation_data.loc[pst.observation_data.obgnme == 'zwells']['obsnme'].tolist()
    zwells = zwells[zwells.name_model.apply(lambda x : x in list_zwells)]
    
    # Read TDEM data
    ztdem = gpd.GeoDataFrame.from_features(records(gis_path + 'tdem_interface' + '.shp', ['Name', 'Name_model']))
    ztdem.Name_model = 'z' + ztdem.Name_model.astype(str)
    ztdem.rename(columns={'Name_model': 'name_model', 'Name': 'name_ori'}, inplace=True)
    
    # Read ERT data
    zert = gpd.GeoDataFrame.from_features(records(gis_path + 'ert_processed_gis' + '.shp', ['agg_name', 'model_name']))
    zert.rename(columns={'model_name': 'name_model', 'agg_name': 'name_ori'}, inplace=True)
    zert = zert[['name_ori', 'name_model', 'geometry']]
    
    # For each observation group, merge dataframe with geopandas information
    df.rename(columns={'name': 'name_model'}, inplace=True)
    hshallow = hshallow.merge(df, on='name_model')
    hdeep = hdeep.merge(df, on='name_model')
    hmuni = hmuni.merge(df, on='name_model')
    hobs = hwells.merge(df, on='name_model') # Dataframe with all h observations
    
    zwells = zwells.merge(df, on='name_model')
    ztdem = ztdem.merge(df, on='name_model')
    zert = zert.merge(df, on='name_model')
    zobs = pd.concat([zwells, ztdem, zert]) # Dataframe with all zeta observations
    
    return(hshallow, hdeep, hmuni, zwells, ztdem, zert, hobs, zobs)


def plot_figA1(pst, df, fcstnme, cmap, gis_path):    
    
    # Export parameters
    fontsze = 11
    figsize = cm2inch(14, 20)
    plt.rcParams['font.size'] = fontsze
    plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Arial']})
    
    # Set up labels
    if 'muni' in fcstnme:
        well_no = fcstnme.split('muni')[1]
        unit_stdev = 'm'
        unit_fcst = 'masl'
        title = 'Interface elevation at municipal well '  + well_no
    
    elif fcstnme == 'fw_storage':
        unit_stdev = 'km$^3$'
        unit_fcst = unit_stdev
        title = 'Freshwater volume (' + r'$V_{\rm f}$'+ ')'
        
    # Load GIS data using gis_path
    wells_muni = gpd.GeoDataFrame.from_features(records(gis_path + 'wells_muni' + '.shp', ['Name', 'Name_model', 'PumpRate']))
    wells_muni.rename(columns={'Name_model': 'name_model', 'Name': 'name_ori'}, inplace=True)
    coastline = gpd.read_file(gis_path + 'polyg_model' + '.shp')
    
    muni_coordinates = get_coordinates(muni_df) # Get coordinates of municipal wells
    xmuni = [coord[0] for name, coord in muni_coordinates.items()]
    ymuni = [coord[1] for name, coord in muni_coordinates.items()]
    
    # Merge dataframe info and geopandas geometry
    hshallow, hdeep, hmuni, zwells, ztdem, zert, hobs, zobs = merge_df_with_gpd(pst, df, gis_path)
    
    figtitle = 'Forecast: ' + title
    cbar_label = '% decrease in prior ' + r'$\sigma$' + '\n when observation used'
        
    # Set up figure
    
    fig, axs = plt.subplots(2, 1, figsize=figsize, sharex=True, sharey=True, subplot_kw=dict(aspect='equal')) 
    fig.suptitle(figtitle)
    axs = axs.ravel()
    
    for i in range(2):
        coastline.plot(ax=axs[i], color='gray', alpha=0.4, edgecolor='k') # Plot coastline
        
        ax_divider = make_axes_locatable(axs[i])
        cax = ax_divider.append_axes("right", size="7%", pad="2%")
        
        # If we're plotting data worth for municipal well forecasts
        
        if 'muni' in fcstnme:
            wells_muni[ wells_muni.name_model.str.endswith(fcstnme[1:]) ].plot(ax=axs[i], 
                       marker='s', markersize=50, facecolors='k', label=None) # Plot municipal wells at which forecast made
            
            idx = [idx for idx, well in enumerate(wells_muni.name_model) if fcstnme[1:] in well][0]
            
            axs[i].annotate(wells_muni.name_model[idx], (xmuni[idx], ymuni[idx]), xytext=(0, 5), 
               textcoords='offset points', ha='center', va='bottom', fontsize=fontsze-5)
        
        # First frame = head observations
        
        if i == 0:
            hobs.plot(ax=axs[i], column=fcstnme, marker='o', markersize=40, edgecolors='k', 
                      linewidths=0.5, cmap=cmap, legend=True, cax=cax, legend_kwds={'label': cbar_label}) # legend=True: tells Geopandas to add the colorbar
            
            # Label observation groups
            
            hshallow.plot(ax=axs[i], marker='o', markersize=40, facecolors='None', 
                          edgecolors='k', linewidths=0.5, label='Shallow wells')
            
            hdeep.plot(ax=axs[i], marker='^', markersize=90, facecolors='None', 
                       edgecolors='white', linewidths=0.75, label='Deep wells')
            
            hmuni.plot(ax=axs[i], marker='s', markersize=40, facecolors='None', 
                       edgecolors='k', linewidths=0.5, label='Municipal wells')
            
            axs[i].set_title('Worth of head observations')
        
        # Second frame = interface observations
        
        elif i == 1:
            zobs.plot(ax=axs[i], column=fcstnme, marker='o', markersize=40, 
                      cmap=cmap, legend=True, cax=cax, legend_kwds={'label': cbar_label}) # markersize=np.abs(zobs.residual)*2
            
            # Label observation groups
            
            zwells.plot(ax=axs[i], marker='^', markersize=90, facecolors='None', 
                        edgecolors='white', linewidths=0.75, label='Deep wells')
            
            ztdem.plot(ax=axs[i], marker='p', markersize=50, facecolors='None', 
                       edgecolors='k', linewidths=0.5, label='TDEM')
            
            zert.plot(ax=axs[i], marker='8', markersize=50, facecolors='None', 
                      edgecolors='gray', linewidths=0.1, label='ERT')
            
            axs[i].set_title('Worth of interface observations')
        
        axs[i].legend(loc='upper left', facecolor='gray', framealpha=0.3)
        
        # Show x and y axis in km instead of m
        
        axs[i].xaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))
        axs[i].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))
        axs[i].set_ylabel('y (km)')
    
    axs[0].set_xlim(300000, 306000)
    axs[0].set_ylim(5267000, 5271000)
    axs[1].set_xlabel('x (km)')
    
#    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))
#    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))    
#    ax.set_xlabel('x (km)')
#    ax.set_ylabel('y (km)')
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    plt.savefig('FigA1_' + fcstnme + '.png', dpi = 300)