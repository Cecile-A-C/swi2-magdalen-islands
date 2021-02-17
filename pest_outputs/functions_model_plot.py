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
import flopy, pickle, copy, pyemu
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter

#---- Load relevant data

# Retrieve data from the bin file generated with QGridder
with open(path_bin_in, 'rb') as handle:
     objects = pickle.load(handle)
     (nrow, ncol, delr, delc, ibound, sea, dem, geology, thk_valley, cx, cy, 
      muni_wells_row_col, muni_pumping, domestic_wells_row_col, ind_old_wells_row_col,
      obs_wells_row_col, obs_head, obs_zeta, tdem_row_col, ert_row_col) = objects

#---- Define functions

def cm2inch(*tupl):
    ''' From a tuple containing centimeters, return a tuple containing inches '''
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def fmt_label(x, pos):
    ''' Function later used to label the colorbar using scientic notation '''
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def plot_fig6a(mf, hk_pp_file, mode):
    
    # Export parameters
    fontsze = 11
    figsze = cm2inch(14,14)
    plt.rcParams['font.size'] = fontsze
    plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Arial']})
        
    # Prepare figure
    hk_arrayP = copy.deepcopy(mf.lpf.hk.array)
    hk_arrayP[hk_arrayP == 1.0e30] = np.nan
    hk_pp = pyemu.pp_utils.pp_file_to_dataframe(hk_pp_file) # Read pilot point file
    
    # Start the plot
    fig = plt.figure(figsize=figsze)
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    fig.subplots_adjust(wspace=0.2, hspace=0.2, left=0.125, right=1, bottom=0, top=1)
    
    mapview = flopy.plot.PlotMapView(model=mf) # Create instance of the PlotMapView class
    
    # Plot hydraulic conductivity grid
    im = mapview.plot_array(hk_arrayP)
    cbar = plt.colorbar(im, cmap = plt.cm.viridis, label = r'$K$' + ' (m/s)', shrink = 0.5,
                 format = ticker.FuncFormatter(fmt_label))
    cbar.ax.tick_params(labelsize=fontsze - 1)
    
    # Plot coastline
    mapview.contour_array(sea, colors='black', levels=np.array([0.5]), linewidths=0.5)
    
    # Plot pilot points
    ax.plot(hk_pp[hk_pp.zone==1].x, hk_pp[hk_pp.zone==1].y, 'k+', 
            label='Pilot points (%1d)' % ( len(hk_pp[hk_pp.zone==1]) ) )
        
    # Other figure parameters
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))    
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title('Hydraulic conductivity (' + r'$K$' + ') map', fontsize = fontsze)
    ax.legend(loc = 'lower right')
    
    plt.tight_layout()
    plt.show()
    plt.savefig('Fig6A_' + mode + '.png', dpi = 600)


def plot_fig6b(mf, zeta, gis_path, mode):
    
    # Export parameters
    fontsze = 11
    figsze = cm2inch(14,14)
    plt.rcParams['font.size'] = fontsze
    plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Arial']})
    cmap = 'Blues'
    
    # Set up labels
    title = 'Interface elevation (' + r'$\zeta$' + ') map'
    label = r'$\zeta$' + ' (masl)'
    isocontours = 10
    cmap = cmap + '_r'
    round_to = 0
    fmt = '%1d'
    
    plt.rcParams['contour.negative_linestyle'] = 'solid' # Negative contour lines are solid and not dashed
    
    # Start plot
    fig = plt.figure(figsize=figsze)
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    fig.subplots_adjust(wspace=0.2, hspace=0.2, left=0.125, right=1, bottom=0, top=1)
    
    mapview = flopy.plot.PlotMapView(model=mf) # Create instance of the PlotMapView class
    
    # Plot model output (interface elevation grid)
    im = mapview.plot_array(zeta, cmap = cmap)
    
    # Properties of the colorbar
    cbar = plt.colorbar(im, label = label, shrink = 0.5, format = fmt)
    cbar.ax.tick_params(labelsize=fontsze - 1)
    
    # Plot contourlines
    mn = round(np.nanmin(zeta), round_to)
    mx = round(np.nanmax(zeta), round_to)
    lvl = np.arange(mn, mx, isocontours) # Isocontours to plot
    cn = mapview.contour_array(zeta, colors='black', linewidths=0.1, levels=lvl) # Plot isocontours
    
    # Label the contourlines
    plt.clabel(cn, colors='k', fmt=fmt, fontsize = fontsze - 2)
    
    # Plot coastline
    mapview.contour_array(sea, colors='black', levels=np.array([0.5]), linewidths=0.5)
    
    # Plot limit of the model (boundary between active & inactive cells)
    mapview.contour_array(ibound, colors='black', levels=np.array([0.5]), linewidths=0.5, linestyles = 'dashed')
    
    # Other figure parameters
    ax.set_title(title, fontsize = fontsze)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.0f')%(x/1e3)))  
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
        
    plt.tight_layout()
    plt.show()
    plt.savefig('Fig6B_' + mode + '.png', dpi = 600)