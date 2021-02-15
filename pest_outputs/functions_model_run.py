# -*- coding: utf-8 -*-
"""
Created on Mon May 25 11:08:02 2020

@author: CECOU50
"""

# Load python libraries
import itertools
import pandas as pd

# Compute the equivalent perpendicular hydraulic conductivity = arithmetic mean of individual hydraulic condutivities
def arithmetic_K(hk_1,thk1,hk_2,thk2):
    hk_arithm = ( hk_1*thk1 + hk_2*thk2 ) / (thk1+thk2)
    return hk_arithm

# Float formatter function to apply when exporting data
FFMT = lambda x: "{0:<20.10E} ".format(float(x))
# <: align left, 20: length of output, 10: # of charact after decimal point, E: sci notation

# String formatter function to apply when exporting data
def SFMT(item):
    try:
        s = "{0:<20s} ".format(item.decode())
    except:
        s = "{0:<20s} ".format(str(item))
    return s

# Define which functions to apply to various exported columns
DIC_FMT = {'name': SFMT, 'value': FFMT, 'ins': SFMT, 'tpl': SFMT}

# Width of values (number of characters)
VAL_START = 12
VAL_CHAR_LEN = 30

# Create unsorted dataframe with index=name and columns=row, col, for specific model cells
def get_row_col_id(dict_row_col):
    rows=[]
    cols=[]
    names=[]
    for elem in dict_row_col.values():
        row, col, corr = itertools.chain(elem[0]) # Unpack tuple
        rows.append(row) # Create list for rows
        cols.append(col) # Create list for columns
    for name in dict_row_col.keys():
        names.append(name) # Create list for names
    df_row_col = pd.DataFrame(data={'row': rows, 'col': cols, 'name': names}) # Create df
#    return(rows, cols, names, df_row_col)
    return df_row_col

# Create sorted df containing the data simulated at specific model cells
def sim_to_df(dict_row_col, model_output, tstp): # (model_output = h_swiP or zetaP)
    df_sim = get_row_col_id(dict_row_col) # Create df with index=name, columns= row, col
    df_sim['value'] = model_output[tstp, 0, df_sim['row'], df_sim['col']] # Add column with model_output for selected timestep
    df_sim = df_sim.sort_values(by=['name'], ascending=True) # Sort data with names in ascending order
    df_sim = df_sim.reset_index(drop=True) # Reset index from 0 after sorting, and do not turn the old index into a column
    return(df_sim)

# Export dataframe to .dat file
def write_df_to_dat(file_path, file_name, df, columns): # (file_name = string, columns=list of strings)
    dest_file = open(file_path + file_name + '.dat','w')
    dest_file.write(df.to_string(col_space=0, columns=columns, 
                                 formatters=DIC_FMT, justify="left", 
                                 header=False, index=False))
    dest_file.close()

# Write parameters to PEST .tpl file
def write_df_to_tpl(file_path, df, columns):
    f_tpl = open(file_path + 'param.tpl','w')
    f_tpl.write("ptf ~\n")
    f_tpl.write(df.to_string(col_space=0, columns=columns, formatters=DIC_FMT, 
                             justify='left', header=False, index=False))

# Write dataframe to PEST .ins file
def write_df_to_ins(file_path, file_name, df, columns):  
    f_ins = open(file_path + file_name + '.ins','w')
    f_ins.write("pif #\n")
    f_ins.write(df.to_string(col_space=0, columns=columns, formatters=DIC_FMT, 
                             justify='left', header=False, index=False))

# In the ERT observations dataframe, remove all rows between two indexes based on the 'name' column
def rmRowsERT(df, col, StartString, EndString): # df=dataframe, col=df.column
    df.reset_index(drop=True, inplace=True)
    StartRow = df[col == StartString].index[0]
    EndRow   = df[col == EndString].index[0] + 1
    df.drop(df.index[StartRow : EndRow], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df
