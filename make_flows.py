# -*- coding: utf-8 -*-
"""
Created on September 7 2020

@author: S. Sterl & A. Devillers
"""

import numpy as np
import pandas as pd
from pandas import ExcelWriter


# %% 1) Set parameters

# [set] year for which to generate hydropower profiles
SelectedYear = 2020

# [preallocate] binary output variable
output_binary = np.zeros(len(index), dtype = np.bool)
index_name_counted = np.zeros(len(index))
                                          
# [calculate] whether or not to write outputs, given the selected year
for n in range(len(index)):
    if ((FirstYear[n] <= SelectedYear)[0]):
        output_binary[n] = True
    else:
        output_binary[n] = False


# %% 2) Correct for reservoir effects

# [preallocate] the parameter f_reg (see Sterl et al. 2020, https://www.nature.com/articles/s41893-020-0539-0; cf. Suppl. Inf.)
f_reg_normal = np.zeros(len(index))
f_reg_dry = np.zeros(len(index))
f_reg_wet = np.zeros(len(index))

# [preallocate] the bias-corrected outflow including reservoir retention calculation
Q_bymonth_normal_RoR = np.full([len(index), months_yr], np.nan)
Q_bymonth_normal_retained = np.full([len(index), months_yr], np.nan)
Q_bymonth_normal_reservoir = np.full([len(index), months_yr], np.nan)
Q_bymonth_dry_RoR = np.full([len(index), months_yr], np.nan)
Q_bymonth_dry_retained = np.full([len(index), months_yr], np.nan)
Q_bymonth_dry_reservoir = np.full([len(index), months_yr], np.nan)
Q_bymonth_wet_RoR = np.full([len(index), months_yr], np.nan)
Q_bymonth_wet_retained = np.full([len(index), months_yr], np.nan)
Q_bymonth_wet_reservoir = np.full([len(index), months_yr], np.nan)


# [assumption] ratio live storage to total storage for reservoirs
f_live = 0.70

# [calculate] correct for reservoir effect (see Sterl et al. 2020)
for n in range(len(index)):
    
    # [loop] for cases where we know both the reservoir size and the yearly average inflow
    if (reservoir_size[n] > 0) and (bias_correction_yearly[n,0] > 0) and (np.mean(Q_bias_corrected_bymonth_normal_natural[n,:]) > 0):
        
        # [calculate] factor f_reg
        f_reg_normal[n] = f_live*reservoir_size[n]*10**6/(365*24*3600*np.mean(Q_bias_corrected_bymonth_normal_natural[n,:]))
        f_reg_wet[n] = f_live*reservoir_size[n]*10**6/(365*24*3600*np.mean(Q_bias_corrected_bymonth_wet_natural[n,:]))
        f_reg_dry[n] = f_live*reservoir_size[n]*10**6/(365*24*3600*np.mean(Q_bias_corrected_bymonth_dry_natural[n,:]))
        
        # [calculate] divide the flow into a "retainable" and a "run-of-river" component
        # [note] the "retainable" flow is assumed to be constant around the year (stable production)
        # alternative approaches, e.g. seasonal complementarity between hydro and VRE, may require other "retainable" flow formulations
        if f_reg_normal[n] >= 1:
            
            Q_bymonth_normal_RoR[n,:] = 0
            Q_bymonth_normal_retained[n,:] = np.mean(Q_bias_corrected_bymonth_normal_natural[n,:])
            
        elif f_reg_normal[n] < 1:
                        
            Q_bymonth_normal_RoR[n,:] = (1 - f_reg_normal[n])*Q_bias_corrected_bymonth_normal_natural[n,:]
            Q_bymonth_normal_retained[n,:] = f_reg_normal[n]*np.mean(Q_bias_corrected_bymonth_normal_natural[n,:])
        
        if f_reg_dry[n] >= 1:
            
            Q_bymonth_dry_RoR[n,:] = 0
            Q_bymonth_dry_retained[n,:] = np.mean(Q_bias_corrected_bymonth_dry_natural[n,:])
            
        elif f_reg_dry[n] < 1:
            
            Q_bymonth_dry_RoR[n,:] = (1 - f_reg_dry[n])*Q_bias_corrected_bymonth_dry_natural[n,:]
            Q_bymonth_dry_retained[n,:] = f_reg_dry[n]*np.mean(Q_bias_corrected_bymonth_dry_natural[n,:])
            
        if f_reg_wet[n] >= 1:
            
            Q_bymonth_wet_RoR[n,:] = 0
            Q_bymonth_wet_retained[n,:] = np.mean(Q_bias_corrected_bymonth_wet_natural[n,:])
            
        elif f_reg_wet[n] < 1:
                        
            Q_bymonth_wet_RoR[n,:] = (1 - f_reg_wet[n])*Q_bias_corrected_bymonth_wet_natural[n,:]
            Q_bymonth_wet_retained[n,:] = f_reg_wet[n]*np.mean(Q_bias_corrected_bymonth_wet_natural[n,:])

# [calculate] sum up retainable and run-of-river component
temp_normal = Q_bymonth_normal_retained + Q_bymonth_normal_RoR
temp_dry = Q_bymonth_dry_retained + Q_bymonth_dry_RoR
temp_wet = Q_bymonth_wet_retained + Q_bymonth_wet_RoR
for n in range(len(index)):
    if np.isnan(temp_normal[n,0]):
        Q_bymonth_normal_reservoir[n,:] = Q_bias_corrected_bymonth_normal_natural[n,:]
        Q_bymonth_dry_reservoir[n,:] = Q_bias_corrected_bymonth_dry_natural[n,:]
        Q_bymonth_wet_reservoir[n,:] = Q_bias_corrected_bymonth_wet_natural[n,:]
    else:
        Q_bymonth_normal_reservoir[n,:] = temp_normal[n,:]
        Q_bymonth_dry_reservoir[n,:] = temp_dry[n,:]
        Q_bymonth_wet_reservoir[n,:] = temp_wet[n,:]

# [calculate] updated seasonality (after including reservoir effects)
for n in range(len(index)):
    if np.mean(Q_bymonth_normal_reservoir[n,:]) > 0:
        seasonality_normal[n,:] = Q_bymonth_normal_reservoir[n,:] / np.mean(Q_bymonth_normal_reservoir[n,:])
        seasonality_dry[n,:] = Q_bymonth_dry_reservoir[n,:] / np.mean(Q_bymonth_normal_reservoir[n,:])
        seasonality_wet[n,:] = Q_bymonth_wet_reservoir[n,:] / np.mean(Q_bymonth_normal_reservoir[n,:])    


print('flows ready (including reservoir retention effects)')


# %% 3) Cascade setup for Blue/White Nile confluence

# MANUAL ADAPTATION #1
# [identify] plants in Nile valley turbining White Nile and Blue Nile flow - Shereik (1st in cascade) & Merowe (1st existing)
# [assume] outflow profile into Nile is sum of GERD (Blue Nile, Ethiopia) and Jebel Aulia (White Nile, Sudan) flow
# [introduce] one-month lag along Nile valley, between GERD and Merowe
 
# [find] relevant dams in Sudan (Nile)
temp_HPP_adapt_1 = np.where(np.array(index_name) == 'Merowe')
temp_HPP_adapt_2 = np.where(np.array(index_name) == 'Shereik')
temp_HPP_adapt_3 = np.where(np.array(index_name) == 'Sabaloka')

# [find] Jebel Aulia in Sudan (White Nile)
temp_source_1 = np.where(np.array(index_name) == 'Jebel Aulia')

# [find] either Roseires or Renaissance, depending on the year (Blue Nile)
temp_source_2 = np.where(np.array(index_name) == 'Renaissance')
if ((FirstYear[temp_source_2] > SelectedYear)[0][0]):
    temp_source_2 = np.where(np.array(index_name) == 'Roseires')


# [change] Sabaloka inflow to sum of Jebel Aulia and Roseires
Q_bymonth_normal_reservoir[temp_HPP_adapt_3, :] = Q_bymonth_normal_reservoir[temp_source_1, :] + Q_bymonth_normal_reservoir[temp_source_2, :]
Q_bymonth_dry_reservoir[temp_HPP_adapt_3, :] = Q_bymonth_dry_reservoir[temp_source_1, :] + Q_bymonth_dry_reservoir[temp_source_2, :]
Q_bymonth_wet_reservoir[temp_HPP_adapt_3, :] = Q_bymonth_wet_reservoir[temp_source_1, :] + Q_bymonth_wet_reservoir[temp_source_2, :]


# [change] Shereik and Merowe inflow to Sabaloka outflow, with one month lag
Q_bymonth_normal_reservoir[temp_HPP_adapt_2, 0] = Q_bymonth_normal_reservoir[temp_HPP_adapt_3, -1]
Q_bymonth_normal_reservoir[temp_HPP_adapt_2, 1:] = Q_bymonth_normal_reservoir[temp_HPP_adapt_3, 0:-1]
Q_bymonth_dry_reservoir[temp_HPP_adapt_2, 0] = Q_bymonth_dry_reservoir[temp_HPP_adapt_3, -1]
Q_bymonth_dry_reservoir[temp_HPP_adapt_2, 1:] = Q_bymonth_dry_reservoir[temp_HPP_adapt_3, 0:-1]
Q_bymonth_wet_reservoir[temp_HPP_adapt_2, 0] = Q_bymonth_wet_reservoir[temp_HPP_adapt_3, -1]
Q_bymonth_wet_reservoir[temp_HPP_adapt_2, 1:] = Q_bymonth_wet_reservoir[temp_HPP_adapt_3, 0:-1]

Q_bymonth_normal_reservoir[temp_HPP_adapt_1, :] = Q_bymonth_normal_reservoir[temp_HPP_adapt_2, :]
Q_bymonth_dry_reservoir[temp_HPP_adapt_1, :] = Q_bymonth_dry_reservoir[temp_HPP_adapt_2, :]
Q_bymonth_wet_reservoir[temp_HPP_adapt_1, :] = Q_bymonth_wet_reservoir[temp_HPP_adapt_2, :]


print('flows ready (including manual changes to set up cascade calculation)')


# %% 4) Correct for all other cascades

# [preallocate] binary index indicating need to include cascade effect
calculate_cascade = np.zeros(len(index))
find_spill_from = np.zeros(len(index))

# [loop] over all hydropower plants
for n in range(len(index)):
    
    # [identify] plants in direct cascade with another, if these already in service in selected year
    if ((SpillFrom[n] != 'nan')[0]) and ((FirstYear[n] <= SelectedYear)[0]):
        
        # [identify] upstream plant in cascade
        find_spill_from[n] = np.where(index_name == SpillFrom[n])[0][0]
        
        temp_found_upstream = 0
        # [check] if upstream plant in cascade already in service in selected year
        while temp_found_upstream < 1:
            # [check] if the upstream plant in cascade is not yet in service BUT there is another plant a bit further upstream in the cascade, take that one instead
            if ((FirstYear[np.int(find_spill_from[n])] > SelectedYear)[0]) and ((SpillFrom[np.int(find_spill_from[n])] != 'nan')[0]):
                find_spill_from[n] = np.where(index_name == SpillFrom[np.int(find_spill_from[n])])[0][0]
            # [check] if the upstream plant in cascase is already in service
            elif ((FirstYear[np.int(find_spill_from[n])] <= SelectedYear)[0]):
                 temp_found_upstream = 1
                 calculate_cascade[n] = 1
            # [check] if the upstream plant in cascade is not yet in service, and there is effectively no cascade yet
            else:
                temp_found_upstream = 1
                
        if calculate_cascade[n] == 1:
            # [loop] for cases where we know the average flow
            if (bias_correction_yearly[n,0] > 0) and (bias_correction_yearly[np.int(find_spill_from[n]),0] > 0):
                
                # [calculate] seasonality based on hydropower dam upstream in cascade
                seasonality_normal[n,:] = Q_bymonth_normal_reservoir[np.int(find_spill_from[n]),:]/np.mean(Q_bymonth_normal_reservoir[np.int(find_spill_from[n]),:])
                seasonality_dry[n,:] = Q_bymonth_dry_reservoir[np.int(find_spill_from[n]),:]/np.mean(Q_bymonth_normal_reservoir[np.int(find_spill_from[n]),:])
                seasonality_wet[n,:] = Q_bymonth_wet_reservoir[np.int(find_spill_from[n]),:]/np.mean(Q_bymonth_normal_reservoir[np.int(find_spill_from[n]),:])
                
                # [calculate] outflow for hydropower dams in cascade
                Q_bymonth_normal_reservoir[n,:] = bias_correction_yearly[n]*seasonality_normal[n,:]
                Q_bymonth_dry_reservoir[n,:] = bias_correction_yearly[n]*seasonality_dry[n,:]
                Q_bymonth_wet_reservoir[n,:] = bias_correction_yearly[n]*seasonality_wet[n,:]
                
            else:
                
                # [copy] seasonality from upstream cascade dams for cases where we do not know the average flow
                seasonality_normal[n,:] = seasonality_normal[np.int(find_spill_from[n]),:]
                seasonality_dry[n,:] = seasonality_dry[np.int(find_spill_from[n]),:]
                seasonality_wet[n,:] = seasonality_wet[np.int(find_spill_from[n]),:]
            

# [note]
# in cascade systems, the first dam in the cascade typically impounds a reservoir to stabilise downstream flow;
# the other dams then do not need a big reservoir anymore
# for this reason, the reservoir calculation is not repeated for hydropower plants in a cascade

print('flows ready (including cascade effects)')


# %% 5) Calculate seasonality of capacity factor

# [preallocate] the modeled flow at which the power production is maximal
utilisation_rate_normal = np.zeros(shape = (len(index), months_yr))
utilisation_rate_dry = np.zeros(shape = (len(index), months_yr))
utilisation_rate_wet = np.zeros(shape = (len(index), months_yr))

# [preallocate] the output bias-corrected by yearly average
CF_bymonth_normal = np.zeros(shape = (len(index), months_yr))
CF_bymonth_dry = np.zeros(shape = (len(index), months_yr))
CF_bymonth_wet = np.zeros(shape = (len(index), months_yr))

# [define] the ratio between design discharge and maximum flow
f_ratio = 0.5

# [write] mean and IQ range of time series of CF, without reservoir taken into account
for n in range(len(index)):
    
    # [loop] over all plants where annual mean discharge is known
    if bias_correction_yearly[n] > 0:
        
        # [calculate] CF for all plants where annual mean discharge and design discharge are known
        if design_discharge[n] > 0:
            utilisation_rate_normal[n,:] = Q_bymonth_normal_reservoir[n,:] / design_discharge[n]
            utilisation_rate_dry[n,:] = Q_bymonth_dry_reservoir[n,:] / design_discharge[n]
            utilisation_rate_wet[n,:] = Q_bymonth_wet_reservoir[n,:] / design_discharge[n]
            
        # [calculate] CF for all plants where annual mean discharge is known, but design discharge is unknown
        else:
            utilisation_rate_normal[n,:] = Q_bymonth_normal_reservoir[n,:] / bias_correction_yearly[n] * availability[n]
            utilisation_rate_dry[n,:] = Q_bymonth_dry_reservoir[n,:] / bias_correction_yearly[n] * availability[n]
            utilisation_rate_wet[n,:] = Q_bymonth_wet_reservoir[n,:] / bias_correction_yearly[n] * availability[n]
    
    # [calculate] CF for all plants where annual mean discharge is unknown
    else:
        utilisation_rate_normal[n,:] = seasonality_normal[n,:] / (f_ratio*np.max(seasonality_normal[n,:]))
        utilisation_rate_dry[n,:] = seasonality_dry[n,:] / (f_ratio*np.max(seasonality_normal[n,:]))
        utilisation_rate_wet[n,:] = seasonality_wet[n,:] / (f_ratio*np.max(seasonality_normal[n,:]))


# [calculate] output time series (CF of hydropower plant)
for n in range(len(index)):
    for m in range(months_yr):
            CF_bymonth_normal[n,m] = min(utilisation_rate_normal[n,m], 1) #*capacity[n]
            CF_bymonth_dry[n,m] = min(utilisation_rate_dry[n,m], 1) #*capacity[n]
            CF_bymonth_wet[n,m] = min(utilisation_rate_wet[n,m], 1) #*capacity[n]
            # [correction] in case the statistical calculations led to "upper" or "lower" being below or above "mean" in certain months
            if CF_bymonth_dry[n,m] > CF_bymonth_normal[n,m]:
                CF_bymonth_dry[n,m] = CF_bymonth_normal[n,m]
            if CF_bymonth_wet[n,m] < CF_bymonth_normal[n,m]:
                CF_bymonth_wet[n,m] = CF_bymonth_normal[n,m]

print('outputs ready')


# %% 6) Write data for visualisation

# [export to Excel] flow time series of bias-corrected inflows
writer : ExcelWriter = pd.ExcelWriter('Results.xlsx')

df_names = pd.DataFrame(np.array(index_name)[output_binary])
df_names.to_excel(writer, 'Results', index = False, startcol = 0, header = ["Name"])

df_countries = pd.DataFrame(np.array(index_country)[output_binary])
df_countries.to_excel(writer, 'Results', index = False, startcol = 1, header = ["Country"])

df_output_normal = pd.DataFrame(CF_bymonth_normal[output_binary,:])
df_output_normal.to_excel(writer, 'Results', index = False, startcol = 2, header = [i+j+k for i,j,k in zip(months_names_short,header_output,header_normal)])

df_output_dry = pd.DataFrame(CF_bymonth_dry[output_binary,:])
df_output_dry.to_excel(writer, 'Results', index = False, startcol = 14, header = [i+j+k for i,j,k in zip(months_names_short,header_output,header_dry)])

df_output_wet = pd.DataFrame(CF_bymonth_wet[output_binary,:])
df_output_wet.to_excel(writer, 'Results', index = False, startcol = 26, header = [i+j+k for i,j,k in zip(months_names_short,header_output,header_wet)])

writer.save()

print('Results for Excel ready (I)')

