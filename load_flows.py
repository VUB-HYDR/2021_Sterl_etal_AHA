# -*- coding: utf-8 -*-
"""
Created on August 13 2020

@author: S. Sterl & A. Devillers
"""

import numpy as np
import pandas as pd
from pandas import ExcelWriter

# NOTE: The SWAT+ data should ONLY be used to infer seasonalities and interannual variabilities, not absolute values of river flow, unless bias-corrected.
# Bias-correction could be done using MONTHLY average flow (which may change the seasonality, but retain the SWAT+ interannual variability for each month)
# or using YEARLY average flow (retaining also the SWAT+ seasonality).

# [load] the relevant columns from the simulation .txt files
# (this takes a while; recommend to save it afterwards to avoid having to repeat this each time)

#channel_mon = np.genfromtxt("SWAT+_channel_mon_EWEMBI_hist.txt", skip_header = 3, usecols = (0,1,2,3,5,7), invalid_raise = False)
#np.save('SWAT+_channel_mon_EWEMBI_hist',channel_mon)

#reservoir_mon = np.genfromtxt("SWAT+_reservoir_mon_EWEMBI_hist.txt", skip_header = 3, usecols = (0,1,2,3,5,29), invalid_raise = False)
#np.save('SWAT+_reservoir_mon_EWEMBI_hist',reservoir_mon)


# %% 0.1) Load data


channel_mon = np.load('SWAT+_channel_mon_EWEMBI_hist.npy')
reservoir_mon = np.load('SWAT+_reservoir_mon_EWEMBI_hist.npy')
print('channels loaded')

# filename of Excel sheet with collected data
filename = 'African Hydropower Atlas v1-0.xlsx'

# [load] country of each hydro plant
index_country = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [0])
index_country = index_country['Country'].tolist()

# [load] names of each hydro plant
index_name = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [1])
index_name = index_name['Unit Name'].tolist()

# [load] channel indices for each hydro plant from prepared Excel file based on GIS data from https://zenodo.org/record/3580663
index = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [2])
index = index.to_numpy(dtype = 'int')

# [load] reservoir indices for each hydro plant from prepared Excel file based on GIS data from https://zenodo.org/record/3580663
index_res = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [3])
index_res = index_res.to_numpy(dtype = 'int')

# [load] design discharge (in m^3/s)
design_discharge = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [5])
design_discharge = design_discharge.to_numpy(dtype = 'float')
design_discharge[np.isnan(design_discharge)] = 0

# [load] yearly data for bias-correction (in m^3/s); if unavailable, use 50% of design discharge if latter available
bias_correction_yearly = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [4])
bias_correction_yearly = bias_correction_yearly.to_numpy(dtype = 'float')
bias_correction_yearly[np.isnan(bias_correction_yearly)] = 0
for n in range(len(bias_correction_yearly)):
    if (bias_correction_yearly[n] == 0) and (design_discharge[n] > 0):
        bias_correction_yearly[n] = (1/2)*design_discharge[n]

# [load] capacity (in MW)
capacity = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [6])
capacity = capacity.to_numpy(dtype = 'float')
capacity[np.isnan(capacity)] = 0

# [load] filling days (in days)
filling_days = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [7])
filling_days = filling_days.to_numpy(dtype = 'float')
filling_days[np.isnan(filling_days)] = 0

# [load] reservoir size (in million m^3)
reservoir_size = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [8])
reservoir_size = reservoir_size.to_numpy(dtype = 'float')
reservoir_size[np.isnan(reservoir_size)] = 0

# [load] average availability (%) [set to 50% in case of "Generic"]
availability = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [9])
availability = availability['Availability'].tolist()
for n in range(len(availability)):
    if availability[n] == 'Generic':
        availability[n] = 0.5
availability = np.array(availability)
availability[np.isnan(availability)] = 0

# [load] river of each hydro plant
index_river_name = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [10])
index_river_name = index_river_name['River Name'].tolist()

# [load] country of each hydro plant
index_basin = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [11])
index_basin = index_basin['River Basin'].tolist()

# [load] first year of hydropower plant
FirstYear = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [12])
FirstYear = FirstYear.to_numpy(dtype = 'float')
FirstYear[np.isnan(FirstYear)] = 0

# [load] variable SpillFrom for cascading analysis
SpillFrom = pd.read_excel (filename, sheet_name = '6 - Inputs code and GIS', usecols = [13])
SpillFrom = SpillFrom.to_numpy(dtype = 'str')

print('input loaded')



# %% 0.2) Extract data

# [extract] uncorrected inflow from SWAT+ simulation
flow_in = channel_mon[:,5]
flow_in_res = reservoir_mon[:,5]
channel_no = channel_mon[:,4]
reservoir_no = reservoir_mon[:,4]

# [define] time period covered by SWAT+ simulations
year_start = 1980
year_end = 2016

# [define] other time-related parameters
years_full = list(range(year_start, year_end + 1))
months_yr = 12
hrs_day = 24
months_names_short = np.array(["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
months_names_short = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
header_output = ["_capacity_factor"]*12
header_flow = ["_flow"]*12
header_normal = ["_normal"]*12
header_dry = ["_dry"]*12
header_wet = ["_wet"]*12
stat_type = ["Lower","Mean","Upper"]
profile = ["Flow","Output"]

# [set] number of days for each month of the year and corresponding hours
days_year = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
positions = np.zeros(shape = (len(days_year) + 1, 1))
for n in range(len(days_year)):
    positions[n+1] = hrs_day*days_year[n] + positions[n]


# [preallocate] the (uncorrected) flow time series and its seasonality (normalised) from the SWAT simulations
flow_extract = np.zeros(shape = (len(index), len(years_full)*months_yr))

# [extract] flow per channel ID
for n in range(len(index)):
    # [write] simulated flow values if flow to be extracted from "channel" dataset
    if (index[n] > 0) and (np.sum(channel_no == np.int(index[n])) > 0) :
        temp = flow_in[channel_no == np.int(index[n])]
        flow_extract[n,:] = temp[0:len(years_full)*months_yr]
    # [write] simulated flow values if flow to be extracted from "reservoir" dataset
    elif (index[n] > 0) and (np.sum(channel_no == np.int(index[n])) == 0) :
        temp = flow_in_res[reservoir_no == np.int(index_res[n])]
        flow_extract[n,:] = temp[0:len(years_full)*months_yr]
    # [write] nan values if simulated flow not available
    else:
        flow_extract[n,:] = np.nan

# [correct] take out negative values for flows based on "reservoirs" (= zero inflow times)        
flow_extract[flow_extract < 0] = 0

# [remove] spin-up years
years_spin = 8
years = years_full[years_spin:len(years_full)]
flow_extract = flow_extract[:,years_spin*months_yr:len(flow_extract)]


# %% 0.3) Bias-correct river flow data for each reservoir location


# [preallocate] the seasonality of the (uncorrected) flow time series
seasonality_extract = np.zeros(shape = (len(index), len(years)*months_yr))


# [preallocate] the corresponding normal and dry/wet range of the flow's seasonality in the simulations
seasonality_normal = np.full([len(index), months_yr], np.nan)
seasonality_dry = np.full([len(index), months_yr], np.nan)
seasonality_wet = np.full([len(index), months_yr], np.nan)


# [preallocate] the yearly average of the (uncorrected) flow time series
yearly_flow_average = np.full([len(index), len(years)],np.nan)
year_of_normal_flow = np.full([len(index)], np.nan)
year_of_dry_flow = np.full([len(index)], np.nan)
year_of_wet_flow = np.full([len(index)], np.nan)
corr_dry = np.ones(shape = len(index))*np.nan
corr_wet = np.ones(shape = len(index))*np.nan


# [preallocate] the corresponding normal and dry/wet bounds of bias-corrected flow
Q_bias_corrected_bymonth_normal_natural = np.full([len(index), months_yr], np.nan)
Q_bias_corrected_bymonth_dry_natural = np.full([len(index), months_yr], np.nan)
Q_bias_corrected_bymonth_wet_natural = np.full([len(index), months_yr], np.nan)


# find years corresponding to median and dry/wet range in terms of simulated flow volume
range_pct_dry = 5
range_pct_wet = 95
range_pct_dry_res = 10
range_pct_wet_res = 90
for n in range(len(index)):
    for y in range(len(years)):
        yearly_flow_average[n,y] = np.mean(flow_extract[n, y*months_yr:(y+1)*months_yr])
    if ~np.isnan(np.mean(yearly_flow_average[n,:])):
        year_of_normal_flow[n] = np.where(yearly_flow_average[n,:] == np.median(yearly_flow_average[n,:]))[0][0]
        # [check] if filling time less than one year, take normal pct range; if higher than one year, take reduced percentile range
        if filling_days[n] < 365:
            year_of_dry_flow[n] = np.where(np.abs(yearly_flow_average[n,:] - np.percentile(yearly_flow_average[n,:],range_pct_dry)) == np.min(np.abs(yearly_flow_average[n,:] - np.percentile(yearly_flow_average[n,:],range_pct_dry))))[0][0]
            year_of_wet_flow[n] = np.where(np.abs(yearly_flow_average[n,:] - np.percentile(yearly_flow_average[n,:],range_pct_wet)) == np.min(np.abs(yearly_flow_average[n,:] - np.percentile(yearly_flow_average[n,:],range_pct_wet))))[0][0]
        else:
            year_of_dry_flow[n] = np.where(np.abs(yearly_flow_average[n,:] - np.percentile(yearly_flow_average[n,:],range_pct_dry_res)) == np.min(np.abs(yearly_flow_average[n,:] - np.percentile(yearly_flow_average[n,:],range_pct_dry_res))))[0][0]
            year_of_wet_flow[n] = np.where(np.abs(yearly_flow_average[n,:] - np.percentile(yearly_flow_average[n,:],range_pct_wet_res)) == np.min(np.abs(yearly_flow_average[n,:] - np.percentile(yearly_flow_average[n,:],range_pct_wet_res))))[0][0]
        
           
# remove extreme years from series (beyond selected dry/wet range)
for n in range(len(index)):
    for y in range(len(years)):
        if ~np.isnan(np.mean(yearly_flow_average[n,:])):
            if yearly_flow_average[n,y] < yearly_flow_average[n,np.int(year_of_dry_flow[n])]:
                flow_extract[n, y*months_yr:(y+1)*months_yr] = np.nan
            if yearly_flow_average[n,y] > yearly_flow_average[n,np.int(year_of_wet_flow[n])]:
                flow_extract[n, y*months_yr:(y+1)*months_yr] = np.nan
                

# [write] normal and dry/wet range of monthly time series of normalised inflow
for n in range(len(index)):
    if ~np.isnan(np.mean(yearly_flow_average[n,:])):
        for m in range(months_yr):
            seasonality_normal[n,m] = np.nanmedian(flow_extract[n, np.arange(m, len(years)*months_yr, months_yr)]) / np.nanmean(flow_extract[n,:])
        # [introduce] corrective factors for dry/wet years
        corr_dry[n] = np.nanmean(flow_extract[n, np.arange(np.int(year_of_dry_flow[n]*months_yr), np.int((year_of_dry_flow[n] + 1)*months_yr), 1)]) / np.nanmean(flow_extract[n,:])
        corr_wet[n] = np.nanmean(flow_extract[n, np.arange(np.int(year_of_wet_flow[n]*months_yr), np.int((year_of_wet_flow[n] + 1)*months_yr), 1)]) / np.nanmean(flow_extract[n,:])


# [cap] the dry and wet ranges in case SWAT+ simulation gives unreasonably extreme results (dry year: less than ca. 10% of mean annual flow; wet year: more than ca. 4x mean annual flow)
for n in range(len(index)):
    if corr_dry[n] < np.nanpercentile(corr_dry, 5):
        corr_dry[n] = np.nanpercentile(corr_dry, 5)
    if corr_wet[n] > np.nanpercentile(corr_wet, 95):
        corr_wet[n] = np.nanpercentile(corr_wet, 95)


# [calculate] dry/wet year seasonality
for n in range(len(index)):
    seasonality_dry[n,:] = seasonality_normal[n,:]*corr_dry[n]
    seasonality_wet[n,:] = seasonality_normal[n,:]*corr_wet[n]


# [calculate] flow time series bias-corrected (in m^3/s)
for n in range(len(index)):
    if bias_correction_yearly[n] > 0:
        Q_bias_corrected_bymonth_normal_natural[n,:] = bias_correction_yearly[n]*seasonality_normal[n,:]
        Q_bias_corrected_bymonth_dry_natural[n,:] = bias_correction_yearly[n]*seasonality_dry[n,:]
        Q_bias_corrected_bymonth_wet_natural[n,:] = bias_correction_yearly[n]*seasonality_wet[n,:]


print('flows ready (initial calculation based on natural flow)')
