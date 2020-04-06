#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 15:02:45 2020

@author: finnroberts
"""

import os
import pandas as pd
import geopandas as gpd
import numpy as np
import Levenshtein as leven

#%%

#####################
## Reading in data ##

## Filepath for local files
fp = '/Users/finnroberts/Minnesota/2020_Spring/GEOG5541/Group/Project-Repo'

## Path to data on github
url = 'https://raw.githubusercontent.com/GEOCOMP-Brucellosis-Project/Project-Repo/master/'

## Read animal data and update columns
animal_data = pd.read_csv(os.path.join(url, 'Data', 'animal_vac_data.csv'))
animal_data.columns = [
        'id', 'unitCode', 'unitType', 'province',
        'county', 'livestock_type', 'time_j', 'time_g',
        'lat' , 'long', 'n_sample', 'n_checked', 'n_infected', 
        'n_rejected', 'n_suspicious'
                      ]

## Create new columns storing month and year of animal testing
animal_data[['month', 'year']] = animal_data['time_g'].str.split('/', expand = True)[[0,2]]

## Sum infection data grouped by county, year, and month
animal_data_grp = animal_data.groupby(['county', 'year', 'month'], as_index = False)[animal_data.columns[10:15]].sum().merge(animal_data)

## Calculate infection rate
animal_data_grp['animal_inf_rate'] = animal_data_grp['n_infected']/animal_data['n_sample']

## Read in Iran data
## This is what doesn't work from github...
#iran_data = gpd.read_file(os.path.join(url, 'Iran_shp', 'iran_admin.shp'))
iran_data = gpd.read_file(fp + '/Iran_shp/iran_admin.shp')

## Subset to relevant columns and update names
iran_data = iran_data[['ADM2_EN','ADM2_FA','ADM1_EN','ADM1_FA','Shape_Leng','Shape_Area','geometry']]
iran_data.columns = ['county_en', 'county_fa', 'province_en', 'province_fa', 'shape_len', 'shape_area', 'geometry']

## Read in human data
human_data = pd.read_csv(os.path.join(url, 'Data', 'Human_Brucellosis_2015-2018_V2.csv')).drop(['Unnamed: 18', 'Unnamed: 19'], axis = 1)

## Left join human data to spatial data
## Gives us a fair amount of missing matches...
# human_sp_data = human_data.merge(iran_data, how = 'outer', left_on = 'County', right_on = 'county_en')

#%%

##########################################
## Identifying Spelling Inconsistencies ##

## Function to match potential misspelled strings
## Takes two pandas series, calculates Levenshtein distance to identify potential matches
## Used in the likely_matches function
def match_names(s1, s2):
    
    ## Calculate Levenshtein distance and ratio
    dists = np.array([leven.distance(name1, name2) for name1 in s1 for name2 in s2])
    ratios = np.array([leven.ratio(name1, name2) for name1 in s1 for name2 in s2])
    
    ## Reshape and convert to df so we can identify which values are for which name combo
    dists_df = pd.DataFrame(data = dists.reshape(len(s1), len(s2)), index = s1, columns = s2)
    ratios_df = pd.DataFrame(data = ratios.reshape(len(s1), len(s2)), index = s1, columns = s2)
    
    ## Get column names where min distance and max ratio occurs
    matches = pd.DataFrame({'name_dist': dists_df.idxmin(axis = 1), 
                            'name_ratio': ratios_df.idxmax(axis = 1), 
                            'dist': dists_df.min(axis = 1), 
                            'ratio': ratios_df.max(axis = 1)})
    
    return(matches)

## Function to return highly probable string matches - the rest will have to be done manually
## Depends on match_names function
def likely_matches(s1, s2):
    
    ## Unique values in each series
    vals1 = s1.unique()
    vals2 = s2.unique()
    
    ## Unique values in series 1 that aren't in series 2
    unique1 = pd.Series(np.setdiff1d(vals1, vals2))
    
    ## Unique values in series 2 that aren't in series 1
    unique2 = pd.Series(np.setdiff1d(vals2, vals1))
    
    ## Create dataframe recording potential name matches
    matched = match_names(unique1, unique2)
    
    ## Add column recording whether distance and ratio identify the same match
    matched['name_match'] = (matched['name_dist'] == matched['name_ratio'])
    
    ## Can be highly confident when nameMatch = True, ratio > .75 - this matches 145 of the 192 that need matches
    matched['matched'] = np.where((matched['name_match'] == True) & (matched['ratio'] > .75), matched['name_dist'], 'NULL')
           
    return(matched)
    
#%%

######################################
## Cleaning our data's county names ##

matched_df = likely_matches(iran_data['county_en'], human_data['County'])

## These names are still unmatched but their proposed names haven't yet been taken
## Visual inspection suggests the name_dist is the right match in all these cases - so we go with that
matched_df.loc[
    ((matched_df['name_dist'].isin(matched_df['matched']) == False) & 
     ((matched_df['name_ratio'].isin(matched_df['matched']) == False))) & 
    (matched_df['matched'] == 'NULL'), 
    'matched'] = matched_df['name_dist']

## We're left with these names, which are definitely wrong (21 names)
## Their proposed matches have already been taken - will have to do by hand
match_by_hand = matched_df.loc[matched_df['matched']=='NULL']

## Will need to match manually here
##
##

## Prepare matched names for joining on human data
matched_pairs = matched_df['matched'].reset_index()
matched_pairs.columns = ['county_shp', 'matched']

## Join matched names on human data to get a column of county names compatible with the shapefile
human_data2 = pd.merge(matched_pairs, human_data, how = 'right', left_on = 'matched', right_on = 'County')

## Wherever county_shp and matched are null, just fill in both with the associated county name
## (These are the values that were the same in the csv and shp from the beginning)
human_data2.loc[(human_data2['county_shp'].isnull()), 'county_shp'] = human_data2['County']

## Join updated human data on shapefile
human_sp_data = human_data2.merge(iran_data, how = 'outer', left_on = 'county_shp', right_on = 'county_en')


## Next:
## Deal with manual matches
## Cross reference final 'matched' column with all the names in the original csv/shp files to see if there are as-yet unmmatched potential names
## Double check some earlier matches - they may not all be correct



#%%

## Also - some places that did merge have different provinces? Double check this once names are updated
# human_sp_data[['County', 'county_en', 'Province', 'province_en']].loc[human_sp_data['county_en'].isnull()]

