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
## This is what doesn't wory from github...
#iran_data = gpd.read_file(os.path.join(url, 'Iran_shp', 'iran_admin.shp'))
iran_data = gpd.read_file(fp + '/Iran_shp/iran_admin.shp')

## Subset to relevant columns and update names
iran_data = iran_data[['ADM2_EN','ADM2_FA','ADM1_EN','ADM1_FA','Shape_Leng','Shape_Area','geometry']]
iran_data.columns = ['county_en', 'county_fa', 'province_en', 'province_fa', 'shape_len', 'shape_area', 'geometry']

## Read in human data
human_data = pd.read_csv(os.path.join(url, 'Data', 'Human_Brucellosis_2015-2018_V2.csv')).drop(['Unnamed: 18', 'Unnamed: 19'], axis = 1)

## Left join human data to spatial data
## Gives us a fair amount of missing matches...
human_sp_data = human_data.merge(iran_data, how = 'outer', left_on = 'County', right_on = 'county_en')

#%%

##########################################
## Identifying Spelling Inconsistencies ##

## Function to match potential misspelled strings
## Takes two pandas series, calculates Levenshtein distance to identify potential matches
def match_names(s1, s2):
    
    ## Testing levenshtein distance method
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

## Investigating our data with this method ##

## Subset to unique county names in each file
shp_counties = iran_data['county_en'].unique()
csv_counties = human_data['County'].unique()

## Unique values in shapefile that aren't in human_data
shp_ct_unique = pd.Series(np.setdiff1d(shp_counties, csv_counties))

## Unique values in csv file that aren't in shapefile
csv_ct_unique = pd.Series(np.setdiff1d(csv_counties, shp_counties))

## Create dataframe recording potential name matches
matched_df = match_names(shp_ct_unique, csv_ct_unique)

## Add column recording whether distance and ratio identify the same match
matched_df['name_match'] = matched_df['name_dist'] == matched_df['name_ratio']

## Can be highly confident when nameMatch = True, ratio > .75 - this matches 145 of the 192 that need matches
matched_df['matched'] = np.where((matched_df['name_match'] == True) & (matched_df['ratio'] > .75), matched_df['name_dist'], 'NULL')
       
## These names are still unmatched but their proposed names haven't yet been taken
## Visual inspection suggests the name_dist is the right match in all these cases - so we go with that
matched_df.loc[
    ((matched_df['name_dist'].isin(matched_df['matched']) == False) & ((matched_df['name_ratio'].isin(matched_df['matched']) == False))) & (matched_df['matched'] == 'NULL'), 
    'matched'] = matched_df['name_dist']

## We're left with these names, which are definitely wrong (21 names)
## Their proposed matches have already been taken - will have to do by hand
match_by_hand = matched_df.loc[matched_df['matched']=='NULL']


## Next:
## Match last 21 names manually somehow?
## Reduce matched_df to a series of name pairs
## in the human_data, map county names to their partners in this matching series just using a join
## Rejoin human_data with spatial data using this updated county name column
## Could also wrap some of the above matching code in another function



#%%

## Also - some places that did merge have different provinces? Double check this once names are updated
human_sp_data[['County', 'county_en', 'Province', 'province_en']].loc[human_sp_data['county_en'].isnull()]
