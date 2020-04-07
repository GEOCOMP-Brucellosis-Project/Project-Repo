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
## The index of the resulting df contains values from the *first* series passed to the function
def likely_matches(s1, s2):
    
    s1 = pd.Series(s1)
    s2 = pd.Series(s2)
    
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
    
    ## Can be highly confident when nameMatch = True, ratio >= .75 - this matches 145 of the 192 that need matches
    matched['matched'] = np.where((matched['name_match'] == True) & (matched['ratio'] >= .75), matched['name_dist'], 'NULL')
           
    return(matched)
    
#%%

######################################
## Cleaning our data's county names ##

## Capitlize strings to improve matching
counties1 = human_data.loc[human_data['County'] != 'Null']['County'].str.capitalize()
counties2 = iran_data['county_en'].str.capitalize()

## Create mapping of likely pairs
matched_df = likely_matches(counties1, counties2)

## These names weren't auto-matched but their proposed matches seem right
## So we match them manually
names = [
    'Ali abad katul', 'Bafgh', 'Bandar qaz', 'Bile Savar', 
    'Dailam', 'Eslam Abad Gharb', 'Gonbad  kavoos', 
    'Ijroud', 'Jovein', 'Mahvalat', 'Menojan', 
    'Neyshabur',  'Orzoieyeh', 'Ray', 'Tehran jonub', 
    'Tehran shomal', 'Tiran o karvan'
    ]

matched_df.loc[matched_df.index.isin(names), 'matched'] = matched_df['name_dist']

## These names were not quite as obvious but seem to have the following matches
match_dict = {
    'Abadeh tashk':'Abadeh',
    'Agh ghala':'Aqqala',
    'Ahvaz e gharb':'Ahvaz',
    'Ahvaz e shargh':'Ahvaz',
    'Gilan qarb':'Gilan-e-gharb',
    'Kharame':'Kherameh',
    'Maraqe':'Maragheh',
    'Mashhad morghab':'Mashhad',
    'Tehran gharb':'Tehran',
    'Tehran shargh':'Tehran',
    'Tehran shomal qarb':'Tehran',
              }

matched_df.loc[matched_df.index.isin(list(match_dict.keys())), 'matched'] = list(match_dict.values())

## So, we're left with these 17 names...
matched_df.loc[(matched_df['matched'] == 'NULL')]

## Prepare matched names for joining on human data
matched_pairs = matched_df['matched'].reset_index()
matched_pairs.columns = ['county_csv', 'matched']

## Capitalize so join works properly
iran_data['county_en'] = iran_data['county_en'].str.capitalize()

## Join matched names on human data to get a column of county names compatible with the shapefile
iran_data2 = pd.merge(matched_pairs, iran_data, how = 'outer', left_on = 'matched', right_on = 'county_en')

## Identify the names that matched precisely and fill them in manually 
## (these are deliberately not captured by the likely_names function to improve its performance)
precise_matches = np.intersect1d(iran_data['county_en'], human_data['County'])
iran_data2.loc[iran_data2['county_en'].isin(precise_matches) & (iran_data2['matched'].isna()), ['county_csv', 'matched']] = iran_data2['county_en']

## Join updated human data on shapefile
human_sp_data = iran_data2.merge(human_data, how = 'outer', left_on = 'county_csv', right_on = 'County')




## TO DO:

## Manually match remaining 17 names
#### Can check provinces as a guide, too. May want to incorporate this into the likely matches function 
#### (i.e. two similar names with the same province are more likely to be matches - at the very least can include province
#### in likely_matches output to aid the user in finding matches manually)

## Deal with animal_data names.

#%%






## Next:

## Match animal names and join onto spatial master dataset (could also just join on spatial data - may not need everything to be in one)
## Merge SES data at province level on master dataset

## Deal with manual matches
## Cross reference final 'matched' column with all the names in the original csv/shp files to see if there are as-yet unmmatched potential names
## Double check some earlier matches - they may not all be correct

test = likely_matches(animal_data['county'], iran_data['county_en']) ## This seems to match way better - maybe function should match both ways to identify good matchs?

test.loc[test['matched'] == 'NULL']



#%%

## Also - some places that did merge have different provinces? Double check this once names are updated
# human_sp_data[['County', 'county_en', 'Province', 'province_en']].loc[human_sp_data['county_en'].isnull()]

