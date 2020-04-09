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

## Fix duplicate provinces
human_data.loc[human_data['Province'] == 'Khorasan jonobi', 'Province'] = 'Khorasan Jonobi'
human_data.loc[human_data['Province'] == 'Khorasan shomali', 'Province'] = 'Khorasan Shomali'

## Left join human data to spatial data
## Gives us a fair amount of missing matches...
# human_sp_data = human_data.merge(iran_data, how = 'outer', left_on = 'County', right_on = 'county_en')

#%%

##########################################
## Identifying Spelling Inconsistencies ##

## Function to match potential misspelled strings
## Takes two pandas series, calculates Levenshtein distance to identify potential matches
## Used in the likely_matches function
def match_names(s1, s2, as_df = True, caps = True, unique = True):
    
    s1 = pd.Series(s1)
    s2 = pd.Series(s2)
    
    if caps == True:
        
        s1 = s1.str.capitalize()
        s2 = s2.str.capitalize()
        
        
    
    ## Unique values in each series
    vals1 = s1.unique()
    vals2 = s2.unique()
    
    ## If unique argument is set to true, only match names that don't have already perfect match
    if unique == True:
        
        ## Unique values in series 1 that aren't in series 2
        vals1 = pd.Series(np.setdiff1d(vals1, vals2))
    
        ## Unique values in series 2 that aren't in series 1
        vals2 = pd.Series(np.setdiff1d(vals2, vals1))
    
    ## Calculate Levenshtein distance and ratio
    dists = np.array([leven.distance(name1, name2) for name1 in vals1 for name2 in vals2])
    ratios = np.array([leven.ratio(name1, name2) for name1 in vals1 for name2 in vals2])
        
    ## Reshape and convert to df so we can identify which values are for which name combo
    dists_df = pd.DataFrame(data = dists.reshape(len(vals1), len(vals2)), index = vals1, columns = vals2)
    ratios_df = pd.DataFrame(data = ratios.reshape(len(vals1), len(vals2)), index = vals1, columns = vals2)
    
    if as_df == True:
        
        ## Get column names where min distance and max ratio occurs
        matches = pd.DataFrame({'name_dist': dists_df.idxmin(axis = 1), 
                                'name_ratio': ratios_df.idxmax(axis = 1), 
                                'dist': dists_df.min(axis = 1), 
                                'ratio': ratios_df.max(axis = 1)})
    
        return(matches)
        
    ## If user wants more detail, we can create a dictionary for each name that contains more graunlar info
    ## This section creates a dictionary with names as keys and dataframes as values
    ## Each dataframe contains all the possible name pairings (as opposed to above, which only supplies best possible values)
    if as_df == False:
        
        ratios_dict = {name1: ratios_df.loc[name1].sort_values(ascending = False) for name1 in vals1}
        dists_dict = {name1: dists_df.loc[name1].sort_values() for name1 in vals1}
        
        comb_dict = {name1: pd.merge(dists_dict[name1], ratios_dict[name1], left_index = True, right_index = True, suffixes=('_dist', '_ratio')) for name1 in vals1}
        
        return(comb_dict)

## Function to return highly probable string matches - the rest will have to be done manually
## Depends on match_names function
## The index of the resulting df contains values from the *first* series passed to the function
def likely_matches(s1, s2, cutoff = 0.75, as_df = True, caps = True, unique = True):
    
    ## Create dataframe recording potential name matches
    matched = match_names(s1, s2, as_df = as_df, caps = caps, unique = unique)
    
    ## Add column recording whether distance and ratio identify the same match
    matched['name_match'] = (matched['name_dist'] == matched['name_ratio'])
    
    ## Can be highly confident when nameMatch = True, ratio >= .75 - this matches 145 of the 192 that need matches
    matched['matched'] = np.where((matched['name_match'] == True) & (matched['ratio'] >= cutoff), matched['name_dist'], 'NULL')
           
    return(matched)
    
#%%

###############################
## Cleaning our data's names ##

## Province Matching
    
provs1 = human_data.loc[human_data['Province'] != 'Null']['Province']
provs2 = iran_data['province_en']

likely_matches(provs1, provs2, caps = False)
match_names(provs1, provs2, as_df = False, caps = False, unique = True)

## Province matchings - this accounts for all discrepancies
match_dict_prov = {
    'West Azerbaijan':'Azarbaijan Gharbi',
    'East Azerbaijan':'Azarbaijan Sharghi',
    'Chaharmahal and Bakhtiari':'Chaharmahal & bakhtiari',
    'Isfahan':'Esfahan',
    'South Khorasan':'Khorasan Jonobi',
    'North Khorasan':'Khorasan Shomali',
    'Razavi Khorasan':'Khorasan Razavi',
    'Kohgiluyeh and Boyer-Ahmad':'Kohgiluyeh & Boyerahmad',
    'Kurdistan':'Kordestan',
    'Sistan and Baluchestan':'Sistan & Bluchestan'
              }

## Invert dictionary
match_dict_prov = {v: k for k, v in match_dict_prov.items()}

## Update province names in human data with dictionary mappings
human_data['Province'] = human_data['Province'].map(match_dict_prov).fillna(human_data['Province'])


#%%

## County Matching

## Remove null values for matching
counties1 = human_data.loc[human_data['County'] != 'Null']['County']
counties2 = iran_data['county_en']

## Create mapping of likely pairs
matched_df = likely_matches(counties1, counties2)
#match_names(counties1, counties2, as_df = False)

matched_df[matched_df['matched'] == 'NULL']
automatched = matched_df[matched_df['matched'] != 'NULL']

match_dict_cty = dict(zip(automatched.index, automatched['matched']))

## Manual matching
match_dict_man = {
    'Ali Abad Katul':'Aliabad', 
    'Bafgh':'Bafq', 
    'Bandar Qaz':'Bandar-e-gaz', 
    'Dailam':'Deylam',
    'Gonbad  kavoos':'Gonbad-e-kavus', 
    'Ijroud':'Eejrud', 
    'Jovein':'Jowayin',
    'Kalale':'Kolaleh',
    'Mahvalat':'Mahvelat', 
    'Menojan':'Manujan', ## auto-matched incorrectly
    'Neyshabur':'Nishapur',  
    'Orzoieyeh':'Arzuiyeh', 
    'Ray':'Rey', 
    'Tehran Jonub':'Tehran', 
    'Tehran Shomal':'Tehran', 
    'Tiran o Karvan':'Tiran-o-korun',
    'Abadeh Tashk':'Abadeh',
    'Agh Ghala':'Aqqala',
    'Ahvaz e gharb':'Ahvaz',
    'Ahvaz e Shargh':'Ahvaz',
    'Gilan Qarb':'Gilan-e-gharb',
    'Kharame':'Kherameh',
    'Maraqe':'Maragheh',
    'Mashhad Morghab':'Mashhad',
    'Tehran Gharb':'Tehran',
    'Tehran Shargh':'Tehran',
    'Tehran Shomal Qarb':'Tehran',
    'Zaveh':'Zave',
    'Bandar Mahshahr':'Mashahr',
    'Qale ganj':'Ghaleye-ganj'
              }

match_dict_cty.update(match_dict_man)




## Map names in dataframe based on dictionary
human_data['County'] = human_data['County'].map(match_dict).fillna(human_data['County'])


## Need to update human_data['County'] with the auto-matched names too before joining
## would be nice to have one big dictinoary of all the mappings for QA afterwords




test = pd.merge(human_data, iran_data, how = 'outer', left_on = 'County', right_on = 'county_en')

test = test[['County','Province','county_en','province_en']]


## QA -  check to see if any counties have multiple provinces after join - that would suggest different provinces in shp vs csv





## But first, should include province info in the matching function since counties with same province are obviously more likely to match
cty_prov_csv = human_data[['County', 'Province']].drop_duplicates('County')
cty_prov_shp = iran_data[['county_en', 'province_en']].drop_duplicates('county_en')

## Little function just to make finding matches by hand easier
## Since provinces are matched, it subsets the possible matches to province to make visual inspection easier
def prov_matcher(name):
    
    cty_prov_csv['County'] = cty_prov_csv['County'].str.capitalize()
    cty_prov_shp['county_en'] = cty_prov_shp['county_en'].str.capitalize()
    
    province = cty_prov_csv.loc[cty_prov_csv['County'] == name]['Province'].iloc[0]
    
    poss_matches = cty_prov_shp[cty_prov_shp['province_en'] == province]

    return(poss_matches)

prov_matcher('Zarghan')

## Dore chagni: Doureh? Dorud?

## Update county names here:











## Don't use this ugly method - delete once matching is working
#matched_df.loc[matched_df.index.isin(list(match_dict.keys())), 'matched'] = list(match_dict.values())

## So, we're ultimately left with these names...
matched_df.loc[(matched_df['matched'] == 'NULL')]

## Records for still unmatched names - helpful for manual matching
unmatched_names = matched_df[matched_df['matched'] == 'NULL'].index
unmatched_dict = {name: matched_dict[name] for name in unmatched_names}




#%%

## JOINING ##

## Prepare matched names for joining on human data
matched_pairs = matched_df['matched'].reset_index()
matched_pairs.columns = ['county_csv', 'matched']

## Capitalize so join works properly
iran_data['county_en'] = iran_data['county_en'].str.capitalize()

## Join matched names on human data to get a column of county names compatible with the shapefile
iran_data2 = pd.merge(matched_pairs, iran_data, how = 'outer', left_on = 'matched', right_on = 'county_en')

## This might be easier to just do in a dictionary:
## (could just assemble a dictionary of all names once matched and then coerce to df)

## Identify the names that matched precisely and fill them in manually 
## (these are deliberately not captured by the likely_names function to improve its performance)
precise_matches = np.intersect1d(iran_data['county_en'], human_data['County'])

iran_data2.loc[iran_data2['county_en'].isin(precise_matches) & 
               (iran_data2['matched'].isna()), 
               ['county_csv', 'matched']] = iran_data2['county_en']

## Join updated human data on shapefile
human_sp_data = iran_data2.merge(human_data, how = 'outer', left_on = 'county_csv', right_on = 'County')


#%%

#######################################
## Cleaning animal data county names ##

counties_ani = animal_data['county'].str.capitalize()

ani_shp_match = likely_matches(counties_ani, counties2)

ani_shp_match.loc[ani_shp_match['matched'] == 'NULL']


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

