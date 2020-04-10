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
    
## Little helper function that creates a dictionary mapping capitalized values
## to original values (so we can go back and forth more easily)
def map_caps(s1):
    
    s1 = pd.Series(s1)
    
    ## Capitalize series values
    s1_caps = s1.str.capitalize().unique()
    
    ## Map original values to capitalized values
    caps_mappings = {name1:name2 for name1, name2 in zip(s1_caps, s1.unique())}
    
    return(caps_mappings)
        
#%%

#########################
## Human Data Cleaning ##

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

## County Matching

## Remove null values for matching
counties1 = human_data.loc[human_data['County'] != 'Null']['County']
counties2 = iran_data['county_en']

## Create mapping of likely pairs
matched_df = likely_matches(counties1, counties2, 0.751)
#match_names(counties1, counties2, as_df = False)

matched_df[matched_df['matched'] == 'NULL']
automatched = matched_df[matched_df['matched'] != 'NULL']

## Revert automatched names back to original form data files and zip matched pairs into dictionary
caps_mappings1 = map_caps(counties1)
caps_mappings2 = map_caps(counties2)

match_dict_cty = dict(zip(automatched.index.map(caps_mappings1), automatched['matched'].map(caps_mappings2)))

## Manual matching
match_dict_man = {
    'Ali Abad Katul':'Aliabad', 
    'Bafgh':'Bafq', 
    'Bandar Qaz':'Bandar-e-Gaz', 
    'Dailam':'Deylam',
    'Gonbad  kavoos':'Gonbad-e-Kavus', 
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
    'Tiran o Karvan':'Tiran-o-Korun',
    'Abadeh Tashk':'Abadeh',
    'Agh Ghala':'Aqqala',
    'Ahvaz e gharb':'Ahvaz',
    'Ahvaz e Shargh':'Ahvaz',
    'Gilan Qarb':'Gilan-e-Gharb',
    'Kharame':'Kherameh',
    'Maraqe':'Maragheh',
    'Tehran Gharb':'Tehran',
    'Tehran Shargh':'Tehran',
    'Tehran Shomal Qarb':'Tehran',
    'Zaveh':'Zave',
    'Bandar Mahshahr':'Mahshahr',
    'Qale ganj':'Ghaleye-Ganj'
              }

match_dict_cty.update(match_dict_man) ## This now has all automatched names and manually matched names

## Add all the perfect matches to this dictionary
## Mapping dictionary now should include all matches
perf_matches = np.intersect1d(iran_data['county_en'], human_data['County'])
match_dict_cty.update(dict(zip(perf_matches, perf_matches)))

## Remove incorrect match for zaboli
del match_dict_cty['zaboli']

## Map names in dataframe based on dictionary
human_data['County'] = human_data['County'].map(match_dict_cty).fillna(human_data['County'])

## Joining ##
human_sp_data = pd.merge(human_data, iran_data, how = 'outer', left_on = 'County', right_on = 'county_en')

## Write mapping dictionary to csv for ease of QA
# pd.DataFrame.from_dict(data=match_dict_cty, orient='index').to_csv(fp + '/human_data_mappings.csv', index_label = ['human_county'], header=['shp_county'])

## QUALITY ASSURANCE NOTES ##

# 'Behbahan' associated with 2 provinces in the human data?
# 'Mashhad Morghab' doesn't go with Mashhhad because provinces don't match. Not sure where it goes

#%%

##########################
## Animal Data Cleaning ##

#provs1 = animal_data['province']
#provs2 = iran_data['province_en']

#likely_matches(provs1, provs2, caps = False)
#match_names(provs1, provs2, as_df = True, caps = False, unique = True)

## Province matchings - this accounts for all discrepancies
match_dict_prov = {
    'West Azerbaijan':'West Azarbayjan',
    'East Azerbaijan':'East Azarbayjan',
    'Chaharmahal and Bakhtiari':'Chaharmahal & bakhtiari',
    'Isfahan':'Esfahan',
    'South Khorasan':'Khorasan Jonobi',
    'North Khorasan':'Khorasan Shomali',
    'Razavi Khorasan':'Khorasan Razavi',
    'Sistan and Baluchestan':'Sistan & Bluchestan',
    'Hamadan':'Hamedan',
    'Kermanshan':'Kermanshah',
    'Kohgiluyeh and Boyer-Ahmad':'Kohgiluyeh and BoyerAhmad',
    'Kurdistan':'Kordestan',
    'Kerman':'South Kerman'
              }

## Invert dictionary
match_dict_prov = {v: k for k, v in match_dict_prov.items()}

## Update province names in human data with dictionary mappings
animal_data['province'] = animal_data['province'].map(match_dict_prov).fillna(animal_data['province'])


#######################################
## Cleaning animal data county names ##

ani_cnties = animal_data['county']

matched_df = likely_matches(ani_cnties, iran_data['county_en'])
#match_names(ani_cnties, counties2, as_df = False)

ani_caps_mappings = map_caps(ani_cnties)

automatched = matched_df[matched_df['matched'] != 'NULL']
#unmatched = matched_df[matched_df['matched'] == 'NULL']

## Start mapping dictionary with the automatched names
match_dict_ani = dict(zip(automatched.index.map(ani_caps_mappings), automatched['matched'].map(caps_mappings2)))

## Manually updated name mappings
match_dict_man_ani = {
    'Aran and Bidgol':'Aran-o-Bidgol',
    'Buin and Miandasht':'Booeino Miyandasht',
    'Deyr':'Dayyer',
    'Haftkel':'Haftgol',
    'Ijrud':'Eejrud',
    'Maneh asd Samalgan':'Maneh-o-Samalqan',
    'Orzueeyeh':'Arzuiyeh',
    'Qaleh Ganj':'Ghaleye-Ganj',
    'Qir and Karzin':'Qir-o-Karzin',
    'Raz and Jargalan':'Razo Jalgelan',
    'Sib and Suran':'Sibo Soran',
    'Tiran and Karvan':'Tiran-o-Korun',
    'Torqebeh and Shandiz(Binalud)':'Torghabe-o-Shandiz',
    'Zaveh':'Zave',
    'Chardavol':'Shirvan-o-Chardavol',
    'Torkaman':'Bandar-e-Torkaman',
    'mahshahr':'Mahshahr'
              }

match_dict_ani.update(match_dict_man_ani)

#unmatched2 = [name for name in unmatched.index.map(ani_caps_mappings) if not name in match_dict_ani.keys()]

## Reamining to be matched:
## Kohgiluyeh and BoyerAhmad: Kohgeluyeh??
## BoyerAhmad also to Kohgeluyeh then?
## Jafarieh and kahak are in Qom province which only has one entry so we could join on this?

## Add identical matches to the dictionary
perf_matches = np.intersect1d(iran_data['county_en'], animal_data['county'])
match_dict_ani.update(dict(zip(perf_matches, perf_matches)))

## Update county names in animal data and do the join
animal_data['county'] = animal_data['county'].map(match_dict_ani).fillna(animal_data['county'])
ani_sp_data = pd.merge(animal_data, iran_data, how = 'outer', left_on = 'county', right_on = 'county_en')

## Write mapping dictionary to csv for ease of QA
# pd.DataFrame.from_dict(data=match_dict_ani, orient='index').to_csv(fp + '/animal_data_mappings.csv', index_label = ['animal_county'], header = ['shp_county'])


#%%

## Also - some places that did merge have different provinces? Double check this once names are updated
# human_sp_data[['County', 'county_en', 'Province', 'province_en']].loc[human_sp_data['county_en'].isnull()]
