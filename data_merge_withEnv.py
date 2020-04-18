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

"""If we all have our folders set up the same (cloned from github) this should work for everyone."""
fp = os.getcwd() #'/Users/finnroberts/Minnesota/2020_Spring/GEOG5541/Group/Project-Repo'

## Path to data on github
#url = 'https://raw.githubusercontent.com/GEOCOMP-Brucellosis-Project/Project-Repo/master/'

## Read animal data and update columns
animal_data = pd.read_csv(os.path.join(fp, 'Data', 'animal_vac_data.csv'))
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

## This accidental escape sequence is problematic later, so deal with manually here
iran_data.loc[iran_data['county_en'] == 'Yasooj\r', 'county_en'] = 'Yasooj'

## Read in human data
human_data = pd.read_csv(os.path.join(fp, 'Data', 'Human_Brucellosis_2015-2018_V3.csv'))

human_data = human_data.rename(columns = {'Urban/Rural/Itinerant/Nomadic':'Pop_setting',
                                          'Prepnancy':'Pregnancy',
                                          'Occuptio':'Occupation',
                                          'Livestock interaction history':'Livestock_int_hist',
                                          'Livestock interaction type':'Livestock_int_type',
                                          'Unpasteurized dairy consumption ':'Unpast_dairy',
                                          'Other family members infection':'Fam_members_inf',
                                          'Outbreak Year':'Outbreak_yr',
                                          'Outbreak Month':'Outbreak_mth',
                                          'Diagnosis Year':'Diagnosis_yr',
                                          'Diagnosis Month':'Diagnosis_mth',
                                          'Livestock vaccination history':'Livestock_vac_hist'})


## Fix duplicate provinces
human_data.loc[human_data['Province'] == 'Khorasan jonobi', 'Province'] = 'Khorasan Jonobi'
human_data.loc[human_data['Province'] == 'Khorasan shomali', 'Province'] = 'Khorasan Shomali'

#%%

##########################################
## Identifying Spelling Inconsistencies ##

## Function to match potential misspelled strings
## Takes two pandas series, calculates Levenshtein distance to identify potential matches
## Used in the likely_matches function
def match_names(s1, s2, as_df = True, caps = True, unique = True):
    
    s1 = pd.Series(s1)
    s2 = pd.Series(s2)
    
    ## Unique values in each series
    vals1 = s1.unique()
    vals2 = s2.unique()
    
    ## If unique argument is set to true, only match names that don't already have perfect match
    if unique == True:
        
        ## Unique values in series 1 that aren't in series 2
        vals1 = pd.Series(np.setdiff1d(vals1, vals2))
    
        ## Unique values in series 2 that aren't in series 1
        vals2 = pd.Series(np.setdiff1d(vals2, vals1))
        
    if caps == True:
        
        ## Capitalize before matching
        vals1 = vals1.str.capitalize()
        vals2 = vals2.str.capitalize()
        
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
    else:
        
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
test = match_names(provs1, provs2, as_df = False, caps = False, unique = True)

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
matched_df = likely_matches(counties1, counties2)
#match_names(counties1, counties2, as_df = False)

#matched_df[matched_df['matched'] == 'NULL']
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
    'Qale ganj':'Ghaleye-Ganj',
    'Sarchahan':'Hajiabad',
    'Kamfirouz':'Marvdasht',
    'Zarghan': 'Shiraz',
    'Beyza':'Sepidan',
    'Dore Chagni':'Doureh',
    'Sepid Dasht':'Khorramabad',
    'Nour Abad':'Mamasani',
    'Aleshtar':'Selseleh',
    'Boyerahmad':'Yasooj',
    'Saduq':'Yazd',
    'Mashhad Morghab':'Khorrambid',
    'Dehdez':'Izeh',
    'zaboli':'Mehrestan',
    'Qaemiyeh':'Kazerun',
    'Samen ol Aemmeh':'Mashhad',
    'kish':'Bandar-Lengeh'
              }

match_dict_cty.update(match_dict_man) ## This now has all automatched names and manually matched names

## Add all the perfect matches to this dictionary
## Mapping dictionary now should include all matches
perf_matches = np.intersect1d(iran_data['county_en'].str.capitalize(), human_data['County'].str.capitalize())

match_dict_cty.update(dict(zip(perf_matches, perf_matches)))

## Map names in dataframe based on dictionary
human_data['County'] = human_data['County'].map(match_dict_cty).fillna(human_data['County'])

## Joining ##
human_sp_data = pd.merge(human_data, iran_data, how = 'outer', left_on = 'County', right_on = 'county_en')

## Write mapping dictionary to csv for ease of QA
# pd.DataFrame.from_dict(data=match_dict_cty, orient='index').to_csv(fp + '/human_data_mappings.csv', index_label = ['human_county'], header=['shp_county'])

## QUALITY ASSURANCE NOTES ##

# 'Behbahan' associated with 2 provinces in the human data?

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
    'mahshahr':'Mahshahr',
    'Jafarieh':'Torbat-e-Jam',
    'Kahak':'Sabzevar',
    'Kohgiluyeh and BoyerAhmad':'Kohgeluyeh',
    'BoyerAhmad':'Yasooj'
              }

match_dict_ani.update(match_dict_man_ani)

#unmatched2 = [name for name in unmatched.index.map(ani_caps_mappings) if not name in match_dict_ani.keys()]

## Add identical matches to the dictionary
perf_matches = np.intersect1d(iran_data['county_en'], animal_data['county'])
match_dict_ani.update(dict(zip(perf_matches, perf_matches)))

## Update county names in animal data and do the join
animal_data['county'] = animal_data['county'].map(match_dict_ani).fillna(animal_data['county'])
ani_sp_data = pd.merge(animal_data, iran_data, how = 'outer', left_on = 'county', right_on = 'county_en')

## Write mapping dictionary to csv for ease of QA
# pd.DataFrame.from_dict(data=match_dict_ani, orient='index').to_csv(fp + '/animal_data_mappings.csv', index_label = ['animal_county'], header = ['shp_county'])

#%%

#######################
## SES Data Cleaning ##

## Read in SES data and rename columns
ses_data = pd.read_csv(os.path.join(fp, 'Data', 'ses_data_clean.csv'))[['province', 'pop', 'hshld_size', 'ses']]
#ses_data.columns = ['province', 'pop', 'hshld_size', 'ses']

## Get province names
ses_provs = ses_data['province'].unique()

## Automatch ses names to spatial data province names
ses_matches = likely_matches(ses_provs, iran_data['province_en'].unique())

## Maps from matched names back to uncapitalized names
ses_caps_map = map_caps(ses_provs)
iran_prov_map = map_caps(iran_data['province_en'].unique())

## Determine names that were successfully matched
automatched = ses_matches[ses_matches['matched'] != 'NULL']

## Create dictionary mapping ses names to spatial data names and update ses_data
match_dict_ses = dict(zip(automatched.index.map(ses_caps_map), automatched['matched'].map(iran_prov_map)))
ses_data['province'] = ses_data['province'].map(match_dict_ses).fillna(ses_data['province'])

## Join data
ses_sp_data = pd.merge(ses_data, iran_data, how = 'outer', left_on = 'province', right_on = 'province_en')

#%%

##########################
##  Env Data Merge      ##

#This may rely on one or two manual edits to the data depending on which source you use
#(some values are float type nan which causes issues)
#so feel free to comment this cell out if it's giving an error.

import jdatetime

def addGregorian(data, yearCol, moCol):
    """
    Adds gregorian month and year columns to a dataframe based on jalali month and year columns.
    Assumes first day of jalali month.
    """
    #Apply function iterates a function over each row in a dataframe (axis=1 specifies row)
    #The function creates a jalali date object for each row from the jalali year and month. 
    #Then, it converts it to a gregorian date and extracts the year or month. 
    #If the jalali year and date are not numeric (e.g. Null as a string) then the year and month get None values. 
    data['year']=data.apply(lambda row: jdatetime.date(int(row[yearCol]), int(row[moCol]), 1).togregorian().strftime('%y') if (row[yearCol].isnumeric() and row[moCol].isnumeric()) else None, axis=1)
    data['month']=data.apply(lambda row: jdatetime.date(int(row[yearCol]), int(row[moCol]), 1).togregorian().strftime('%m') if (row[yearCol].isnumeric() and row[moCol].isnumeric()) else None, axis=1)

def addEnvData(data, envDF, yearCol, moCol):
    """
    Gets environmental variables from a dataframe containing them for each row in a dataframe based on county and date.
    """
    #Apply function iterates a function over each row in a dataframe (axis=1 specifies row)
    #envDF.loc[row.County] sorts the environmental data to the correct county. 
    #The next bracket selects a column by creating a datestring in the correct format. .zfill is required to make sure months are in '04' format instead of '4'
    data['mean_2m_air_temperature'] = data.apply(lambda row: envDF.loc[row.County]['mean_2m_air_temperature_20'+row[yearCol]+row[moCol].zfill(2)] if (row[yearCol]!=None and row.County in envDF.index) else None, axis=1)
    data['mean_total_precipitation'] = data.apply(lambda row: envDF.loc[row.County]['total_precipitation_20'+row[yearCol]+row[moCol].zfill(2)] if (row[yearCol]!=None and row.County in envDF.index) else None, axis=1)
    data['mean_ndvi'] = data.apply(lambda row: envDF.loc[row.County]['mean_20'+row[yearCol]+row[moCol].zfill(2)] if (row[yearCol]!=None and row.County in envDF.index) else None, axis=1)
    data['mean_elevation'] = data.apply(lambda row: envDF.loc[row.County]['mean_elevation'] if (row.County in envDF.index) else None, axis=1)

## Read Environmental Data
envFP=os.path.join(fp, 'Data', 'allParams.csv')
envData = pd.read_csv(envFP, index_col='ADM2_EN')


## Clean up human data, add date column
human_sp_data.loc[pd.isna(human_sp_data['Outbreak_yr']), 'Outbreak_yr']='Null'
addGregorian(human_sp_data, 'Outbreak_yr', 'Outbreak_mth')

## Clean up animal data
ani_sp_data.rename(columns={"county": "County"}, inplace=True)

#addEnvData(human_sp_data, envData, 'year', 'month')
addEnvData(ani_sp_data, envData, 'year', 'month')

#%%
'''
## Write files
human_sp_data = gpd.GeoDataFrame(human_sp_data, crs = 'EPSG:4326', geometry = 'geometry')
human_sp_data.to_file(os.path.join(fp, 'human_shp', 'human_data_clean.shp'))

ani_sp_data = gpd.GeoDataFrame(ani_sp_data, crs = 'EPSG:4326', geometry = 'geometry')
ani_sp_data.to_file(os.path.join(fp, 'animal_shp', 'animal_data_clean.shp'))

ses_sp_data = gpd.GeoDataFrame(ses_sp_data, crs = 'EPSG:4326', geometry = 'geometry')
ses_sp_data.to_file(os.path.join(fp, 'ses_shp', 'ses_data_clean.shp'))
'''

