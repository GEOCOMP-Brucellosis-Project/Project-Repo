#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:04:02 2020

@author: finnroberts
"""

import pandas as pd

## Read data from GitHub repository
url = 'https://raw.githubusercontent.com/GEOCOMP-Brucellosis-Project/Project-Repo/master/Data/'
animal_data = pd.read_csv(url + 'animal_vac_data.csv')

## Update column names
animal_data.columns = [
        'id', 'unitCode', 'unitType', 'province',
        'county', 'livestock_type', 'time_j', 'time_g',
        'lat' , 'long', 'n_sample', 'n_checked', 'n_infected', 
        'n_rejected', 'n_suspicious'
                      ]

## Calculate infection rate
animal_data['animal_inf_rate'] = animal_data['n_infected']/animal_data['n_sample']

## Create new column storing year of animal testing
animal_data['year'] = '20' + animal_data['time_g'].str.slice(-2)

## Read in human data
human_data = pd.read_csv(url + 'Human_Brucellosis_09_11.csv')

## Merge data on ID key - I don't think this is right. Need to discuss proper merge key
merged_data = animal_data.merge(human_data, how='outer', left_on=['id'], right_on=['ID'])













#%% This section is just for random checks/calculations - not intended for inclusion in final code

## Some checks on the sparsity of our data
test = merged_data[merged_data['province'].isin(['Zanjan', 'Kordestan', 'Hamadan', 'Kermanshah'])]

## We have 12 observations with BOTH animal infection and human data
## We have 32 total observations with an animal infection
test[test['n_infected']>0].count()
