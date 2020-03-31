#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:04:02 2020

@author: finnroberts
"""

import os
import sys
import pandas as pd
import geopandas as gpd
import fiona
import folium
import json
import pyproj
import shapely
import rtree

## Filepath for local files
fp = '/Users/finnroberts/Minnesota/2020_Spring/GEOG5541/Group/Project-Repo'

## Read animal data and update columns
animal_data = pd.read_csv(fp + '/Data/animal_vac_data.csv')
animal_data.columns = [
        'id', 'unitCode', 'unitType', 'province',
        'county', 'livestock_type', 'time_j', 'time_g',
        'lat' , 'long', 'n_sample', 'n_checked', 'n_infected', 
        'n_rejected', 'n_suspicious'
                      ]

## Create new columns storing month and year of animal testing
animal_data[['month', 'year']] = animal_data['time_g'].str.split('/', expand = True)[[0,2]]

## Read in human data
human_data = pd.read_csv(fp + '/Data/Human_Brucellosis_09_11.csv')

## This appears to be at province level but the github folder says county level?
iran_data = gpd.read_file(fp + '/Iran_shp/irn_admbnda_adm1_unhcr_20190514.shp')

## Subset columns
iran_data = iran_data[['ADM1_EN','Shape_Leng','Shape_Area','geometry']]

## Define projection for UTM 39N/WGS84 coordinates
utm39Proj = pyproj.Proj("+proj=utm +zone=39N +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## Project UTM 39N centroids to lat/lon
human_data['centroid_lon'] = utm39Proj(human_data['X_Centroid'].values, 
                                       human_data['Y_Centroid'].values, inverse = True)[1]

human_data['centroid_lat'] = utm39Proj(human_data['X_Centroid'].values, 
                                       human_data['Y_Centroid'].values, inverse = True)[0]

## Convert human data to geodataframe
human_data_geo = gpd.GeoDataFrame(human_data, 
                                  geometry = gpd.points_from_xy(human_data['centroid_lat'], human_data['centroid_lon']))

human_data_geo.crs = 'EPSG:4326'

## Spatial join
joined_data = gpd.sjoin(iran_data, human_data_geo)

## Aggregate animal data by county and year
animal_data_grp = animal_data.groupby(['county', 'year', 'month'], as_index = False)[animal_data.columns[10:14]].sum().merge(animal_data)

## Calculate infection rate
animal_data_grp['animal_inf_rate'] = animal_data_grp['n_infected']/animal_data['n_sample']

## Join human and animal data at county level
animal_data_grp.merge(joined_data, left_on = 'county', right_on = 'ADM1_EN')


#### Should spend some more time confirming all these joins are working as expected.
#### Need to update iran shapefile with something that has county names. Just provinces so far
#### Do some QA on animal data -- some instances have negative values for n_rejected...(?)













#%% Basic Visualization

'''
## This section creates a folium map of iran with animal barn sites
## as points. Just to get a sense of what we're dealing with.

## Read in shapefile
fh = fiona.open(fp + '/Iran_shp/irn_admbnda_adm1_unhcr_20190514.shp', 'r')

## Create base map of Iran
m = folium.Map(location=[33, 53], zoom_start=6)

## Add province(?) polygons to map
for f in fh:
    
    j = json.dumps(f)
    folium.features.GeoJson(j).add_to(m)

## Get animal barn site coordinates
points = animal_data[['lat','long']]

## Iterate through coordinates and plot on map. (Takes a bit - lots of points)
for point in range(len(points)):
    
    lat = points.iloc[point]['lat']
    long = points.iloc[point]['long']
    
    folium.CircleMarker(location=[lat, long],
                        radius=1,
                        color='red').add_to(m)

## Save html output
# m.save(fp + '/Iran_test.html')

fh.close()
'''

