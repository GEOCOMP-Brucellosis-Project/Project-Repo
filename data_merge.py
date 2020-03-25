#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:04:02 2020

@author: finnroberts
"""

import pandas as pd
import geopandas as gpd
import fiona
import folium
import json
import pyproj

## Read data from GitHub repository
url = 'https://raw.githubusercontent.com/GEOCOMP-Brucellosis-Project/Project-Repo/master/'

## Filepath for local files
fp = '/Users/finnroberts/Minnesota/2020_Spring/GEOG5541/Group/Project-Repo'

animal_data = pd.read_csv(url + '/Data/animal_vac_data.csv')

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
human_data = pd.read_csv(url + '/Data/Human_Brucellosis_09_11.csv')

#%% Basic Visualization

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
m.save(fp + '/Iran_test.html')

fh.close()

#%% Merging

## NEXT:
## Merge human data on county shapefile by centroid location to get county names
## Aggregate animal data to county level and monthly period
## Join human and animal data at county level

## For some reason cannot get this to work using .shp straight from github
iran_data = gpd.read_file(fp + '/Iran_shp/irn_admbnda_adm1_unhcr_20190514.shp')

## Subset columns
iran_data = iran_data[['ADM1_EN','Shape_Leng','Shape_Area','geometry']]

## Get centroids of each feature in shapefile
iran_data['centroid'] = iran_data.centroid

## Define projection for UTM 39N/WGS84 coordinates
utm39Proj = pyproj.Proj("+proj=utm +zone=39N, +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## Project UTM 39N centroids to lat/lon
human_data['centroid_lat'] = utm39Proj(human_data['X_Centroid'].values, human_data['Y_Centroid'].values, inverse = True)[0]
human_data['centroid_lon'] = utm39Proj(human_data['X_Centroid'].values, human_data['Y_Centroid'].values, inverse = True)[1]




