# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:18:07 2020

@author: Tanner
"""
###Imports
import ee
import pandas
ee.Initialize()


###Main
dailyWeather = ee.ImageCollection("ECMWF/ERA5/DAILY")
iran = ee.FeatureCollection("users/tannerjohnson56/iran_admin")
iran = iran.select("ADM2_EN")
animal_path=r"C:\Users\tjons_000\Google Drive\#Geocomputing\project\BRUCELLOSIS 2008-18.xlsx"
animal_data=pandas.read_excel(animal_path)



startDate = '1997-06-11'
endDate = '1997-06-12'
day=dailyWeather.filterDate(startDate, endDate).first().select("mean_2m_air_temperature","total_precipitation" )

#Gets the average 
withEnv=day.reduceRegions(iran, ee.Reducer.mean())



features = withEnv.getInfo()['features']
cols=[]
for feature in features:
    cols.append(feature['properties'])
envParams=pandas.DataFrame(cols)

# about = withEnv.getInfo()
# features=about['features']
# cols=[]
# for feature in features:
#     cols.append(feature['properties'])
# envParams=pandas.DataFrame(cols)
