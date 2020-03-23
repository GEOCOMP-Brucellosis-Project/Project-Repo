# -*- coding: utf-8 -*-
"""
Bovine Brucellosis Project Code

To-do:
Write code to read-in animal and human data
Write code to calculate a brucellosis data metric that we want to use (can be just a placeholder for now! maybe infected/total bovine?)
Create dataframe with brucellosis data metric AND temperature+precipitation data
Create function to plot the data with Matplotlib!


Back-burner:
Create function to get mean monthly temperature
Create function to get total monthly precipitation



#ideas!
Create buffer around each animal data point for more accurate env. cond
"""


###Imports
import ee
import pandas
ee.Initialize()


###Main
dailyWeather = ee.ImageCollection("ECMWF/ERA5/DAILY") #Loads the ERA5 Daily collection
iran = ee.FeatureCollection("users/tannerjohnson56/iran_admin") #Loads the iran boundaries asset
iran = iran.select("ADM2_EN") #keeps only the county name field, which is administrative level 2
animal_path=r"C:\Users\tjons_000\Google Drive\#Geocomputing\project\BRUCELLOSIS 2008-18.xlsx"
animal_data=pandas.read_excel(animal_path) #loads the animal data as a pandas data frame



startDate = '1997-06-11'
endDate = '1997-06-12'

#dailyWeather is an image COLLECTION with many bands
#After the next line, we will have a single IMAGE with only two bands

#.filterDate() gets images in the collection that fall within the date range
#.first() gets the first image in an image collection and returns an image
#.select() returns the image with only the selected bands
day=dailyWeather.filterDate(startDate, endDate)\
.first()\
.select("mean_2m_air_temperature","total_precipitation" )




#Gets the average of the band values in each feature in the iran asset. See: https://developers.google.com/earth-engine/reducers_reduce_region
withEnv=day.reduceRegions(iran, ee.Reducer.mean()) 

#.getInfo() is used in ee to actually store the values locally since values are normally stored, manipulated, and mapped in the cloud. 
features = withEnv.getInfo()['features'] 

#Now, features is a messy list of dictionaries. We need to extract the "properties" from each feature. 
cols=[] #blank variable to store properties
for feature in features:
    cols.append(feature['properties']) #adds the properties list to cols

envParams=pandas.DataFrame(cols) #creates a new pandas dataframe with the data in cols


#one issue so far: unmatched county names
s1=set(animal_data.County)
s2=set(envParams.ADM2_EN)
matched = s1.intersection(s2) 
unmatched = s1.symmetric_difference(s2)
len(matched)
len(unmatched)
