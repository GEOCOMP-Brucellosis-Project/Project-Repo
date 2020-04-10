# -*- coding: utf-8 -*-
"""
Bovine Brucellosis Project Code

To-do:
NDVI
Elevation/Terrain


#ideas!
Create buffer around each animal data point for more accurate env. cond
"""


###Imports
import ee
import pandas
ee.Initialize()


######### Functions #############
def imageParams(image, shapefile):
    """
    Takes an image and a shapefile, returns a dataframe with the mean band values summarized by shapefile feature. 
    """
    #Gets the average of the band values in each feature in the iran asset. See: https://developers.google.com/earth-engine/reducers_reduce_region
    envFC=image.reduceRegions(shapefile, ee.Reducer.mean(), image.projection().nominalScale(), image.projection())
    envDict=envFC.getInfo()['features']

    #Now, features is a messy list of dictionaries. We need to extract the "properties" from each feature. 
    cols=[] #blank variable to store properties
    for feature in envDict:
        cols.append(feature['properties']) #adds the properties list to cols
    
    return pandas.DataFrame(cols) #creates a new pandas dataframe with the data in cols
    


def getYearlyParams(collection, shapefile, startYear, endYear, first=True):
    """
    Gets the selected yearly parameters as a dataframe for an image collection.
    """
    yearlyParams=pandas.DataFrame(columns=['ADM2_EN'])
    for year in range(startYear, endYear+1):
        print(year)
        for month in range(1, 13):
            print(month, end=' ')
            if(first):
                image=collection.filter(ee.Filter.calendarRange(year,year,'year')).filter(ee.Filter.calendarRange(month, month,'month'))\
                .first()
            else:
                image=collection.filter(ee.Filter.calendarRange(year,year,'year')).filter(ee.Filter.calendarRange(month, month,'month'))\
                .mean()
                
            
            yearlyParams=yearlyParams.merge(imageParams(image, shapefile),\
                                            how='outer',\
                                            on='ADM2_EN',\
                                            suffixes=('', '_'+str(year)+str(month).zfill(2)))
    print("done")
    return yearlyParams




############# Main ######################
iran = ee.FeatureCollection("users/tannerjohnson56/iran_admin").select("ADM2_EN")  #Loads the iran boundaries asset with only county name field kept

dem=ee.Image("JAXA/ALOS/AW3D30/V2_2").select('AVE_DSM')
ndvi=ee.ImageCollection("NOAA/CDR/AVHRR/NDVI/V5").select('NDVI')
ndvi2=ee.ImageCollection("MODIS/006/MOD13A2").select('NDVI')
moWeather = ee.ImageCollection("ECMWF/ERA5/MONTHLY").select("mean_2m_air_temperature","total_precipitation") #Loads the ERA5 monthly collection


weatherParams=getYearlyParams(moWeather, iran, 2008, 2018)
ndviParams=getYearlyParams(ndvi, iran, 2008, 2018, False)
elevParams=imageParams(dem, iran)

allParams=weatherParams.merge(ndviParams, on='ADM2_EN').merge(elevParams, on='ADM2_EN')
allParams.to_csv(r'allParams.csv', index=False)