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
import calendar
import ee
import pandas
ee.Initialize()


######### Functions #############
def imageParams(image, shapefile):
    """
    Takes an image and a shapefile, returns a dataframe with the parameters for each county for each month. 
    """
    #Gets the average of the band values in each feature in the iran asset. See: https://developers.google.com/earth-engine/reducers_reduce_region
    envFC=image.reduceRegions(shapefile, ee.Reducer.mean())
    envDict=envFC.getInfo()['features']

    #Now, features is a messy list of dictionaries. We need to extract the "properties" from each feature. 
    cols=[] #blank variable to store properties
    for feature in envDict:
        cols.append(feature['properties']) #adds the properties list to cols
    
    return pandas.DataFrame(cols) #creates a new pandas dataframe with the data in cols
    


def getYearlyParams(collection, shapefile, startYear, endYear):
    """
    Gets the selected yearly parameters as a dataframe for an image collection.
    """
    yearlyParams=pandas.DataFrame(columns=['ADM2_EN'])
    for year in range(startYear, endYear+1):
        print(year)
        for month in range(1, 13):
            print(month)
            image=collection.filter(ee.Filter.calendarRange(year,year,'year'))\
                .filter(ee.Filter.calendarRange(month, month,'month'))\
                .first()
            
            yearlyParams=yearlyParams.merge(imageParams(image, shapefile),\
                                            how='outer',\
                                            on='ADM2_EN',\
                                            suffixes=('', '_'+str(year)+str(month).zfill(2)))
    print("done")
    return yearlyParams




############# Main ######################

moWeather = ee.ImageCollection("ECMWF/ERA5/MONTHLY").select("mean_2m_air_temperature","total_precipitation") #Loads the ERA5 monthly collection
iran = ee.FeatureCollection("users/tannerjohnson56/iran_admin").select("ADM2_EN") #Loads the iran boundaries asset with only county name field kept


yearlyParams=getYearlyParams(moWeather, iran, 2008, 2018)
yearlyParams.to_csv(r'yearlyParams.csv',index = False)















#################OLD CODE


#old fashioned way
# image1=moWeather.filter(ee.Filter.calendarRange(2010, 2010, 'year'))\
#     .filter(ee.Filter.calendarRange(1, 1, 'month')).first()

# image2=moWeather.filter(ee.Filter.calendarRange(2010, 2010, 'year'))\
#     .filter(ee.Filter.calendarRange(2, 2, 'month')).first()
    

# db0=pandas.DataFrame(columns=['ADM2_EN'])
# db1=imageParams(image1, iran)
# db2=imageParams(image2, iran)

# db3=db0.merge(db1, how='outer', on='ADM2_EN', suffixes=('', '_201001'))
# db3=db1.merge(db2, how='outer', on='ADM2_EN', suffixes=('', '_201002'))






"""withEnv=day.reduceRegions(iran, ee.Reducer.mean()) 

#.getInfo() is used in ee to actually store the values locally since values are normally stored, manipulated, and mapped in the cloud. 
features = withEnv.getInfo()['features'] 

#Now, features is a messy list of dictionaries. We need to extract the "properties" from each feature. 
cols=[] #blank variable to store properties
for feature in features:
    cols.append(feature['properties']) #adds the properties list to cols

envParams=pandas.DataFrame(cols) #creates a new pandas dataframe with the data in cols
"""


#envParams.to_csv(r'envParams.csv',index = False)



# #one issue so far: unmatched county names
# s1=set(animal_data.County)
# s2=set(envParams.ADM2_EN)
# matched = s1.intersection(s2) 
# unmatched = s1.symmetric_difference(s2)
# len(matched)
# len(unmatched)
