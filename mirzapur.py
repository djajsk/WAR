import ee, datetime
import pandas as pd
import numpy as np
import datetime as dt

ee.Authenticate()
ee.Initialize()

geometry = ee.Geometry.Polygon([[[82.62738797515874,25.231463730353102],
[82.62459847778325,25.23068730938747],
[82.62142274230962,25.229600311706324],
[82.61979195922856,25.229017987522806],
[82.61914822906499,25.22847448243552],
[82.61631581634526,25.227736864502994],
[82.61562917083745,25.22614514740447],
[82.61717412323003,25.22583456602935],
[82.6192769750977,25.225329869603126],
[82.62228104919438,25.22470870189572],
[82.62567136138921,25.224126354288742],
[82.62858960479741,25.22432047046754],
[82.63133618682866,25.2240875310158],
[82.63215157836919,25.223078121569724],
[82.63386819213872,25.222029879816937],
[82.63627145141606,25.219234524313098],
[82.63755891174321,25.217953298240086],
[82.6409063085938,25.217759171900433],
[82.6434812292481,25.217953298240086],
[82.64652821868901,25.217759171900433],
[82.64978978485112,25.21810859908878],
[82.65275094360356,25.218419200191455],
[82.65433881134038,25.21865215049807],
[82.65815827697759,25.219933369211127],
[82.65721413940435,25.220088667532146],
[82.65781495422368,25.222573413705287],
[82.65781495422368,25.223660474195594],
[82.6564845785523,25.224863994119975],
[82.65549752563481,25.22529104671435],
[82.64794442504888,25.228047441019903],
[82.64468285888677,25.228823878841926],
[82.64120671600347,25.229522668643],
[82.63987634033208,25.229949704877846],
[82.6379880651856,25.23018263310113],
[82.6368722662354,25.230609667018218],
[82.63588521331792,25.230764951707126],
[82.63506982177739,25.231230804583987],
[82.63352486938481,25.23142490942258],
[82.63099286407476,25.23220132567812],
[82.6284608587647,25.23255071137552],
[82.62738797515874,25.231463730353102]]])

def get_turbidity(start_date, end_date):

# Selecting the satellite and AOI  
# Sentinel 2A
# copernicus/s2_sr 
    sentinel = ee.ImageCollection("COPERNICUS/S2_SR").\
               filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).\
               filterDate(start_date, end_date)
    AOI = geometry

    sentinel_AOI = sentinel.filterBounds(AOI)

#calculate NDTI
    def calculate_NDTI(image):
        ndti = image.normalizedDifference(['B4', 'B3']).rename('NDTI')
        return image.addBands(ndti)
    ndti = sentinel_AOI.map(calculate_NDTI)

# Mean NDTI
    def calculate_mean_NDTI(image):
        image = ee.Image(image)
        mean = image.reduceRegion(reducer = ee.Reducer.mean().setOutputs(['NDTI']),
                                geometry = AOI,
                                scale = image.projection().nominalScale().getInfo(),
                                maxPixels = 100000,
                                bestEffort = True);
        return mean.get('NDTI').getInfo()
        print("mean NDTI calculated")
    
# NDTI mean collection 
    Images_ndti = ndti.select('NDTI').toList(ndti.size())
    ndti_coll = []
    for i in range(Images_ndti.length().getInfo()):
        image = ee.Image(Images_ndti.get(i-1))
        temp_ndti = calculate_mean_NDTI(image)
        ndti_coll.append(temp_ndti)

# Dates Collection
    dates = np.array(ndti.aggregate_array("system:time_start").getInfo())
    day = [datetime.datetime.fromtimestamp(i/1000).strftime('%Y-%m-%d') for i in (dates)]

# Dataframe for Turtbidity

    df = pd.DataFrame(ndti_coll, index = day, columns = ['Turbidity'])
    df.index = pd.to_datetime(df.index, format="%Y/%m/%d")
    df.sort_index(ascending = True, inplace = True)

    return df

#------------------------------------------------------------------------------------------------------------------------------

def get_Salanity(start_date, end_date):

# Selecting the satellite and AOI  
# Sentinel 2A
# copernicus/s2_sr 
    sentinel = ee.ImageCollection("COPERNICUS/S2_SR").\
               filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).\
               filterDate(start_date, end_date)
    AOI = geometry

    sentinel_AOI = sentinel.filterBounds(AOI)

#calculate NDSI
    def calculate_NDSI(image):
        ndsi = image.normalizedDifference(['B11', 'B12']).rename('NDSI')
        return image.addBands(ndsi)
    ndsi = sentinel_AOI.map(calculate_NDSI)

# Mean NDSI
    def calculate_mean_NDSI(image):
        image = ee.Image(image)
        mean = image.reduceRegion(reducer = ee.Reducer.mean().setOutputs(['NDSI']),
                                geometry = AOI,
                                scale = image.projection().nominalScale().getInfo(),
                                maxPixels = 100000,
                                bestEffort = True);
        return mean.get('NDSI').getInfo()
        
# NDSI Mean Collection
    Images_ndsi = ndsi.select('NDSI').toList(ndsi.size())
    ndsi_coll = []
    for i in range(Images_ndsi.length().getInfo()):
        image = ee.Image(Images_ndsi.get(i-1))
        temp_ndsi = calculate_mean_NDSI(image)
        ndsi_coll.append(temp_ndsi)

# Dates Collection
    dates = np.array(ndsi.aggregate_array("system:time_start").getInfo())
    day = [datetime.datetime.fromtimestamp(i/1000).strftime('%Y-%m-%d') for i in (dates)]

# Dataframe for Salinity

    df = pd.DataFrame(ndsi_coll, index = day, columns = ['Salinity'])
    df.index = pd.to_datetime(df.index, format="%Y/%m/%d")
    df.sort_index(ascending = True, inplace = True)

    return df

#------------------------------------------------------------------------------------------------------------------------------

def get_DO(start_date, end_date):

# Selecting the satellite and AOI  
# Sentinel 2A
# copernicus/s2_sr 
    sentinel = ee.ImageCollection("COPERNICUS/S2_SR").\
               filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).\
               filterDate(start_date, end_date)
    AOI = geometry

    sentinel_AOI = sentinel.filterBounds(AOI)

# calculate Dissolved Oxygen
    def calculate_DO(image):
        do = image.normalizedDifference(['B4', 'B5']).rename('DO')
        return image.addBands(do)
    do = sentinel_AOI.map(calculate_DO)

# Mean DO

    def calcualte_mean_DO(image):
        image = ee.Image(image)
        mean = image.reduceRegion(reducer = ee.Reducer.mean().setOutputs(['DO']),
                            geometry = AOI,
                            scale = image.projection().nominalScale().getInfo(),
                            maxPixels = 100000,
                            bestEffort = True);
        return mean.get('DO').getInfo()
        
# DO Mean Collection

    Images_do = do.select('DO').toList(do.size())
    do_coll = []
    for i in range(Images_do.length().getInfo()):
        image = ee.Image(Images_do.get(i-1))
        temp_do = calcualte_mean_DO(image)
        do_coll.append(temp_do)

# Dates Collection
    dates = np.array(do.aggregate_array("system:time_start").getInfo())
    day = [datetime.datetime.fromtimestamp(i/1000).strftime('%Y-%m-%d') for i in (dates)]

# Dataframe for DO

    df = pd.DataFrame(do_coll, index = day, columns = ['Dissolved Oxygen'])
    df.index = pd.to_datetime(df.index, format="%Y/%m/%d")
    df.sort_index(ascending = True, inplace = True)

    return df

#------------------------------------------------------------------------------------------------------------------------------

def get_pH(start_date, end_date):

# Selecting the satellite and AOI  
# Sentinel 2A
# copernicus/s2_sr 
    sentinel = ee.ImageCollection("COPERNICUS/S2_SR").\
               filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20)).\
               filterDate(start_date, end_date)
    AOI = geometry

    sentinel_AOI = sentinel.filterBounds(AOI)

# calculate pH

    def calculate_pH(image):
        ph = ee.Image(8.339).subtract(ee.Image(0.827).multiply(image.select('B1').divide(image.select('B8')))).rename('PH')
        return image.addBands(ph)
    pH = sentinel_AOI.map(calculate_pH)

# Mean pH

    def calculate_mean_pH(image):
        image = ee.Image(image)
        mean = image.reduceRegion(reducer = ee.Reducer.mean().setOutputs(['PH']),
                                geometry = AOI,
                                scale = image.projection().nominalScale().getInfo(),
                                maxPixels = 100000,
                                bestEffort = True);
        return mean.get('PH').getInfo()
        
# pH Mean Collection

    Images_ph = pH.select('PH').toList(pH.size())
    ph_coll= []
    for i in range(Images_ph.length().getInfo()):
        image = ee.Image(Images_ph.get(i-1))
        temp_ph = calculate_mean_pH(image)
        ph_coll.append(temp_ph)

# Dates Collection

    dates = np.array(pH.aggregate_array("system:time_start").getInfo())
    day = [datetime.datetime.fromtimestamp(i/1000).strftime('%Y-%m-%d') for i in (dates)]

# Dataframe for pH

    df = pd.DataFrame(ph_coll, index = day, columns = ['pH'])
    df.index = pd.to_datetime(df.index, format="%Y/%m/%d")
    df.sort_index(ascending = True, inplace = True)

    return df

#------------------------------------------------------------------------------------------------------------------------------

def get_DOM(start_date, end_date):

# Selecting the satellite and AOI  
# Sentinel 3 OLCI
# Copernicus Sentinel-3 OLCI 
    sentinel3 = ee.ImageCollection("COPERNICUS/S3/OLCI").\
              filterDate(start_date, end_date)
    AOI = geometry

    sentinel3_AOI = sentinel3.filterBounds(AOI)

# calculate Dissolved Organic Matter (DOM) from Sentinel-3 OLCI

    def calculate_DM(image):
        rgb = image.select(['Oa08_radiance', 'Oa06_radiance', 'Oa04_radiance'])\
              .multiply(ee.Image([0.00876539, 0.0123538, 0.0115198]))
        dm = rgb.select('Oa08_radiance').divide(rgb.select('Oa04_radiance')).rename('dom')
        return image.addBands(dm)
    dm = sentinel3_AOI.map(calculate_DM)

# Mean DOM

    def calcualte_mean_DM(image):
        image = ee.Image(image)
        mean = image.reduceRegion(reducer = ee.Reducer.mean().setOutputs(['dom']),
                            geometry = AOI,
                            scale = image.projection().nominalScale().getInfo(),
                            maxPixels = 100000,
                            bestEffort = True);
        return mean.get('dom').getInfo()
    
# DOM Mean Collection

    Images_dm = dm.select('dom').toList(dm.size())
    dm_coll= []
    for i in range(Images_dm.length().getInfo()):
        image = ee.Image(Images_dm.get(i-1))
        temp_dm = calcualte_mean_DM(image)
        dm_coll.append(temp_dm)

# Dates Collection

    dates = np.array(dm.aggregate_array("system:time_start").getInfo())
    day = [datetime.datetime.fromtimestamp(i/1000).strftime('%Y-%m-%d') for i in (dates)]

# Dataframe for DOM

    df = pd.DataFrame(dm_coll, index = day, columns = ['dom'])
    df.index = pd.to_datetime(df.index, format="%Y/%m/%d")
    df.sort_index(ascending = True, inplace = True)

    return df

#------------------------------------------------------------------------------------------------------------------------------

def get_sm(start_date, end_date):

# Selecting the satellite and AOI  
# Sentinel 3 OLCI
# Copernicus Sentinel-3 OLCI 
    sentinel3 = ee.ImageCollection("COPERNICUS/S3/OLCI").\
              filterDate(start_date, end_date)
    AOI = geometry

    sentinel3_AOI = sentinel3.filterBounds(AOI)

# calculate suspended matter

    def calculate_SM(image):
        rgb = image.select(['Oa08_radiance', 'Oa06_radiance', 'Oa04_radiance'])\
              .multiply(ee.Image([0.00876539, 0.0123538, 0.0115198]))
        suspended_matter = rgb.select('Oa08_radiance').divide(rgb.select('Oa06_radiance')).rename('suspended_matter')
        return image.addBands(suspended_matter)
    sm = sentinel3_AOI.map(calculate_SM)

# Mean of Suspended Matter

    def meanSM(image):
        image = ee.Image(image)
        mean = image.reduceRegion(reducer = ee.Reducer.mean().setOutputs(['suspended_matter']),
                           geometry = AOI,
                           scale = image.projection().nominalScale().getInfo(),
                           maxPixels = 100000,
                           bestEffort = True);
        return mean.get('suspended_matter').getInfo()
            
# Suspended-Matter Mean Collection

    Images_sm = sm.select('suspended_matter').toList(sm.size())
    sm_coll= []
    for i in range(Images_sm.length().getInfo()):
        image = ee.Image(Images_sm.get(i-1))
        temp_sm = meanSM(image)
        sm_coll.append(temp_sm)

# Dates Collection

    dates = np.array(sm.aggregate_array("system:time_start").getInfo())
    day = [datetime.datetime.fromtimestamp(i/1000).strftime('%Y-%m-%d') for i in (dates)]

# Dataframe for suspended matter

    df = pd.DataFrame(sm_coll, index = day, columns = ['Suspended Matter'])
    df.index = pd.to_datetime(df.index, format="%Y/%m/%d")
    df.sort_index(ascending = True, inplace = True)

    return df


#------------------------------------------------------------------------------------------------------------------------------

turbidity = get_turbidity('2018-01-01', '2022-12-31')
print("Turbidity Collection Completed")

DOM = get_DOM('2018-01-01', '2022-12-31')
print("DOM Collection Completed")

suspended_matter = get_sm('2018-01-01', '2022-12-31')
print("Suspended Matter Collection Completed")

pH = get_pH('2018-01-01', '2022-12-31')
print("pH Collection Completed")

dissolved_oxygen = get_DO('2018-01-01', '2022-12-31')
print("Dissolved Oxygen Collection Completed")

# Exporting the dataframes to csv files

turbidity.to_csv('turbidity_mirzapur.csv', index=True)
pH.to_csv('ph_mirzapur.csv', index=True)
DOM.to_csv('dom_mirzapur.csv', index=True)
dissolved_oxygen.to_csv('dissolved_oxygen_mirzapur.csv', index=True)
suspended_matter.to_csv('suspended_matter_mirzapur.csv', index=True)

