#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 10:43:00 2020

@author: Mitchell

Code to predict dryest period of year in areas in Northern Bangladesh
"""

import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import ee
import matplotlib.dates as mdates
import geopandas as gp
ee.Initialize()

#-----Assets-----
#nation = ee.FeatureCollection("users/mlt2177/Adm0")
#divisions = ee.FeatureCollection("users/mlt2177/Adm1")
districts = ee.FeatureCollection("users/mlt2177/Adm2")
#upazilas = ee.FeatureCollection("users/mlt2177/Adm3")
#PRECIP = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD")
#district_li = districts.aggregate_array("NAME_2").getInfo()
#regions_of_interest = ["Sirajgonj", 'Kurigram', 'Gaibanda','Jamalpur', 'Sylhet']
sirajganj = districts.filterMetadata("NAME_2","equals","Sirajganj").first()
sylhet = districts.filterMetadata("NAME_2","equals","Sylhet").first()
#-----------------
SoilMoisture = ee.ImageCollection("NASA_USDA/HSL/soil_moisture")
start = '2018-10-01';
end = '2019-04-01';

scale = 50
#conversion of earth engine dates to python dates
def to_dt(date_str):
    return dt.date.fromisoformat(date_str[0:10])



def create_sm_ts(geometry):

    def vals_and_dates(image):
        val = image.reduceRegion(ee.Reducer.mean(), geometry, 100, maxPixels = 10e11)\
                    .values().get(0)
        date_str = ee.Date(image.get('system:time_start')).format()
        return image.set("Date",date_str).set("Mean",val)
    # image collections
    start, end = '2018-10-01', '2019-04-01'
    sm_coll = SoilMoisture.filterDate(start, end).filterBounds(geometry)\
                          .select('ssm')\
                          .map(vals_and_dates)
    s_vals = sm_coll.aggregate_array('Mean').getInfo()
    s_dates = sm_coll.aggregate_array('Date').getInfo()
    index = pd.to_datetime(s_dates , format = '%Y-%m-%dT%H:%M:%S')
    ts = pd.Series(s_vals, name = 'Soil Moisture', index = index)
    return ts

def find_min(sm_ts,geometry, baseline_len = 45):
    smooth = False
    if smooth:
        sm_ts = sm_ts.rolling(window=3).median()
    daily_ts = sm_ts.resample('D').mean().interpolate()
#    get asc/dsc info
    def add_dt_str(image):
        date = ee.Date(image.get("system:time_start"))
        return image.set("Date", date.format()).copyProperties(image)
    s1_coll = ee.ImageCollection("COPERNICUS/S1_GRD").filterDate(start, end)\
                .filterBounds(geometry).map(add_dt_str)
    s1_dates = s1_coll.aggregate_array('Date').getInfo()
    s1_orbit = pd.Series(s1_coll.aggregate_array('orbitProperties_pass').getInfo())
    s1_orbit.index = pd.to_datetime(s1_dates,format = '%Y-%m-%dT%H:%M:%S')
    

    bl_start_index = 0
    bl_end_index = baseline_len - 1
    
#    start with very high number
    lowest_mean = 10e10
    for i in range(daily_ts.values.size - baseline_len - 1):
        study_arr = daily_ts.values[i: i + baseline_len - 1]
        mean = study_arr.mean()
#        check if this mean is lower
        if mean < lowest_mean:
#            check if there are ascending and descending images
            test_start = daily_ts.index.values[i]
            test_end = daily_ts.index.values[i + baseline_len - 1]
            test_orbits = s1_orbit[test_start: test_end].values
            if ('ASCENDING' in test_orbits) and ('DESCENDING' in test_orbits):
                lowest_mean = mean
                bl_start_index = i
                bl_end_index = i + baseline_len - 1
    date_range = (daily_ts.index.values[bl_start_index].astype(str)[0:10],
                    daily_ts.index.values[ bl_end_index].astype(str)[0:10])
    return date_range
    

#plt.plot(ts)
#plt.show()
#plt.plot(ts.rolling(window=3).median())
#plt.show()
#
#
#
#tup = find_min(ts, district)



def calc_single_baseline(study_geometry, baseline_len = 45):

    ts = create_sm_ts(study_geometry)
    start, end = find_min(ts, study_geometry, baseline_len)
    return start, end



#sw175_5 = ee.Geometry.Point([91.68,   24.6277]).buffer(15*1000)
#print(calc_single_baseline(sw175_5))





def make_paper_plot():
    sirajganj_ssm = create_sm_ts(sirajganj.geometry())
    siraj_bs, siraj_be = calc_single_baseline(sirajganj.geometry())
    # "4th weeek of Jan to 3rd week of Feb"
    siraj_rs, siraj_re = pd.Timestamp(2019, 1, 4*7) , pd.Timestamp(2019, 2, 4*7)
    
    sylhet_ssm = create_sm_ts(sylhet.geometry())
    sylhet_bs,sylhet_be = calc_single_baseline(sylhet.geometry())
    #"Mid novermber to mid december"
    sylhet_rs, sylhet_re = pd.Timestamp(2018, 11, 15) , pd.Timestamp(2018, 12, 16)
    
    plt.rcParams["font.family"] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'
    
    f,  axs = plt.subplots(2,2, figsize =(8,8), sharey=True)
    ((noneAxis, mapax), (ax1, ax2)) = axs
    gs = axs[0,0].get_gridspec()
    for ax in axs[0]:
        ax.remove()
    axbig = f.add_subplot(gs[0,:])
    f.subplots_adjust(hspace= -0.2, wspace=0.025)
    axbig.axis('off')
    
    districts = gp.read_file("/Users/Mitchell/Documents/LDEOinternship/Shapefiles_Final_30_7_2020/Districts_Shapefile_7_15_2020/district.shp")
    siraj_sylhet = districts[(districts['NAME_2'] == 'Sylhet' ) | (districts['NAME_2'] == 'Sirajganj' )]
    districts.boundary.plot(ax = axbig, color = '#7e7e7e', alpha = 0.5)
    districts.plot(ax = axbig, color = '#a9a9a9', alpha = 0.5)
    siraj_sylhet = siraj_sylhet.plot(ax = axbig, color = 'purple', alpha = 0.6)
    
    #siraj arrow 
    color = '#444444'
    axbig.arrow(89.5,24.4,-2,-4.5,fc=color,ec=color,clip_on=False, width = 0.06, length_includes_head = True, color = 'black',alpha = 0.6)
    #sylhet arrow 
    axbig.arrow(92,24.8,2,-4.9,fc=color,ec=color,clip_on=False, width = 0.06, length_includes_head = True, color = 'black',alpha = 0.6)
    #axbig.arrow(0.5,0.4,0.1,0.1,fc='b',ec='b',clip_on=False)
    
    
    
    ax1.plot(sirajganj_ssm)
    #ax1.plot(sirajganj_bl_ssm, color = 'red')
    #ax1.legend(['Soil Moisture Time Series','Sirajganj Baseline Period'], fontsize = 11)
    ax1.set_title('A. Sirajganj', fontsize = 14)
    ax1.set_xlabel('Date', fontsize = 12)
    ax1.set_ylabel('SMOS Soil Moisture (mm)', fontsize = 12)
    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=2)) 
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m')) 
    ax1.axvspan(pd.Timestamp(siraj_bs), pd.Timestamp(siraj_be), alpha=0.2, color='red', hatch = "\\" )
    ax1.axvspan(siraj_rs, siraj_re, alpha=0.2, color='blue', hatch = '/')
    
    ax2.plot(sylhet_ssm)
    ax2.set_title('B. Sylhet', fontsize = 14)
    ax2.set_xlabel('Date', fontsize =12)
    ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=2)) 
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m')) 
    
    ax2.axvspan(pd.Timestamp(sylhet_bs), pd.Timestamp(sylhet_be), alpha=0.2, color='red', hatch = "\\" )
    ax2.axvspan(sylhet_rs, sylhet_re, alpha=0.2, color='blue', hatch = '/')
    
    lgd = f.legend(['soil moisture','baseline period', 'rice transplantation period'],
            bbox_to_anchor = (0.85,0.65),
            loc = 2)
    
    plt.tight_layout()
    plt.savefig('/Users/Mitchell/Documents/LDEOinternship/AguImages/UpdatedBaselineComparisonWMap.png',dpi = 700,
                 bbox_extra_artists=(lgd,), bbox_inches='tight', pad_inches = 0.1)
    
    plt.show()
#make_paper_plot()




#study_upazilas = upazilas.filter(ee.Filter.And(ee.Filter.Or(
#    ee.Filter.eq('NAME_3', 'Bhuapur'),ee.Filter.eq('NAME_3', 'Nagarpur'),
#    ee.Filter.eq('NAME_3', 'Tangail Sadar'),ee.Filter.eq('NAME_3', 'Dhunat'),
#    ee.Filter.eq('NAME_3', 'Sariakandi'),ee.Filter.eq('NAME_3', 'Bera'),
#    ee.Filter.eq('NAME_3', 'Kazipur'),ee.Filter.eq('NAME_3', 'Belkuchi'),
#    ee.Filter.eq('NAME_3', 'Chauhali'),ee.Filter.eq('NAME_3', 'Shahjadpur'),
#    ee.Filter.eq('NAME_3', 'Sirajganj Sadar'))))
#
#upa_list = study_upazilas.toList(15).getInfo()
#
##for upa in upa_list:
##    feature = ee.Feature(upa)
##    geom= feature.geometry()
##    start,end = calc_single_baseline(geom)
##    print(feature.get('NAME_3').getInfo(), start, end)
    
    
    
    
    
    
    
    
#    --------- OLD VERSION: ---------

# takes a study geometry a and returns precipitation and soilmoisture ts
#def process_coll(study_geom):
#
#    def vals_and_dates(image):
#        val = image.reduceRegion(ee.Reducer.mean(), study_geom, scale, maxPixels = 10e11)\
#                    .values().get(0)
#        date_str = ee.Date(image.get('system:time_start')).format()
#        return image.set("Date",date_str).set("Mean",val)
#    # image collections
#    precip_coll = PRECIP.select("precipitation").filterDate(start, end) \
#                  .filterBounds(study_geom)\
#                  .map(vals_and_dates)
#    sm_coll = SoilMoisture.filterDate(start, end).filterBounds(study_geom)\
#                          .select('ssm')\
#                          .map(vals_and_dates)
#
#    p_vals = precip_coll.aggregate_array('Mean').getInfo()
#    p_dates = precip_coll.aggregate_array('Date').getInfo()
#    s_vals = sm_coll.aggregate_array('Mean').getInfo()
#    s_dates = sm_coll.aggregate_array('Date').getInfo()
#
#    return p_vals, p_dates, s_vals, s_dates



#interpolates time series to provide daily data for both soil moisture and precipitation
#returns data frame with soil moisture and precipitation time series
#def interpolate_data(p_vals, p_dates, s_vals, s_dates):
#
#    precip = pd.DataFrame.from_dict({"Date" : p_dates , "precipitation": p_vals})
#    precip["datetime"] = pd.to_datetime(precip["Date"].apply(to_dt))
#    precip.index = precip["datetime"]
#    del precip["datetime"]
#    del precip["Date"]
#    precip_interpol = precip.resample('D').mean()
#    precip_interpol['precipitation'] = precip_interpol['precipitation'].interpolate()
#    precip_interpol.head(365)
#
##  adds soil moisture to
#    sm = pd.DataFrame.from_dict({"Date" : s_dates , "ssm": s_vals})
#    sm["datetime"] = pd.to_datetime(sm["Date"].apply(to_dt))
#    sm["ssm"] = sm["ssm"]
##    print(1)
##    plt.plot(sm["datetime"], sm["ssm"])
##    plt.show()
##    print(2)
#    sm.index = sm["datetime"]
#    
#    del sm["datetime"]
#    del sm["Date"]
#    sm_interpol = sm.resample('D').mean()
#    sm_interpol['ssm'] = sm_interpol['ssm'].interpolate()
#    sm_interpol.head(365)
##    start and end date for each time series such that they are the same length
#    start_date = max(sm_interpol.index.values[0] , precip_interpol.index.values[0])
#    end_date = min(sm_interpol.index.values[-1] , precip_interpol.index.values[-1])
##    slices each by start/end date to provide equal length time series
#    precip_interpol_sliced = precip_interpol.loc[start_date:end_date]
#    sm_interpol_sliced = sm_interpol.loc[start_date:end_date]
##    puts precip and soil moisture data frame together
#    final_df = precip_interpol_sliced.copy()
#    final_df["ssm"] = sm_interpol_sliced["ssm"]
#    return final_df



#calculates baseline from data frame
##returns tuple of strings: (start, end)
#def calc_baseline(df, geom, baseline_len = 45, plot = False):
##    normalizes soil moisture and precipitation
##    df["SMNorm"] = df['ssm'] / df['ssm'].values.max()
##    df["PrecipNorm"] = df["precipitation"] / df["precipitation"].values.max()
##    df["1/precip*sm"] = 1 / ((0.5* df["PrecipNorm"].values) * df["SMNorm"].values)
##    df["1/sm"] = 1 / (df["SMNorm"].rolling(window=4).median())# + df["PrecipNorm"])
##    df["1/sm+precip"] = 1 / (df["SMNorm"].rolling(window=5).median() + df["PrecipNorm"])
#
##    df.to_csv("validation/soilmoisturets{}.csv".format("Jamalpur"))
#    
##    make s1 stack to check ascending and descending
#    def add_dt_str(image):
#        date = ee.Date(image.get("system:time_start"))
#        return image.set("Date", date.format()).copyProperties(image)
#
#    s1_coll = ee.ImageCollection("COPERNICUS/S1_GRD").filterDate(start, end)\
#                .filterBounds(geom).map(add_dt_str)
#    s1_dates = s1_coll.aggregate_array('Date').getInfo()
#    global s1_orbit
#    s1_orbit = pd.Series(s1_coll.aggregate_array('orbitProperties_pass').getInfo())
#    s1_orbit.index = pd.to_datetime(s1_dates,format = '%Y-%m-%dT%H:%M:%S')
#    
#
#    bl_start_index = 0
#    bl_end_index = baseline_len - 1
#    study_vals = df['ssm'].rolling(window=3).median()
#    
##    plt.plot(study_vals)
#    
#    highest_mean = 10e10
#    for i in range(study_vals.values.size - baseline_len):
#        study_arr = study_vals.values[i: i + baseline_len - 1]
#        mean = study_arr.mean()
##        print(mean)
#        if mean < highest_mean:
##            check if there are ascending and descending images
#            test_start = df.index.values[i]
#            test_end = df.index.values[i + baseline_len - 1]
#            test_orbits = s1_orbit[test_start: test_end].values
#            if ('ASCENDING' in test_orbits) and ('DESCENDING' in test_orbits):
#
#                highest_mean = mean
#                bl_start_index = i
#                bl_end_index = i + baseline_len - 1
##    print(std)
#    final_start = df.index.values[bl_start_index].astype(str)[0:10]
#    final_end = df.index.values[bl_end_index].astype(str)[0:10]
#    
#    global sirajganj_ssm
#    global sirajganj_bl_ssm
#    sirajganj_ssm = study_vals
#    sirajganj_bl_ssm = study_vals[final_start:final_end]
#    
#
#
#    return final_start, final_end, df


#calculates baseline for single geometry
#params- district: district name as string
#        baseline_len = 45 desired length of baseline period


#gets all baselines for list of region names as strings (ensure proper spelling)
#def get_all_baselines(regions_li, baseline_len = 45):
#    bl_dict = {}
#    for district in regions_li:
#        print(district)
#        p_v, p_d, s_v, s_d = process_coll(district)
#        df = interpolate_data(p_v, p_d, s_v, s_d)
#        start, end = calc_baseline(df, baseline_len = baseline_len)
#        bl_dict[district]  = [start,end]
#    return bl_dict


#dis = ee.FeatureCollection('users/mlt2177/Adm2')
#district = dis.filterMetadata('NAME_2','equals','Sirajganj').first().geometry()
#start, end = calc_single_baseline(district)
#print(start,end)










#study_dis = ['Sirajganj','Sylhet','Natore','Jamalpur']
#for district_name in study_dis:
#    district = dis.filterMetadata('NAME_2','equals',district_name).first().geometry()
#    start, end = calc_single_baseline(district)
#    print(district_name, start,end)

#bl_dict = get_all_baselines(regions_of_interest, baseline_len = 45)


#2020_flood_dis = ['Kurigram', 'Gaibanda', 'Bogra', 'Sirajgonj', 'Rajbari', 'Tangaif', 'Mankiganj']
    
    
#testing results
#    sylhet: 2019-01-03 2019-02-16 F1 = 0.901, bias = 1.05
#            2019-02-03 2019-03-19
    
    

    
    
    
