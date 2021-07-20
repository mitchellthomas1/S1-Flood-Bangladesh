#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 14:48:18 2020

@author: Mitchell

"""

import ee
import sys
import pandas as pd
import numpy as np
import zscore
from baselinecreation import calc_single_baseline
ee.Initialize()

#path for saving csv locally, change this based on your file system
path = 'timeseriesoutput'

# Earth engine bangaldesh admin data
country = ee.FeatureCollection("users/mlt2177/Adm0old")
divisions = ee.FeatureCollection("users/mlt2177/Adm1")
districts = ee.FeatureCollection("users/mlt2177/Adm2")
upazilas = ee.FeatureCollection("users/mlt2177/Adm3")


#returns: dataframe with year's data
#saves dataframe to csv with filename: "{}-{}-{}-scale{}-bline{}{}.csv"
def make_s1_ts(study_area, area_descrip,targdates, setScale, basedates = None, baseline_len = 45, retry_times = 5,  include_med = False):
    ''' Function to make a sentinel-1 time series of flooded area over either a 
    Bangladesh admin district, buffer zone around a given point, or around any GEE geometry

    Parameters:

    study_area  and area_descrip:

    ***** 3 modes for geometries to study **********
    ----- for BGD admin districts (default mode): ---------
          study_area (str): name of districts/division/upazila
          area_descrip (str): adminstrative unit, choose 'country', 'division', 'district', 'upazila'
    ------ to evaluate buffer around point -------
         study_area (tuple) = ([lon, lat], radius) : tuple with list of coordinate points 
            (list of floats) and radius (scalar value) in meters
         area_descrip = choose 'buffer' 
    ------- to evaluate custom geometry ----------
         study_area = some ee.Geometry object
         area_descrip = some label/description for the object

    targdates (tuple of str): tuple of start and end dates both in 'YYYY-MM-DD' format

    setScale (int or float): scale of reduction for earth engine aggregations. This will affect the 
        speed in which operations are run. Recommended 10 for most accurate data 
        and higher values (ex: 100 or 200) for fast performance.

    basedates (tuple of str, optional) : pass manual baseline period (basestart, baseend) to override automatic baseline creation
    
    retry_time (int, optional) : times to retry while getting an error

    
    Returns:
    dataframe with full data
    Continuously saves dataframe to csv in path with filename: "{}-{}-{}-scale{}-bline{}{}.csv"

        




    '''

    #    Threholds for algorithm:
    # Z-score thresholds
    zvv_thd = -2
    zvh_thd = -2
    # threshold for raw VH
    rawvh_th = -22

#    smoothing radius and threshold in meters
    smoothing_radius = 30
    smoothing_threshold = 0.6


    #    Shapefile assets
    targstart, targend = targdates

    area_descrips = {"division": (divisions, 'NAME_1'),
                  "district": (districts, "NAME_2"),
                  "upazila": (upazilas, 'NAME_3'),
                  "country": (country, 'NAME_0')}



    def make_collections(study_area, area_descrip,targdates, setScale, basedates, baseline_len):

        if area_descrip == 'buffer':
            buff_rad = int(study_area[1]) * 1000
            buff_coord = list(study_area[0])
            print(buff_coord)
            study_feature = ee.Feature(ee.Geometry.Point(buff_coord).buffer(buff_rad))
            study_area = str(buff_coord[0])+'-'+str(buff_coord[1]) +'-'+str(int(buff_rad/1000))+'km'
        elif type(study_area) == ee.geometry.Geometry:
            study_feature = ee.Feature(study_area)
            study_area = 'passedgeom'
        else:
            #   create region in Earth engine
            study_feature  = area_descrips[area_descrip][0]\
                                .filter(ee.Filter.eq(area_descrips[area_descrip][1], study_area)).first()
    #    area of region in m^2
        total_area = study_feature.geometry().area(10).getInfo()

        targstart, targend = targdates[0], targdates[1]
    #    set baseline to override if basedates = True, otherwise automatically create baseline
        if basedates:
            basestart, baseend = basedates[0], basedates[1]
        else:
            basestart, baseend = calc_single_baseline(study_feature.geometry(), baseline_len)
            print("Baseline Automatically Chosen: ({})-({})".format(basestart, baseend))


        print(study_area, targstart, targend, "scale = ", setScale)
        # provides a month of images before the start of the target period in order to avoid mosaicking issues
        targbefore = ee.Date(targstart).advance(-1,'month')
        # Define collection filters (VV+VH, IW mode)
        filters = [
          ee.Filter.listContains("transmitterReceiverPolarisation", "VV"),
          ee.Filter.listContains("transmitterReceiverPolarisation", "VH"),
          ee.Filter.equals("instrumentMode", "IW"),
          ]
        # function to remove image border inaccuracy
        def correctBorders(image):
          maskedImage = (image.updateMask(image.select('VV').gt(-45)))
          return maskedImage
        # Load S1 collection with filters
        s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filter(filters).filterBounds(study_feature.geometry())
         # create and merge baseline and target period stacks
        s1_combined = (s1.filterDate(basestart,baseend)
                    .merge(s1.filterDate(targbefore,ee.Date(targend).advance(1,'day')))
                    .map(correctBorders))
    #    create list of dates and times where images occur
        timeList = s1_combined.filterDate(targstart, ee.Date(targend).advance(1,'day')) \
                                    .sort('system:time_start').aggregate_array('system:time_start')

    #     converts system:time_start to date
        def to_date(time):
            return ee.Date(time).format()
        date_li_py = timeList.map(to_date).getInfo()
        time_li_py = timeList.getInfo()
        orbit_list = s1_combined.filterDate(targstart, ee.Date(targend).advance(1,'day')) \
                                    .sort('system:time_start').aggregate_array('orbitProperties_pass').getInfo()
    #    get rid of several images on a day:
        prev_date = ""
        remove_indices = []
        for i, date_time in enumerate(date_li_py):
            date = date_time[:10]
            if date == prev_date:
                remove_indices.append(i-1)
            prev_date = date
#       remove those indices from time and orbit lists
        time_li_py = [i for j, i in enumerate(time_li_py) if j not in remove_indices]
        orbit_list = [i for j, i in enumerate(orbit_list) if j not in remove_indices]
        # filter s1 by ascending and descending
        s1_asc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'ASCENDING'))
        s1_dsc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'DESCENDING'))
        # Compute Z-scores per orbital direction and merge into zscore
        z_iwasc = zscore.calc_zscore(s1_asc, basestart, baseend, 'IW', 'ASCENDING')
        z_iwdsc = zscore.calc_zscore(s1_dsc, basestart, baseend, 'IW', 'DESCENDING')
        z_coll = ee.ImageCollection(z_iwasc.merge(z_iwdsc)).sort('system:time_start')

        return z_coll, time_li_py, study_feature, total_area, study_area, area_descrip,  s1_combined, orbit_list




    # main GEE function to loop and create a time series
    def create_data_point(time_start, s1_combined, z_coll, orbit):

#          powArea = ee.Number(POWImage.multiply(ee.Image.pixelArea()).reduceRegion(reducer = ee.Reducer.sum(), geometry = study_feature.geometry(),
#                                              scale = setScale, maxPixels =  1e15).values().get(0)).divide(total_area).getInfo()
#          print('pow:',powArea)
          date_val = ee.Date(time_start)
          date_start = date_val.advance(1,'day')
          # create custom z-collection based on desired study date
          z_time_adjust = z_coll.filterDate(date_val.advance(-1,'month'), date_start).sort('system:time_start')

          s1_image = s1_combined.filterDate(date_val, date_start).first()
          s1_proj = s1_image.select('VH').projection()

          z_image = z_time_adjust.mosaic().setDefaultProjection(s1_proj, None, 10)

          s1_study = s1_combined.filterDate(date_val.advance(-1,'month'), date_start).mosaic().setDefaultProjection(s1_proj, None, 10)

          vvflag = z_image.select('VV').lte(zvv_thd)
          vhflag = z_image.select('VH').lte(zvh_thd)

          rawvhflag = s1_study.select('VH').lte(rawvh_th)

          flood_class = (vvflag).add(vhflag).add(rawvhflag).rename(['flood_sum'])
          bool_map = flood_class.gte(1)

          smooth = bool_map.focal_mean(radius = smoothing_radius, units = 'meters').gt(smoothing_threshold).setDefaultProjection(s1_proj, None, 10)

          boolArea= smooth.multiply(ee.Image.pixelArea())

          zonal_stats = boolArea.reduceRegion(reducer = ee.Reducer.sum(), geometry = study_feature.geometry(),
                                              scale = setScale, maxPixels =  1e15)

          tries = 0
          passed3 = False
          while tries < 5 and passed3 == False:
              tries +=1
              try:
    #                  print(zonal_stats.getInfo())
#                  print(zonal_stats.getInfo())
                  sum_val = zonal_stats.get('flood_sum').getInfo()
                  passed3 = True

              except Exception as inst3:
                  print("3rd inner exception")
                  print("type:", type(inst3))    # the exception instance
                  print("instance:" , inst3)   #example of exception I was seeing: User memory limit exceeded.
                  pass
          if passed3 == False:
              raise(ee.ee_exception.EEException)

          frac_flood = sum_val / total_area
          return date_val.format().getInfo(), sum_val, frac_flood


    z_coll, time_li_py, study_feature, total_area, new_name, new_type,  s1_coll, orbit_list = make_collections(study_area, area_descrip,targdates, setScale, basedates, baseline_len)
    dates = []
    vals = []
    ratios = []
    csv_filename = path + "/{}-{}-{}-{}-scale{}.csv".format(new_name,new_type,
                                     targstart.replace('-',''),targend.replace('-',''), setScale)

    for i, time_start in enumerate(time_li_py):
        passed1 = False
        while passed1 == False:
            try:
                print(i +1 , "/" , len(time_li_py))
                passed2 = False
                j = 0
                tries2 = 3
                while j < tries2 and passed2 == False:
                    j +=1
                    try:
                        date, val, ratio = create_data_point(time_start, s1_coll, z_coll, orbit_list[i])
                        passed2 = True
                        passed1 = True
                        break
                    except ee.ee_exception.EEException as inst2:
                        print("2nd inner exception")
                        print("type:", type(inst2))    # the exception instance and its type
                        print("instance:" , inst2)
                if passed2 == False:
                    raise(ee.ee_exception.EEException)


            except ee.ee_exception.EEException as inst1:
                print("1st outer exception")
                print("info:", sys.exc_info()[0])
                print("type:", type(inst1))    # the exception instance  type
                print("args:", inst1.args)     # arguments stored in .args
                print("instance:" , inst1)
#                redefine z
                created_new = False
                while created_new == False:
                    try:
                        ee.Initialize()
                        z, f1, f2, f3, f4, f5, f6 = make_collections(study_area,
                                                area_descrip, targdates, setScale, basedates, baseline_len)
                        created_new = True
                    except ee.ee_exception.EEException:
                        print('create_new_failed')
                        pass



        print("date: ", date, "  fraction: ", ratio)
        dates.append(date)
        vals.append(val)
        ratios.append(ratio)
#        add orbital direction
        orbits= ee.ImageCollection("COPERNICUS/S1_GRD").filterDate(targstart, targend).filterBounds(study_feature.geometry())\
                            .aggregate_array('orbitProperties_pass').getInfo()
        pd_dict = {"Dates": dates, "Flooded Area (m^2)": vals, "Fractional flooded Area" : ratios,'orbitalDirection':orbits[:i+1]}
        df = pd.DataFrame.from_dict(pd_dict, orient = "columns")
        df.to_csv( csv_filename, index = False)

    return df