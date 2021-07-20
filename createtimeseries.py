#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 14:48:18 2020

@author: Mitchell

Creates Sentinel-1 flooded area time series csv
Set parameters at bottom of code
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


#function to make a sentinel-1 time series of flooded area
#through variation of DeVries algorithm


#Parameters:

# ***** 3 modes for geometries to study **********
#----- for BGD admin districts (default mode): ---------
#       area_name: name of districts/division/upazila
#       area_type: pass 'division', 'district', or 'upazila'
#------ to evaluate buffer around point -------
#      area_name (tuple) = ([lon, lat], radius)
#      area_type = 'buffer'
#------- to evaluate custom geometry ----------
#      area_name = some ee.Geometry object
#      area_type = some label/description for the object


#targdates: pass tuple with (targstart, targend)
#setScale: (int) scale of reduction (recommended = 10 for best results)
# basedates (optional) : pass manual baseline period (basestart, baseend) to override automatic baseline creation
#retry_time (optional) : times to retry while getting an error

#returns: dataframe with year's data
#saves dataframe to csv with filename: "{}-{}-{}-scale{}-bline{}{}.csv"
def make_s1_ts(area_name, area_type,targdates, setScale, basedates = None, baseline_len = 45, retry_times = 5,  include_med = False):

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
    country = ee.FeatureCollection("users/mlt2177/Adm0old")
    divisions = ee.FeatureCollection("users/mlt2177/Adm1")
    districts = ee.FeatureCollection("users/mlt2177/Adm2")
    upazilas = ee.FeatureCollection("users/mlt2177/Adm3")
    area_types = {"division": (divisions, 'NAME_1'),
                  "district": (districts, "NAME_2"),
                  "upazila": (upazilas, 'NAME_3'),
                  "country": (country, 'NAME_0')}



    def make_collections(area_name, area_type,targdates, setScale, basedates, baseline_len):

        if area_type == 'buffer':
            buff_rad = int(area_name[1]) * 1000
            buff_coord = list(area_name[0])
            print(buff_coord)
            study_feature = ee.Feature(ee.Geometry.Point(buff_coord).buffer(buff_rad))
            area_name = str(buff_coord[0])+'-'+str(buff_coord[1]) +'-'+str(int(buff_rad/1000))+'km'
        elif type(area_name) == ee.geometry.Geometry:
            study_feature = ee.Feature(area_name)
            area_name = 'passedgeom'
        else:
            #   create region in Earth engine
            study_feature  = area_types[area_type][0]\
                                .filter(ee.Filter.eq(area_types[area_type][1], area_name)).first()
    #    area of region in m^2
        total_area = study_feature.geometry().area(10).getInfo()

        targstart, targend = targdates[0], targdates[1]
    #    set baseline to override if basedates = True, otherwise automatically create baseline
        if basedates:
            basestart, baseend = basedates[0], basedates[1]
        else:
            basestart, baseend = calc_single_baseline(study_feature.geometry(), baseline_len)
            print("Baseline Automatically Chosen: ({})-({})".format(basestart, baseend))


        print(area_name, targstart, targend, "scale = ", setScale)
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

        return z_coll, time_li_py, study_feature, total_area, area_name, area_type,  s1_combined, orbit_list




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


#          bool_map = ( (z_image.select('VV').lte(zvv_thd)).add(z_image.select('VH').lte(zvh_thd)).add(POWImage)).gte(1).rename(['flood_sum'])


        # create boolean image
#          image_sum = floods.eq(3).add(floods.eq(13))
#
#          if include_pow:
#              image_sum = image_sum.add(POWImage)
#          if include_med:
#              s1_study = s1_combined.filterDate(date_val.advance(-1,'month'), date_start).mosaic()
#              rawVVth, rawVHth =  -14, -14
#              lowConf = floods.gte(1).lte(2).add(floods.gte(11).lte(12))
#              belowVVRaw = s1_study.select('VV').lt(rawVVth)
#              belowVHRaw = s1_study.select('VH').lt(rawVHth)
#              lowConfFloodMap = lowConf.add(belowVVRaw).add(belowVHRaw).eq(3)
#              image_sum = image_sum.add(lowConfFloodMap)

#          bool_map = image_sum.gte(1)
          smooth = bool_map.focal_mean(radius = smoothing_radius, units = 'meters').gt(smoothing_threshold).setDefaultProjection(s1_proj, None, 10)
#          print('scale', bool_smooth.projection().nominalScale().getInfo())

          #multiply the pixel by its area
          boolArea= smooth.multiply(ee.Image.pixelArea())
          # reduce image into single variable

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


    z_coll, time_li_py, study_feature, total_area, new_name, new_type,  s1_coll, orbit_list = make_collections(area_name, area_type,targdates, setScale, basedates, baseline_len)
    dates = []
    vals = []
    ratios = []
    csv_filename = path + "/{}-{}-{}-{}-scale{}-test2-4-27.csv".format(new_name,new_type,
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
                        z, f1, f2, f3, f4, f5, f6 = make_collections(area_name,
                                                area_type, targdates, setScale, basedates, baseline_len)
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




#main function to run code and save csv
# pass a_name, a_type, targstart, targend,  scale as
# command line args or system will ask for them
if __name__ == '__main__':

#    if len(sys.argv) != 6:
#        a_name = input("Name of Region: ")
#        a_type = input("Type of Region (division, district, upazila): ")
#        targstart = input("Target Start: ")
#        targend = input("Target End: ")
#        scale = int(input("scale: "))
#    else:
#        filename, a_name, a_type, targstart, targend,   scale = sys.argv
#
##    a_name, a_type, targstart, targend,   scale = 'Sunamganj','district','2017-05-01','2017-06-01',20
#    df = make_s1_ts(a_name, a_type,(targstart,targend), int(scale), retry_times = 100)

#    targstart, targend = '2017-01-01', '2020-12-31'
    targstart, targend = '2017-01-01', '2020-12-31'
    upazilas = ee.FeatureCollection("users/mlt2177/Adm3")
    upa_names = upazilas.aggregate_array('NAME_3').getInfo()
    dis_names = upazilas.aggregate_array('NAME_2').getInfo()
    div_names = upazilas.aggregate_array('NAME_1').getInfo()
    
    
    filelog_file = 'datasetV2/filelog.csv'
    try:
        log_df = pd.read_csv(filelog_file)
    except FileNotFoundError:
        
        log_df = pd.DataFrame.from_dict({'Index': list(range(len(upa_names))) ,
                                'Upazila': upa_names,'District':dis_names,
                                'Division': div_names,'Started': np.zeros(len(upa_names), dtype = bool),
                                'Completed': np.zeros( len(upa_names) , dtype = bool)})
        log_df.to_csv(filelog_file)
        

    
    todo = log_df['Index'][ ~ log_df['Started'].astype(bool)]
    del log_df
    
    while len(todo) > 1:
        print(todo)
        i = todo.values[0]
#        update log
        log_df  = pd.read_csv(filelog_file)
        log_df.index = log_df['Index']
        log_df.loc[i, 'Started'] = True
        log_df.to_csv(filelog_file)
        del log_df
        
#        run time series
        upa, dis, div = upa_names[i], dis_names[i], div_names[i]
        name_str = div + '_' + dis + '_' + upa
        study_upa = upazilas.filterMetadata('NAME_3','equals',upa) \
                    .filterMetadata('NAME_2','equals',dis).filterMetadata('NAME_1','equals',div) \
                    .first().geometry()
    
        make_s1_ts(study_upa, name_str, (targstart,targend), 10)
        
        
        # once finished, update log
        #        update log
        log_df  = pd.read_csv(filelog_file)
        log_df.index = log_df['Index']
        log_df.loc[i, 'Completed'] = True
        log_df.to_csv(filelog_file)
        del log_df



#    IIASA_coord = pd.read_csv('/Users/Mitchell/Documents/LDEOinternship/S1Dataset/IIASAFiles/CoordinatesForIIASA.csv')
#    coordinates = IIASA_coord[['Longitude','Latitude']]
#    coordinates.index = IIASA_coord['Community Name']
#    start_index, end_index = None, None
#    for name_key in coordinates.index.values[start_index : end_index]:
#        lon,lat = coordinates.loc[name_key]['Longitude'], coordinates.loc[name_key]['Latitude']
#        targstart,targend = '2020-01-01','2021-01-01'
#        df = make_s1_ts(([lon, lat],3), 'buffer',(targstart,targend), int(10), retry_times = 100)
#


#    gauges_df = pd.read_csv('/Users/Mitchell/Documents/LDEOinternship/gauge_paper_correlation_names.csv')
#    coordinates = gauges_df[['Longitude','Latitude']]
#    coordinates.index = gauges_df['gauge_name']
#    start_index, end_index = 30,None
#    for name_key in coordinates.index.values[start_index : end_index]:
#        lon,lat = coordinates.loc[name_key]['Longitude'], coordinates.loc[name_key]['Latitude']
#        targstart,targend = '2017-01-01','2020-12-31'
#        df = make_s1_ts(([lon, lat],15), 'buffer',(targstart,targend), int(10), retry_times = 100)
    # select_gauge = 'SW228.5'
    # gauges_df = pd.read_csv('gauge_paper_correlation_names.csv')
    # coordinates = gauges_df[['Longitude','Latitude']]
    # coordinates.index = gauges_df['gauge_name']
    # select_coordinates = coordinates.loc[select_gauge]
    #
    # lon,lat = select_coordinates['Longitude'], select_coordinates['Latitude']
    # targstart,targend = '2020-07-01','2020-10-31'
    # df = make_s1_ts(([lon, lat],15), 'buffer',(targstart,targend), int(100), retry_times = 100)



#    siraj_wide_roi = ee.FeatureCollection('users/mlt2177/sirajRoiGeoms').filterMetadata('NAME','equals','WideTestPixel') \
#                            .first().geometry()
#    df = make_s1_ts(siraj_wide_roi, 'SirajWidePixel', ('2019-09-10','2020-11-01'), int(20), retry_times = 100, include_pow = True)
