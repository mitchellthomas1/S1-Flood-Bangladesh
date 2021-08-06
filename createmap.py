

import ee
import zscore
from baselinecreation import calc_single_baseline



### ------TEST PARAMETERS-------###
zvv_thd = -2.0
zvh_thd = -2.0
rawvh_thd = -22
smoothing_radius = 30 #(meters)
smoothing_threshold = 0.6



def make_collections(geom, targstart, targend, basestart, baseend):

     # Load S1 collection with filters
    filters = [
        ee.Filter.listContains("transmitterReceiverPolarisation", "VV"),
        ee.Filter.listContains("transmitterReceiverPolarisation", "VH"),
        ee.Filter.equals("instrumentMode", "IW"),
    ]
    s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filter(filters).filterBounds(geom)

    if int(s1.filterDate(targstart, targend).size().getInfo()) == 0:
        raise('''Error: no dates within the chosen study period. 
        Pick a new start (inclusive) and end (exclusive) date''')
    

    # function to correct border inaccuracies
    def correctBorders(image):
        maskedImage = (image.updateMask(image.select('VV').gt(-45)))
        return maskedImage

        # create and merge baseline and target period stacks
    s1_combined = (s1.filterDate(basestart,baseend)
                .merge(s1.filterDate(targstart, targend))
                .map(correctBorders)).sort('system:time_start')
    #    create list of dates and times where images occur

    # filter s1 by ascending and descending
    s1_asc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'ASCENDING'))
    s1_dsc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'DESCENDING'))
    # Compute Z-scores per orbital direction and merge into zscore
    z_iwasc = zscore.calc_zscore(s1_asc, basestart, baseend, 'IW', 'ASCENDING')
    z_iwdsc = zscore.calc_zscore(s1_dsc, basestart, baseend, 'IW', 'DESCENDING')

    z_collection = ee.ImageCollection(z_iwasc.merge(z_iwdsc)).sort('system:time_start')

    return z_collection, s1_combined

def classify_image(z_image, s1_image):
    '''
    Function to actually classify flooded area based on thresholds set at top of script

    Parameters:
    z_image (ee.Image): Image with z-score for flood day

    s1_image (ee.Image): Image with raw S1 Values


    Returns: a flood image
    '''

    vvflag = z_image.select('VV').lte(zvv_thd)
    vhflag = z_image.select('VH').lte(zvh_thd)

    rawvhflag = s1_image.select('VH').lte(rawvh_thd)

    flood_class = (vvflag).add(vhflag).add(rawvhflag).rename(['flood_sum'])
    bool_map = flood_class.gte(1)

    smooth = bool_map.focal_mean(radius = smoothing_radius, units = 'meters').gt(smoothing_threshold)

    return smooth


def make_map_image(roi, start, end, scale, export_folder , descrip_str = '', baseline_len = 45, clip = False):

    '''

    Parameters:
    roi (ee.Geometry): roi over which to make map and export tiff

    start, end (str, ee.Date): start(inclusive) and end(exclusive) of study period. This range will be mosaicked so
         that the latest pixel in the period is chosen

    export_folder (str): folder in Google Drive in which to save tif

    
    '''

    # calculate baseline
    basestart, baseend  = calc_single_baseline(roi, baseline_len)

    #make z and raw s1 collections
    z_collection, s1_collection = make_collections(roi, start, end, basestart, baseend)

    # make z and s1 mosaicks
    z_image, s1_image = z_collection.filterDate(start, end).mosaic(), s1_collection.filterDate(start,end).mosaic()

    # classify a flood image
    flood_image = classify_image(z_image, s1_image)

#     Export to drive
    task = ee.batch.Export.image.toDrive(image=flood_image,  # an ee.Image object.
                                     region=roi,  # an ee.Geometry object.
                                     description= 'S1FloodMap' + descrip_str,
                                     folder= export_folder,
                                     fileNamePrefix= 'S1FloodMap' + descrip_str,
                                     scale=scale,
                                     crs='EPSG:4326')
    
    # Export to Asset
#    task = ee.batch.Export.image.toAsset(image=flood_image,  # an ee.Image object.
#                                     region=roi,  # an ee.Geometry object.
#                                     description= 'S1FloodMap' + descrip_str,
#                                     assetId= 'users/mlt2177/S1FloodMap' + descrip_str,
#                                     scale=scale,
#                                     crs='EPSG:4326')
        

    task.start()

    return task

       

