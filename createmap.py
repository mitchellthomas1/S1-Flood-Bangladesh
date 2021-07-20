

import ee
import zscore

def make_map_image(roi, start, end, baseline_len = 45, clip = False):
    '''

    Parameters:
    roi (ee.Geometry): roi over which to make map

    
    '''

        # calculate baseline
        basestart, baseend  = calc_single_baseline(roi, baseline_len)


        # Load S1 collection with filters
        filters = [
          ee.Filter.listContains("transmitterReceiverPolarisation", "VV"),
          ee.Filter.listContains("transmitterReceiverPolarisation", "VH"),
          ee.Filter.equals("instrumentMode", "IW"),
        ]
        s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filter(filters).filterBounds(study_feature.geometry())
        
        if int(s1.filterDate(start, end).size().getInfo()) == 0:
            raise('''Error: no dates within the chosen study period. 
            Pick a new start (inclusive) and end (exclusive) date''')
        

        # function to correct border inaccuracies
        def correctBorders(image):
          maskedImage = (image.updateMask(image.select('VV').gt(-45)))
          return maskedImage

         # create and merge baseline and target period stacks
        s1_combined = (s1.filterDate(basestart,baseend)
                    .merge(s1.filterDate(start, end))
                    .map(correctBorders))
        #    create list of dates and times where images occur

        # filter s1 by ascending and descending
        s1_asc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'ASCENDING'))
        s1_dsc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'DESCENDING'))
        # Compute Z-scores per orbital direction and merge into zscore
        z_iwasc = zscore.calc_zscore(s1_asc, basestart, baseend, 'IW', 'ASCENDING')
        z_iwdsc = zscore.calc_zscore(s1_dsc, basestart, baseend, 'IW', 'DESCENDING')
        z_coll = ee.ImageCollection(z_iwasc.merge(z_iwdsc)).sort('system:time_start')

