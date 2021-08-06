import ee
ee.Initialize()

from createmap import make_map_image

#ROI being studied
sirajganj_roi = ee.Geometry.Polygon(
        [[[89.61086842059495, 24.519158572887243],
          [89.61086842059495, 24.30344374797097],
          [89.78596302508714, 24.30344374797097],
          [89.78596302508714, 24.519158572887243]]], None, False)


# Start (inclusive) and end (exclusive) of study period. Latest pixel in range is taken over ROI
start = '2019-07-12'
end   =	'2019-07-22'

#scale of exporting. Recommend: 10. Higher values will be faster.
scale = 10

#Description for filename and task name
descrip_str = 'SirajganjJuly2019'

#Folder in GEE to export to
export_folder = 'FloodMapping'

 
task = make_map_image(sirajganj_roi, start, end, scale, export_folder )

# print status until completion
status = task.status()['state']
print('Task: ', status)
while status != 'COMPLETED':
    if task.status()['state'] != status:
        status = task.status()['state']
        print('Task: ',status)
