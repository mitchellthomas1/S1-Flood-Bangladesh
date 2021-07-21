import ee
from createtimeseries import make_s1_ts

# Any ee geometry object. Here an example is shown over the Brahmaputra river near Sirajganj District
custom_roi =  ee.Geometry.Polygon(
        [[[89.60290221563474, 24.5689522895152],
          [89.60290221563474, 24.278860653015048],
          [89.90777282110349, 24.278860653015048],
          [89.90777282110349, 24.5689522895152]]], None, False)
# A name/description for the area
admin_level = 'sirajganj_custom_roi'

start = '2017-01-01'
end = '2021-01-01'

#reduction scale  (meters)
scale = 100

flood_df = make_s1_ts(custom_roi, description, (start,end), scale)