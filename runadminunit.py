import ee
from createtimeseries import make_s1_ts

# name of adminstrative unit in Bangladesh
admin_name  = 'Jamalpur'
# admin level
admin_level = 'district'

start = '2017-01-01'
end = '2021-01-01'

#reduction scale  (meters)
scale = 100

flood_df = make_s1_ts(admin_name, admin_level,(start,end), scale)