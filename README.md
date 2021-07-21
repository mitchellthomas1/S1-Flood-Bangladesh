# S1-Flood-Bangladesh
Data creation tools for running Sentinel-1 flood algorithm for index insurnace applications described in Thomas et al. (In Review).

An account on the GEE is required to use S1-Flood-Bangladesh. To sign up for an account, go to https://earthengine.google.com.


## Instructions
S1-Flood-Bangladesh offers several tools to implement the algorithm for Sentinel-1 based flood mapping described in Thomas et al. (In Review). 
The repo allows a user to create and export a single flood map or a flooded area time series using the algorithm.

### Make a flood map
The script createmap.py houses tools to create a single flood map. The function make_map_image() in this script creates and exports this flood map to Google Drive. 

make_map_image() takes an ROI (ee.Geometry object), a start (inclusive) and end (exclusive) date for the study period (latest pixel is taken within this range),
a scale for exporting (10 is recommended but a larger value will export more quickly), and a folder to export to in Google Drive.

The script makemapexample shows this function being used over a flood in Sirajganj district, 2019, in Bangladesh.

### Make a flood time series
The goal of this algorithm is ultimately to provide a reliable and consistent surface water time series. 

The script createtimeseries.py houses the code to make this time series. Due to Earth Engine limitations, this script is designed
to run mostly client side and make many successive calls to the server side. Generally this goes aginst Earth Engine recommendations,
but to make hundreds of flood maps at 10 meter resolution we have found this method to be more reliable.

The function make_s1_ts() is housed in this script and creates and saves the time series. It can either run the time series over an adminstrative district in Bangladesh,
a buffer zone around a certain point, or a custom EE geometry. An example in an administrative district is given in runadminunit.py and an example on a custom geometry is
shown in runcustomgeometry.py



## Examples
Example codes are mentioned above but summarized again here. This codes can be edited based on your area of interest to create data most easily.

1. Make and save a flood map over a custom EE geometry: makemapexample.py

2. Make a S1 flood time series over an adminstrative unit in Bangladesh: runadminunit.py

3. Make a S1 flood time series over a custom EE geometry: runcustomgeometry.py




