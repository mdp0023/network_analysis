# this script contains functions to collect resource data for AOI
import logging
import warnings
import geopandas as gpd
import pandas as pd
import contextily as cx
import numpy as np
import matplotlib.pyplot as plt
import osmnx as ox
import plotly.express as px

# CAPUTRE THE SHAPELY DEPRECIAITION WARNINGS
logging.captureWarnings(True)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
######################################
# TODOS
# TODO: Ideally I will need a building footprints, center points, and parcels for each resource 
# TODO: For now: using place function -> adjust to taking shapefile
# TODO: Examine hospital layer some more - which facilities should we be pulling?
# TODO: Might be able to pull all footprint data from OSM, but unsure
# TODO: EMS Stations (ambulance stations) are severily lacking
# TODO: Gas stations/convenient stores - lots of cross over? Might need to do some investigating here
# TODO: Should I (how should I) consider a location where multiple resources can be accessed at 


######################################
# NOTES
# NOTE: Shapely is set to change dramatically so might require edits in the future
# NOTE: Because OSM layer might have points or footprints, I am going to standardize and just keep centroid of each location
# NOTE: WRT hospitals, documentaiton of settings is going to be KEY - what gets included, what doesn'y, why/why not, etc. 
#   FURTHERMORE: Lots of double checking extracted layers need to be done before running algorithms

########################################################################
# VARIABLES
place = 'Austin, Texas, USA'
aoi_crs = 32614

########################################################################
# extract and reproject supermarket locations
tags = {'shop':'supermarket'}
supermarkets = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
supermarket_points = supermarkets.set_geometry(supermarkets.centroid)
supermarket_points['Resource'] = 'Supermarket'

# Extract and reproject convenience store locations
tags = {'shop':'convenience'}
convenience = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
convenience_points = convenience.set_geometry(convenience.centroid)
convenience_points['Resource'] = 'Convenience Store'

# Extract and convert hospital locations
# has information regarding if their is an emergency room, and potentially the number of beds
# lots of settings that can be pulled - hospital grounds, individual buildings across a campus, etc.
# To avoid future depreciation issues, going to use healthcare key
# only extracting those with emergency rooms
tags = {'healthcare': 'hospital'}
hospitals = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
hospitals = hospitals[hospitals['emergency']=='yes']
hospitals_points = hospitals.set_geometry(hospitals.centroid)
hospitals_points['Resource'] = 'Hospitals (ERs)'

# Extract and reproject pharmacy locations
# can potentially make a distinction between one that can/cannot sell prescription drugs - but mostly missing
tags={'amenity':'pharmacy'}
pharmacy = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
pharmacy_points = pharmacy.set_geometry(pharmacy.centroid)
pharmacy_points['Resource'] = 'Pharmacies'

# Extract and reproject fire station locations
tags = {'amenity':'fire_station'}
fire_station = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
fire_station_points = fire_station.set_geometry(fire_station.centroid)
fire_station_points['Resource'] = 'Fire Stations'

# Extract and reproject ambulance station locations
tags = {'emergency': 'ambulance_station'}
ambulance_station = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
ambulance_station_points = ambulance_station.set_geometry(
    ambulance_station.centroid)
ambulance_station_points['Resource'] = 'Ambulance Stations'

# Extract and reproject police station locations
tags={'amenity':'police'}
police_station = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
police_station_points = police_station.set_geometry(police_station.centroid)
police_station_points['Resource'] = 'Police Stations'

# Extract and reproject gas station locations
tags = {'amenity':'fuel'}
gas_station = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
gas_station_points = gas_station.set_geometry(gas_station.centroid)
gas_station_points['Resource'] = 'Gas Stations'

# Merge data into singular dataframe
resource_points = pd.concat(
    [supermarket_points, convenience_points, hospitals_points, 
    pharmacy_points, fire_station_points, ambulance_station_points, 
    police_station_points, gas_station_points])
#####################################################
# Plot
fig, ax = plt.subplots(figsize=[12,8])
ax.axis('off')
supermarket_points.plot('Resource', categorical=True, ax=ax, legend=True, 
                    legend_kwds={'loc': 'center left', 'bbox_to_anchor': (1, 0.5)})
plt.show()



# # I don't think I like plotly and will be using folium for interactive mapping
# # Create interactive plotly plot
# fig = px.scatter_mapbox(resource_points, 
#                         lon=resource_points.to_crs(epsg=4326).geometry.x, 
#                         lat=resource_points.to_crs(epsg=4326).geometry.y, 
#                         color='resource')
# fig.update_layout(mapbox_style="open-street-map")
# fig.update_layout(autosize=True)
# fig.update_geos(fitbounds="locations")

