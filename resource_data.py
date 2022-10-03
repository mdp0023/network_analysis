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
# TODO: Ideally I will need a [building footprints], center points, and parcels for each resource 
# TODO: Examine hospital layer some more - which facilities should we be pulling?
# TODO: Might be able to pull all footprint data from OSM, but unsure
# TODO: EMS Stations (ambulance stations) are severily lacking, at least in the Austin area
# TODO: Gas stations/convenient stores - lots of cross over? Might need to do some investigating here
# TODO: Should I (how should I) consider a location where multiple resources can be accessed at 
# TODO: explore if resources have other attributes I should consider
    # e.g., hospitals has attribute of if ER or not, might want to consider including that


######################################
# NOTES
# NOTE: Shapely is set to change dramatically so might require edits in the future
# NOTE: Because OSM layer might have points or footprints, I am going to standardize and just keep centroid of each location
# NOTE: WRT hospitals, documentaiton of settings is going to be KEY - what gets included, what doesn't, why/why not, etc. 
#   FURTHERMORE: Lots of double checking extracted layers need to be done before running algorithms

########################################################################

########################################################################
def resource_centroids(boundary, resources = 'default', custom_r=None, disused_filter = True, save=None, cache=False):
    """
    A function to pull resource centroids Open Street Map.
    
    Utilizes functions already created in OSMNx. This function should streamline the data pulling process ny combining the resource that have identified as being critical.

    REMEMBER: Input boundary shapefile should be projected, and code should correctly manage projections.

    :param boundary: Filepath to .shp file of boundary of interest. Should have projected coordinates of study area
    :type boundary: string
    :param resources: 'default' or 'custom'. Default will use default list of critical resources, custom will be user input list of dictionaries
    :type resources: string 
    :param custom_r: list of key-value pairs refering to what resource sshould be extracted from openstreetmap
    :type custom_r: list of key-value pairs
    :param disused_filter: If True, will filter out resources that have an input for disused column, signifiying that it is most likely closed
    :type disused_filter: bool
    :param save: If not None, save to file path
    :type save: string
    :param cache: sets OSMNx cache setting - if cache folder exists, and need new data, must be set to False. Default is False, set to True if repeatedly pulling same information, it will increase speed.
    :type cache: bool
    
    :return: **resource_points**, geopandas shapefile of centroids of resource locations within boundary
    :rtype: geopandas shapefile

    """

    # set cache to setting. Set to True to save time, set to False if need new data and folder exists
    ox.settings.use_cache = cache

    #open boundary file, save crs, and reproject to EPSG:4326 (WGS84)- required for accurate pull from openstreetmap
    geo = gpd.read_file(boundary)
    crs_proj = geo.crs
    geo=geo.to_crs("EPSG:4326")
    
    # determine list of key value pairs
    if resources == 'default':
        resource_list = [{'shop': 'supermarket'}, 
                         {'amenity': 'fuel'}, {'shop': 'convenience'},
                         {'healthcare': 'hospital'}, {'amenity':'pharmacy'},
                         {'amenity': 'fire_station'}, {'emergency': 'ambulance_station'},
                         {'amenity':'police'}]

    elif resources == 'custom':
        resource_list = custom_r

    # create emtpy dataframe to store data as it is pulled
    resource_points=pd.DataFrame()

    # iterate through key value pairs
    for key_value in resource_list:
        data = ox.geometries_from_polygon(geo['geometry'][0], tags=key_value).to_crs(crs_proj)
        data_points = data.set_geometry(data.centroid)
        data_points['Resource'] = list(key_value.values())[0]
    
    # Merge data into singular dataframe
        resource_points = pd.concat([resource_points,data_points])

    # drop columns that has a list in - for some reason throwing an error
    for col in resource_points.columns:
        if any(isinstance(val, list) for val in resource_points[col]):
            #print(f'Column: {col}, has a list in it')
            resource_points = resource_points.drop(labels=col, axis='columns')

    # filter out disused locations
    if disused_filter is True:
        # find columns in output that begin with disused 
        # Do this incase it is not disused_sh (shop) like for supermarkets
        filter_col = [col for col in resource_points if col.startswith('disused')]
        for col in filter_col:
            resource_points = resource_points[resource_points[col].isnull()]

    # save file if necessary
    if save is not None:
        resource_points.to_file(save)

    return resource_points


# VARIABLES
place = 'Austin, Texas, USA'
aoi_crs = 32614
folder = '/home/mdp0023/Desktop/external/Data/Network_Data'
SC_boundary_file = f'{folder}/Shoal_Creek/SC_Boundary/SC_Boundary.shp'
AN_boundary_file = f'{folder}/Austin_North/AN_Boundary/AN_Boundary.shp'
BPA_boundary_file = f'{folder}/Beaumont_Port_Arthur/BPA_Boundary/BPA_Boundary.shp'

SC_output = f'{folder}/Shoal_Creek/SC_Resource_Parcel_Shapefiles/SC_Supermarket_Parcel_Points.shp'
AN_output = f'{folder}/Austin_North/AN_Resource_Parcel_Shapefiles/AN_Supermarket_Parcel_Points.shp'
BPA_output = f'{folder}/Beaumont_Port_Arthur/BPA_Resource_Parcel_Shapefiles/BPA_Supermarket_Parcel_Points.shp'

custom_r = [{'shop': 'supermarket'}]

resource_points = resource_centroids(boundary = BPA_boundary_file,
                                    resources = 'custom',
                                    custom_r = custom_r,
                                    disused_filter = True,
                                     save=BPA_output)


####################################################################################
#PLOT

# fig, ax = plt.subplots(figsize=[12, 8])
# ax.axis('off')
# resource_points.plot(column='Resource', categorical=True, ax=ax, legend=True,
#                      legend_kwds={'loc': 'center left', 'bbox_to_anchor': (1, 0.5)})
# SC_boundary.boundary.plot(ax=ax)
# plt.show()


# # has information regarding if their is an emergency room, and potentially the number of beds
# # only extracting those with emergency rooms
# tags = {'healthcare': 'hospital'}
# hospitals = ox.geometries_from_place(place, tags=tags).to_crs(epsg=aoi_crs)
# hospitals = hospitals[hospitals['emergency']=='yes']
# hospitals_points = hospitals.set_geometry(hospitals.centroid)
# hospitals_points['Resource'] = 'Hospitals (ERs)'


# # I don't think I like plotly and will be using folium for interactive mapping
# # Create interactive plotly plot
# fig = px.scatter_mapbox(resource_points, 
#                         lon=resource_points.to_crs(epsg=4326).geometry.x, 
#                         lat=resource_points.to_crs(epsg=4326).geometry.y, 
#                         color='resource')
# fig.update_layout(mapbox_style="open-street-map")
# fig.update_layout(autosize=True)
# fig.update_geos(fitbounds="locations")

