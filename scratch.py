import network_exploration_stuff as mynet
import logging
import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
import geopandas as gpd
import rasterio as rio
import numpy as np
import random
import pandas as pd
pd.set_option('display.max_rows', None)

# CAPUTRE THE SHAPELY DEPRECIAITION WARNINGS
logging.captureWarnings(True)

# VARIABLES USED #############################################################
# file path
path = "/home/mdp0023/Documents/Codes_Projects/\
network_analysis/Network_Testing_Data"
image_path = "/home/mdp0023/Documents/Codes_Projects/\
network_analysis/Poster_Graphics"
inset_path = "/home/mdp0023/Documents/Codes_Projects/network_analysis/bboxes"
# AOI without buffer
aoi_area = f'{path}/Neighborhood_Network_AOI.shp'
# AOI with buffer
aoi_buffer = f'{path}/Neighborhood_Network_AOI_Buf_1km.shp'
# centroids of res parcels and food marts
res_points_loc = f'{path}/Residential_Parcels_Points_Network_AOI.shp'
food_points_loc = f'{path}/Food_Marts_Points_Network_AOI.shp'
# shapefiles of res parcels and food marts
res_parcels = f'{path}/Residential_Parcels_Network_AOI.shp'
food_parcels = f'{path}/Food_Marts_Network_AOI.shp'
# path for inundation
inundation = f'{path}/Network_Inun.tif'
raster = rio.open(inundation)

# LOADING WORK ###############################################################
#G = mynet.shape_2_graph(source=aoi_buffer)
G = mynet.read_graph_from_disk(
    path='/home/mdp0023/Documents/Codes_Projects/network_analysis/Network_Testing_Data', name='AOI_Graph_Inundated')

G = mynet.read_graph_from_disk(
    path='/home/mdp0023/Documents/Codes_Projects/network_analysis/Network_Testing_Data', name='AOI_Graph')

G = mynet.rename(G=G)


# LOAD OTHER DATA ############################################################
# shapefile centroids of residental plots
res_points = gpd.read_file(res_points_loc)
# shapefile of res parcels
res_locs = gpd.read_file(res_parcels)
# shapefile centroids of 3 foodmart plots
food_points = gpd.read_file(food_points_loc)
# shapefile of food mart parcels
food_locs = gpd.read_file(food_parcels)
# shapefile of area of interest
aoi_area = gpd.read_file(aoi_area)

# ensure proper projection
G = ox.projection.project_graph(G, to_crs=32614)
res_locs = res_locs.to_crs(epsg=32614)
food_locs = food_locs.to_crs(epsg=32614)
res_points = res_points.to_crs(epsg=32614)
food_points = food_points.to_crs(epsg=32614)
aoi_area = aoi_area.to_crs(epsg=32614)


# decomp, sink_insights, res_locs, dest_locs = mynet.flow_decomposition(G=G, 
#                                     res_points=res_points, 
#                                     dest_points=food_points, 
#                                     res_locs=res_locs,
#                                     dest_locs=food_locs,
#                                     G_demand='inundation_demand', 
#                                     G_capacity='inundation_capacity',
#                                      G_weight='inundation_travel_time')

decomp, sink_insights, res_locs, dest_locs = mynet.flow_decomposition(G=G,
                                            res_points=res_points,
                                            dest_points=food_points,
                                            res_locs=res_locs,
                                            dest_locs=food_locs,
                                            G_demand='demand',
                                            G_capacity='capacity',
                                            G_weight='travel_time')
                                     




res_locs.plot(column='cost_of_flow', legend=True)
plt.show()
