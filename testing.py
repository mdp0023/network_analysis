# CREATED ON JANUARY 25TH TO TEST FUNCTIONS AS I FIX BUGS, ORGANIZE DOCUMENTATION, ETC. 
# MAJOR CODE CLEANING, MAKING SURE EVERYTHING IS DOCUMENTED AND WORKING PROPERLY BEFORE MAJOR TA BUG FIXXES

# setting variables
import network_exploration_stuff as mynet
import logging
import sys
import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
import random
import pandas as pd
import rasterio
pd.set_option('display.max_rows', None)

# CAPUTRE THE SHAPELY DEPRECIAITION WARNINGS
logging.captureWarnings(True)

# don't truncate output of numpy array
np.set_printoptions(threshold=sys.maxsize)

# for testing purposes, should set use_cahce=False
ox.config(use_cache=True)

# AOI VARIABLES
# Folder paths
f_path = '/home/mdp0023/Desktop/external/Data/Network_Data/AOI_Testing'
f_boundaries = f'{f_path}/AOI_Boundary'
f_food_shp = f'{f_path}/AOI_Food_Mart_Shapefiles'
f_graphs = f'{f_path}/AOI_Graphs'
f_parcel_access_shp = f'{f_path}/AOI_Parcel_Access_Shapefiles'
f_res_shp = f'{f_path}/AOI_Residental_Parcel_Shapefiles'

#file paths
# boundary
aoi_area = f'{f_boundaries}/Neighborhood_Network_AOI.shp'
# boundary with buffer
aoi_buffer = f'{f_boundaries}/Neighborhood_Network_AOI_Buf_1km.shp'
# Centroids of res parcels and food marts
res_points_loc = f'{f_res_shp}/Residential_Parcels_Points_Network_AOI.shp'
food_points_loc = f'{f_food_shp}/Food_Marts_Points_Network_AOI.shp'
# Shapefiles of res parcels and food marts
res_parcels = f'{f_res_shp}/Residential_Parcels_Network_AOI.shp'
food_parcels = f'{f_food_shp}/Food_Marts_Network_AOI.shp'
inundation_raster = f'{f_path}/AOI_Inundation.tif'
# raster = rio.open(inundation_raster)

# do crs_check
# EPSG:32614 -> WGS84 UTM zone 14N projection
# EPSG:26914 -> NAD83 UTM zone 14N projection
# EPSG:4326  -> WGS84 
mynet.crs_check(source=aoi_area, shapefiles=[
                aoi_buffer, res_points_loc, food_points_loc, res_parcels, food_parcels], crs=26914)


# read in files 
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


# FUNCTION TESTING

# shape_2_graph
network = mynet.shape_2_graph(aoi_buffer)

# save_2_disk
mynet.save_2_disk(G=network, path=f_graphs, name='AOI_Graph')

# read_graph_from_disk
network = mynet.read_graph_from_disk(path=f_graphs, name='AOI_Graph')

# print list of all the edge attributes available to us
print(list(list(network.edges(data=True))[0][-1].keys()))

# parallel_edges
parallel_edges, self_loop_edges = mynet.parallel_edges(network)

# nearest_nodes
output1 = mynet.nearest_nodes(G=network, res_points=res_points, dest_points=food_points)

# nearest_nodes_vertices
output2 = mynet.nearest_nodes_vertices(network, res_points, food_locs)

# random_shortest_path
paths = mynet.random_shortest_path(network, res_points, food_points, plot=False)

# min_cost_flow_parcels
output3 = mynet.min_cost_flow_parcels(network, res_points, food_points, food_locs, dest_method='multiple')

# max_flow_parcels
output4 = mynet.max_flow_parcels(network, res_points, food_points, food_locs, dest_method='multiple')

# plot_aoi -> lots can be added to this function
mynet.plot_aoi(network, res_locs, food_locs)

# summary_function -> lots can be added to this function in the future
output5 = mynet.summary_function(network)

# inundate_network
# inundated_net = mynet.inundate_network(network, f_graphs, inundation_raster)
# print(list(list(inundated_net.edges(data=True))[0][-1].keys()))

# flow_decomposition
output6 = mynet.flow_decomposition(network, res_points, food_points, res_locs, food_locs, dest_method='multiple')


