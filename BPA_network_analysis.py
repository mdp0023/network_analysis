# this is the code where I will do scratch work for the Beaumont Port Arthur Network

# Packages
import network_exploration_stuff as mynet
import matplotlib.pyplot as plt
import geopandas as gpd
import rasterio as rio
import networkx as nx
import pandas as pd
import numpy as np
import logging

# ignore shapely deprectiation warnings
logging.captureWarnings(True)

############################################################################################
# INPUT VARIABLES
# Folder path
folder_path = '/home/mdp0023/Desktop/external/Data/Network_Data/Beaumont_Port_Arthur'
# bounadry and boundary with buffer
aoi_area = f'{folder_path}/BPA_Boundary/BPA_Boundary.shp'
aoi_buffer = f'{folder_path}/BPA_Boundary/BPA_Boundary_3km.shp'
# Centroids of res parcels and grocery stores
res_parcel_loc = f'{folder_path}/BPA_Residential_Parcel_Shapefiles/BPA_Residential_Parcels.shp'
res_points_loc = f'{folder_path}/BPA_Residential_Parcel_Shapefiles/BPA_Residential_Parcels_Points.shp'
food_parcel_loc = f'{folder_path}/BPA_Resource_Parcel_Shapefiles/BPA_Supermarket_Parcel.shp'
food_points_loc = f'{folder_path}/BPA_Resource_Parcel_Shapefiles/BPA_Supermarket_Parcel_Points.shp'

res_parcels = gpd.read_file(res_parcel_loc)
res_points = gpd.read_file(res_points_loc)
food_parcels = gpd.read_file(food_parcel_loc)
food_points = gpd.read_file(food_points_loc)

############################################################################################

G = mynet.shape_2_graph(aoi_buffer)

# available edge and node attributes
edge_attributes = list(list(G.edges(data=True))[0][-1].keys())
node_attributes = list(list(G.nodes(data=True))[0][-1].keys())
print(edge_attributes)

# get unique speed values across network (check for major errors)
unique_speeds = np.unique(
    list(nx.get_edge_attributes(G, 'speed_kph').values()))
percent_edge_w_speed = len(list(nx.get_edge_attributes(
    G, 'speed_kph').values()))/G.number_of_edges()

# get list of unique road classifications
unique_rtypes = list(set(list(nx.get_edge_attributes(G, 'highway').values())))
percent_edge_w_rtype = len(list(nx.get_edge_attributes(
    G, 'highway').values()))/G.number_of_edges()
print(unique_rtypes)

# get unique capacity values (check for errors)
unique_caps = np.unique(list(nx.get_edge_attributes(G, 'capacity').values()))
percent_edge_w_cap = len(list(nx.get_edge_attributes(
    G, 'capacity').values()))/G.number_of_edges()
print(unique_caps)

# get unique number of lanes (check for errors)
num_lanes = np.unique(list(nx.get_edge_attributes(G, 'lanes').values()))
percent_edge_w_lane = len(
    list(nx.get_edge_attributes(G, 'lanes').values()))/G.number_of_edges()
print(num_lanes)

# see percentage of edges that have respected stats
print(f"percentage edges with lane info : {percent_edge_w_lane}")
print(f"percentage edges with speed info: {percent_edge_w_speed}")
print(f"percentage edges with rtype info: {percent_edge_w_rtype}")
print(f"percentage edges with cap info  : {percent_edge_w_cap}")

print(edge_attributes)
print(len(G.nodes()))
print(len(G.edges()))

mynet.plot_aoi(G=G, res_parcels=res_parcels, resource_parcels=food_parcels)
plt.show()
