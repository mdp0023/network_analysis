# this is the code where I will do scratch work for the Austin North Network

# Packages
import network_exploration_stuff as mynet
import rasterio as rio
import networkx as nx
import pandas as pd
import numpy as np
import logging

# ignore shapely deprectiation warnings
logging.captureWarnings(True)

# AOI VARIABLES
# Folder path
folder_path = '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North'
# bounadry and boundary with buffer
aoi_area = f'{folder_path}/AN_Boundary/AN_Boundary.shp'
aoi_buffer = f'{folder_path}/AN_Boundary/AN_Boundary_3km.shp'
# Centroids of res parcels and grocery stores
res_points_loc = f'{folder_path}/Austin_North/AN_Residential_Parcel_Shapefiles/AN_Residential_Parcels_Points.shp.shp'
food_points_loc = f'{folder_path}/N_Resource_Parcel_Shapefiles/AN_Supermarket_Parcel_Points.shp.shp'
# Inundation raster
# inundation_raster = f'{folder_path}/AN_Inundation/AN_Memorial_Day_Compound.tif'
# raster = rio.open(inundation_raster)

G = mynet.shape_2_graph(aoi_buffer)

# available edge and node attributes
edge_attributes = list(list(G.edges(data=True))[0][-1].keys())
node_attributes = list(list(G.nodes(data=True))[0][-1].keys())
print(edge_attributes)

# get unique speed values across network (check for major errors)
unique_speeds = np.unique(list(nx.get_edge_attributes(G, 'speed_kph').values()))
percent_edge_w_speed = len(list(nx.get_edge_attributes(G, 'speed_kph').values()))/G.number_of_edges()

# get list of unique road classifications
unique_rtypes = list(set(list(nx.get_edge_attributes(G, 'highway').values())))
percent_edge_w_rtype = len(list(nx.get_edge_attributes(G, 'highway').values()))/G.number_of_edges()
print(unique_rtypes)

# get unique capacity values (check for errors)
unique_caps = np.unique(list(nx.get_edge_attributes(G, 'capacity').values()))
percent_edge_w_cap = len(list(nx.get_edge_attributes(G, 'capacity').values()))/G.number_of_edges()
print(unique_caps)

# get unique number of lanes (check for errors)
num_lanes = np.unique(list(nx.get_edge_attributes(G, 'lanes').values()))
percent_edge_w_lane = len(list(nx.get_edge_attributes(G, 'lanes').values()))/G.number_of_edges()
print(num_lanes)

# see percentage of edges that have respected stats
print(f"percentage edges with lane info : {percent_edge_w_lane}")
print(f"percentage edges with speed info: {percent_edge_w_speed}")
print(f"percentage edges with rtype info: {percent_edge_w_rtype}")
print(f"percentage edges with cap info  : {percent_edge_w_cap}")

