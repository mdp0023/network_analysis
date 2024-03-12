# this is the code where I will do scratch work for the Beaumont Port Arthur Network

# Packages
import network_analysis_base as mynet
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
folder_path = '/media/mdp0023/extradrive1/Data/Network_Data/Beaumont_Port_Arthur'
# bounadry and boundary with buffer
aoi_area = f'{folder_path}/BPA_Boundary/BPA_Boundary.shp'
aoi_buffer = f'{folder_path}/BPA_Boundary/BPA_Boundary_3km.shp'
# Centroids of res parcels and grocery stores
res_parcel_loc = f'{folder_path}/BPA_Residential_Parcel_Shapefiles/BPA_Residential_Parcels.shp'
res_points_loc = f'{folder_path}/BPA_Residential_Parcel_Shapefiles/BPA_Residential_Parcels_Points.shp'
food_parcel_loc = f'{folder_path}/BPA_Resource_Parcel_Shapefiles/BPA_Supermarket_Parcel.shp'
food_points_loc = f'{folder_path}/BPA_Resource_Parcel_Shapefiles/BPA_Supermarket_Parcel_Points.shp'

res_parcels = gpd.read_file(res_parcel_loc)
# res_points = gpd.read_file(res_points_loc)
# food_parcels = gpd.read_file(food_parcel_loc)
# food_points = gpd.read_file(food_points_loc)

############################################################################################

G = mynet.shape_2_graph(aoi_buffer)
G = mynet.rename(G=G)

# # available edge and node attributes
# edge_attributes = list(list(G.edges(data=True))[0][-1].keys())
# node_attributes = list(list(G.nodes(data=True))[0][-1].keys())
# print(edge_attributes)

# # get unique speed values across network (check for major errors)
# unique_speeds = np.unique(
#     list(nx.get_edge_attributes(G, 'speed_kph').values()))
# percent_edge_w_speed = len(list(nx.get_edge_attributes(
#     G, 'speed_kph').values()))/G.number_of_edges()

# # get list of unique road classifications
# unique_rtypes = list(set(list(nx.get_edge_attributes(G, 'highway').values())))
# percent_edge_w_rtype = len(list(nx.get_edge_attributes(
#     G, 'highway').values()))/G.number_of_edges()
# print(unique_rtypes)

# # get unique capacity values (check for errors)
# unique_caps = np.unique(list(nx.get_edge_attributes(G, 'capacity').values()))
# percent_edge_w_cap = len(list(nx.get_edge_attributes(
#     G, 'capacity').values()))/G.number_of_edges()
# print(unique_caps)

# # get unique number of lanes (check for errors)
# num_lanes = np.unique(list(nx.get_edge_attributes(G, 'lanes').values()))
# percent_edge_w_lane = len(
#     list(nx.get_edge_attributes(G, 'lanes').values()))/G.number_of_edges()
# print(num_lanes)

# # see percentage of edges that have respected stats
# print(f"percentage edges with lane info : {percent_edge_w_lane}")
# print(f"percentage edges with speed info: {percent_edge_w_speed}")
# print(f"percentage edges with rtype info: {percent_edge_w_rtype}")
# print(f"percentage edges with cap info  : {percent_edge_w_cap}")

# print(edge_attributes)
# print(len(G.nodes()))
# print(len(G.edges()))


# output = mynet.nearest_nodes(G=G,
#                                 res_points=res_points,
#                                 dest_points=food_points,
#                                 G_demand='demand')
                                
# G = output[0] 
# unique_origin_nodes = output[1] 
# unique_dest_nodes = output[2] 
# positive_demand = output[3] 
# shared_nodes = output[4] 
# res_points = output[5] 
# dest_points = output[6]

# output = mynet.max_flow_parcels(G=G, 
#                                 res_points=res_points, 
#                                 dest_points=food_points,
#                                 G_capacity='capacity', 
#                                 G_weight='travel_time', 
#                                 G_demand='demand',
#                                 dest_method='multiple', 
#                                 dest_parcels=food_parcels, 
#                                 ignore_capacity=False)

# flow_dictionary = output[0]
# cost_of_flow = output[1]
# max_flow = output[2]
# access = output[3]

# print(f'level of access: {access}')
# print(f'Maximum amount of flow: {max_flow}')


# #relate flow dictionary back to DRY graph for plotting purposes
# G_map = nx.DiGraph(G)
# # relate
# for edge in G_map.edges:
#     values = {(edge[0], edge[1]): {'dry_flow':
#                                    flow_dictionary[edge[0]][edge[1]]}}
#     nx.set_edge_attributes(G=G_map, values=values)
# G_map = nx.MultiDiGraph(G_map)


mynet.plot_aoi(G=G, 
               res_parcels=res_parcels, 
               scalebar=True)
plt.show()
