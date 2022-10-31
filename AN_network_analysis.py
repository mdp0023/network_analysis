# this is the code where I will do scratch work for the Austin North Network

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
folder_path = '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North'
# bounadry and boundary with buffer
aoi_area = f'{folder_path}/AN_Boundary/AN_Boundary.shp'
aoi_buffer = f'{folder_path}/AN_Boundary/AN_Boundary_3km.shp'
# Centroids of res parcels and grocery stores
res_parcel_loc = f'{folder_path}/AN_Residential_Parcel_Shapefiles/AN_Residential_Parcels.shp'
res_points_loc = f'{folder_path}/AN_Residential_Parcel_Shapefiles/AN_Residential_Parcels_Points.shp'
food_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_Supermarket_Parcel.shp'
food_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_Supermarket_Parcel_Points.shp'

res_parcels = gpd.read_file(res_parcel_loc)
res_points = gpd.read_file(res_points_loc)
food_parcels = gpd.read_file(food_parcel_loc)
food_points = gpd.read_file(food_points_loc)
# Inundation raster
# inundation_raster = f'{folder_path}/AN_Inundation/AN_Memorial_Day_Compound.tif'
# raster = rio.open(inundation_raster)
############################################################################################

# # download graph
# G = mynet.shape_2_graph(aoi_buffer)
# G = mynet.rename(G=G)

# # save graph
# mynet.save_2_disk(G=G, path=f'{folder_path}/AN_Graphs', name='AN_Graph')

# load saved graph
G = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph')

G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels = mynet.nearest_nodes_vertices(G=G,
                                                                                                                                       res_points=res_points, dest_parcels=food_parcels, G_demand='demand')



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

edge_attributes = list(list(G.edges(data=True))[0][-1].keys())
print(edge_attributes)
print(f'number of nodes: {len(G.nodes())}')
print(f'number of edges: {len(G.edges())}')

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

# # need to convert back to multidigraph to plot properly
# G_map = nx.MultiDiGraph(G_map)

# flow decomposition test
output = mynet.flow_decomposition(G=G, 
                                        res_points=res_points,
                                        dest_points=food_points, 
                                        res_parcels=res_parcels, 
                                        dest_parcels=food_parcels, 
                                        G_demand='demand', 
                                        G_capacity='capacity', 
                                        G_weight='travel_time',
                                        dest_method='multiple')


decomposed_paths=output[0]
sink_insights=output[1]
res_parcels=output[2]
dest_parcels=output[3]
print(dest_parcels)
dest_parcels.to_file('/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/test_decomp')
res_parcels.plot(column='cost_of_flow', legend=True)

# mynet.plot_aoi(G=G_map, res_parcels=res_parcels,
#                resource_parcels=food_parcels, edge_width='dry_flow')
plt.show()


