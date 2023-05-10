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
import time
import cProfile
import pstats
from pstats import SortKey
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
f_other_shp2 = f'{f_path}/AOI_Other_Resource_Shapefiles'
f_other_shp = f'{f_path}/AOI_Other_Resource2_Shapefiles'
f_graphs = f'{f_path}/AOI_Graphs'
f_parcel_access_shp = f'{f_path}/AOI_Parcel_Access_Shapefiles'
f_res_shp = f'{f_path}/AOI_Residental_Parcel_Shapefiles'

#file paths
# boundary
aoi_area = f'{f_boundaries}/Neighborhood_Network_AOI.shp'
# boundary with buffer
aoi_buffer = f'{f_boundaries}/Neighborhood_Network_AOI_Buf_1km.shp'
# Centroids of res parcels, food marts, and 2nd resource
res_points_loc = f'{f_res_shp}/Residential_Parcels_Points_Network_AOI.shp'
food_points_loc = f'{f_food_shp}/Food_Marts_Points_Network_AOI.shp'
other_points_loc = f'{f_other_shp}/Other_Resource2_Points_Network_AOI.shp'
other_points_loc2 = f'{f_other_shp2}/Other_Resource_Points_Network_AOI_edited.shp'
# Shapefiles of res parcels, food marts, and 2nd resource
res_parcels = f'{f_res_shp}/Residential_Parcels_Network_AOI.shp'
food_parcels = f'{f_food_shp}/Food_Marts_Network_AOI.shp'
other_parcels = f'{f_other_shp}/Other_Resource2_Network_AOI.shp'
other_parcels2 = f'{f_other_shp2}/Other_Resource_Network_AOI_edited.shp'
# Inundation raster
inundation_raster = f'{f_path}/AOI_Inundation.tif'
# raster = rio.open(inundation_raster)

# do crs_check
# EPSG:32614 -> WGS84 UTM zone 14N projection
# EPSG:26914 -> NAD83 UTM zone 14N projection
# EPSG:4326  -> WGS84 
mynet.crs_check(source=aoi_area, shapefiles=[
                aoi_buffer, res_points_loc, food_points_loc, res_parcels, food_parcels], crs=26914)


# read in files 
# shapefile centroids of residental plots
res_points = gpd.read_file(res_points_loc)
# shapefile of res parcels
res_locs = gpd.read_file(res_parcels)
# shapefile centroids of 3 foodmart plots
food_points = gpd.read_file(food_points_loc)
# shapefile of food mart parcels
food_locs = gpd.read_file(food_parcels)
# shapefile centroids of other resource
other_points = gpd.read_file(other_points_loc)
# shapefile of other resource parcels
other_locs = gpd.read_file(other_parcels)
# shapefile of area of interest
aoi_area = gpd.read_file(aoi_area)

# shapefile centroids of other resource
other_points2 = gpd.read_file(other_points_loc2)
# shapefile of other resource parcels
other_locs2 = gpd.read_file(other_parcels2)

# FUNCTION TESTING

# shape_2_graph
# network = mynet.shape_2_graph(aoi_buffer)

# save_2_disk
# mynet.save_2_disk(G=network, path=f_graphs, name='AOI_Graph')

# read_graph_from_disk
network = mynet.read_graph_from_disk(path=f_graphs, name='AOI_Graph_Inundated')
# ox.io.save_graph_geopackage(network, filepath='/home/mdp0023/Desktop/external/Data/Network_Data/AOI_Testing/test')
# print list of all the edge attributes available to us
print(list(list(network.edges(data=True))[0][-1].keys()))

# # parallel_edges
# parallel_edges, self_loop_edges = mynet.parallel_edges(network)

# # nearest_nodes, single destination
# output1a = mynet.nearest_nodes(G=network, res_points=res_points, dest_points=food_points)
# # nearest_nodes, multiple destinations
# output1b = mynet.nearest_nodes(G=network, res_points=res_points, dest_points=[food_points,other_points])

# print(output1a[1])
# print(output1b[1])
# print(output1a[2])
# print(output1b[2])
# print(output1a[3])
# print(output1b[3])
# print(output1a[4])
# print(output1b[4])

# nearest_nodes_vertices, single destination
# output2a = mynet.nearest_nodes_vertices(network, res_points, food_locs, food_points)
# # nearest_nodes_vertices, multiple destinations
# output2b = mynet.nearest_nodes_vertices(network, res_points, [food_locs, other_locs], [food_points, other_points])


# # # # random_shortest_path
# # paths = mynet.random_shortest_path(network, res_points, food_points, plot=False)

# MIN_COST_FLOW_PARCELS
# single destination file, dest_method='single'
# output3a = mynet.min_cost_flow_parcels(network, res_points, food_points, food_locs, dest_method='single')
# # multiple destination files, dest_method='single'
output3b = mynet.min_cost_flow_parcels(network, res_points, [food_points, other_points], [food_locs, other_locs], dest_method='single')
# # single destination file, dest_method='multiple'
output3c = mynet.min_cost_flow_parcels(network, res_points, food_points, food_locs, dest_method='multiple')
# # # multiple destination files, dest_method='multiple'
output3d = mynet.min_cost_flow_parcels(network, res_points, [food_points, other_points], [food_locs, other_locs], dest_method='multiple')
# print(output3b[1])
# print(output3c[1])
# print(output3d[1])
# # # # relate output3a/b/c/d to network to plot
# for i in output3a[0]:
#     for j in output3a[0][i]:
#         nx.set_edge_attributes(network, {(i,j,0): {'test_flow1':output3a[0][i][j]}})
# # # for i in output3b[0]:
# # #     for j in output3b[0][i]:
# # #         nx.set_edge_attributes(network, {(i,j,0): {'test_flow2':output3b[0][i][j]}})
# # # for i in output3c[0]:
# # #     for j in output3c[0][i]:
# # #         nx.set_edge_attributes(network, {(i,j,0): {'test_flow3':output3c[0][i][j]}})
# # # for i in output3d[0]:
# # #     for j in output3d[0][i]:
# # #         nx.set_edge_attributes(network, {(i,j,0): {'test_flow4':output3d[0][i][j]}})      

#MAX_FLOW_PARCELS
# single destination file, dest_method='single'
# output4a = mynet.max_flow_parcels(network, res_points, food_points, food_locs, dest_method='single')
# # multiple destination files, dest_method='single'
# output4b = mynet.max_flow_parcels(network, res_points, [food_points, other_points], [food_locs, other_locs], dest_method='single')
# # single destination file, dest_method='mulitple'
# output4c = mynet.max_flow_parcels(network, res_points, food_points, food_locs, dest_method='multiple')
# # multiple destination files, dest_method='multiple'
# output4d = mynet.max_flow_parcels(network, res_points, [food_points, other_points], [food_locs, other_locs], dest_method='multiple')


# plot_aoi -> lots can be added to this function
# mynet.plot_aoi(network, res_locs, food_locs, edge_width='test_flow1')
# # mynet.plot_aoi(network, res_locs, [food_locs, other_locs], edge_width='test_flow2')
# # mynet.plot_aoi(network, res_locs, food_locs, edge_width='test_flow3')
# mynet.plot_aoi(network, res_locs, [food_locs, other_locs], edge_width='test_flow4')
# plt.show()

# summary_function -> lots can be added to this function in the future
# output5 = mynet.summary_function(network)

# # inundate_network
# inundated_net = mynet.inundate_network(network, f_graphs, inundation_raster)
# print(list(list(inundated_net.edges(data=True))[0][-1].keys()))

# FLOW_DECOMPOSITION
# single destination file, dest_method='single'
# output6a = mynet.flow_decomposition(network, 
#                                     res_points, 
#                                     food_points, 
#                                     res_locs, 
#                                     food_locs, 
#                                     dest_method='single')
# multiple destination files, dest_method='single
# output6b = mynet.flow_decomposition(network, 
#                                     res_points, 
#                                     [food_points, other_points], 
#                                     res_locs, 
#                                     [food_locs, other_locs], 
#                                     dest_method='single')
# # # single destination file, dest_method='multiple'
# output6c = mynet.flow_decomposition(network,res_points_c
#                                     res_points,
#                                     food_points,
#                                     res_locs,
#                                     food_locs,
#                                     dest_method='multiple')
# # # multiple destination files, dest_method='multiple
# # output6d = mynet.flow_decomposition(network,
# #                                     res_points,
# #                                     [food_points, other_points],
# #                                     res_locs,
# #                                     [food_locs, other_locs],
# #                                     dest_method='multiple')


# # # traffic_assignment

output7a = mynet.traffic_assignment(network, 
                                    res_points,
                                    food_points,
                                    food_locs,
                                    # G_weight='travel_time',
                                    G_weight='inundation_travel_time_agr',
                                    G_capacity='inundation_capacity_agr',
                                   dest_method='multiple',
                                   termination_criteria=['iter',3], 
                                   algorithm='path_based',
                                   method='MSA',
                                   sparse_array=True)
tstta = [round(num, 6) for num in output7a[2]]
sptta = [round(num, 6) for num in output7a[3]]
rga = [round(num, 6) for num in output7a[4]]
print(rga)
for u, v, data in network.edges(data=True):
    data["capacity"] *= 0.005

print('newrun')
time1=time.time()
output7b = mynet.traffic_assignment(network,
                                    res_points,
                                    [food_points, other_points, other_points2],
                                    [food_locs, other_locs, other_locs2],
                                    dest_method='multiple',
                                    termination_criteria=['iter', 10],
                                    algorithm='path_based',
                                    method='MSA',
                                    sparse_array=True,
                                    # G_weight='inundation_travel_time_agr',
                                    # G_capacity='inundation_capacity_agr'
                                    )
time2=time.time()
print(time2-time1)



# output6c = mynet.flow_decomposition(output7a[0],
#                                     res_points,
#                                     food_points,
#                                     res_locs,
#                                     food_locs,
#                                     dest_method='multiple',
#                                     G_weight='Weight_Array_Iter',
#                                     G_capacity='inundation_capacity_agr')

# output6c[2].to_file('/home/mdp0023/Desktop/external/Data/Network_Data/AOI_Testing/decomp_res_test.shp')




# tsttb = [round(num, 6) for num in output7b[2]]
# spttb = [round(num, 6) for num in output7b[3]]
# rgb = [round(num, 6) for num in output7b[4]]

# # print(tstta)
# # print(sptta)
# # print(rga)

# print(tsttb)
# print(spttb)
# print(rgb)


# # print(list(list(output7[0].edges(data=True))[0][-1].keys()))

# # mynet.plot_aoi(output7b[0], res_locs, [food_locs, other_locs, other_locs2], edge_width='TA_Flow')
# # plt.show()
