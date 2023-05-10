# sample_methodology_figures
# this code is to create the 4 small inset figures for the methodology components slide of EGU
# four figures are as follows:
    # 1. res parcels and colored resource parcels
    # 2. Inundation overlay
    # 3. classified road depths
    # 4. final traffic assignment (road utilization and decomposition)

# FIGURES WERE SAVED HERE:
# /home/mdp0023/Desktop/external/Data/Network_Data/AOI_Testing/example_figs

import rasterio as rio
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import osmnx as ox
import networkx as nx
import sys
import logging
import network_exploration_stuff as mynet

# AOI VARIABLES
# Folder paths
f_path = '/home/mdp0023/Desktop/external/Data/Network_Data/AOI_Testing'
f_boundaries = f'{f_path}/AOI_Boundary'
f_graphs = f'{f_path}/AOI_Graphs'
f_parcel_access_shp = f'{f_path}/AOI_Parcel_Access_Shapefiles'
f_res_shp = f'{f_path}/AOI_Residental_Parcel_Shapefiles'
f_resources=f'{f_path}/resource_parcel_samples'
# boundary
aoi_area = f'{f_boundaries}/Neighborhood_Network_AOI.shp'
# boundary with buffer
aoi_buffer = f'{f_boundaries}/Neighborhood_Network_AOI_Buf_1km.shp'
# Centroids of res parcels, and 3 resources
res = f'{f_res_shp}/Residential_Parcels_Points_Network_AOI.shp'
resource1 = f'{f_resources}/resource1_points.shp'
resource2 = f'{f_resources}/resource2_points.shp'
resource3 = f'{f_resources}/resource3_points.shp'
# Shapefiles of res parcels, and 3 resources
res_par = f'{f_res_shp}/Residential_Parcels_Network_AOI.shp'
resource1_par = f'{f_resources}/resource1.shp'
resource2_par = f'{f_resources}/resource2.shp'
resource3_par = f'{f_resources}/resource3.shp'
# Inundation raster
inundation_raster = f'{f_path}/AOI_Inundation.tif'
# bounding box
bbox=f'{f_resources}/bbox.shp'
raster = rio.open(inundation_raster)

# load shapefiles
res=gpd.read_file(res)
resource1=gpd.read_file(resource1)
resource2=gpd.read_file(resource2)
resource3=gpd.read_file(resource3)
res_par = gpd.read_file(res_par)
resource1_par = gpd.read_file(resource1_par)
resource2_par = gpd.read_file(resource2_par)
resource3_par = gpd.read_file(resource3_par)
bbox = gpd.read_file(bbox)

# load networks
network = mynet.read_graph_from_disk(path=f_graphs, name='AOI_Graph')
inundated_net = mynet.read_graph_from_disk(path=f_graphs, name='AOI_Graph_Inundated')
print(list(list(inundated_net.edges(data=True))[0][-1].keys()))

# using inundated_network, solve TA problem - use conservative inundation impact
output = mynet.traffic_assignment(inundated_net,
                                    res,
                                    [resource1, resource2, resource3],
                                    [resource1_par, resource2_par, resource3_par],
                                    dest_method='multiple',
                                    termination_criteria=['iter', 5],
                                    algorithm='path_based',
                                    method='MSA',
                                    sparse_array=True,
                                  G_weight='inundation_travel_time_con',
                                  G_capacity='inundation_capacity_con')

# decompose the TA assignment flow
flow_decomp = mynet.flow_decomposition(output[0],
                                    res,
                                       [resource1, resource2, resource3],
                                    res_par,
                                       [resource1_par, resource2_par, resource3_par],
                                    dest_method='multiple')
# save parcels to determine loss access
# flow_decomp[2].to_file(
#     '/home/mdp0023/Desktop/external/Data/Network_Data/AOI_Testing/resource_parcel_samples/loss_access_parcels.shp')
# figure 1: res and resource parcels
mynet.plot_aoi(network, 
               res_par, 
               [resource1_par, resource2_par, resource3_par],
               background_edges=network,
               bbox=bbox)

# figure 2: res and resource parcels with inundation
mynet.plot_aoi(network,
               res_par,
               [resource1_par, resource2_par, resource3_par],
               background_edges=network,
               bbox=bbox,
               inundation=raster)

# figure 3: classified road depths
mynet.plot_aoi(inundated_net,
               res_par,
               [resource1_par, resource2_par, resource3_par],
               background_edges=network,
               bbox=bbox,
               edge_color='max_inundation_mm',
               edge_width_weight=3)

# load lost access parcels
lost_access_parcels=gpd.read_file(f'{f_resources}/lost_access_parcels_only.shp')
# figure 4: routed vehicles and accessibility
mynet.plot_aoi(output[0],
               res_par,
               [resource1_par, resource2_par, resource3_par],
               background_edges=network,
               bbox=bbox,
               edge_width='TA_Flow',
               loss_access_parcels=lost_access_parcels)



plt.show()




