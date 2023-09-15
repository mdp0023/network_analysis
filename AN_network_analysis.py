# this is the code where I will do scratch work for the Austin North Network

# Packages
import network_exploration_stuff as mynet
import matplotlib.pyplot as plt
import osmnx as ox
import geopandas as gpd
import rasterio as rio
import networkx as nx
import pandas as pd
import numpy as np
import logging
import time
import math
import copy
import os

# ignore shapely deprectiation warnings
logging.captureWarnings(True)

############################################################################################
# INPUT VARIABLES
# Folder path
folder_path = '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North'
inundation_path = '/home/mdp0023/Documents/Codes_Projects/TEST_FOLDER/IO/combined_inundation'

# bounadry and boundary with buffer
aoi_area = f'{folder_path}/AN_Boundary/AN_Boundary.shp'
aoi_buffer = f'{folder_path}/AN_Boundary/AN_Boundary_3km.shp'

# background water
bg_water = f'/home/mdp0023/Desktop/external/Data/Inundation/Water_Boundaries/Austin_OSM_water_boundaries.shp'

# Centroids and parcels
res_parcel_loc = f'{folder_path}/AN_Residential_Parcel_Shapefiles/AN_Residential_Parcels.shp'
res_points_loc = f'{folder_path}/AN_Residential_Parcel_Shapefiles/AN_Residential_Parcels_Points.shp'
# supermarkets
food_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_Supermarket_Parcel.shp'
food_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_Supermarket_Parcel_Points.shp'
# emergency rooms
er_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_hospital_parcels.shp'
er_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_hospital_points.shp'
# pharmacies
pharm_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_pharmacy_parcels.shp'
pharm_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_pharmacy_points.shp'
# police stations
police_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_police_parcels.shp'
police_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_police_points.shp'
# convenience stores
conv_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_convenience_parcels.shp'
conv_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_convenience_points.shp'
# fire stations
fire_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_fire_station_parcels.shp'
fire_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_fire_station_points.shp'
# ems stations
ems_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_ambulance_station_parcels.shp'
ems_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_ambulance_station_points.shp'
# fuel stations
fuel_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_fuel_parcels.shp'
fuel_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_fuel_points.shp'

# read shapefiles
res_parcels = gpd.read_file(res_parcel_loc)
res_points = gpd.read_file(res_points_loc)

food_parcels = gpd.read_file(food_parcel_loc)
food_points = gpd.read_file(food_points_loc)

er_parcels = gpd.read_file(er_parcel_loc)
er_points = gpd.read_file(er_points_loc)

pharm_parcels = gpd.read_file(pharm_parcel_loc)
pharm_points = gpd.read_file(pharm_points_loc)

police_parcels = gpd.read_file(police_parcel_loc)
police_points = gpd.read_file(police_points_loc)

conv_parcels = gpd.read_file(conv_parcel_loc)
conv_points = gpd.read_file(conv_points_loc)

fire_parcels = gpd.read_file(fire_parcel_loc)
fire_points = gpd.read_file(fire_points_loc)

ems_parcels = gpd.read_file(ems_parcel_loc)
ems_points = gpd.read_file(ems_points_loc)

fuel_parcels = gpd.read_file(fuel_parcel_loc)
fuel_points = gpd.read_file(fuel_points_loc)


all_resource_parcels=[food_parcels, er_parcels, pharm_parcels, police_parcels,
                      conv_parcels, fire_parcels, ems_parcels, fuel_parcels]
all_resource_points=[food_points, er_points, pharm_points, police_points,
                      conv_points, fire_points, ems_points, fuel_points]

aoi_shape = gpd.read_file(aoi_area)

# Inundation raster
bg_water = gpd.read_file(bg_water)
# inundation_raster = f'{folder_path}/AN_Inundation/AN_Memorial_Day_Compound.tif'
# raster = rio.open(inundation_raster)
###############################################################################cond#############

# download graph
# ox.config(use_cache=True)
# G = mynet.shape_2_graph(aoi_buffer)
# G = mynet.rename(G=G)

# save graph
# mynet.save_2_disk(G=G, path=f'{folder_path}/AN_Graphs', name='AN_Graph_w_bridges')

# open saved graph
G = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_w_bridges')
# G_inun = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_Inundated')
# # G_TA_flood = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_Traffic_Assignment_AGR_Inundation')
# G_TA_no_flood = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_Traffic_Assignment_No_Inundation')
# G = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_2015052518_inundation_TA')

print(list(list(G.edges(data=True))[0][-1].keys()))



# Inundate network
# G_inun = mynet.inundate_network(G=G, 
#                                 path=f'{folder_path}/AN_Graphs', 
#                                 inundation=f'{folder_path}/AN_Inundation/AN_Memorial_Day_Compound.tif')

# available edge and node attributes
# edge_attributes = list(list(G.edges(data=True))[0][-1].keys())
# node_attributes = list(list(G.nodes(data=True))[0][-1].keys())
# print(edge_attributes)

# get unique speed values across network (check for major errors)
# unique_speeds = np.unique(list(nx.get_edge_attributes(G, 'speed_kph').values()))
# percent_edge_w_speed = len(list(nx.get_edge_attributes(G, 'speed_kph').values()))/G.number_of_edges()
# plt.hist(list(nx.get_edge_attributes(G, 'speed_kph').values()))
# plt.show()

# # get list of unique road classifications
# unique_rtypes = list(set(list(nx.get_edge_attributes(G, 'highway').values())))
# percent_edge_w_rtype = len(list(nx.get_edge_attributes(G, 'highway').values()))/G.number_of_edges()
# print(unique_rtypes)

# # get unique capacity values (check for errors)
# unique_caps = np.unique(list(nx.get_edge_attributes(G, 'capacity').values()))
# percent_edge_w_cap = len(list(nx.get_edge_attributes(G_inun, 'capacity').values()))/G.number_of_edges()
# print(unique_caps)

# # get unique number of lanes (check for errors)
# num_lanes = np.unique(list(nx.get_edge_attributes(G, 'lanes').values()))
# percent_edge_w_lane = len(list(nx.get_edge_attributes(G, 'lanes').values()))/G.number_of_edges()
# print(num_lanes)

# # see percentage of edges that have respected stats
# print(f"percentage edges with lane info : {percent_edge_w_lane}")
# print(f"percentage edges with speed info: {percent_edge_w_speed}")
# print(f"percentage edges with rtype info: {percent_edge_w_rtype}")
# print(f"percentage edges with cap info  : {percent_edge_w_cap}")

# edge_attributes = list(list(G.edges(data=True))[0][-1].keys())
# print(edge_attributes)
# print(f'number of nodes: {len(G.nodes())}')
# print(f'number of edges: {len(G.edges())}')

# output = mynet.max_flow_parcels(G=G, 
#                                 res_points=res_points, 
#                                 dest_points=food_points,
#                                 G_capacity='capacity', 
#                                 G_weight='travel_time', 
#                                 G_demand='demand',
#                                 dest_method='multiple', 
#                                 dest_parcels=food_parcels, 
#                                 ignore_capacity=False)

# output = mynet.max_flow_parcels(G=G_inun,
#                                 res_points=res_points,
#                                 dest_points=food_points,
#                                 G_capacity='inundation_capacity_agr',
#                                 G_weight='inundation_travel_time_agr',
#                                 G_demand='demand',
#                                 dest_method='multiple',
#                                 dest_parcels=food_parcels,
#                                 ignore_capacity=False)

# flow_dictionary = output[0]
# cost_of_flow = output[1]
# max_flow = output[2]
# access = output[3]

# print(f'level of access: {access}')
# # print(f'cost of flow: {cost_of_flow}') 
# print(f'Maximum amount of flow: {max_flow}')

# #relate flow dictionary back to DRY graph for plotting purposes
# G_map = nx.DiGraph(G)
# # relate
# for edge in G_map.edges:
#     values = {(edge[0], edge[1]): {'wet_flow':
#                                    flow_dictionary[edge[0]][edge[1]]}}
#     nx.set_edge_attributes(G=G_map, values=values)

# # need to convert back to multidigraph to plot properly
# G_map = nx.MultiDiGraph(G_map)



########################################################################
#FLOW DECOMPOSITION TEST
# output = mynet.flow_decomposition(G=G_inun, 
#                                         res_points=res_points,
#                                         dest_points=food_points, 
#                                         res_parcels=res_parcels, 
#                                         dest_parcels=food_parcels, 
#                                         G_demand='demand', 
#                                         G_capacity='inundation_capacity_agr', 
#                                         G_weight='inundation_travel_time_agr',
#                                         dest_method='multiple')

# output = mynet.flow_decomposition(G=G, 
#                                         res_points=res_points,
#                                         dest_points=food_points, 
#                                         res_parcels=res_parcels, 
#                                         dest_parcels=food_parcels, 
#                                         G_demand='demand', 
#                                         G_capacity='capacity', 
#                                         G_weight='travel_time',
#                                         dest_method='multiple')

  
# decomposed_paths=output[0]
# sink_insights=output[1]
# res_parcels=output[2]
# dest_parcels=output[3]

# res_parcels.plot(column='cost_of_flow', legend=True)

# costs = sorted(res_parcels['cost_of_flow'].tolist())
# costs = [x for x in costs if math.isnan(x)==False]
# plt.hist(costs)
# plt.show()


##########################################################################
# PLOTTING

# mynet.plot_aoi(G=G_map, 
#                 res_parcels=res_parcels,
#                 resource_parcels=food_parcels, 
#                 edge_width='wet_flow',
#                 decomp_flow=True,
#                 loss_access_parcels=res_parcels.loc[res_parcels['cost_of_flow'].isna()])
# plt.show()


#####


# # inundate networks for each time step of flood estimate FLUVIAL FLOODING
# rasters = ['geoflood_merge_2015052521',
#            'geoflood_merge_2015052522',
#            'geoflood_merge_2015052523',
#            'geoflood_merge_2015052600',
#            'geoflood_merge_2015052601',
#            'geoflood_merge_2015052602', 
#            'geoflood_merge_2015052603',
#            'geoflood_merge_2015052604',
#            'geoflood_merge_2015052605',
#            'geoflood_merge_2015052606',
#            'geoflood_merge_2015052607',
#            'geoflood_merge_2015052608',
#            'geoflood_merge_2015052609',
#            'geoflood_merge_2015052610',
#            'geoflood_merge_2015052611',
#            'geoflood_merge_2015052612',
#            'geoflood_merge_2015052613',
#            'geoflood_merge_2015052614',
#            'geoflood_merge_2015052615',
#            'geoflood_merge_2015052616',
#            'geoflood_merge_2015052617'
#            ]
# for inundation_raster in rasters:
#     mynet.inundate_network(G,
#                             f'{folder_path}/AN_Graphs',
#                             f'/home/mdp0023/Documents/Codes_Projects/TEST_FOLDER/IO/combined_inundation_copy/{inundation_raster}',
#                             name=f'AN_Graph_{inundation_raster[15:]}_inundation_fluvial')
#     print(f'{inundation_raster} calculated')

# inundate networks for each time step of flood estimate PLUVIAL FLOODING
# rasters = ['fsm_merge_2015052603',
#            'fsm_merge_2015052604',
#            'fsm_merge_2015052605',
#            'fsm_merge_2015052606',
#            'fsm_merge_2015052607',
#            'fsm_merge_2015052608',
#            'fsm_merge_2015052609',
#            'fsm_merge_2015052610',
#            'fsm_merge_2015052611',
#            'fsm_merge_2015052612',
#            'fsm_merge_2015052613',
#            'fsm_merge_2015052614',
#            'fsm_merge_2015052615',
#            'fsm_merge_2015052616',
#            'fsm_merge_2015052617']
# for inundation_raster in rasters:
#     mynet.inundate_network(G,
#                            f'{folder_path}/AN_Graphs',
#                            f'/home/mdp0023/Documents/Codes_Projects/TEST_FOLDER/IO/combined_inundation_copy/{inundation_raster}',
#                            name=f'AN_Graph_{inundation_raster[10:]}_inundation_pluvial')
#     print(f'{inundation_raster} calculated')

# # inundate networks for each time step of flood estimate COMPOUND FLOODING
# rasters = ['2015052603_inundation.tif',
#            '2015052604_inundation.tif',
#            '2015052605_inundation.tif',
#            '2015052606_inundation.tif',
#            '2015052607_inundation.tif',
#            '2015052608_inundation.tif',
#            '2015052609_inundation.tif',
#            '2015052610_inundation.tif',
#            '2015052611_inundation.tif',
#            '2015052612_inundation.tif',
#            '2015052613_inundation.tif',
#            '2015052614_inundation.tif',
#            '2015052615_inundation.tif',
#            '2015052616_inundation.tif',
#            '2015052617_inundation.tif']
# for inundation_raster in rasters:
#     mynet.inundate_network(G,
#                            f'{folder_path}/AN_Graphs',
#                            f'/home/mdp0023/Documents/Codes_Projects/TEST_FOLDER/IO/combined_inundation_copy/{inundation_raster}',
#                            name=f'AN_Graph_{inundation_raster[0:10]}_inundation')
#     print(f'{inundation_raster} calculated')


##########################################################################
# MIN_COST_FLOW_PARCELS TEST
# flow_dictionary, cost_of_flow = mynet.min_cost_flow_parcels(G=G, 
#                                                                 res_points=res_points, 
#                                                                 dest_points=food_points, 
#                                                                 dest_parcels=food_parcels, 
#                                                                 G_demand='demand', 
#                                                                 G_capacity='capacity', 
#                                                                 G_weight='travel_time', 
#                                                                 dest_method='multiple')
# print(cost_of_flow)

##########################################################################
# TRAFFIC ASSIGNMENT AND CORRESPONDING SPTT
# G_output, AEC_list, TSTT_list, SPTT_list, RG_list, iter = mynet.traffic_assignment(G=G,
#                                                                                    res_points=res_points,
#                                                                                    dest_points=all_resource_points,
#                                                                                    dest_parcels=all_resource_parcels,
#                                                                                    G_capacity='capacity',
#                                                                                    G_weight='travel_time',
#                                                                                    algorithm='path_based',
#                                                                                    method='MSA',
#                                                                                    link_performance='BPR',
#                                                                                    termination_criteria=[
#                                                                                        'iter', 10],
#                                                                                    dest_method='multiple',
#                                                                                    sparse_array=True)

# iterate through graphs and run TA assignment for all resources

tstamps =   ['2015052604',
           '2015052605',
           '2015052606',
           '2015052607',
           '2015052608',
           '2015052609',
           '2015052610',
           '2015052611',
           '2015052612',
           '2015052613',
           '2015052614',
           '2015052615',
           '2015052616',
           '2015052617']


# for tstamp in tstamps:
#     G = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name=f'AN_Graph_{tstamp}_inundation')


#     time1=time.time()
#     G_output, AEC_list, TSTT_list, SPTT_list, RG_list, iter = mynet.traffic_assignment(G=G,
#                                                                                         res_points=res_points, 
#                                                                                         dest_points=all_resource_points,
#                                                                                         dest_parcels=all_resource_parcels,  
#                                                                                         G_capacity='inundation_capacity_agr',
#                                                                                         G_weight='inundation_travel_time_agr',
#                                                                                         algorithm='path_based', 
#                                                                                         method='MSA',
#                                                                                         link_performance='BPR',
#                                                                                         termination_criteria=['AEC',100],
#                                                                                         dest_method='multiple',
#                                                                                         sparse_array=True)
#     time2=time.time()
#     # save graph
#     mynet.save_2_disk(G=G_output, path=f'{folder_path}/AN_Graphs', name=f'AN_Graph_{tstamp}_inundation_TA')
    
#     # save outputss
#     list_dict={'AECs' : AEC_list,
#             'TSTTs' : TSTT_list,
#             'SPTTs' : SPTT_list,
#             'RGs' : RG_list,
#             'iters' : iter}
#     outputdf = pd.DataFrame(list_dict)
#     outputdf.to_csv(f'{folder_path}/AN_Graphs/TA_Summary_Files/{tstamp}_inundation.csv', index=False)

# FLOW DECOMPOSITION
# for tstamp in tstamps:
#     G = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name=f'AN_Graph_{tstamp}_inundation_TA')

#     res_points1=copy.copy(res_points)
#     res_parcels1=copy.copy(res_parcels)
#     all_resource_points1=copy.copy(all_resource_points)
#     all_resource_parcels1 = copy.copy(all_resource_parcels)
#     output = mynet.flow_decomposition(G=G,
#                                         res_points=res_points1,
#                                         res_parcels=res_parcels1,
#                                         dest_points=all_resource_points1,
#                                         dest_parcels=all_resource_parcels1,
#                                         G_demand='demand',
#                                       G_capacity='inundation_capacity_agr',
#                                         G_weight='Weight_Array_Iter',
#                                         dest_method='multiple')
#     output[2].to_file(f'/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs/Flow_decomp/{tstamp}_inundation_res_parcel_flow_decomp.shp')

#     print(f'{tstamp} calculated')


# print(f'AEC_list: {AEC_list}')
# print(f'TSTT_list: {TSTT_list}')
# print(f'SPTT_list: {SPTT_list}')
# print(f'RG_list: {RG_list}')
# print(f'iter: {iter}')

##########################################################################
# # NO FLOOD TRAVEL TIME 
# G = mynet.read_graph_from_disk(
#     path=f'{folder_path}/AN_Graphs', name=f'AN_Graph_2015052516_inundation_TA_pre_inundation')

# output = mynet.flow_decomposition(G=G,
#                                   res_points=res_points,
#                                   res_parcels=res_parcels,
#                                   dest_points=all_resource_points,
#                                   dest_parcels=all_resource_parcels,
#                                   G_demand='demand',
#                                   G_capacity='capacity',
#                                   G_weight='Weight_Array_Iter',
#                                   dest_method='multiple')
# output[2].to_file(
#     f'/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs/Flow_decomp/AN_Graph_2015052516_No_Flood_res_parcel_flow_decomp.shp')

# decomposed_paths=output[0]
# sink_insights=output[1]
# res_parcels=output[2]
# dest_parcels=output[3]


# mynet.plot_aoi(G=G_TA_no_flood,
#                 background_edges=G,
#                 bbox=aoi_shape,
#                 res_parcels=res_parcels,
#                 resource_parcels=food_parcels,
#                 edge_width='TA_Flow',
#                 decomp_flow=True,
#                 # loss_access_parcels=res_parcels.loc[res_parcels['cost_of_flow'].isna()],
#                 scalebar=True)


# plt.show()


############################################################################################
# INDIVIDUAL REDUNDNACY CALCULATION
# times=['2015052516',
#            '2015052517',
#            '2015052518',
#            '2015052519',
#            '2015052520',
#            '2015052521',
#            '2015052522',
#            '2015052523',
#            '2015052600',
#            '2015052601',
#            '2015052602',
#            '2015052603',
#            '2015052617']
# for t in times:

#     G = mynet.read_graph_from_disk(
#         path=f'{folder_path}/AN_Graphs', name=f'AN_Graph_{t}_inundation_TA')

#     time1=time.time()
#     output8 = mynet.redundancy_metric(G,
#                                     res_points=res_points,
#                                     res_parcels=res_parcels,
#                                     dest_points=all_resource_points,
#                                     dest_parcels=all_resource_parcels,
#                                     inc=1500,    # additional for every 15 seconds 
#                                     max_length=90000,  # set to 15 min X 60 for seconds (X100 for current iteration of TA)
#                                     G_capacity='inundation_capacity_agr',
#                                     G_weight='Weight_Array_Iter')

#     # save to pull into qgis
#     output8.to_file(
#         f'/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Redundancy/individual_redundancy/individual_redundancy_{t}.shp')

#     time2=time.time()
#     print(time2-time1)
#     print(f'{t} individual redundancy calculated')

# NETWORK REDUNDNACY CALCULATION
name=['_', '_fluvial_', '_pluvial_']
times = ['2015052603']
for idx, ftype in enumerate(['', '_fluvial','_pluvial']):
    for t in times:

        G = mynet.read_graph_from_disk(
            path=f'{folder_path}/AN_Graphs', name=f'AN_Graph_{t}_inundation{ftype}')
        
        time1=time.time()
        output9 = mynet.network_redundancy_metric(G=G,
                                                res_points=res_points,
                                                res_parcels=res_parcels,
                                                dest_points=all_resource_points,
                                                dest_parcels=all_resource_parcels,
                                                G_capacity='inundation_capacity_agr', # edges with 0 capacity removed
                                                G_weight='travel_time')        # I don't think TA weight needs to be used, weight shouldn't play role in this calc
                                                # G_weight='Weight_Array_Iter')

        # save to pull into qgis
        output9.to_file(
            f'/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Redundancy/network_redundancy/network_redundancy{name[idx]}{t}_n2.shp')
        time2=time.time()
        print(time2-time1)
        print(f'{ftype} {t} individual redundancy calculated')



##############################################################################
# # PLOT INUNDATED NETWORK BY DEPTH ON ROADWAY
# G_temp = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_2015052516_inundation')
# mynet.plot_aoi(G=G_temp,
#                res_parcels=res_parcels,
#                resource_parcels=food_parcels,
#                rotation=-90,
#                scalebar=True,
#                bbox=aoi_shape,
#                bg_water=bg_water)


# mynet.plot_aoi(G=G_temp,
#                 res_parcels=res_parcels,
#                 bbox=aoi_shape,
#                 resource_parcels=food_parcels,
#                 edge_color='max_inundation_mm', 
#                 scalebar=True)

# plt.show()


# See roads in each category
# road_depths_0_1=[]
# road_depths_1_15=[]
# road_depths_15_30=[]
# road_depths_30_60=[]
# road_depths_60_over = []
# for x in G_inun.edges(data=True):
#     # each x has from node, to node, and data variable
#     data = x[2]
#     if data['max_inundation_mm'] < 1:
#         road_depths_0_1.append(data['max_inundation_mm'])
#     elif data['max_inundation_mm'] < 15:
#         road_depths_1_15.append(data['max_inundation_mm'])
#     elif data['max_inundation_mm'] < 30: 
#         road_depths_15_30.append(data['max_inundation_mm'])
#     elif data['max_inundation_mm'] < 60:
#         road_depths_30_60.append(data['max_inundation_mm'])
#     elif data['max_inundation_mm'] >= 60:
#         road_depths_60_over.append(data['max_inundation_mm'])


# print(len(road_depths_0_1))
# print(len(road_depths_1_15))
# print(len(road_depths_15_30))
# print(len(road_depths_30_60))
# print(len(road_depths_60_over))


