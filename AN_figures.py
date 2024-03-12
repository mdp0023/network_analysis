import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning)

import matplotlib.patches as mpatches
import os
import math
import time
import logging
import numpy as np
import pandas as pd
import networkx as nx
import rasterio as rio
import geopandas as gpd
import osmnx as ox
import matplotlib.pyplot as plt
import network_analysis_base as mynet
import matplotlib as mpl
import rasterio as rio

# this is the code where I will do scratch work for the Austin North Network FIGURE MAKING

# ignore shapely deprectiation warnings
logging.captureWarnings(True)

############################################################################################
# INPUT VARIABLES
# Folder path
folder_path = '/media/mdp0023/extradrive1/Data/Network_Data/Austin_North'
inundation_path = '/home/mdp0023/Documents/Codes_Projects/TEST_FOLDER/IO/combined_inundation'

# bounadry and boundary with buffer
aoi_area = f'{folder_path}/AN_Boundary/AN_Boundary.shp'
aoi_buffer = f'{folder_path}/AN_Boundary/AN_Boundary_3km.shp'

# background water
bg_water = f'/media/mdp0023/extradrive1/Data/Inundation/Water_Boundaries/Austin_OSM_water_boundaries.shp'

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


all_resource_parcels = [food_parcels, er_parcels, pharm_parcels, police_parcels,
                        conv_parcels, fire_parcels, ems_parcels, fuel_parcels]
all_resource_points = [food_points, er_points, pharm_points, police_points,
                       conv_points, fire_points, ems_points, fuel_points]

aoi_shape = gpd.read_file(aoi_area)

# Inundation raster
bg_water = gpd.read_file(bg_water)

# svi
bg_svi = gpd.read_file(
    '/home/mdp0023/Documents/Codes_Projects/SVI_Code/Travis_County/SVI_Shapefiles/Travis_county_svi_2015_selected.shp')

#network_redun
net_redun = gpd.read_file('/media/mdp0023/extradrive1/Data/Network_Data/Austin_North/AN_Redundancy/network_redundancy_n2/network_redundancy_2015052522.shp')

# household redun
house_redun = gpd.read_file(
    '/media/mdp0023/extradrive1/Data/Network_Data/Austin_North/AN_Redundancy/individual_redundancy/individual_redundancy_2015052522.shp')
# set font family
plt.rcParams['font.family'] = "ubuntu"


# figure size converstion
mm = 1/25.4
###############################################################################################
# plot PEAK inundation on roads
# G_temp = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_2015052522_inundation')
# figure,ax=mynet.plot_aoi(G=G_temp,
#                res_parcels=res_parcels,
#                resource_parcels=all_resource_parcels,
#                rotation=-90,
#                scalebar=True,
#                bbox=aoi_shape,
#                bg_water=bg_water,
#                edge_color='max_inundation_mm')

# figure.savefig('/home/mdp0023/Documents/Codes_Projects/Network_Analysis/preliminary_results/test.pdf',bbox_inches='tight')
# plt.show()

###############################################################################################
# # Plot of Traffic Patterns
# G_temp = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_2015052522_inundation_TA')
# figure,ax=mynet.plot_aoi(G=G_temp,
#                res_parcels=res_parcels,
#                resource_parcels=all_resource_parcels,
#                rotation=-90,
#                scalebar=True,
#                bbox=aoi_shape,
#                bg_water=bg_water,
#                edge_width='TA_Flow')

# figure.savefig('/home/mdp0023/Documents/Codes_Projects/Network_Analysis/preliminary_results/test.pdf',bbox_inches='tight')
# plt.show()

# ###############################################################################################
# # # Plot study area map and svi
# G_temp = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_2015052522_inundation_TA')
# figure,ax=mynet.plot_aoi(G=G_temp,
#                res_parcels=res_parcels,
#                resource_parcels=all_resource_parcels,
#                rotation=-90,
#                scalebar=True,
#                bbox=aoi_shape,
#                bg_water=bg_water)


# # set bbox to localized area
# total_bounds = aoi_shape.total_bounds
# # convert bounds to N, S, E, W
# total_bounds = (total_bounds[3],
#                 total_bounds[1],
#                 total_bounds[2],
#                 total_bounds[0])
# center_point=((total_bounds[2]-total_bounds[3])/2 + total_bounds[3],
#                 (total_bounds[0]-total_bounds[1])/2 + total_bounds[1])

# # # add svi as background
# # bg_svi_rot = gpd.GeoSeries(g for g in bg_svi['geometry'])
# # bg_svi_rot = bg_svi_rot.rotate(-90, origin=center_point)


# # envgdf = gpd.GeoDataFrame(bg_svi, geometry=bg_svi_rot)
# # # envgdf = envgdf.rename(columns={0: 'geometry_new'}).set_geometry('geometry_new')


# # # bg_svi = bg_svi.merge(envgdf, on='index')
# # # bg_svi = bg_svi.set_geometry('geometry_new')

# # envgdf.plot(ax=ax, column='Rank', cmap='bwr_r', zorder=100, alpha=0.75)
# # envgdf.boundary.plot(ax=ax,  zorder=101, alpha=0.75, edgecolor='black')


# figure.savefig('/home/mdp0023/Documents/Codes_Projects/Network_Analysis/preliminary_results/test.pdf',bbox_inches='tight')
# plt.show()


###########################################################################################################
# Austria Network
# graph = ox.graph.graph_from_bbox(east=16.416622+.075,west=16.41662-0.075,north=48.233427+0.025,south=48.233427-0.025, network_type='drive')
# G_gdf_edges = ox.graph_to_gdfs(G=graph, nodes=False)

# bg_water = gpd.read_file('/home/mdp0023/Desktop/external/Data/Inundation/Vienna_water/vienna_waterways.shp')
# figure, ax = mynet.plot_aoi(G=graph,
#                res_parcels=res_parcels,
#                resource_parcels=all_resource_parcels,
#                bg_water=bg_water)

# figure.savefig('/home/mdp0023/Documents/Codes_Projects/Network_Analysis/preliminary_results/test.pdf',bbox_inches='tight')
# plt.show()


#################################################################################################################
# # SVI, SS, and ES
# # read in data
# G_temp = mynet.read_graph_from_disk(
#     path=f'{folder_path}/AN_Graphs', name='AN_Graph_2015052522_inundation_TA')
# # create the figure 
# fig, ax = plt.subplots(1, 3, figsize=(180*mm, 70*mm))

# # Plot study area network info and remove axis information 
# trues = [True, False, False]
# for n in [0, 1, 2]:
#     figure, axe = mynet.plot_aoi(G=G_temp,
#                 res_parcels=res_parcels,
#                 resource_parcels=all_resource_parcels,
#                 rotation=-90,
#                 scalebar=trues[n],
#                 bbox=aoi_shape,
#                 bg_water=bg_water,
#                 default=[fig, ax[n]])
    
#     # Hide X and Y axes label marks
#     ax[n].xaxis.set_tick_params(labelbottom=False)
#     ax[n].yaxis.set_tick_params(labelleft=False)
#     # Hide X and Y axes tick marks
#     ax[n].set_xticks([])
#     ax[n].set_yticks([])

# # DETERMINE PIVOT POINT
# # set bbox to localized area
# total_bounds = aoi_shape.total_bounds
# # convert bounds to N, S, E, W
# total_bounds = (total_bounds[3],
#                 total_bounds[1],
#                 total_bounds[2],
#                 total_bounds[0])
# center_point=((total_bounds[2]-total_bounds[3])/2 + total_bounds[3],
#                 (total_bounds[0]-total_bounds[1])/2 + total_bounds[1])

# # get geometry transformation data
# bg_svi_rot = gpd.GeoSeries(g for g in bg_svi['geometry'])
# bg_svi_rot = bg_svi_rot.rotate(-90, origin=center_point)

# # open svi file
# envgdf = gpd.GeoDataFrame(bg_svi, geometry=bg_svi_rot)

# # plot svi and  border
# envgdf.plot(ax=ax[0], column='Rank', cmap='coolwarm_r',
#             zorder=100, alpha=0.75, scheme='quantiles', k=10)
# envgdf.boundary.plot(ax=ax[0],  zorder=101, alpha=0.75, edgecolor='black',linewidth=0.5)

# # plot ss and border
# envgdf.plot(ax=ax[1], column='Factor 1', cmap='coolwarm', zorder=100, alpha=0.75, scheme='quantiles',k=10)
# envgdf.boundary.plot(ax=ax[1],  zorder=101, alpha=0.75, edgecolor='black',linewidth=0.5)

# # plot es and border
# envgdf.plot(ax=ax[2], column='Factor 2', cmap='coolwarm_r', zorder=100, alpha=0.75, scheme='quantiles',k=10)
# envgdf.boundary.plot(ax=ax[2],  zorder=101, alpha=0.75, edgecolor='black',linewidth=0.5)

# # ax[0].set_xlabel('SVI', style ='italic',fontsize=12, weight='bold')
# # ax[1].set_xlabel('Social Status', style ='italic',fontsize=12, weight='bold')
# # ax[2].set_xlabel('Economic Status', style ='italic',fontsize=12, weight='bold')

# ax[0].set_title('SVI', style ='italic',fontsize=12, weight='bold')
# ax[1].set_title('Social Status', style ='italic',fontsize=12, weight='bold')
# ax[2].set_title('Economic Status', style ='italic',fontsize=12, weight='bold')


# # add colorbar
# n_bins=10
# # Custom color ramp
# cmap = mpl.cm.coolwarm
# # using quantile method, determine number of bins
# bounds = np.linspace(0, 1, n_bins+1)
# np.delete(bounds, 0)
# bins = []
# for bound in bounds:
#     bins.append(envgdf["SVI_scaled"].quantile(bound))
# norm = mpl.colors.BoundaryNorm(bins, cmap.N)
# bins_r = [round(elem, 2) for elem in bins]

# cax = fig.add_axes([0.05,0.175,0.90,0.075])
# cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#                     cax=cax, orientation='horizontal',
#                     label="Indicator Percentile")
# cbar.ax.set_xticklabels(['0','10','20','30','40','50','60','70','80','90','100'])
# plt.xticks(bins)

# idx1=0

# for letter in ['A','B','C']:
# # add a,b,c, etc.
#     ax[idx1].annotate(letter, 
#                 xy=(0.05,0.85), 
#                 xycoords='axes fraction',
#                 xytext=(0.05, 0.85),
#                 textcoords='axes fraction',
#                 va='bottom', 
#                 ha='center',
#                 fontsize=10,
#                 bbox=dict(boxstyle="Round",
#                         fc="lightgray", ec="darkgray", lw=0.5))
#     idx1+=1



# plt.subplots_adjust(top= 1.0,
#                     bottom=0.2,
#                     left=0.025,
#                     right=0.975,
#                     hspace=1.0,
#                     wspace=0.05)


# ##########################################################################
# # methods figure

# # read in data
# G_temp = mynet.read_graph_from_disk(
#     path=f'{folder_path}/AN_Graphs', name='AN_Graph_2015052522_inundation_TA')
# G_temp_inun = mynet.read_graph_from_disk(path=f'{folder_path}/AN_Graphs', name='AN_Graph_2015052522_inundation')
# # create the figure
# fig, ax = plt.subplots(3, 2, figsize=(180*mm, 210*mm))

# # DETERMINE PIVOT POINT 
# # set bbox to localized area
# total_bounds = aoi_shape.total_bounds
# # convert bounds to N, S, E, W
# total_bounds = (total_bounds[3],
#                 total_bounds[1],
#                 total_bounds[2],
#                 total_bounds[0])
# center_point = ((total_bounds[2]-total_bounds[3])/2 + total_bounds[3],
#                 (total_bounds[0]-total_bounds[1])/2 + total_bounds[1])

# # get geometry transformation data for parcel data
# net_redun_rot = gpd.GeoSeries(g for g in net_redun['geometry'])
# net_redun_rot = net_redun_rot.rotate(-90, origin=center_point)

# # get geometry transfomration data for SVI data
# svi_rot = gpd.GeoSeries(g for g in bg_svi['geometry'])
# svi_rot = svi_rot.rotate(-90, origin=center_point)


# ####################################################################################################################################################
# # # [0][0] inundation
# figure, axe = mynet.plot_aoi(G=G_temp_inun,
#                              res_parcels=res_parcels,
#                              resource_parcels=all_resource_parcels,
#                              rotation=-90,
#                              scalebar=True,
#                              bbox=aoi_shape,
#                              bg_water=bg_water,
#                              edge_color='max_inundation_mm',
#                              default=[fig, ax[0][0]])
# # Hide X and Y axes label marks
# ax[0][0].xaxis.set_tick_params(labelbottom=False)
# ax[0][0].yaxis.set_tick_params(labelleft=False)
# # Hide X and Y axes tick marks
# ax[0][0].set_xticks([])
# ax[0][0].set_yticks([])

# ####################################################################################################################################################
# # [0][1] traffic assignment
# figure, axe = mynet.plot_aoi(G=G_temp,
#                              res_parcels=res_parcels,
#                             resource_parcels=all_resource_parcels,
#                             rotation=-90,
#                             scalebar=False,
#                             bbox=aoi_shape,
#                             bg_water=bg_water,
#                             edge_width='TA_Flow',
#                              default=[fig, ax[0][1]])
# # Hide X and Y axes label marks
# ax[0][1].xaxis.set_tick_params(labelbottom=False)
# ax[0][1].yaxis.set_tick_params(labelleft=False)
# # Hide X and Y axes tick marks
# ax[0][1].set_xticks([])
# ax[0][1].set_yticks([])


# ####################################################################################################################################################
# # [1][0] newtork redundancy
# figure, axe = mynet.plot_aoi(G=G_temp,
#                              res_parcels=res_parcels,
#                             resource_parcels=all_resource_parcels,
#                             rotation=-90,
#                             scalebar=False,
#                             bbox=aoi_shape,
#                             bg_water=bg_water,
#                              default=[fig, ax[1][0]])

# # open file to plot
# envgdf = gpd.GeoDataFrame(net_redun, geometry=net_redun_rot)

# # plot network redundancy
# envgdf.plot(ax=ax[1][0], column='net_redun', cmap='coolwarm_r',
#             zorder=100, alpha=0.95, scheme='quantiles',linewidth=0)

# # Hide X and Y axes label marks
# ax[1][0].xaxis.set_tick_params(labelbottom=False)
# ax[1][0].yaxis.set_tick_params(labelleft=False)
# # Hide X and Y axes tick marks
# ax[1][0].set_xticks([])
# ax[1][0].set_yticks([])


# ####################################################################################################################################################
# # [2][0] Household Redundancy
# figure, axe = mynet.plot_aoi(G=G_temp,
#                              res_parcels=res_parcels,
#                              resource_parcels=all_resource_parcels,
#                              rotation=-90,
#                              scalebar=False,
#                              bbox=aoi_shape,
#                              bg_water=bg_water,
#                              default=[fig, ax[2][0]])

# # open file to plot
# envgdf = gpd.GeoDataFrame(house_redun, geometry=net_redun_rot)

# # plot svi and  border
# envgdf.plot(ax=ax[2][0], column='redundancy', cmap='coolwarm_r',
#             zorder=100, alpha=0.95, scheme='quantiles',linewidth=0)

# # Hide X and Y axes label marks
# ax[2][0].xaxis.set_tick_params(labelbottom=False)
# ax[2][0].yaxis.set_tick_params(labelleft=False)
# # Hide X and Y axes tick marks
# ax[2][0].set_xticks([])
# ax[2][0].set_yticks([])


# ####################################################################################################################################################
# # [1][1] newtork reliability
# figure, axe = mynet.plot_aoi(G=G_temp,
#                              res_parcels=res_parcels,
#                              resource_parcels=all_resource_parcels,
#                              rotation=-90,
#                              scalebar=False,
#                              bbox=aoi_shape,
#                              bg_water=bg_water,
#                              default=[fig, ax[1][1]])

# # perfrom network reliability calcuation at appropraite time
# fpath = '/media/mdp0023/extradrive1/Data/Network_Data/Austin_North/AN_Graphs'
# prefix='AN_Graph'
# time = 2015052522
# ftype = 'inundation'

# network = mynet.read_graph_from_disk(path=fpath,
#                                         name=f'{prefix}_{time}_{ftype}')
# gdf_edges = ox.graph_to_gdfs(G=network, nodes=False)

# # spatial join
# sjoined_data = gpd.sjoin(left_df=bg_svi,
#                             right_df=gdf_edges,
#                             how='left')

# # count the number of roads within each block group and relate back to svi geodataframe
# count_dict = sjoined_data['GEOID'].value_counts().to_dict()
# bg_svi["count"] = bg_svi["GEOID"].apply(lambda x: count_dict.get(x))

# # count the number of roads within each block group with 0 capacity under agressive flood relationship
# subset_df = sjoined_data.loc[sjoined_data['inundation_capacity_agr'] == 0]
# count_dict = subset_df['GEOID'].value_counts().to_dict()
# bg_svi["agr_no_cap"] = bg_svi["GEOID"].apply(lambda x: count_dict.get(x))
# bg_svi['agr_no_cap'] = bg_svi['agr_no_cap'].fillna(0)

# # rotate
# envgdf = gpd.GeoDataFrame(bg_svi, geometry=svi_rot)

# # plot svi and  border
# envgdf.plot(ax=ax[1][1], column='agr_no_cap', cmap='coolwarm',
#             zorder=100, alpha=0.75, scheme='quantiles',linewidth=0)
# envgdf.boundary.plot(ax=ax[1][1],  zorder=101,alpha=0.75, edgecolor='black', linewidth=0.5)

# # Hide X and Y axes label marks
# ax[1][1].xaxis.set_tick_params(labelbottom=False)
# ax[1][1].yaxis.set_tick_params(labelleft=False)
# # Hide X and Y axes tick marks
# ax[1][1].set_xticks([])
# ax[1][1].set_yticks([])


# ####################################################################################################################################################
# # [2][1] household reliability
# figure, axe = mynet.plot_aoi(G=G_temp,
#                              res_parcels=res_parcels,
#                              resource_parcels=all_resource_parcels,
#                              rotation=-90,
#                              scalebar=False,
#                              bbox=aoi_shape,
#                              bg_water=bg_water,
#                              default=[fig, ax[2][1]])

# # determine household reliability at the appropriate time 
# time = '2015052522'
# cost_atrs = ['cost_of_fl',
#             'cost_of__1',
#             'cost_of__2',
#             'cost_of__3',
#             'cost_of__4',
#             'cost_of__5',
#             'cost_of__6',
#             'cost_of__7']

# weight = 4
# load_fpath = '/media/mdp0023/extradrive1/Data/Network_Data/Austin_North/AN_Graphs/Flow_decomp'
# time = 2015052522
# extension = '_inundation_res_parcel_flow_decomp.shp'

# res_parcels = gpd.read_file(f'{load_fpath}/{time}{extension}')
# # get geometry transformation datadata
# res_parcels_rot = gpd.GeoSeries(g for g in res_parcels['geometry'])
# res_parcels_rot = res_parcels_rot.rotate(-90, origin=center_point)

# # Reproject cost of access values into a "feasible set"
# # removing the impact of outliers by setting them all equal to the 3rd quantile plus 3*the IQR
# for cost_atr in cost_atrs:
#     # # IQR METHOD to mask impact of outliers
#     costs = sorted(res_parcels[cost_atr].tolist())
#     costs = [x for x in costs if math.isnan(x) == False]
#     q1, q3, = np.percentile(costs, [25, 75])
#     iqr = q3-q1
#     upper_bound = q3+(3*iqr)
#     res_parcels.loc[res_parcels[cost_atr] >= upper_bound, [cost_atr]] = upper_bound

# # fill column NaNs with appropriately weighted values
# max_values = res_parcels[cost_atrs].max().to_dict()
# for key in max_values:
#     max_values[key]*=weight
# res_parcels.fillna(value=max_values, inplace=True)

# # sum across the columns
# res_parcels['aggregate'] = res_parcels[cost_atrs].sum(axis=1)

# # rotate
# envgdf = gpd.GeoDataFrame(res_parcels, geometry=res_parcels_rot)
# # plot appropriate layer
# envgdf.plot(ax=ax[2][1], column='aggregate', cmap='coolwarm',
#             zorder=100, alpha=0.95, scheme='quantiles', linewidth=0)

# # Hide X and Y axes label marks
# ax[2][1].xaxis.set_tick_params(labelbottom=False)
# ax[2][1].yaxis.set_tick_params(labelleft=False)
# # Hide X and Y axes tick marks
# ax[2][1].set_xticks([])
# ax[2][1].set_yticks([])

# # adjust titles and plot annotations
# ax[0][0].set_title('Inundated Roads', style ='italic',fontsize=12, weight='bold', pad=-2)
# ax[0][1].set_title('Traffic Assignment', style='italic',fontsize=12, weight='bold', pad=-2)
# ax[1][0].set_title('Network Redundancy', style='italic', fontsize=12, weight='bold', pad=-2)
# ax[1][1].set_title('Network Reliability', style='italic', fontsize=12, weight='bold', pad=-2)
# ax[2][0].set_title('Household Redundancy', style='italic', fontsize=12, weight='bold', pad=-2)
# ax[2][1].set_title('Household Reliability', style='italic', fontsize=12, weight='bold', pad=-2)

# idx1=0
# idx2=0
# for letter in ['A','B','C','D','E','F']:
# # add a,b,c, etc.
#     ax[idx1][idx2].annotate(letter, 
#                 xy=(0.05,0.90), 
#                 xycoords='axes fraction',
#                 xytext=(0.05, 0.90),
#                 textcoords='axes fraction',
#                 va='bottom', 
#                 ha='center',
#                 fontsize=10,
#                 bbox=dict(boxstyle="Round",
#                         fc="lightgray", ec="darkgray", lw=0.5))
#     idx2+=1
#     if idx2 == 2:
#         idx1+=1
#         idx2=0


# ####################################################################################################################################################
# # add legend information

# #INUNDATION BAR
# # extract the values to plot colors to
# vals = pd.Series(nx.get_edge_attributes(G_temp_inun, 'max_inundation_mm'))
# # set the bounds of the color
# bounds=[0,10,50,100,150,200,250,300,350,400,500,600,vals.dropna().max()]
# #define the colormap
# cmap = plt.cm.Blues  

# # extract all colors from the .jet map
# cmaplist = [cmap(i) for i in range(cmap.N)]
# # force the first color entry to be clear and only show underlying roads
# cmaplist[0] = (.5, .5, .5, 0.0)

# # create the new cmap
# cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

# bounds=[0,10,150,300,600,vals.dropna().max()]
# cmap = (mpl.colors.ListedColormap(['lightgray', 'paleturquoise','deepskyblue', 'royalblue', 'mediumblue']))

# # normalize the colors and create mapping color object
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
# cm_map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
# color_series=vals.map(cm_map.to_rgba)

# cax = fig.add_axes([0.05, 0.075, 0.325, 0.025])
# cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#                     cax=cax, orientation='horizontal',
#                     label='Max Inundation on Road (cm)')
# cbar.ax.set_xticklabels(['0', '1', '15', '30', '60', '+60'])

# # NUM CARS BAR
# x0=1.25
# x1=0.8
# y0=-0.14
# y1=-0.14
# arrow = mpatches.FancyArrowPatch((x0, y0), (x1, y1),
#                                  facecolor='black', 
#                                  edgecolor='black', 
#                                  linewidth=0, 
#                                  transform=ax[2][0].transAxes,
#                                  arrowstyle='wedge,tail_width=0.8,shrink_factor=0.5',
#                                  mutation_scale=25,)
# arrow.set_clip_on(False)
# ax[2][0].add_patch(arrow)
# xy = (1, -0.35)
# xy_text = (1, -0.35)
# ax[2][0].annotate('Num. of vehicles \non road link',
#             xy=xy,
#             xycoords='axes fraction',
#             xytext=xy_text,
#             textcoords='axes fraction',
#             va='bottom',
#             ha='center',
#             fontsize=10,
#             verticalalignment='center')

# # COOLWARM BAR
# n_bins=5
# # Custom color ramp
# cmap = mpl.cm.coolwarm
# # using quantile method, determine number of bins
# bounds = np.linspace(0, 1, n_bins+1)
# np.delete(bounds, 0)
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
# cax = fig.add_axes([0.62,0.075,0.325,0.025])
# cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#                     cax=cax, orientation='horizontal',
#                     label="Metric Percentile")
# cbar.ax.set_xticklabels(['0','20','40','60','80','100'])


# # add space for legends
# fig.subplots_adjust(top=0.975,
#                     bottom=0.115,
#                     left=0.025,
#                     right=0.975,
#                     hspace=0.1,
#                     wspace=0.0)




plt.show()
