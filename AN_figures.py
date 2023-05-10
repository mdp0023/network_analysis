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
import network_exploration_stuff as mynet

# this is the code where I will do scratch work for the Austin North Network FIGURE MAKING

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
