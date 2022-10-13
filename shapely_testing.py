# I am using this to test some different python methods to find multiple nearest nodes 
# I don't want to just do a buffer around centroid or a buffer around the building for a variety of reasons:
    # too big or too small depending on buffer size - highly dependent on rural vs. urban
    # if too large could cross inundated area - not good
    # if too small, won't have any advantage

# THEREFORE: Idea is to find the corners of parcel boundaries (or subset of simplified verticies) and find the nearest points 
    # I can simplify boundaries less to include more search areas, and can also include the centroid 

import network_exploration_stuff as mynet
import matplotlib.pyplot as plt
import geopandas as gpd
import networkx as nx
import pandas as pd
import osmnx as ox
import numpy as np
from shapely.geometry import Polygon
import logging

# ignore shapely deprectiation warnings
logging.captureWarnings(True)


# INPUT VARIABLES
# Folder path
folder_path = '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North'

# Centroids of res parcels and grocery stores
res_parcel_loc = f'{folder_path}/AN_Residential_Parcel_Shapefiles/AN_Residential_Parcels.shp'
res_points_loc = f'{folder_path}/AN_Residential_Parcel_Shapefiles/AN_Residential_Parcels_Points.shp'
food_parcel_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_Supermarket_Parcel.shp'
food_points_loc = f'{folder_path}/AN_Resource_Parcel_Shapefiles/AN_Supermarket_Parcel_Points.shp'

# read files
food_parcels = gpd.read_file(food_parcel_loc)
food_points = gpd.read_file(food_points_loc)
res_parcels = gpd.read_file(res_parcel_loc)
res_points = gpd.read_file(res_points_loc)


# bounadry and boundary with buffer
aoi_area = f'{folder_path}/AN_Boundary/AN_Boundary.shp'
aoi_buffer = f'{folder_path}/AN_Boundary/AN_Boundary_3km.shp'




######################################################################################


G = mynet.shape_2_graph(aoi_buffer)
G = mynet.rename(G=G)

# # simplify food_parcels
food_parcels_simp = food_parcels.simplify(3, preserve_topology=False)

# account for multipolygons by creating convex hulls of each parcel
convex_hull = food_parcels_simp.convex_hull
food_parcels1 = gpd.GeoDataFrame({'geometry': convex_hull})


coords = [list(food_parcels1.geometry.exterior[row_id].coords) for row_id in range(food_parcels1.shape[0])]

# for coord in coords:
#     print(len(coord))
#     print(coord)
#     x = [i for i,j in coord]
#     y = [j for i,j in coord]
#     print(x)
#     print(y)

destinations_int = []
destinations_str = []
for idx, coord in enumerate(coords):
    longs = [i for i, j in coord]
    lats = [j for i, j in coord]
    dests = ox.distance.nearest_nodes(G, longs, lats)
    destinations_int.append(dests)
    destinations_str.append(' '.join(str(x) for x in dests))



