import network_exploration_stuff as mynet

import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
import geopandas as gpd
import rasterio as rio
import numpy as np


import network_exploration_stuff as mynet

# VARIABLES USED #############################################################
# file path
path = "/home/mdp0023/Documents/Codes_Projects/\
network_analysis/Network_Testing_Data"
image_path = "/home/mdp0023/Documents/Codes_Projects/\
network_analysis/Poster_Graphics"
inset_path = "/home/mdp0023/Documents/Codes_Projects/network_analysis/bboxes"
# AOI without buffer
aoi_area = f'{path}/Neighborhood_Network_AOI.shp'
# AOI with buffer
aoi_buffer = f'{path}/Neighborhood_Network_AOI_Buf_1km.shp'
# centroids of res parcels and food marts
res_points_loc = f'{path}/Residential_Parcels_Points_Network_AOI.shp'
food_points_loc = f'{path}/Food_Marts_Points_Network_AOI.shp'
# shapefiles of res parcels and food marts
res_parcels = f'{path}/Residential_Parcels_Network_AOI.shp'
food_parcels = f'{path}/Food_Marts_Network_AOI.shp'
# path for inundation
inundation = f'{path}/Network_Inun.tif'
raster = rio.open(inundation)

# LOADING WORK ###############################################################
G = mynet.shape_2_graph(source=aoi_buffer)
G = mynet.rename(G=G)


# LOAD OTHER DATA ############################################################
# shapefile centroids of residental plots
res_points = gpd.read_file(res_points_loc)
# shapefile of res parcels
res_locs = gpd.read_file(res_parcels)
# shapefile centroids of 3 foodmart plots
food_points = gpd.read_file(food_points_loc)
# shapefile of food mart parcels
food_locs = gpd.read_file(food_parcels)
# shapefile of area of interest
aoi_area = gpd.read_file(aoi_area)

# ensure proper projection
G = ox.projection.project_graph(G, to_crs=32614)
res_locs = res_locs.to_crs(epsg=32614)
food_locs = food_locs.to_crs(epsg=32614)
res_points = res_points.to_crs(epsg=32614)
food_points = food_points.to_crs(epsg=32614)
aoi_area = aoi_area.to_crs(epsg=32614)



flow_dict, flow_cost = mynet.min_cost_flow_parcels(G=G,
                                                   res_points=res_points,
                                                   dest_points=food_points)


################ Develop flow decomposition algorithm function
# Be able to take in some sort of flow dictionary and decompose it into a list of paths




def flow_decomposition(G=G,
                       res_points=res_points,
                       dest_points=food_points, G_demand='demand', G_capacity='capacity', G_weight='travel_time'):

    # COPIED FROM MIN COST FLOW - CREATES THE NETWORK WE NEED TO MANIPULATE

    # Travel times must be whole numbers -  round values if not whole numbers
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])
    # find the nearest node for each parcel (origins and destinations)
    # create empty lists for nodes
    origins = []
    destinations = []
    # append origins to have nearest intersection of all res_points
    for x in list(range(0, len(res_points))):
        res_point = res_points.iloc[[x]]
        res_lon = res_point['geometry'].x.iloc[0]
        res_lat = res_point['geometry'].y.iloc[0]
        origin = ox.distance.nearest_nodes(G, res_lon, res_lat)
        origins.append(origin)
    # append destinations to have nearest intersections of all dest_points
    for x in list(range(0, len(dest_points))):
        des_point = dest_points.iloc[[x]]
        des_lon = des_point['geometry'].x.iloc[0]
        des_lat = des_point['geometry'].y.iloc[0]
        destination = ox.distance.nearest_nodes(G, des_lon, des_lat)
        destinations.append(destination)
    # find the number of unique origin nodes: count*-1 = demand
    unique_origin_nodes = np.unique(origins)
    unique_origin_nodes_counts = {}
    for unique_node in unique_origin_nodes:
        count = np.count_nonzero(origins == unique_node)
        unique_origin_nodes_counts[unique_node] = {G_demand: -1*count}
    # determine the unique destination nodes
    unique_dest_nodes = np.unique(destinations)
    # set demand at these unique dest nodes to 0 (SEE TODO)
    for x in unique_dest_nodes:
        unique_origin_nodes_counts[x] = {G_demand: 0}
        unique_origin_nodes = np.delete(
            unique_origin_nodes, np.where(unique_origin_nodes==x))
    # Graph must be in digraph format: convert multidigraph to digraph
    # TODO: Need to see if this is having an impact on the output
    G = nx.DiGraph(G)
    # add source information (the negative demand value prev. calculated)
    nx.set_node_attributes(G, unique_origin_nodes_counts)
    # Calculate the positive demand: the sink demand
    demand = 0

    for x in unique_origin_nodes_counts:
        demand += unique_origin_nodes_counts[x][G_demand]
    positive_demand = demand*-1

    # create artificial node with demand calc.
    #    All sinks will go to this demand in order to balance equation
    G.add_nodes_from([(99999999, {G_demand: positive_demand})])
    # add edges from sinks to this artificial node

    for x in unique_dest_nodes:
        G.add_edge(x, 99999999)
        G[x][99999999][G_weight]=0
    flow_dictionary = nx.min_cost_flow(G, demand=G_demand, weight=G_weight, capacity=G_capacity)

    cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
    print(f"Cost of flow:{cost_of_flow}")

    # calculate flow across the entire network (i.e., each instance of a unit of flow)
    total_flow = 0
    for x in flow_dictionary.values():
        total_flow += sum(x.values())
    print(f"Total flow on network:{total_flow}")


    decomposed_paths ={}
    for origin in unique_origin_nodes:
        path =  ox.distance.shortest_path(G,
                                  origin,
                                         99999999,
                                         weight=G_weight)
        # TODO: I believe this should be EITHER: source amount, or minimum flow along an edge along path, whichever is SMALLEST
        '''
        This algorithm needs improvement. The issue is that multiple decompositions are feasible, so need to find the justifiable option.
        I BELIEVE, what should be done, is a "greedy" algorithm, successively looking at the shortest path from a source to sink
        This really will only change the result when flows are split - flow from a single source is going in different directions

        Also important to note that due to nature of network, there should never be any cycles therefore it should always decompose to series of paths

        Here is how this needs to work:

        1. find the shortest paths from a super source to a super sink
            Super source connected to actual sources, capacity on edge is original source
            Might consider finding shorest paths along an intermediate network that only contains pathways with flows on them from flow dictinoary
        2. Based on flow dictionary from min-cost-flow algorithm, determine amount of flow to reduce by
            What this means: reduce flow along these edge paths by the smallest edge flow
        3. Go back to step 1, find new shorest path and recalculate
            
        '''
        flow_to_eliminate = list(unique_origin_nodes_counts[origin].values())[0]
        cost = nx.path_weight(G, path, weight=G_weight)
        decomposed_paths[f'{origin}']={'sink':path[-2], 'path': path, 'flow': flow_to_eliminate*-1, 'Cost Per Flow': cost}



    output='test_output'
    return output

test_function = flow_decomposition()
print(test_function)
