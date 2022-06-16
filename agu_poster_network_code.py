import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
import geopandas as gpd
import rasterio as rio


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
# G = shape_2_graph(source=aoi_buffer)
# save_2_disk(G=G, path=path)
G = mynet.read_graph_from_disk(path=path, name='AOI_Graph')
G = mynet.rename(G=G)
inundated_G = mynet.read_graph_from_disk(path=path, name='AOI_Graph_Inundated')
inundated_G = mynet.rename(G=inundated_G)

# parallel_edges(G=G)

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

# #RANDOM SHORTEST PATH #######################################################
# routes = mynet.random_shortest_path(G=G,
#                               res_points=res_points,
#                               dest_points=food_points,
#                               plot=False)


# MIN COST FLOW OF PARCELS ####################################################
flow_dict, flow_cost = mynet.min_cost_flow_parcels(G=G,
                                             res_points=res_points,
                                             dest_points=food_points)


# MAX FLOW OF PARCELS ########################################################
#   therefore we can calculate demand
# If we do not know if everyone can reach a destination by either being cut off
#   or some other capacity issue, we must calculate the max flow possible
# flow_dict, flow_cost, max_flow, access = mynet.max_flow_parcels(G=G,
#                                                           res_points=res_points,
#                                                           dest_points=food_points)
# print(flow_cost)

# PLOT BUILDING PARCELS AND ROAD NETWORK TOGETHER ############################
# create instance of figure and axis
# figure = plot_aoi(G=G, res_parcels=res_locs, resource_parcels=food_locs)


# SUMMARY FUNCTION ###########################################################
# I want to be able to print a convenient summary table of the types of edges
# and nodes that are in the data
# num_edges, num_nodes, fig = summary_function(G=G)


# INUNDATED NETWORK ########################################################
#create a network with increased travel times and decreased capacities
# inundated_G = mynet.inundate_network(G=G,
#                                CRS=32614,
#                                path=path,
#                                inundation=inundation)

inundated_G = mynet.read_graph_from_disk(path=path, name='AOI_Graph_Inundated')
inundated_G = ox.projection.project_graph(inundated_G, to_crs=32614)

#INCREASE IN COST FROM INUNDATION ###########################################
# cost and flow in dry network
dry_flow_dic, dry_cost_of_flow, dry_max_flow, dry_access = mynet.max_flow_parcels(
    G=G,
    res_points=res_points,
    dest_points=food_points,
    G_capacity='capacity',
    G_weight='travel_time')
print(f"Dry cost of flow: {dry_cost_of_flow}")
print(f"Dry maximum flow: {dry_max_flow}")

# inundated flow
wet_flow_dic, wet_cost_of_flow, wet_max_flow, wet_access = mynet.max_flow_parcels(
    G=inundated_G,
    res_points=res_points,
    dest_points=food_points,
    G_capacity='inundation_capacity',
    G_weight='inundation_travel_time')
print(f"Wet cost of flow: {wet_cost_of_flow}")
print(f"Wet maximum flow: {wet_max_flow}")

##


# # DETERMINE WHAT RESIDENTS LOSE ACCESS TO RESOURCE #########################
# # determine locations of all of the dest_points by iterating
# dest_points = food_points
# dest_locs = []
# lst = list(range(0, len(dest_points)))
# for loc in lst:
#     dest_point = dest_points.iloc[[loc]]
#     dest_lon = dest_point['geometry'].x.iloc[0]
#     dest_lat = dest_point['geometry'].y.iloc[0]
#     dest_locs.append((dest_lon, dest_lat))

# # determine the nearest node for all destinations
# destinations = []
# for dest_loc in dest_locs:
#     destination = ox.distance.nearest_nodes(G, dest_loc[0], dest_loc[1])
#     destinations.append(destination)

# # create a new graph that removes any link that has a capacity of 0
# subgraph = nx.MultiDiGraph([(u, v, d) for u, v, d in inundated_G.edges(
#     data=True) if d['inundation_capacity'] > 0])
# subgraph.add_nodes_from(inundated_G.nodes())

# # create a blank dataframe of residents w/o access to resource
# data = []
# # iterate through res_points
# for idx, res_point in res_points.iterrows():
#     # extract lat and lon

#     res_lon = res_point['geometry'].x
#     res_lat = res_point['geometry'].y
#     res_loc = res_lon, res_lat
#     # find the nearest node for this origin
#     origin = ox.distance.nearest_nodes(G, res_loc[0], res_loc[1])

#     # Find  path from origin to all destinations (if it exists)

#     routes = []
#     for destination in destinations:
#         routes.append(nx.has_path(G=subgraph,
#                                   source=origin,
#                                   target=destination))
#     routes = [not a for a in routes]
#     if all(routes) is True:
#         data.append(res_point)

# no_access_pt = gpd.GeoDataFrame(data)
# no_access_propids = no_access_pt['PROP_ID'].unique()
# no_access_res_parcels = res_locs[res_locs['PROP_ID'].isin(no_access_propids)]
# access_parcels = res_locs[~res_locs['PROP_ID'].isin(no_access_propids)]
# no_access_res_parcels.to_file(f'{path}/parcels_wo_access.shp')
# access_parcels.to_file(f'{path}/parcels_w_access.shp')

# ##############################################################################
# CREATE SOME MAPS TO VISUALIZE IMPACT OF INUNDATION ####################
# first need to open clean copies of graphs that do not have artifical nodes
G_map = mynet.read_graph_from_disk(path=path, name='AOI_Graph')
G_map = ox.projection.project_graph(G_map, to_crs=32614)
G_map = mynet.rename(G=G_map)
inundated_G_map = mynet.read_graph_from_disk(path=path, name='AOI_Graph_Inundated')

# relate flow dictionary back to INUNDATED graph for plotting purposes
inundated_G_map = nx.DiGraph(inundated_G_map)
# relate
for edge in inundated_G_map.edges:
    values = {(edge[0], edge[1]): {'inundation_flow':
                                   wet_flow_dic[edge[0]][edge[1]]}}
    nx.set_edge_attributes(G=inundated_G_map, values=values)
inundated_G_map = nx.MultiDiGraph(inundated_G_map)

#relate flow dictionary back to DRY graph for plotting purposes
G_map = nx.DiGraph(G_map)
# relate
for edge in G_map.edges:
    values = {(edge[0], edge[1]): {'dry_flow':
                                   dry_flow_dic[edge[0]][edge[1]]}}
    nx.set_edge_attributes(G=G_map, values=values)
G_map = nx.MultiDiGraph(G_map)


# create some inundated plots
# inundated flow
# plot = mynet.plot_aoi(G=inundated_G_map,
#                       res_parcels=res_locs,
#                       resource_parcels=food_locs,
#                       edge_width='inundation_flow',
#                       bbox=aoi_area,
#                       loss_access_parcels=no_access_res_parcels,
#                       insets=[f"{inset_path}/bbox2.shp",
#                               f"{inset_path}/bbox3.shp",
#                               f"{inset_path}/bbox4.shp"],
#                       scalebar=True,
#                       save_loc=f"{image_path}/inundated_flow.pdf")

# non-inundated flow plot for paper
plot = mynet.plot_aoi(G=G_map,
                      res_parcels=res_locs,
                      resource_parcels=food_locs,
                      edge_width='dry_flow',
                      bbox=aoi_area,
                      scalebar=True,
                      save_loc=f"{image_path}/non_inundated_flow.pdf")
plt.show()
# # inundation map
# plot = mynet.plot_aoi(G=G_map,
#                       res_parcels=res_locs,
#                       resource_parcels=food_locs,
#                       bbox=aoi_area,
#                       inundation=raster,
#                       insets=[f"{inset_path}/bbox1.shp"],
#                       scalebar=True,
#                       save_loc=f"{image_path}/inundation.pdf")

# # inset maps
# # inset 1
# # plot = plot_aoi(G=G_map,
# #                 res_parcels=res_locs,
# #                 resource_parcels=food_locs,
# #                 bbox=gpd.read_file(f"{inset_path}/bbox1.shp"),
# #                 inundation=raster,
# #                 save_loc=f"{image_path}/inset1.pdf")
# # # inset 2
# # plot = plot_aoi(G=inundated_G_map,
# #                 res_parcels=res_locs,
# #                 resource_parcels=food_locs,
# #                 edge_width='inundation_flow',
# #                 bbox=gpd.read_file(f"{inset_path}/bbox2.shp"),
# #                 loss_access_parcels=no_access_res_parcels,
# #                 save_loc=f"{image_path}/inset2.pdf")
# # # inset 3
# # plot = plot_aoi(G=inundated_G_map,
# #                 res_parcels=res_locs,
# #                 resource_parcels=food_locs,
# #                 edge_width='inundation_flow',
# #                 bbox=gpd.read_file(f"{inset_path}/bbox3.shp"),
# #                 loss_access_parcels=no_access_res_parcels,
# #                 save_loc=f"{image_path}/inset3.pdf")
# # # inset 4
# # plot = plot_aoi(G=inundated_G_map,
# #                 res_parcels=res_locs,
# #                 resource_parcels=food_locs,
# #                 edge_width='inundation_flow',
# #                 bbox=gpd.read_file(f"{inset_path}/bbox4.shp"),
# #                 loss_access_parcels=no_access_res_parcels,
# #                 save_loc=f"{image_path}/inset4.pdf")
# plt.show()

# #############################################################################
# #############################################################################

# # OTHER FUNCTIONS TO TRY OUT AT SOME POINT
# # REMOVE PERCENTAGE OF LINKS RANDOMLY #########################################
# # See how number of people access resources changes based on change in access

# # WHAT LINK(S) HAVE GREATEST IMPACT ON COST ###################################
# # of the links with flow, the removal of what link(s) has biggest impact
# #   on cost
# # print(flow_dict)
# # links_w_flow = []
# # for u, vs in flow_dict.items():
# #     for v in vs:
# #         if flow_dict[u][v] > 0:
# #             links_w_flow.append([u, v])
# # print(links_w_flow)
# # print(len(links_w_flow))
