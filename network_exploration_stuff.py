# network exploration stuff
# using a number of websites for background information
# https://automating-gis-processes.github.io/2017/lessons/L7/network-analysis.html
# https://networkx.org/documentation/stable/tutorial.html
# https://geoffboeing.com/2016/11/osmnx-python-street-networks/

import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import numpy as np
from scipy.stats import iqr
import rasterio as rio
from rasterio.plot import show
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap
# import momepy as mpy
from rasterstats import zonal_stats
from matplotlib_scalebar.scalebar import ScaleBar
import random
import copy

# set some panda display options
pd.set_option('display.max_columns', None)

# todo: Is renaming the nodes from 1 to N really necessary?
# TODO: Nearest node - some of the closest origin nodes are the destination
#   nodes
#   Solution right now is that after sources are marked, dest nodes are set to
#   have a 0 demand (justificaiton: walking distance?)
# TODO: Add road capacities in min cost flow calculation
# TODO: Need to use different method to determine random node variable name
#   in the min cost flow function
# TODO: Min cost flow and max cost flow functions have some repeat code
#   I could condense a lot of the work
# TODO: travel disruption - right now using arbitrary disruption
#   for example - flood classificaiton 1 is set to reduce speed by 75%
#                 flood classificaiton 2 is set to reduce speed by 95%
#                 flood classificaiton 3 and higher reduces capacity to 0
# TODO: change inundate network function to be more broad
#   see above TODO, just needs to have more inputs
#   this includes flood classificaitons and impacts on speed, etc.
# TODO: Use kwargs methodology to create new edges with appropriate attributes
#   see # TODO: KWARGS below in code for example
# TODO: Move CRS check into functions
# TODO: move digraph conversion into function
#   There are a lot of cross overs between digraph and multigraph, need to
#   examine this in more detail
# TODO: In plot aoi - allow it to accept a list of dest parcels
#   i.e., food and gas locations
# TODO: When it comes to plotting, doing the different algorithms adds nodes
#   and edges. Therefore, have to plot the un-manipulated networks
#   this also means that they should be loaded seperately
# TODO: read_graph_from_disk - problem with preserving dtypes so I hardcoded a
#   temporary fix - right now it goes through the attributes and determines if
#    it could be converted to a float - if it can it does
#   this involves manually skipping oneway as well (i think b/c its bool)
#   this was needed b/c  inundation_G load not preserve datatypes
# TODO: rework methodology for determing what residents lose access to resource
#   right now doing brute force calculation to determine which house one by one
#   should be able to relate maximum flow results to this but ran out of time
# TODO: convert parcel identificaiton (access vs. no access) to function



# Functions to examine:
    #https://osmnx.readthedocs.io/en/stable/osmnx.html?highlight=interpolate#osmnx.utils_geo.interpolate_points 
    #https://stackoverflow.com/questions/64104884/osmnx-project-point-to-street-segments 
    #https://github.com/gboeing/osmnx/issues/269 

# Possible explorations?
# Survey informed model - using survey/census data to determine mvmnt patterns?
#   Maybe stochastic methods of movements
#   e.g., 50% of the people in this area only ride public transportation

# Assumptions that need to be documented
# resources ae connected only 1 of the nearest intersections
#   realistically, there are multiple intersections for many resources

# Snapping parcels to any intersections may snap to the wrong ones
#   if a place is close to freeway interchange, might snap there instead of
#   actual road it neighbors
#   Possible solution: when snapping to roads, ignore roads above a specific
#   classification (e.g., freeways)


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

# THE FUNCTIONS #############################################################


def shape_2_graph(source=''):
    '''
    Extracts road network from shapefile.

    This function takes an input path to a .shp file (assuming the associated .shp, .shx, .dbf, etc. are all in the same folder)
    and extracts/returns the OpenStreetMap road network. This function automatically pulls the speed and travel time information,
    which are necessary components in further calculations. 

    **Important Notes**
    
    Ensure that the shapefile area has a sufficent buffer to capture enough of the road network.
    For example, if the AOI is a block group, residents might require driving outside bounadry to get to areas within
    and excluding roadnetworks, espcially those that border the AOI, can lead to errors.

    .. TODO::

        Only potential improvement is to include function to add a buffer to AOI instead of requiring that the user automatically do this,
        but this is not the most meaningful improvement


    :param source: Path to .shp file of the area of interest
    :type source: str
    :return: **G**, the resultant road network from OpenStreetMap
    :rtype: networkx.MultiDiGraph
    
    '''

    G = ox.graph_from_polygon(
        gpd.read_file(source)['geometry'][0],
        network_type='drive')

    ox.add_edge_speeds(G)
    ox.add_edge_travel_times(G)

    return(G)


def save_2_disk(G='', path='', name='AOI_Graph'):
    '''
    Save network as graphxml file.

    Once saved, the graph can be reloaded in for future analysis. The User avoids having to redownload 
    the same graph multiple times.
    
    :param G: Input graph to be saved
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param path: Path to folder where graph will be saved
    :type path: str
    :param name: Name the graph will be saved as. *Default='AOI_Graph*
    :type name: str
    

    '''
    ox.io.save_graphml(
        G, f'{path}/{name}')


def read_graph_from_disk(path='', name='AOI_Graph'):
    '''
    Open a saved graphxml file.

    This works in conjunction with :func:`save_2_disk`, opening a saved graphxml file.

    :param path: Path to folder where graphxml file is saved.
    :type path: str
    :param name: Name of the graphxml file. *Default='AOI_Graph'*
    :type name: str
    :return: **G**, the opened graph network.
    :rtype: networkx.Graph [Multi, MultiDi, Di]
    '''
    
    G = ox.io.load_graphml(
        f'{path}/{name}')
    # Having an issue with preserving datatypes so going to try and convert any attribute to float that can be a float
    nodes, edges = ox.graph_to_gdfs(G)
    list_floats = {}
    for edge in edges.keys():
        try:
            edges[edge].astype(float)
            list_floats[edge] = float
        except ValueError:
            pass
        except TypeError:
            pass
    list_floats.pop('oneway')

    G = ox.io.load_graphml(
        f'{path}/{name}', edge_dtypes=list_floats)
    return(G)


def rename(G=''):
    '''
    Rename nodes from 1 to N, the total number of nodes.

    This function is useful after extracting network using OSMnx, and makes referencing nodes/paths easier.
    It serves no specific purpose but I think it is still useful.

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :return: **G**, renamed graph network
    :rtype: networkx.Graph [Multi, MultiDi, Di] 
    '''
   
    # Number of nodes
    num_nodes = len(G)
    # list from 1 to [num_nodes]
    lst = list(range(1, num_nodes+1))
    # Create a dictionary (value for each key is None) from sorted G nodes
    mapping = dict.fromkeys(sorted(G))
    # iterate through this dictionary (with enumerate to have count var.)
    for idx, key in enumerate(mapping):
        # set each key value (previously None) to new value
        mapping[key] = lst[idx]
    # nx.relabel_nodes relabels the nodes in G based on mapping dict
    G = nx.relabel_nodes(G, mapping)
    return(G)


def parallel_edges(G=''):
    '''
    Find parallel edges in graph network.

    A lot of network algorithms in NetworkX require that networks be DiGraphs, or directional graphs.
    When extracting road networks with OSMnx, they are often returned as MultiDiGraphs, meaning they could potentially have 
    parallel edges. This function determines if parallel edges exists, and converting between MultiDiGraphs to DiGraphs could
    lead to some sort of error. 

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :returns: 
        - **Bool**, *True* if parallel edges exist, otherwise *False*
        - **parallel_edges**, list of parallel edges, *[[u,v,num_edges],...]*

    :rtype: tuple

    '''
   
    # creat empty list of parallel edges
    parallel_edges = []
    # Number of nodes
    num_nodes = len(G)
    # list from 1 to [num_nodes]
    lst = list(range(1, num_nodes+1))
    # iterate through all possible combinations of (u,v)
    for u in lst:
        for v in lst:
            # extract number of edges between nodes u and v
            num_edges = G.number_of_edges(u=u, v=v)
            # if the number of edges is greater than 2, we have parallel edges
            if num_edges >= 2:
                # append list, with the two nodes and number of edges
                parallel_edge = [u, v, num_edges]
                parallel_edges.append(parallel_edge)
    if len(parallel_edges) >= 1:
        return True, parallel_edges
    if len(parallel_edges) < 1:
        return False, [0]


def nearest_nodes(G='', res_points='', dest_points='', G_demand='demand'):
    '''
    Return Graph with demand attribute based on snapping sources/sinks to nearest nodes (intersections).
    
    The purpose of this function is take a series of input locations and determine the nearest nodes (intersections) within a Graph network. It takes this information to create source and sink information. 

    **Important Note**: The output graph has an attribute labeled G_demand, which shows the number of closest residential parcels to each unique intersection. Because these are considered sources, they have a negative demand. Sink locations, or the intersections that are closest to the the dest_points, will not have a demand value (G_demand == 0) because we do not know how much flow is going to that node until after a min-cost flow algorithm is run. I.e., to route properly, we end up creating an artifical sink and then decompose the flow to determine how much flow is going to each destination.

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the residential parcels are to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param G_demand: name of attribuet in G refering to node demand, *Default='demand*
    :type G_demand: string
    :returns: 
        - **G**, *networkx.DiGraph* with G_demand attribute showing source values (negative demand)
        - **unique_origin_nodes**, *lst* of unique origin nodes
        - **unique_dest_nodes**, *lst* of unqiue destination nodes
        - **positive_demand**, *int* of the total positive demand across the network. **Does not include the residents that share a closest intersection with a resource parcel.**
        - **Shared_nodes**, *lst* of nodes that res and destination parcels share same nearest 
        - **res_points**, *geopandas.GeoDataFrame*, residential points with appended attribute of 'nearest_node'
        - **dest_points**, *geopandas.GeoDataFrame*, destination/sink points with appended attribute of 'nearest_node'
    :rtype: tuple


    '''

    # create empty lists for nodes of origins and destinations
    origins = []
    destinations = []
    # append origins to have nearest node (intersection) of all res_points
    for x in list(range(0, len(res_points))):
        res_point = res_points.iloc[[x]]
        res_lon = res_point['geometry'].x.iloc[0]
        res_lat = res_point['geometry'].y.iloc[0]
        origin = ox.distance.nearest_nodes(G, res_lon, res_lat)
        origins.append(origin)
        # add new column to res_points dataframe with nearest node
        res_points.loc[res_points.index[x],'nearest_node'] = int(origin)
        
    # append destinations to have nearest node (intersections) of all dest_points
    for x in list(range(0, len(dest_points))):
        des_point = dest_points.iloc[[x]]
        des_lon = des_point['geometry'].x.iloc[0]
        des_lat = des_point['geometry'].y.iloc[0]
        destination = ox.distance.nearest_nodes(G, des_lon, des_lat)
        destinations.append(destination)
        # add new column to des_points dataframe with nearest node
        dest_points.loc[res_points.index[x], 'nearest_node'] = int(destination)

    # Creat the demand attribute and set all equal to 0
    #TODO: make sure this has no impact on other functions
        # it's not necessary - but is meant to keep values for every attribute instead of having missing data/holes in attributes
        # THIS IS IMPACTING G ELSEWHERE AS WELL - MIGHT BE MAJOR ERROR
    nx.set_node_attributes(G, values=0, name=G_demand)

    # Create list of unique origins and destinations
    unique_origin_nodes = np.unique(origins)
    unique_dest_nodes = np.unique(destinations)
    # Based on unique origins, determine negative demands (source values)
    unique_origin_nodes_counts = {}
    for unique_node in unique_origin_nodes:
        count = np.count_nonzero(origins == unique_node)
        unique_origin_nodes_counts[unique_node] = {G_demand: -1*count}
    # set demand at  unique dest nodes to 0 (SEE TODO) - if a source and sink share an intersection, then we assume it is within walking distance and is always accessible
    shared_nodes = []
    for x in unique_dest_nodes:
        unique_origin_nodes_counts[x] = {G_demand: 0}
        # Remove from the unique origin nodes because no longer an origin
        unique_origin_nodes = np.delete(
            unique_origin_nodes, np.where(unique_origin_nodes==x))
        # Create a list of shared nodes
        shared_nodes.append(x)

    # Convert graph to digraph format
    # TODO: Need to see if this is having an impact on the output
    G = nx.DiGraph(G)
    # add source information (the negative demand value prev. calculated)
    nx.set_node_attributes(G, unique_origin_nodes_counts)

    # Calculate the positive demand: the sink demand
    demand = 0
    for x in unique_origin_nodes_counts:
        demand += unique_origin_nodes_counts[x][G_demand]
    positive_demand = demand*-1

    return(G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points)


def random_shortest_path(G='', res_points='', dest_points='', plot=False):
    '''
    Shortest path between a random residential parcel and all of a given resource.

    This function takes a set of residential points, chooses one at random, and then 
    routes it to every destination point (i.e., some given resource). This function is capable of producing
    a figure that also shows how that residential parcel would be routed to each of the destination points.
    The OSMnx function used in this function is able to accept MultiDiGraphs, unlike some others used later on.
    
    **Important Notes** 
    
    This function finds the shortest path from the closest node (i.e., intersection) of a residential parcel and the 
    closest node (i.e., intersection) for each destination. **It is worth considering changing this methodology**, adding nodes
    and edges from the points themselves, but this requires further exploration.


    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the random residential parcel is to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param plot: Whether or not a plot should be created, *Default=False*
    :type plot: bool    
    :returns: 
        - **routes**, list of shortest paths, each path is a list of node IDs.
        - **figure**, matplotlib figure object *(optional)*
    :rtype: list or tuple

    '''
  
    # select a random residential parcel, extract lon and lat variable
    res_point = res_points.sample()
    res_lon = res_point['geometry'].x.iloc[0]
    res_lat = res_point['geometry'].y.iloc[0]
    res_loc = res_lon, res_lat

    # determine locations of all of the dest_points by iterating
    dest_locs = []
    lst = list(range(0, len(dest_points)))
    for loc in lst:
        dest_point = dest_points.iloc[[loc]]
        dest_lon = dest_point['geometry'].x.iloc[0]
        dest_lat = dest_point['geometry'].y.iloc[0]
        dest_locs.append((dest_lon, dest_lat))

    # find the nearest node to the origin & all the destinations
    origin = ox.distance.nearest_nodes(G, res_loc[0], res_loc[1])
    destinations = []
    for dest_loc in dest_locs:
        destination = ox.distance.nearest_nodes(G, dest_loc[0], dest_loc[1])
        destinations.append(destination)

    # find shortest path from origin to all destinations
    # Cost (WEIGHT) is selected as travel time
    routes = []
    for destination in destinations:
        routes.append(ox.distance.shortest_path(G,
                                                origin,
                                                destination,
                                                weight='travel_time'))

    if plot is True:
        fig, ax = ox.plot_graph_routes(G, routes, route_colors=['r', 'b', 'g'],
                                       route_linewidth=6,
                                       orig_dest_size=0,
                                       bgcolor='k')
        return routes, fig
    else:
        return routes


def min_cost_flow_parcels(G='', res_points='', dest_points='', G_demand='demand', G_capacity='capacity', G_weight='travel_time'):
    '''
    Find the minimum cost flow of all residential parcels to a given resource.

    Given a set of sources (residential parcels), the minimum cost flow of all sources to sinks (given resource type or destinations)
    is determined. Sources send demand, and are therefore given a negative value. Sinks recieve demand, and are therefore given a positive 
    value. This notation is derived from the NetworkX package.

    Due to snapping methodology (*see important notes below*), source demand is the sum of the nearest parcels. If the nearest source and sink 
    share the same node (i.e., a residential parcel and a grocery store share the same intersection), demand is set to zero. 
    Justification for this is if you share the same nearest intersection then the resource is within walking distince.

    Because we do not know how much flow goes to each resource destination, each destination point is connected to the same artifical sink,
    which is given the total positive demand.

    If this function returns an error, than every residential parcel cannot reach the resource.

    **Important Notes**
    
    Similar to :func:`random_shortest_path`, this function finds the shortest path from the closest node 
    (i.e., intersection) of residential parcels and the closest node (i.e., intersection) for each destination. 
    **It is worth considering changing this methodology**, adding nodes and edges from the points themselves, 
    but this requires further exploration.

    This function also does not fact in limits to a resource. In other words, there is no limit on the number of residential parcels
    that can go to any single resource location. **This feature also needs further exploration.**

    :param G: Graph network. Will be converted to DiGraph in function
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the random residential parcel is to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param G_demand: name of attribute in G refering to node demands, *Default='demand'* 
    :type G_demand: string
    :param G_capacity: name of attribute in G refering to edge capacities, *Default='capacity'* 
    :type G_capacity: string
    :param G_weight: name of attribute in G refering to edge weights, *Default='travel_time'* 
    :type G_weight: string  
    :returns: 
        - **flow_dictionary**, dictionary of dictionaries keyed by nodes for edge flows
        - **cost_of_flow**, integer, total cost of all flow. If weight/cost is travel time, this is the total time for everyone to reach the resource (seconds).
    :rtype: tuple
    :Raises: nx.NetworkXUnfeasible if all demand cannot be satisfied, i.e., all residential parcels cannot reach the resource.


    ''' 
    # Travel times must be whole numbers -  round values if not whole numbers
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])
    
    # Find the source and sink demands and append to graph G
    G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points = nearest_nodes(
        G=G, res_points=res_points, dest_points=dest_points, G_demand=G_demand)


    # create artificial sink node with calculated total demand
    #    All sinks will go to this demand in order to balance equation
    G.add_nodes_from([(99999999, {G_demand: positive_demand})])
    # add edges from sinks to this artificial node
    for x in unique_dest_nodes:
        kwargs = {f"{G_weight}": 0}
        G.add_edge(x, 99999999, **kwargs)


    # run the min_cost_flow function to retrieve FlowDict
    try:
        flow_dictionary = nx.min_cost_flow(
            G, demand=G_demand, weight=G_weight, capacity=G_capacity)
        # run the cost_of_flow function to retrieve total cost
        cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
        return flow_dictionary, cost_of_flow
    except nx.NetworkXUnfeasible:
        return None, None


def max_flow_parcels(G='', res_points='', dest_points='', G_capacity='capacity', G_weight='travel_time', G_demand='demand'):
    '''
    Find maximum flow with minimum cost of a network between residential parcels and a given resource.
    
    Unlike :func:`min_cost_flow_parcels`, this function calculates the maximum flow from all of the residential parcels to a given resource.
    Meaning, this function is independent of whether or not every residential parcel can actually reach a resource.
    This function uses the same snapping methodology as :func:`min_cost_flow_parcels`, see there for comments on how this works.
    
    This function creates both an artifical source and sink between residential parcels and resouce locations respectively. Furthermore,
    the capacity that connects the artifical source to residential parcel intersections (see :func:`min_cost_flow_parcels` for discussion on method)
    is equal to the number of residential parcels that are close to that intersection. This way the maximum flow ends up being the maximum number
    of households.


    :param G: Graph network. Will be converted to DiGraph in function
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the residential parcels are to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param G_capacity: name of attribute in G refering to edge capacities, *Default='capacity'* 
    :type G_capacity: string
    :param G_weight: name of attribute in G refering to edge weights, *Default='travel_time'* 
    :type G_weight: string  
    :returns: 
        - **flow_dictionary**, dictionary of dictionaries keyed by nodes for edge flows
        - **cost_of_flow**, integer, total cost of all flow. If weight/cost is travel time, this is the total time for everyone to reach the resource (seconds).
        - **max_flow**, integer, the maximum flow (i.e., households) that can access the resource
        - **access**, string, either 'Complete' or'Partial', representing if every residential parcel can/cannot access the resource
    :rtype: tuple

    '''

    # Travel times must be whole numbers - just round values
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])

    # Find the source and sink demands and append to graph G
    G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points = nearest_nodes(
        G=G, res_points=res_points, dest_points=dest_points, G_demand=G_demand)

    # add artifical source node
    G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
    # add edge from artifical source node to real source nodes with 0 weight and capacity equal to demand
    sums=0
    for unique_node in unique_origin_nodes:
        sums -= G.nodes[unique_node][G_demand]
        kwargs = {f"{G_weight}": 0, f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}  # TODO: KWARGS
        G.add_edge(0, unique_node, **kwargs)
        # since we added an artificial source node, all original source nodes must have a zero demand
        G.nodes[unique_node][G_demand]=0
        
    # create artificial sink node 
    G.add_nodes_from([(99999999, {G_demand: positive_demand})])
    # add edges from sinks to this artificial node
    for x in unique_dest_nodes:
        kwargs={f"{G_weight}":0}
        G.add_edge(x, 99999999, **kwargs)



    # run the max_flow_min_cost function to retrieve FlowDict
    flow_dictionary = nx.max_flow_min_cost(
        G=G, s=0, t=99999999, weight=G_weight, capacity=G_capacity)
    # run the cost_of_flow function to retrieve total cost
    cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
    max_flow, max_dict = nx.maximum_flow(G,
                                         _s=0,
                                         _t=99999999,
                                         capacity=G_capacity)
    if max_flow == sums:
        access = 'Complete'
    else:
        access = 'Partial'
    return flow_dictionary, cost_of_flow, max_flow, access


def plot_aoi(G='', res_parcels='', 
                    resource_parcels='', 
                    edge_width=None, 
                    bbox=None, 
                    loss_access_parcels=None, 
                    scalebar=False,
                    inundation=None, 
                    insets=None, 
                    save_loc=None):
    '''
    Create a plot with commonly used features.

    This function has a variety of parameters to create quick figures that I commonly look at. This can be actively changed to add more/less features.
    **Returns final figure, use plt.show() to see plot.**


    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_parcels: The residential parcels to be plotted. Default color is tan.
    :type res_parcels: geopandas.GeoDataframe
    :param resource_parcels: The resource parcels to be plotted. Default color is blue.
    :type resource_parcels: geopandas.GeoDataframe
    :param edge_width: *Optional*, attribute of *G* to scale road edges by (*Default=None*).
    :type edge_width: string
    :param bbox: set boundary of figures. If *None*, bbox set to max boundary of *G*. 
    :type bbox: geopandas.GeoDataframe
    :param lose_access_parcels: The residential parcels that can no longer access a resource (i.e., due to flooding) (*Default=None*).
    :type lose_access_parcels: geopandas.GeoDataframe
    :param scalebar: If *True*, plots a scalebar (*Default=False*).
    :type scalebar: bool
    :param inundation: Raster of inundation to plot over figure. Use rasterio.open
    :type inundadtion: rasterio Object
    :param insets: File locations of areas to create inset box outlines (*Default=None*). 
    :type insets: list of strings
    :param save_loc: Location to save figure (*Default=None*).
    :type save_loc: string

    :return: **fig**, the produced figure
    :rtype: matplotlib figure

    '''
    
    
    # convert graph to gdf of edges and nodes
    G_gdf_edges = ox.graph_to_gdfs(G=G, nodes=False)
    G_gdf_nodes = ox.graph_to_gdfs(G=G, edges=False)
    # determing bounding box for plotting purposes
    if bbox is None:
        # set bbox to area of graph
        total_bounds = G_gdf_edges.total_bounds
        # convert bounds to N, S, E, W
        total_bounds = (total_bounds[3],
                        total_bounds[1],
                        total_bounds[2],
                        total_bounds[0])
    else:
        # set bbox to localized area
        total_bounds = bbox.total_bounds
        # convert bounds to N, S, E, W
        total_bounds = (total_bounds[3],
                        total_bounds[1],
                        total_bounds[2],
                        total_bounds[0])

    fig, ax = plt.subplots(facecolor='white', figsize=(12, 12))
    ax.axis('off')

    # plot roads, residential parcels, and resource parcels
    res_parcels.plot(ax=ax, color='antiquewhite', edgecolor='tan')
    resource_parcels.plot(ax=ax, color='cornflowerblue', edgecolor='royalblue')

    #TODO: ADD FEATURE TO HIGHLIGHT NODES - USEFUL IN ORDER TO TRANSLATE RESULTS TO GRAPH QUICKLY WHILE TESTING 
    # G_gdf_nodes.loc[[57],'geometry'].plot(ax=ax,color='red')
    # G_gdf_nodes.loc[[56], 'geometry'].plot(ax=ax, color='green')
    #### END TEST


    # option to plot loss of access parcels
    if loss_access_parcels is None:
        pass
    else:
        loss_access_parcels.plot(ax=ax, color='firebrick', edgecolor='darkred')

    # plot background light gray roads
    ox.plot.plot_graph(G,
                       ax=ax,
                       node_size=0,
                       edge_color='lightgray',
                       show=False,
                       edge_linewidth=1,
                       bbox=total_bounds)

    if edge_width is None:
        pass
    else:

        flow_dict = nx.get_edge_attributes(G, edge_width)
        flows = []
        for x in flow_dict:
            flows.append(flow_dict[x])
        unique_flows = np.unique(flows).tolist()
        unique_flows.pop(0)

        old_min = min(unique_flows)
        old_max = max(unique_flows)
        new_min = 1
        new_max = 7

        for flow in unique_flows:
            new_value = (((flow - old_min) * (new_max - new_min)) /
                         (old_max - old_min)) + new_min

            # the following is selection of graph nodes and edges natively
            # having trouble plotting this for unknown reason though
            # selected_edges = [(u, v, e)
            #                  for u, v, e in G.edges(data=True) if e[edge_width] == flow]

            # select edges of gdf based on flow value
            selected_edges = G_gdf_edges[G_gdf_edges[edge_width] == flow]
            sub = ox.graph_from_gdfs(gdf_edges=selected_edges,
                                     gdf_nodes=G_gdf_nodes)

            ox.plot.plot_graph(sub,
                               ax=ax,
                               node_size=0,
                               edge_color='black',
                               show=False,
                               edge_linewidth=new_value,
                               bbox=total_bounds)

    # optional plot inundation raster
    if inundation is None:
        pass
    else:
        # Choose colormap
        cmap = pl.cm.Blues
        # Get the colormap colors
        my_cmap = cmap(np.arange(cmap.N))
        # Set alpha
        my_cmap[:, -1] = np.linspace(0.5, 1, cmap.N)
        my_cmap[0][3] = 0
        # Create new colormap
        my_cmap = ListedColormap(my_cmap)
        show(raster, ax=ax, cmap=my_cmap, zorder=100)

    # optional plot inset boundaries
    if insets is None:
        pass
    else:
        for inset in insets:
            inset = gpd.read_file(inset)
            inset.boundary.plot(ax=ax,
                                edgecolor='forestgreen',
                                linewidth=3,
                                zorder=1000)

    # optional add scale bar
    if scalebar is False:
        pass
    else:
        ax.add_artist(ScaleBar(1,
                               frameon=True,
                               box_color='lightgray',
                               location='lower left'))
    fig.tight_layout()

    if save_loc is None:
        pass
    else:
        plt.savefig(save_loc)
    return(fig)


def summary_function(G=''):
    '''
    Return quick summary information regarding a network.

    This could be expanded upon in the future to be able to quickly access various summary statistics.

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :returns: 
        - **num_edges**, number of edges in the network
        - **num_nodes**, number of nodes in the network
        - **figure**, multiple histograms of network including types of roads, speeds, lengths, travel_times
    :rtype: tuple

    '''
    
    num_edges = len(G.edges)
    num_nodes = len(G.nodes)

    road_types = {}
    speeds = []
    lengths = []
    travel_times = []

    # for each edge in the graph network
    for x in G.edges(data=True):
        # each x has from node, to node, and data variable
        data = x[2]
        # obtain unique types of roads
        # first extract type of road and number of those roads
        if data['highway'] not in road_types:
            road_types[data['highway']] = 1
        else:
            road_types[data['highway']] += 1
        # create a list of speeds (kph)  lengths (meters) and travel time (sec)
        speeds.append(data["speed_kph"])
        lengths.append(data["length"])
        travel_times.append(data["travel_time"])

    # plot histograms of types of roads, speeds, lengths, travel_times
    # for travel times , num bins calculated using Freeman-Diaconis rules

    def fd_rule(values):
        """
        Calculates the number of bins in histogram using freeman-diaconis rules
        """
        num = round(max(values)-min(values) /
                    (2*iqr(values)*len(values)**(-1/3)))
        return(num)

    n_bins_travel_times = fd_rule(travel_times)
    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2)
    ax0.bar(list(road_types.keys()), road_types.values())
    ax1.hist(speeds, 8)
    ax2.hist(lengths, 8)
    ax3.hist(travel_times, n_bins_travel_times)
    # set titles of different plots
    ax0.set_title('Road Types')
    ax1.set_title('Sppeds (kmh)')
    ax2.set_title('Length of Segment (meters)')
    ax3.set_title('Travel Time (secs)')
    # rotate road axis labels
    ax0.tick_params(axis='x', labelrotation=45, labelsize=6)
    fig.tight_layout()

    return(num_edges, num_nodes, fig)


def inundate_network(G='', CRS=32614, path='', inundation=''):
    '''
    Create a new graph network based on the imapct of an inundation layer.

    Currently - An inundation layer is overlayed on a graph network. The maximum depth that intersects each road segment (within a given 
    distance based on the type of road, **this needs to become a parameter and should be customizable**). That maximum depth is then equated 
    to a decrease in speed and/or capacity on that entire road segment (**also needs attention** because water on a road segment might only impact
    poritions of that segment). 

    This function adds attributes including:

    - max_inundation
    - inundation_class
    - inundation_capacity
    - inundation_speed_kph
    - inundation_travel_time
    
    **Important Note** - There are a number of methodologies in which this can be accomplished, and others need to be explored further. The 
    documenation listed must be updated regularly whenever changes are made to reflect current options of this function. Just a few things that need
    to be addressed include:

    - customizable road widths
    - customizable impacts - based on the literature
    - impact on partial segments compared to entire segment


    :param G: The graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param CRS: Coordinate reference system to conduct overlay in. *Default=32614, eqivalent to UTM14N*
    :type CRS: int
    :param path: File path to save output graph network to.
    :type path: string
    :param inundation: File path to inundation raster .tif
    :type inundation: string
    

    :return: **inundated_G**, the impacted graph network.
    :rtype: networkx.Graph [Multi, MultiDi, Di]
    
    '''
    
    nodes, edges = ox.graph_to_gdfs(G)

    # need to set a CRS for the edges to plot properly (UNITS FOR BUFFER)
    edges = edges.to_crs(epsg=CRS)  # WGS84/UTM 14N
    nodes = nodes.to_crs(epsg=CRS)
    # B/C map is in meters, buffer will also be in meters
    # Each road has different number of lanes and therefore diferent buffer
    # Set width variables from intuition:
    # average lane width is 3.7-meter in the US, for now use generic # of lanes
    lane_width = 3.7  # meters, standard for Interstate highway system
    motorway_w = lane_width*6/2
    primary_w = lane_width*4/2
    secondary_w = lane_width*2/2
    tertiary_w = lane_width*2/2
    residential_w = lane_width*2/2
    motorway_link_w = lane_width*1/2
    primary_link_w = lane_width*1/2
    secondary_link_w = lane_width*1/2
    # set conditions for relate
    conditions = [
        edges['highway'] == 'motorway',
        edges['highway'] == 'primary',
        edges['highway'] == 'secondary',
        edges['highway'] == 'tertiary',
        edges['highway'] == 'residential',
        edges['highway'] == 'motorway_link',
        edges['highway'] == 'primary_link',
        edges['highway'] == 'secondary_link',
    ]
    # set choices for relate
    choices = [motorway_w,
               primary_w,
               secondary_w,
               tertiary_w,
               residential_w,
               motorway_link_w,
               primary_link_w,
               secondary_link_w]
    # relate road type to lane/ buffer width
    edges['buffer_width'] = np.select(conditions, choices, default=lane_width)
    # buffer roads with flat end caps: buffer geometry is set as new attribute
    edges['buffer_geometry'] = edges.buffer(distance=edges['buffer_width'],
                                            cap_style=2)
    # zonal stats returns a dictionary, we only want to keep the value
    # This is function that relates maximum flood depth to each road segment
    edges['max_inundation'] = [x['max'] for x in
                               zonal_stats(edges['buffer_geometry'],
                                           inundation,
                                           stats='max')]
    # set flood classification attribute/column
    # TODO FOR NOW: Use paper depth classifications, needs to be changed to an input
    flood_classes = [0, 1, 2, 3, 4, 5, 6]
    bounds = [0.01, 0.15, 0.29, 0.49, 0.91, 1.07]

    for idx, f_class in enumerate(flood_classes):
        if idx == 0:
            edges.loc[edges['max_inundation'] < bounds[0],
                      'inundation_class'] = f_class
        elif idx == len(flood_classes)-1:
            edges.loc[edges['max_inundation'] > bounds[-1],
                      'inundation_class'] = f_class
        else:
            edges.loc[(edges['max_inundation'] > bounds[idx-1]) &
                      (edges['max_inundation'] <= bounds[idx]),
                      'inundation_class'] = f_class

    # Relate inundation class to decrease in speed
    # create new inundated_capacity, inundation_speed, inundation_tavel_time
    edges['inundation_capacity'] = 999999
    edges.loc[edges['inundation_class'] >= 3,
              'inundation_capacity'] = 0
    # if inundation class is 0, speed remains the same
    edges.loc[edges['inundation_class'] == 0,
              'inundation_speed_kph'] = edges['speed_kph']

    # if inundation  class is 1, speed reduced to 25% maximum
    edges.loc[edges['inundation_class'] == 1,
              'inundation_speed_kph'] = edges['speed_kph']*0.25

    # if inundation  class is 2, speed reduced to 5% maximum
    edges.loc[edges['inundation_class'] == 2,
              'inundation_speed_kph'] = edges['speed_kph']*0.05

    # if inundation is anything above 2, set speed to 0
    edges.loc[edges['inundation_class'] >= 3,
              'inundation_speed_kph'] = 0

    # calculate inundation travel time, replace inf with 0
    edges['inundation_travel_time'] = edges['length']/(edges[
        'inundation_speed_kph']/60/60*1000)
    edges.loc[edges['inundation_travel_time'] == np.inf,
              'inundation_travel_time'] = 0
    # save edited graph
    new_graph = ox.graph_from_gdfs(nodes, edges)
    save_2_disk(G=new_graph, path=path, name='AOI_Graph_Inundated')
    return (new_graph)


def flow_decomposition(G='', res_points='', dest_points='',res_locs='',dest_locs='', G_demand='demand', G_capacity='capacity', G_weight='travel_time'):
    '''
    Given a network, G, decompose the maximum flow (w/ minimum cost) into individual flow paths.

    This function serves the role of taking a given network, calculating its maximum flow (with minimum cost), and decomposing the resultant network solution into a feasible series of individual flow pathways. A number of flow decomposition solutions can exist for any given network, *especially* under circumstance where a flow from the same sources travels along different paths. The decomposed information including path, sink, and cost of the path are related back to the residential and destination parcel boundaries for mapping purposes.

    **Important Note** - This function is set up with the intention that all of the flow from one source may not go to the same sink or travel along the same pathway. This is not the case with how the maximum flow is currently calculated, but under more realisitic travel cost functions this will likely change.

    **Furthermore,** this flow decomposition algorithm is just one possible solution. Flows are decomposed by looking at the shortest unit path (taking into account edge weight/costs) from the aritical source to sink. This "greedy" algorithm produces individual flow paths that will have a higher variance, giving preference to paths that are shorter. 
    
    If multiple divering flow paths exist (i.e., same source but different paths/sinks), the function becomes **stochastic** and assigns a residential parcel (who has that source) a random sink, path, and cost that are representative of the distribution of sinks-paths-costs from that source (the sink-path-costs are all appropriately selected in that the cost reflects the path to that specific sink). 


    :param G: The graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the residential parcels are to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param res_locs: Polygons of all residential parcels
    :type res_locs: geopandas.GeoDataFrame
    :param dest_locs: Polygons of all of a specific destination/resource
    :type dest_locs: geopandas.GeodataFrame
    :param G_demand: name of attribuet in G refering to node demand, *Default='demand*
    :type G_demand: string
    :param G_capacity: name of attribute in G refering to edge capacities, *Default='capacity'* 
    :type G_capacity: string
    :param G_weight: name of attribute in G refering to edge weights, *Default='travel_time'* 
    :type G_weight: string      

    :returns:
        - **decomposed_paths**, dictionary of dictionaries keyed by unique source nodes
        - **sink_insights**, dictionary of dictionaries containing information about each sink
        - **res_locs**, geopandas.GeoDataFrame, residential parcels with flow information
        - **dest_locs**, geopandas.GeoDataFrame, resource/destination parcels with flow information
    :rtype: tuple



    '''
    
    # Run max_flow_parcels to get flow dictionary
    flow_dict, flow_cost, max_flow, access = max_flow_parcels(G=G, res_points=res_points, dest_points=dest_points, G_capacity=G_capacity, G_weight=G_weight, G_demand=G_demand)
    
    # Travel times must be whole numbers -  round values if not whole numbers
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])

    #Run nearest_nodes to get unique origins and destinations
    G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points = nearest_nodes(
        G=G, res_points=res_points, dest_points=dest_points, G_demand=G_demand)
    
    # Create artificial source and sink
    # add artifical source node
    G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
    # add edge from artifical source node to real source nodes with:
    #   edge weight = 0
    #   node demand = 0
    #   edge capacity = source demand -> the flow from that original source node
    sums = 0
    for unique_node in unique_origin_nodes:
        sums -= G.nodes[unique_node][G_demand]
        # TODO: KWARGS
        kwargs = {f"{G_weight}": 0,
                  f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}
        G.add_edge(0, unique_node, **kwargs)
        G.nodes[unique_node][G_demand] = 0

    # create artificial sink node with 0 weight
    G.add_nodes_from([(99999999, {G_demand: positive_demand})])
    for x in unique_dest_nodes:
        kwargs={f"{G_weight}":0}
        G.add_edge(x, 99999999, **kwargs)


    # Create intermediate graph from max_flow_parcels dictionary, consisting of only edges that have a flod going across them
    # Set the allowable flow of each edge in intermediate graph to the flow going across that edge in original max_flow_parcels solution
    G_inter = nx.DiGraph()
    for i in flow_dict:
        for j in flow_dict[i]:
            if flow_dict[i][j] >0:
                kwargs = {f"{G_weight}": G[i][j][G_weight],
                        'allowable_flow': flow_dict[i][j]}
                G_inter.add_edge(i,j,**kwargs)
    
    # creat empty dictionaries for decomposed path and sink insights
    decomposed_paths = {}
    sink_insights = {}
    
    # Begin greedy algorithm for flow decomposition
    # loop through all shortest paths from artifical source to artificial sink
    #TODO: CURRENTLY-> set up as greedy algorithm, will produce highest variance in final flows if flow paths from common sources split
    while len(G_inter.edges)>0:
        # Find the shortest path from the artifical source to artificial sink
        path = ox.distance.shortest_path(G_inter,0,99999999,weight=G_weight)
        # create empty list of flows -> allowable flow through each edge along shortest path
        flows = []
        # length of the shortest path
        len_path=len(path)
        # calculate the cost per unit flow along this path
        cost = nx.path_weight(G_inter, path, weight=G_weight)
        # for each edge in the path, append flows with the allowable flow
        for idx, x in enumerate(path):
            if idx==len_path-1:
                pass
            else:
                flows.append(G_inter[path[idx]][path[idx+1]]['allowable_flow'])
        # limiting flow is the minimum allowable flow along the edges in the path
        limiting_flow=min(flows)
        # for each edge along the path, subtract the limiting flow from allowable flow
        # if the allowable flow is reduced to zero, remove the edge
        for idx, x in enumerate(path):
            if idx == len_path-1:
                pass
            else:
                G_inter[path[idx]][path[idx+1]]['allowable_flow'] -= limiting_flow
                if G_inter[path[idx]][path[idx+1]]['allowable_flow'] == 0:
                    G_inter.remove_edge(path[idx],path[idx+1])
        
        # append the decomposed path dictionary with the necessary information
        # the key is the origin
        source = path[1]
        sink = path[-2]
        if source in decomposed_paths.keys():
            decomposed_paths[source]['Sink'].append(sink)
            decomposed_paths[source]['Path'].append(path[1:-1])
            decomposed_paths[source]['Flow'].append(limiting_flow)
            decomposed_paths[source]['Cost Per Flow'].append(cost)
            decomposed_paths[source]['Number of Unique Paths'] += 1
            decomposed_paths[source]['Total Flow'] += limiting_flow  
        else:
            decomposed_paths[source] = {'Source': [source], 'Sink': [sink], 'Path': [path[1:-1]],'Flow': [limiting_flow], 'Cost Per Flow': [cost],'Number of Unique Paths':1, 'Total Flow': limiting_flow}

        # Add sink insights, keyed by the sink   
        if sink in sink_insights.keys():
            sink_insights[sink]['Number of unique paths'] +=1
            sink_insights[sink]['Total flow w/o walking']+= limiting_flow
            sink_insights[sink]['Total flow w/ walking']+= limiting_flow
        else:
            sink_insights[sink] = {'Number of unique paths': 1, 'Total flow w/o walking': limiting_flow, 'Total flow w/ walking': limiting_flow + len(
            res_points.loc[res_points['nearest_node'] == sink])}
    # End flow decomposition algorithm

    # Begin relating decomposition results to outputs 

    # Relating origin results:
    for node in unique_origin_nodes:
        # Create empty lists of possible sinks, paths, and cost_per_flow
        sinks=[]
        paths=[]
        cost_per_flow=[]
        # extract decomposition information -> if decomp_info doesn't exist, that means unique origin node has NO path to sink
        try:
            decomp_info = decomposed_paths[node]
        except: 
            missing_paths = True
        else:
            decomp_info = decomposed_paths[node]
            missing_paths=False
        # for each path from a unique origin to any sink, need to append lists of possible attributes
        if missing_paths is True:
            # If no decomp_info, source node is not serviceable -> unreachable
            res_points.loc[res_points['nearest_node'] == node,['service']] = 'no'
            res_points.loc[res_points['nearest_node']
                           == node, ['cost_of_flow']] = np.NaN

        else:
            # If multiple paths from one source to sink exist, need to create lists that contain a representative distribution of possible paths
            for idx, pathway in enumerate(decomp_info['Sink']):
                flow=decomp_info['Flow'][idx]
                sinks.extend([decomp_info['Sink'][idx]] for x in range(flow))
                path = decomp_info['Path'][idx]
                path= ' '.join(str(e) for e in path)
                paths.extend(path for x in range(flow))
                cost_per_flow.extend([decomp_info['Cost Per Flow'][idx]] for x in range(flow))
            
            # Create subset of res_points that have the source node as their nearest_node
            subset = res_points.loc[res_points['nearest_node'] == node]

            # for each of the res_points in the subset, randomly select an appropriate sink, path, and cost
            for index, res in subset.iterrows():
                i=random.choice(range(len(sinks)))
                sink_pop = sinks.pop(i)
                path_pop = paths.pop(i)
                cost_pop = cost_per_flow.pop(i)

                # Append the res_points geodataframe with the appropriate information include:
                    # The sink that resident goes to
                    # the path that resident takes 
                    # the cost of that path
                    # that it is NaN walkable -> see shared nodes
                    # that it is serviceable -> can reach destination
                res_points.loc[res_points['PROP_ID']==res['PROP_ID'],['sink']] = sink_pop
                res_points.loc[res_points['PROP_ID']==res['PROP_ID'],['path']] = [path_pop]
                res_points.loc[res_points['PROP_ID'] ==
                            res['PROP_ID'],['cost_of_flow']] = cost_pop
                res_points.loc[res_points['PROP_ID'] ==
                            res['PROP_ID'],['walkable?']] = np.NaN
                res_points.loc[res_points['PROP_ID'] ==
                            res['PROP_ID'],['service']] = 'yes'
    
    # For shared_nodes, we consider these res_points walkable for sharing the same nearest node
    for node in shared_nodes:
        res_points.loc[res_points['nearest_node'] == node,['walkable?']] = 'yes'
        res_points.loc[res_points['nearest_node'] == node,['service']] = 'yes'
        res_points.loc[res_points['nearest_node'] == node,['cost_of_flow']] = 0
        res_points.loc[res_points['nearest_node'] == node,['sink']] = node
        res_points.loc[res_points['nearest_node'] == node,['path']] = 0
        # Append the dest_points geodataframe with walking nodes to accurately represent flow totals
        dest_points.loc[dest_points['nearest_node'] == node, ['total_flow_w_walking']] = len(res_points.loc[res_points['nearest_node'] == node])

    # Append dest_nodes with appropriate flow information
    for node in unique_dest_nodes:
        try: 
            total_flow = sink_insights[node]['Total flow w/o walking']
        except:
            no_flow=True
        else: 
            total_flow = sink_insights[node]['Total flow w/o walking']
            no_flow=False
        if no_flow is True:
            dest_points.loc[dest_points['nearest_node'] == node, ['total_flow_wo_walking']] = 0
        else:
            dest_points.loc[dest_points['nearest_node'] == node, ['total_flow_wo_walking']] = total_flow
            dest_points.loc[dest_points['nearest_node'] == node, ['total_flow_w_walking']] += total_flow

    # Sync the res_points with res_locs data
    res_locs = res_locs.merge(res_points[['PROP_ID','nearest_node','service','sink','path','cost_of_flow','walkable?']], on='PROP_ID')
    # Sync the dest_points with dest_locs data
    dest_locs = dest_locs.merge(dest_points[['PROP_ID', 'nearest_node','total_flow_wo_walking','total_flow_w_walking']], on='PROP_ID')

    # return the decomposed paths dict of dict, the sink_insights dict of dict, and the res_locs/dest_locs geodataframes
    return decomposed_paths, sink_insights, res_locs, dest_locs




def shortestPath(origin, capacity_array, weight_array):
    """
    This method finds the shortest path from origin to all other nodes
    in the network.  
    """
    # Need the number of nodes for calculations
    numnodes = np.shape(capacity_array)[0]

    # Set up backnode and cost lists to return
    backnode = [-1] * numnodes
    costlabel = [np.inf] * numnodes

    # I am going to be implementing Dijkstra's Algorithm

    # Step 1: Initialize all labels Li=infinity origin Ls=0
    # costlabel already initialized, need to change cost label for origin
    costlabel[origin] = 0
    # Also creating f_label for final labels
    f_label = []

    # Step 2: Initialize SEL with origin
    sel = [origin]
    x = 0
    # Create termination loop, if sel empty, terminate
    while len(sel) > 0:
        # Step 3: Choose a node i from SEL with minimum L
        # doing this by creating a list of cost labels for each node in SEL
        labels = []
        for node in sel:
            labels.append(costlabel[node])
        # node we want is popped from SEL (the index of min in label)
        node_i = sel.pop(labels.index(min(labels)))

        # Step 4: for all arc (i,j) repeat following
        for node_j, matrix_value in enumerate(capacity_array[node_i]):
            # if matrix_value is equal to -1, no link exist
            if matrix_value == -1:
                pass
            # also need to skip nodes in f, because already finalized
            elif node_j in f_label:
                pass
            else:
                # Step 4a: calculate costlabel temp
                l_temp = costlabel[node_i] + weight_array[node_i][node_j]
                # Step 4b: if new path is shorter, update cost and backnode
                if l_temp < costlabel[node_j]:
                    costlabel[node_j] = l_temp
                    backnode[node_j] = node_i
                # Step 4c: add j tos can list if it is not already there
                if node_j not in sel:
                    sel.append(node_j)
                # Step 5: add i to f_label
                f_label.append(node_i)
    backnode=np.asarray(backnode)
    costlabel=np.asarray(costlabel)
    # backnodes equal to list indicating previous node in shortest path
    # costlabel is list of costs associated with each node]
    return (backnode, costlabel)

def pathTo(backnode, costlabel, origin, destination):
    cost = costlabel[destination]
    path = [destination]
    idx = destination
    while idx > 0:
        idx = backnode[idx]
        path.append(idx)
    path.append(origin)
    path.reverse()
    path.pop(0)
    path.pop(0)
    path=np.asarray(path)
    cost=np.asarray(cost)
    return(path, cost)

#I DON'T THINK I NEED maxFlow OR mod_search
def maxFlow(source, sink, numnodes, capacity_array, weight_array):
    """
    This method finds a maximum flow in the network.  Return a tuple
    containing the total amount of flow shipped, and the flow on each link
    in a list-of-lists (which takes the same form as the adjacency matrix).
    These are already created for you at the start of the method, fill them
    in with your implementation of minimum cost flow.
    """

    # Set up network flows matrix; flow[i][j] should have the
    # value of x_ij in the max-flow algorithm.
    totalFlow = 0
    flow = list()
    for i in range(numnodes):
        flow.append([0] * numnodes)

    # I am using the augmenting path method
    # Step 1: initialize the flow (done above with total flow and flow variables)

    # Step 2: identify unidriected path pi connecting s and t
    # this is my modified search method (see method for specifics)
    reachable, backnode = mod_search(source=source, flow=flow, numnodes=numnodes, capacity_array=capacity_array)
    # Based on criteria in mod_search, if sink is in reachable,
    #   than an augmented path exists
    while sink in reachable:
        # Determine the augmented path using the backnodes variable
        path = [sink]
        idx = sink
        while idx != source:
            idx = backnode[idx]
            path.append(idx)
        path.reverse()
        # Determine which nodes in path are forward and reverse arcs
        forward_arcs = []
        reverse_arcs = []
        for index, i in enumerate(path):
            if index+1 == len(path):
                pass
            else:
                j = path[index+1]
                if capacity_array[i][j] != -1:
                    forward_arcs.append([i, j])
                else:
                    reverse_arcs.append([j, i])

        # Step 3: Calculate residual capacity, and add to total flow
        residual_capacity_list = []
        for arc in forward_arcs:
            i = arc[0]
            j = arc[1]
            residual_capacity_list.append(capacity_array[i][j] - flow[i][j])

        for arc in reverse_arcs:
            i = arc[1]
            j = arc[0]
            residual_capacity_list.append(capacity_array[i][j] - flow[i][j])

        residual_capacity = min(residual_capacity_list)
        totalFlow += residual_capacity

        # Step 4: add/subtract residual_capacity from arcs
        for arc in forward_arcs:
            i = arc[0]
            j = arc[1]
            flow[i][j] += residual_capacity

        for arc in reverse_arcs:
            i = arc[1]
            j = arc[0]
            flow[i][j] -= residual_capacity

        # recalculate the reachable nodes to determine if should terminate or not
        reachable, backnode = mod_search(
            source=source, flow=flow, numnodes=numnodes, capacity_array=capacity_array)
    flow = np.asarray(flow)
    return (totalFlow, flow)

def mod_search(source, flow, numnodes, capacity_array):
# creates a list of reachable nodes from source based on augmenting path rules
# a reachable node is either forward or reverse ([i][j] or [j][i] and capacity greater than flow)
    reachable = [source]
    backnode = [-1]*numnodes
    SEL = [source]
    while len(SEL) > 0:
        i = SEL[0]
        SEL = SEL[1:]
        for j in range(0, numnodes):
            # selection criteria for augmenting path:
            # foward path with available capacity
            # reverse path with available capacity
            if ((capacity_array[i][j] != -1 and capacity_array[i][j] > flow[i][j]) or (capacity_array[j][i] != -1 and capacity_array[j][i] < flow[j][i])) and j not in reachable:
                reachable.append(j)
                SEL.append(j)
                backnode[j] = i
    return reachable, backnode







def traffic_assignment(G='',
                        res_points='', 
                        dest_points='', 
                        G_demand='demand', 
                        G_capacity='capacity', 
                        G_weight='travel_time', 
                        algorithm='link_based', 
                        method='CFW',
                        link_performance='BPR'):
    """
    This is my attempt at a traffic assignment problem, with options for user equilibirum (UE), and system optimal (SO) solutions,
    Ideally, I would also like to add different algorithms (link-based, path-based, and bush-based) to test their operations, 
    To scale this computationally, I am thinking of writing this all in terms of numpy arrays - which may be faster
    TODO: look into rewriting other code to also put in terms of numpy arrays? Not sure on this if necessary, feasible, or what
            It's definitely feasible it's probably just not necessary - NetworkX is probably optimized for what we would need it for 

    Plan is to also include options within each type of algorithm (e.g., MSA and Frank-Wolfe for link-based algorithms)
    Also plan on using different types of termination criteria (e.g., AEC option, number of iterations option for testing purposes, etc. )

    Algorithm argument: Whether to use link, path, or bush based algorithm
    method argument: if the algorithm is multiple options (e.g., MSA versus CFW for link based algorithms). For now, only link based algorithms have different methods.
        These are the algorithms that will can be utilized:
            - Link Based
                - Method of Successive Averages (MSA)
                - Bisection
                - Newton's Method
                - Conjugate Frank-Wolfe
            
            - Path Based
                - Gradient Projection
            
            -Bush Based 
                - Algorithm B 
                - (Maybe later) Origin-based assignment
                - (Maybe later) Linear user cost equilibrium

    Link_performance argument: determines which link peformance equation to be utilized (e.g., BPR).
        Right now, only focused on BPR but could expand to include the following:
            - Davidson's delay model
            - Akcelik function
            - Conical delay model

    """


    # Travel times must be whole numbers -  round values if not whole numbers
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])

    # Convert graph to digraph format
    # TODO: Need to see if this is having an impact on the output
    G = nx.DiGraph(G)

    # Set variables that will be used as constants throughout algorithm 
    super_sink = 99999999
    super_origin = 0
    artificial_weight = 999999999999
    artificial_capacity = 999999999999


    # Theoretical capacity for each edge based on how many lanes there are
    # TODO: FOR NOW, SET AS A CONSTANT - I believe typically 2300 vehicles/hour per lane is a standard - need to research tho
    # This should also most likely be done before this section anyways
    nx.set_edge_attributes(G, values=2300, name=G_capacity)
    
    # Set BPR alpha and beta constants
    # TODO: FOR NOW, SET AS CONSTANT - MIGHT BE A FUNCTION OF THE TYPE OF ROAD - THIS SHOULD ALSO LIKELY BE DONE ELSEWHERE
        # HOWEVER: as an array/matrix, increases the computational time greatly
    nx.set_edge_attributes(G, values=0.15, name='alpha')
    nx.set_edge_attributes(G, values=4, name='beta')

    # PREPROCESSING STEPS: 
    # a.   Determine sources/sinks, 
    # b.   set aritificial components of network, 
    # c.   create OD Matrix
    # d.   create list of ordered nodes
    # e.   Convert capacity, weight, alpha, and beta arrays

    # a. Determine sources and sinks
    G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points = nearest_nodes(G=G, res_points=res_points, dest_points=dest_points, G_demand='demand')

    # b_1. Add an artifical 0 node to maintain other nodes positions
    #   Doesn't actually need attributes since not connected, but filling for continuity sake
    G.add_node(super_origin, **{G_weight: artificial_weight, G_capacity: 0})

    # b_2. Add artifical edge from each destination node to the artifical sink with zero cost and maximum capacity
    # This is to allow minimum cost routing to whatever resource - destination of OD pair can change based on cost
    for idx, sink in enumerate(unique_dest_nodes):
        G.add_edge(sink, super_sink, **{G_weight:0, G_capacity:artificial_capacity})

    # c. Create OD_matrix and add artifical edge from each origin to super_sink with extreme weight and capacity to capture loss of service 
    #    2 reasons to add origin->super_sink edge: 
    #   (1): Can build in elastic demand cost if people are priced out of accessign resource, or can set to arbitrarily high value to
    #   (2): By having artifical route, can always send max flow w/o considering if accessible - if path goes through artifical link, than no access when cost is arbitrarily high
    OD_matrix = np.empty([len(unique_origin_nodes),3])
    for idx,x in enumerate(unique_origin_nodes):
        OD_matrix[idx][0] = x                            # origin nodes
        OD_matrix[idx][1] = G.nodes[x]['demand']*-1      # flow associated with the origin node
        OD_matrix[idx][2] = len(G.nodes())-1             # destination node (when nodes reordered for numpy array, it has a value of length-1 (b/c of artifical origin 0 ))
        G.add_edge(x, super_sink, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

    # d. Sort nodes in proper order
    #   This should in theory then preserve saying node '2' in array is node '2' in network g, especially with adding the artifical node, 0
    nodes_in_order = sorted(G.nodes())

    # convert to capacity and weight arrays
    capacity_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight=G_capacity, nonedge=-1)    # array of theoretical link capacities
    weight_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight=G_weight, nonedge=-1)        # array of free flow weights (costs)
    alpha_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='alpha', nonedge=-1) 
    beta_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='beta', nonedge=-1)
    # TODO: Check the timing on this -> I believe if I actually used a full array (above) instead of constant (below) the compute time would be significantly longer
    alpha_array = 0.15
    beta_array=4
    #weight_array_star = np.copy(weight_array)  # copy of free flow weights, what I am going to be changin with each iteration
    

    # set variables
    dest_id = len(G.nodes())-1   # destination node
    numnodes=len(G.nodes())      # number of nodes 
    sum_d = positive_demand      # sum of the demand


    ##################################################################################################################################
    # EXAMPLE ARRAYS THAT WERE USED FOR TESTING PURPOSES


    # Going to create an test network for building this function based on steve boyles textbox
    # G = nx.DiGraph()
    # G.add_nodes_from([0,1,2,3,4,5,6])
    # G.add_edges_from([(1,3),(1,5),(5,6),(6,3),(2,5),(6,4),(2,4)])
    # nx.set_edge_attributes(G,10,'Weight')
    # nx.set_edge_attributes(G,15000,'Capacity')
    # # OD Pair: 1-3: 5000, 2-4: 10000
    # OD_matrix = np.array([[1,5000,3],
    #                     [2,10000,4]])
    # # BPR function = 10 + x/100

    # nodes_in_order = sorted(G.nodes())
    # capacity_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='Capacity', nonedge=-1)
    # weight_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='Weight', nonedge=-1)

    # # set variables for link algorithm
    # dest_id = len(G.nodes())-1
    # numnodes = len(G.nodes())
    # sum_d = 15000

    # This is another network example specifically for bush-based algorithms FROM TEXTBOOK
    # G = nx.DiGraph()
    # G.add_nodes_from([0,1,2,3,4,5,6,7,8,9])
    # G.add_edges_from([(1,2),(2,3),(4,1),(4,2),(4,5),(5,2),(5,6),(6,3),(7,4),(7,8),(8,5),(8,9),(9,6)])
    # attrs = {(1,2):{'Weight':4},(2,3):{'Weight':2},(4,1):{'Weight':4},(4,2):{'Weight':40},(4,5):{'Weight':2},(5,2):{'Weight':2},(5,6):{'Weight':2},
    #                 (6,3):{'Weight':2},(7,4):{'Weight':2},(7,8):{'Weight':2},(8,5):{'Weight':2},(8,9):{'Weight':2},(9,6):{'Weight':2}}
    # nx.set_edge_attributes(G,attrs)
    # nx.set_edge_attributes(G,15000,'Capacity')
    # # OD Pair: 1-3: 5000, 2-4: 10000
    # OD_matrix = np.array([[7,10,3]])
    # # BPR function depends on the link

    # nodes_in_order = sorted(G.nodes())
    # capacity_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='Capacity', nonedge=-1)
    # weight_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='Weight', nonedge=-1)
    # # set variables for link algorithm
    # sum_d = 10

    # This is another network example specifically for bush-based algorithms FROM NOTES
    G = nx.DiGraph()
    G.add_nodes_from([0,1,2,3,4,5,6,7,8,9])
    G.add_edges_from([(1,2),(1,4),(2,3),(2,5),(3,6),(4,5),(4,7),(5,6),(5,8),(6,9),(7,8),(8,9)])
    attrs = {(1,2):{'Weight':3},(1,4):{'Weight':3},(2,3):{'Weight':3},(2,5):{'Weight':5},(3,6):{'Weight':3},(4,5):{'Weight':5},(4,7):{'Weight':3},
                    (5,6):{'Weight':5},(5,8):{'Weight':5},(6,9):{'Weight':3},(7,8):{'Weight':3},(8,9):{'Weight':5}}
    nx.set_edge_attributes(G,attrs)
    nx.set_edge_attributes(G,15000,'Capacity')
    # OD Pair: 1-3: 5000, 2-4: 10000
    OD_matrix = np.array([[4,1000,9],
                          [1,1000,9]])
    # BPR function depends on the link

    nodes_in_order = sorted(G.nodes())
    capacity_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='Capacity', nonedge=-1)
    weight_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='Weight', nonedge=-1)
    # set variables for link algorithm
    sum_d = 2000

    ################################################################################################################################


    # WRITE FUNCTIONS FOR CALCULATING BPR FUNCTIONS FOR WEIGHT ARRAY AND WEIGHT ARRAY DERIVATIVE
    # I think for now it would be easiest to assume that every road has the same equation?
    # weight_array variable should remain the free flow time array
        # weight_array_itter (or weight_array_iter, need to fix spelling for different algorithms) will be the changing weight array with each iteration

    # LINK PERFORMANCE HELPER FUNCTIONS 
    # FIXME: new helper functions go here
    # So far, only link_performance_function and link_performance derivative require equation inputs 

    def link_performance_function(capacity_array,flow_array,weight_array,eq=link_performance,alpha_array=0,beta_array=0):
        """ Function returning the weight_array_iter array based on user selected equation
        
        Universival arguments:
            capacity_array == node-node array of link capacities
            flow_array     == node-node array of link flows
            weight_array   == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
            TODO: right now I am not sure if array of values will work, need to double check
            alpha_array == either node-node array of alpha values for each link OR single value 
            beta_array  == either node-node array of beta values for each link OR single value     

        returns: weight_array_iter
        """
        
        if eq=="BPR":
            weight_array_iter = weight_array*(1+alpha_array*(flow_array/capacity_array)**beta_array)
        # elif eq == "Davidsons":
            # pass
        
        return weight_array_iter


    def link_performance_derivative(capacity_array, flow_array, weight_array, eq=link_performance, alpha_array=0, beta_array=0):
        """ Function returning the derivative array based on user selected equation
        
        Universival arguments:
            capacity_array == node-node array of link capacities
            flow_array     == node-node array of link flows
            weight_array   == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
            TODO: right now I am not sure if array of values will work, need to double check
            alpha_array == either node-node array of alpha values for each link OR single value 
            beta_array  == either node-node array of beta values for each link OR single value     

        returns: link_performance_derivative
        """

        if eq == "BPR":
            link_performance_derivative = (
                beta_array*weight_array*alpha_array/capacity_array)*flow_array**(beta_array-1)
        # elif eq == "Davidsons":
            # pass

        return link_performance_derivative


    def bisection_zeta(lambda_val, flow_array_star, capacity_array, flow_array, weight_array, eq=link_performance, alpha_array=0, beta_array=0):
        """ Function returning the zeta value when using the bisection link-based method
        
        Universival arguments:
            lambda_val      == lambda value to shift flow
            flow_array_star == updated flow array
            capacity_array  == node-node array of link capacities
            flow_array      == node-node array of link flows
            weight_array    == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
            TODO: right now I am not sure if array of values will work, need to double check
            alpha_array == either node-node array of alpha values for each link OR single value 
            beta_array  == either node-node array of beta values for each link OR single value     

        returns: bisection zeta value
        """
        x_hat = lambda_val*flow_array_star+(1-lambda_val)*flow_array

        x_hat_link_performance = link_performance_function(capacity_array=capacity_array, 
                                                           flow_array=x_hat,
                                                           weight_array=weight_array,
                                                            eq=eq, 
                                                            alpha_array=alpha_array, 
                                                            beta_array=beta_array)

        zeta = np.sum(x_hat_link_performance*(flow_array_star-flow_array))

        return(zeta)
        

    def newton_f(lambda_val, flow_array_star, capacity_array, flow_array, weight_array, eq=link_performance, alpha_array=0, beta_array=0):
        """ Function returning the f and f_prime values when using the newton's and CFW link-based methods
        
        Universival arguments:
            lambda_val      == lambda value to shift flow
            flow_array_star == updated flow array
            capacity_array  == node-node array of link capacities
            flow_array      == node-node array of link flows
            weight_array    == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
            TODO: right now I am not sure if array of values will work, need to double check
            alpha_array == either node-node array of alpha values for each link OR single value 
            beta_array  == either node-node array of beta values for each link OR single value     

        returns: f and f_prime
        """
        x_hat = lambda_val*flow_array_star+(1-lambda_val)*flow_array

        x_hat_link_performance = link_performance_function(capacity_array=capacity_array,
                                                           flow_array=x_hat,
                                                           weight_array=weight_array,
                                                           eq=eq,
                                                           alpha_array=alpha_array,
                                                           beta_array=beta_array)

        x_hat_link_performance_derivative = link_performance_derivative(capacity_array=capacity_array, 
                                                                        flow_array=x_hat,
                                                                        weight_array=weight_array, 
                                                                        eq=eq, 
                                                                        alpha_array=alpha_array, 
                                                                        beta_array=beta_array)

        f = np.sum(x_hat_link_performance*(flow_array_star-flow_array))
        f_prime = np.sum(x_hat_link_performance_derivative*(flow_array_star-flow_array)**2)

        return f, f_prime


    def cfw_alpha(flow_array_star_old, flow_array_star, capacity_array, flow_array, weight_array, eq=link_performance, alpha_array=0, beta_array=0):
        """ Function returning the numerator and denominator of alpha value for CFW link-based method

        Universival arguments:
            flow_array_star_old == previous updated flow array
            flow_array_star     == updated flow array (the AON, or All-Or-Nothing assignment)
            capacity_array      == node-node array of link capacities
            flow_array          == node-node array of link flows
            weight_array        == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
            TODO: right now I am not sure if array of values will work, need to double check
            alpha_array == either node-node array of alpha values for each link OR single value 
            beta_array  == either node-node array of beta values for each link OR single value     

        returns: f and f_prime
        """

        deriv = link_performance_derivative(capacity_array=capacity_array, 
                                            flow_array=flow_array, 
                                            weight_array=weight_array, 
                                            eq=eq, 
                                            alpha_array=alpha_array, 
                                            beta_array=beta_array)

        #eq (6.3) in boyles textbook page 176
        item1 = flow_array_star_old - flow_array
        item2 = flow_array_star - flow_array
        item3 = flow_array_star - flow_array_star_old
        

        num = np.sum(deriv*item1*item2)
        den = np.sum(deriv*item1*item3)
        return num, den



    # ALGORITHMS
    # LINK-BASED ALGORITHM
    if algorithm == "link_based":

        # INITIALIZE LINK BASED ALGORITHM
        # a.    Create empty flow array - same size as the capacity_array (b/c node-node relationship)
        # b.    Create initial feasible flow array with all-or-nothing minimum cost flow assignment (free flow travel cost)

        # a. Create empty flow array 
        flow_array = np.zeros_like(capacity_array)

        # b. Create initial feasible flow array
        for x in OD_matrix:
            origin = np.int(x[0])
            flow = x[1]
            destination = np.int(x[2])
            
            # calculate shortest path from an origin to all other locations
            backnode, costlabel = shortestPath(origin=origin, capacity_array=capacity_array, weight_array=weight_array)
            # determine path and cost from origin to super sink
            path, cost = pathTo(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)

            # update the flow array by iterating through the shortest paths just determined
            i=0
            while i < len(path)-1:
                flow_array[path[i]][path[i+1]] += flow
                i+=1


        # TERMINATION CRITERIA
        # Could do # of iterations, or AEC <= a value, or RG <= value
        # TODO: Fix termination criteria loop
        iterations = 3
        iter=0

        # LINK BASED ALGORITHM: Within loop:
        # 1.    Recalculate weight_array (travel cost times) using BPR function and initial flow_array
        # 2.    Create empty kappa array to hold OD shortest path costs, and empty flow_array_star to hold new flows
        # 3.    Calculate shortest paths with new weight_array
        # 4.    Determine cost of each path and fill kappa array
        # 5.    Fill flow_array_star with new flows
        # 6.    Calculate lambda 
        # 7.    shift flow, and repeat

        # variables I am interested in keeping track off:
        AEC_list = []
        TSTT_list = []
        SPTT_list = []
        RG_list = []

        # CFW marker -> CFW algorithm has different 1st iteration, use this marker to get around 
        CFW_marker = 0

        while iter <= iterations:
            # 1. recalculate weight_array
            # FIXME: EQUATION WRITTEN: Replace with BPR function
            weight_array_iter = flow_array/100 +10         # Example network BPR
            # weight_array_iter = weight_array*(1+alpha_array*(flow_array/capacity_array)**beta_array)      # Real network BPR
            
            # 2. Create empty kappa array and flow_array_star
            kappa=np.zeros([1,np.shape(OD_matrix)[0]])
            flow_array_star = np.zeros_like(capacity_array)
            # For shortest path travel time calculation (termination criteria)
            SPTT=0      
            for idx, x in enumerate(OD_matrix):
                origin = np.int(x[0])
                destination = np.int(x[2])
                flow = x[1]
                # 3. Calculate shortest paths
                backnode, costlabel = shortestPath(origin=origin, capacity_array=capacity_array, weight_array=weight_array_iter)
                path, cost = pathTo(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
                # 4. Fill shorest path cost array, kappa
                kappa[0][idx]=cost      
                # For shortest path travel time calculation (termination criteria)
                SPTT += cost*flow
                # 5. Update the flow array
                i = 0
                while i < len(path)-1:
                    flow_array_star[path[i]][path[i+1]] += flow
                    i += 1

            # Calculate termination variables
            TSTT = np.sum(np.multiply(flow_array, weight_array_iter))
            AEC = (TSTT-SPTT)/sum_d
            RG = TSTT/SPTT - 1

            # 6. Calculate Lambda using the method given 

            if method == 'MSA':
                # MSA
                lambda_val = 1/(iter+2)
                    
            elif method == 'Bisection': 
                # Frank-Wolfe: Bisection
                # termination criteria for bisection method -> when high and low lambda_val are within a given distance of each other
                term_criteria = 1/32
                # initialize lambda_ lo, and hi values
                lambda_lo = 0
                lambda_hi = 1

                # set up termination criteria:
                while (lambda_hi - lambda_lo) >= term_criteria:
                    # find bisection of lo and hi values
                    lambda_val = (lambda_hi-lambda_lo)/2 + lambda_lo
                    # Calculate x_hat, and subsequently zeta
                    x_hat = lambda_val*flow_array_star+(1-lambda_val)*flow_array
                    zeta = np.sum((x_hat/100 + 10)*(flow_array_star-flow_array))   # Example network BPR
                    # FIXME: EQUATION WRITTEN: Check BPR math below and remove above
                    # zeta = np.sum(((((lambda_val*flow_array_star+(1-lambda_val)*flow_array)/capacity_array)**beta_array*alpha_array+1)*weight_array)*(flow_array_star-flow_array)) 

                    # determine shift
                    if zeta > 0:
                        lambda_hi = lambda_val
                    
                    elif zeta < 0:
                        lambda_lo = lambda_val

                    elif zeta == 0:
                        lambda_hi = lambda_val
                        lambda_lo = lambda_val
            
            elif method == 'Newtons Method':
                # Frank-Wolfe: Newton's Method
                # initialize lambda - just going to choose 0.5
                lambda_val = 0.5
                # calculate x_hat
                x_hat = lambda_val*flow_array_star + (1-lambda_val)*flow_array
                # calculate f and f_prime
                f = np.sum((x_hat/100 + 10)*(flow_array_star-flow_array))  #  example network BPR
                f_prime = np.sum((1/100)*(flow_array_star-flow_array)**2)      # example network BPR
                # FIXME: EQUATION WRITTEN: Program more complex BPR function
                # Update Lambda
                lambda_val -= f/f_prime
                if lambda_val > 1:
                    lambda_val = 1
                elif lambda_val <0:
                    lambda_val = 0
                #TODO: see PPT, add criteria about hitting the same endpoint twice in a row to terminate

            elif method == 'CFW':
                # Frank-wolfe: Conjugate Frank-Wolfe
                # REFER TO PAGE 177 in BOYLES TEXTBOOT
                # first iteration: use Newton's method only
                # subsequent interations, determine new x*star based on alpha calcuation 
                # TODO: add functionality for real BPR function
                epsilon=0.01
                if CFW_marker > 0:
                    # solve for alpha, eq (6.3) in boyles textbook, page 176
                    item1 = flow_array_star_old - flow_array
                    item2 = flow_array_star - flow_array
                    item3 = flow_array_star - flow_array_star_old
                    # FIXME : EQUATION WRITTEN: BPR function
                    num = np.sum((1/100)*item1*item2)
                    den = np.sum((1/100)*item1*item3)
                    if den == 0:
                        alpha=0
                    else:
                        alpha = num/den
                        if alpha >= 1:
                            alpha = 1 - epsilon
                        if alpha < 0:
                            alpha = 0
                    flow_array_star = alpha*flow_array_star_old+(1-alpha)*flow_array_star

                # Calculate Lambda using new flow_array_star estimation utilizing Newton's Method
                # initialize lambda - just going to choose 0.5
                lambda_val = 0.5
                # calculate x_hat
                x_hat = lambda_val*flow_array_star + (1-lambda_val)*flow_array
                # calculate f and f_prime
                # example network BPR
                f = np.sum((x_hat/100 + 10)*(flow_array_star-flow_array))
                # example network BPR
                f_prime = np.sum((1/100)*(flow_array_star-flow_array)**2)
                # FIXME: EQUATION WRITTEN: Program more complex BPR function (SEE ABOVE EQUATIONS)
                # Update Lambda
                lambda_val -= f/f_prime
                if lambda_val > 1:
                    lambda_val = 1
                elif lambda_val < 0:
                    lambda_val = 0

                # flow_array_star_old value for the next iteration
                flow_array_star_old = np.copy(flow_array_star)
                # set marker to pass for subsequent interations
                CFW_marker = 1

            # 7. Shift flows
            flow_array = lambda_val*flow_array_star + (1-lambda_val)*flow_array

            # Fill termination variable lists
            AEC_list.append(AEC)
            TSTT_list.append(TSTT)
            SPTT_list.append(SPTT)
            RG_list.append(RG)
            print(RG)
            iter +=1 

    
    elif algorithm == 'path_based':
        # Specifically using Newton's Method of Gradient Projection (not to be confused with Projected Gradient)

        # Termination criteria variables I am interested in keeping track off:
        AEC_list = []
        TSTT_list = []
        SPTT_list = []
        RG_list = []

        # Initialize the flow and weight arrays
        flow_array = np.zeros_like(capacity_array)
        # FIXME : EQUATION WRITTEN: Example BPR : Update to real - I don't think this needs to be here? 
        weight_array_iter = flow_array/100 + 10

        # Initialize nested list to store OD path information 
        # TODO: IDK if this should be a nested list, nested dictionary, or pandas dataframe
            # for now, using a nested list because it's the easiest to create
        paths_array = []
        for OD_pair in OD_matrix:
            # paths_array[x][0]: [O, D]
            # paths_array[x][1]: [list of arrays that are paths for OD-pair]
            # paths_array[x][2]: [list of flows (associated with paths)]
            # paths_array[x][3]: [cost per unit flow (associated with paths)]

            origin = OD_pair[0]
            flow = OD_pair[1]
            destination = OD_pair[2]
            # Calculate shortest path for OD pair
            backnode, costlabel = shortestPath(origin=origin, capacity_array=capacity_array, weight_array=weight_array_iter)
            path, cost = pathTo(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)

            # update the flow array in order to update the weight array
            i = 0
            while i < len(path)-1:
                flow_array[path[i]][path[i+1]] += flow
                i += 1
            # FIXME : EQUATION WRITTEN: Example BPR : Update to real
            weight_array_iter = flow_array/100 + 10             

            # calculate new cost of shortest path
            cost=0
            i = 0
            while i < len(path)-1:
                cost += weight_array_iter[path[i]][path[i+1]]
                i += 1

            # append the paths_array
            paths_array.append([[origin,destination],[list(path)],[flow],[float(cost)]])

        # calculate initial termination criteria variables
        SPTT = 0
        for idx, x in enumerate(OD_matrix):
            origin = np.int(x[0])
            destination = np.int(x[2])
            flow = x[1]
            # Calculate shortest paths
            backnode, costlabel = shortestPath(
                origin=origin, capacity_array=capacity_array, weight_array=weight_array_iter)
            path, cost = pathTo(
                backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
            # For shortest path travel time calculation (termination criteria)
            SPTT += cost*flow

        # termination criteria
        TSTT = np.sum(np.multiply(flow_array, weight_array_iter))
        AEC = (TSTT-SPTT)/sum_d
        RG = TSTT/SPTT - 1

        # begin iterations
        iterations = 2
        iter=0

        while iter <= iterations:

            for OD_pair in paths_array:
                origin = OD_pair[0][0]
                destination = OD_pair[0][1]
                paths = OD_pair[1]
                flows = OD_pair[2]
                costs = OD_pair[3]
                # Find the new shortest path 
                backnode, costlabel = shortestPath( origin=origin, capacity_array=capacity_array, weight_array=weight_array_iter)
                path_hat, cost_hat = pathTo(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)

                path_hat = list(path_hat)

                # if path already in, don't append, just move it, and associated information to the end
                # TODO : Need to fix this to address if path_hat is only path in paths -  I think for now it should still work fine?
                    # can add an if else statement: if len paths == 1 [do this] else enter for loop
                
                if path_hat in paths:
                    index = paths.index(path_hat)
                    paths.append(paths.pop(index))
                    flows.append(flows.pop(index))
                    costs.append(costs.pop(index))
             
                # if not, append the new path and perform path functions
                else:
                    paths.append(path_hat)
                    # set to 0 to act as place holder for upcoming calculations
                    flows.append(0)
                    costs.append(0)

                # perform functions on each path except the shortest path (in the last index of list)
                for idx, path in enumerate(paths[:-1]):

                    all_links=[]
                    # extract links from path_hat (new shortest path)
                    i = 0
                    while i < len(path_hat)-1:
                        all_links.append([path_hat[i],path_hat[i+1]])
                        i += 1
                    # extract links from current path iteration
                    i = 0
                    while i < len(path)-1:
                        all_links.append([path[i],path[i+1]])
                        i += 1
                    # determine the unique node-node pairs (i.e., the unique links across all paths)
                    unique_links = [list(x) for x in set(tuple(x) for x in all_links)]

                    # determine cost difference (numerator of delta_h equation)
                    num = costs[idx]-cost_hat

                    # determine the sum of derivatives (denominator of delta_h equation)
                    den = 0
                    for link in unique_links:
                        # FIXME: EQUATION WRITTEN: update to derivative of real BPR function - currently set as sample network equation
                        den += 1/100
                    
                    # calculate delta_h
                    delta_h = num/den

                    # adjusted path flow value
                    adj = min(flows[idx], delta_h)

                    # update paths_array flow values
                    flows[idx] -= adj
                    flows[-1] += adj

                    # TODO: add if statement if flow is reduced to 0, need to remove pathway b/c no longer an active path
                    # Need to make sure this function properly if it occurs
                    # if flow reduced to 0 remove from list of active paths
                    if flows[idx] == 0:
                        paths.pop(idx)
                        flows.pop(idx)
                        costs.pop(idx)

                    # for the path under examination and path_hat, need to update flow array to calc new weight array
                    i = 0
                    while i < len(path)-1:
                        flow_array[path[i]][path[i+1]] -= adj
                        i += 1
                    i = 0
                    while i < len(path_hat)-1:
                        flow_array[path_hat[i]][path_hat[i+1]] += adj
                        i += 1
                    # FIXME : EQUATION WRITTEN: BPR: Update to real
                    weight_array_iter = flow_array/100 + 10             

                    # WITH NEW WEIGHT ARRAY, NEW COSTS NEED TO BE CALCULATED AND ASSOCIATED WITH EACH PATH UNDER EXAMINATION
                    # THIS SHOULD FIX THE EXTRA LOOP ERROR
                    for OD_pair in paths_array:
                        paths = OD_pair[1]
                        costs = OD_pair[3]
                        for idx, path in enumerate(paths):
                            cost= 0
                            i = 0
                            while i < len(path)-1:
                                cost += weight_array_iter[path[i]][path[i+1]]
                                i += 1
                                costs[idx] = cost

            # calculate  termination criteria variables
            SPTT = 0
            for idx, x in enumerate(OD_matrix):
                origin = np.int(x[0])
                destination = np.int(x[2])
                flow = x[1]
                # Calculate shortest paths
                backnode, costlabel = shortestPath(
                    origin=origin, capacity_array=capacity_array, weight_array=weight_array_iter)
                path, cost = pathTo(
                    backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
                # For shortest path travel time calculation (termination criteria)
                SPTT += cost*flow

            TSTT = np.sum(np.multiply(flow_array, weight_array_iter))
            AEC = (TSTT-SPTT)/sum_d
            RG = TSTT/SPTT - 1
            print(RG)
                        
            iter+=1

        #print(paths_array)
        # pass ---> I don't think this pass needs to be here


    elif algorithm == 'bush_based':
        

        # Initializing
        num_nodes = weight_array.shape[0]
        node_list = [num for num in range(0, num_nodes)]


        def topo_order_function(bush, origin, weight_array, num_nodes = num_nodes, node_list=node_list):
            """Function to calculate topographic order on a bush."""
            # variables
            visited = [False] * (num_nodes)

            # Find shortest path from origin to all locations - need to determine what nodes are unreachable 
            backnode, costlabel = shortestPath(origin, capacity_array, weight_array)
            # where backnode == -1, unreachable, and won't be considered in topological ordering (including the origin)
            unreachable = [i for i,j in enumerate(backnode) if j == -1]
            # creat visited labels 
            visited = [False]* (num_nodes)
            # set the unreachable nodes and origin nodes to True in order to skip them
            for val in unreachable:
                visited[val]=True
            
            # create adjacency table
            adj_table = {}
            for node in node_list:
                if visited[node] is True:
                    pass
                else:
                    adj_table[node] = [
                        i for i, j in enumerate(bush[node]) if j > 0]

            # conduct topological ordering from the origin node
            topo_order = []

            # Helper function for topo-sort
            def topo_sort_utils(node, visited, topo_order):
                """ Function to call to sort """
                # Mark current node as visited
                visited[node] = True
                # recur for all the nodes adjacent to this node
                for i in adj_table[node]:
                    if visited[i] is False:
                        topo_sort_utils(i, visited, topo_order)

                # append node to stack
                topo_order.append(node)

            for node in node_list:
                if visited[node] is False:
                    topo_sort_utils(node, visited, topo_order)

            topo_order.append(origin)
            topo_order = topo_order[::-1]

            return topo_order


        # Convert network into a set of bushes for each OD pair
        #topo_orders = []
        bushes = []
        bush_flows = []

        for pair in OD_matrix:
            origin = pair[0]
            flow = pair[1]
            destination = pair[2]

            # Create an empty  bush that is the same shape as weight_array
            bush = np.full_like(weight_array, -1)

            # Find shortest path from origin to all locations - need to determine what nodes are unreachable 
            backnode, costlabel = shortestPath(origin, capacity_array, weight_array)
            # where backnode == -1, unreachable, and won't be considered in topological ordering (including the origin)
            unreachable = [i for i,j in enumerate(backnode) if j == -1]
            # creat visited labels 
            visited = [False]* (num_nodes)
            # set the unreachable nodes and origin nodes to True in order to skip them
            for val in unreachable:
                visited[val]=True
            
            # create the reachable set - include the origin for Prim's algorithm to work properly
            reachable=[origin]
            for i in node_list:
                if i in unreachable:
                    pass
                else:
                    reachable.append(i)
            
            # Use Prim's algorithm to create minimum spanning tree from origin - append the bush
            prims_selected = [False]*(len(reachable))
            num_edges = 1
            # origin should be first number, so set equal to true
            prims_selected[0]=True
            # Begin Prim's algorithm
            while num_edges < (num_nodes - 1):
                minimum=999999999999
                x=0
                y=0
                for i in reachable:
                    if prims_selected[reachable.index(i)]:
                        for j in reachable:
                            if ((weight_array[i][j] > 0) and (not prims_selected[reachable.index(j)])):
                                if minimum > weight_array[i][j]:
                                    minimum = weight_array[i][j]
                                    x=i
                                    y=j
                #print(str(x) + "-" + str(y) + ":" + str(weight_array[x][y]))
                bush[x][y] = 1
                if y == 0:
                    pass
                else:
                    prims_selected[reachable.index(y)] = True
                num_edges+=1

            # Set initial flow for each bush - flow along shortest path
            bush_flow = np.zeros_like(weight_array)
            backnode, costlabel = shortestPath(origin, bush, weight_array)
            path, cost = pathTo(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
            i = 0
            while i < len(path)-1:
                bush_flow[path[i]][path[i+1]] += flow
                i += 1

            # append the trackers with topo orders and bushes
            #topo_orders.append(topo_order)
            bushes.append(bush)
            bush_flows.append(bush_flow)



        # sum flow on bushes for initial network flow
        flow_array = np.zeros_like(weight_array)
        for bush in bush_flows:
            flow_array += bush
        
        # calculate initial network travel times
        # FIXME : EQUATION WRITTEN: Change to fit BPR function THIS IS A MUST
        division_array = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                        [0, 0, 200, 0, 200, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 200, 0, 100, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 200, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 100, 0, 200, 0, 0],
                                        [0, 0, 0, 0, 0, 0, 100, 0, 100, 0],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 200],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 200, 0],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 100],
                                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        addition_array = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 3, 0, 3, 0, 0, 0, 0, 0],
                                [0, 0, 0, 3, 0, 5, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 3, 0, 0, 0],
                                [0, 0, 0, 0, 0, 5, 0, 3, 0, 0],
                                [0, 0, 0, 0, 0, 0, 5, 0, 5, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 3],
                                [0, 0, 0, 0, 0, 0, 0, 0, 3, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 5],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    
        # TODO change to BPR function
        weight_array_itter = (flow_array/division_array)**2+addition_array
        weight_array_itter = np.nan_to_num(weight_array_itter)
    

        # Begin the loops of, for each bush,:
            # 1. find if shortcuts exist (i.e., find new shortest path from O-D, if it uses links not on bush, add them)
                # 1a. find shortest path
                # 1b. add shortcuts to bush
            # 2. Calculate L and U labels in forward topological order
            # 3. conduct divergence node determination loop to shift flows using Newton's method
            # 4. After all Newton's adjustments for this bush, calculate new network flow and adjust travel times
            # repeat for each bush

        # ITER ETC. ETC. GOES HERE
        iter=0
        while iter <= 10:

            for idx, bush in enumerate(bushes):
                # set variables
                origin = OD_matrix[idx][0]
                flow = OD_matrix[idx][1]
                destination = OD_matrix[idx][2]
                #topo_order = topo_orders[idx]
                bush_flow = bush_flows[idx]

                # 1. FIND SHORTCUTS AND ADD TO BUSH
                # find shortest path from origin to destination
                backnode, costlabel = shortestPath(origin, capacity_array, weight_array_itter)
                path, cost = pathTo(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
                # append short cuts on bush if they exist
                id = 0
                while id < len(path)-1:
                    i=path[id]
                    j=path[id+1]
                    if bush[i][j] == -1:
                        bush[i][j] = 1
                    id += 1

                # 2. CALCULATE L AND U LABELS IN FORWARD TOPOLOGICAL ORDER
                # conduct topological ordering from the origin node
                topo_order = topo_order_function(bush, origin=origin, weight_array=weight_array_itter, num_nodes=num_nodes, node_list=node_list)
                # Set L and U variables
                L_link = np.full_like(weight_array, np.inf)
                L_node = np.full(len(topo_order), np.inf)
                U_link = np.full_like(weight_array, -np.inf)
                U_node = np.full(len(topo_order), -np.inf)
                # L and U at r (origin) are equal to 0
                L_node[0] = 0
                U_node[0] = 0

                # determine labels in forward topological ordering
                id=0
                while id <= len(topo_order)-1:
                    
                    # for first in topological order
                    i = topo_order[id]
                    if id==0:
                        # TODO: I think this can be optimized to only search for values in topo_order in position higher than what we are currently are at - the whole benefit of using topo order also check on this feature below in this same loop - CAN MAKE IT BETTER
                        for j in topo_order:
                            if bush[i][j] == 1:
                                L_link[i][j] = weight_array_itter[i][j]
                                U_link[i][j] = weight_array_itter[i][j]
                    
                    # for last in topological order
                    elif id == len(topo_order)-1:
                        l_val = L_node[id]
                        u_val = U_node[id]
                        for h in topo_order:
                            if bush[h][i] == 1:
                                l_val = min(l_val, L_link[h][i])
                                u_val = max(l_val, U_link[h][i])
                        L_node[id] = l_val
                        U_node[id] = u_val
                    
                    # for all others
                    else:
                        j=topo_order[id+1]
                        l_val = L_node[id]
                        u_val = U_node[id]
                        for h in topo_order:
                            if bush[h][i] == 1:
                                l_val = min(l_val, L_link[h][i])
                                u_val = max(u_val, U_link[h][i])
                        L_node[id] = l_val
                        U_node[id] = u_val
                        for j in topo_order:
                            if bush[i][j] == 1:
                                L_link[i][j] = min(L_link[i][j], L_node[id]+weight_array_itter[i][j])
                                U_link[i][j] = max(U_link[i][j], U_node[id]+weight_array_itter[i][j])
                            
                    id += 1
                # remove infs for clarity
                U_link[np.isinf(U_link)] = 0
                L_link[np.isinf(L_link)] = 0



                # DIVERGENCE NODE LOOP
                loop = True
                while loop is True:
                    # determine topoligcally last node in the bush
                    i = topo_order[-1]
                    # determine the longest and shortest paths from U and L Labels
                    max_val = -np.inf
                    min_val = np.inf
                    max_path = [i]
                    min_path = [i]
                    for h in topo_order:
                        # TODO: same todo as last loop, i think this can be optimized due to topological nature of list
                        if bush[h][i] == 1:
                            if max_val < U_link[h][i]:
                                max_val = U_link[h][i]
                                next_node_max = h
                            if min_val > L_link[h][i]:
                                min_val = L_link[h][i]
                                next_node_min = h
                    max_path.append(next_node_max)
                    min_path.append(next_node_min)

                    # loop for max path
                    while max_path[-1] != origin:
                        max_val = -np.inf
                        for h in topo_order:
                            if bush[h][next_node_max] == 1: 
                                if max_val < U_link[h][next_node_max]:
                                    max_val = U_link[h][next_node_max]
                                    next_node_holder = h
                        max_path.append(next_node_holder)
                        next_node_max = next_node_holder

                    # loop for min path
                    while min_path[-1] != origin:
                        min_val = np.inf
                        for h in topo_order:
                            if bush[h][next_node_min] == 1: 
                                if min_val > L_link[h][next_node_min]:
                                    min_val = L_link[h][next_node_min]
                                    next_node_holder = h

                        min_path.append(next_node_holder)
                        next_node_min = next_node_holder

                    # reverse order of lists for better comprehension
                    min_path = min_path[::-1]
                    max_path = max_path[::-1]

                    # if max and min path are the same, there is only one path and therefore we do not need to shift flow
                    if min_path == max_path:
                        loop=False
                    
                    # determine divergence node "a" from min and max paths
                    # Divergence node is the last node common to both
                    else:
                        for val in min_path[-2::-1]:
                            if val in max_path[-2::-1]:
                                a = val
                                break

                        # TODO: I don't know if this below alteration of i and a is valid or fixing the negative delta h problem
                        # need to determine if set of a to i only has a single path: if so, a and I need to shift back common nodes
                        separation = min_path.index(i) - min_path.index(a) 
                        while separation == 1:
                            i = copy.deepcopy(a)
                            new_index_min = min_path.index(a)-1
                            new_index_max = max_path.index(a)-1
                            for val in min_path[new_index_min::-1]:
                                if val in max_path[new_index_max::-1]:
                                    a=val
                            separation = min_path.index(i) - min_path.index(a)


                    
                        # extract sigma_l and sigma_u, the list of links from a to i
                        sigma_U = max_path[max_path.index(a):(max_path.index(i)+1)]
                        sigma_L = min_path[min_path.index(a):(min_path.index(i)+1)]

                        # calculate delta h, equation 6.61 in the textbook
                        numerator = (U_node[topo_order.index(i)] - U_node[topo_order.index(a)]) - (L_node[topo_order.index(i)] - L_node[topo_order.index(a)])
                    
                        # TODO: need to optimize the calculation of denominator based on if links have the same or different BPR
                        # FOR NOW: Just calcualting the whole matrix due to the sample data being used
                        # FIXME: EQUATION WRITTEN: Fix BPR function
                        derivative_division_array = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                                    [0, 0, 20000, 0, 20000, 0, 0, 0, 0, 0],
                                                    [0, 0, 0, 20000, 0, 5000, 0, 0, 0, 0],
                                                    [0, 0, 0, 0, 0, 0, 20000, 0, 0, 0],
                                                    [0, 0, 0, 0, 0, 5000, 0, 20000, 0, 0],
                                                    [0, 0, 0, 0, 0, 0, 5000, 0, 5000, 0],
                                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 20000],
                                                    [0, 0, 0, 0, 0, 0, 0, 0, 20000, 0],
                                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 5000],
                                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])


                        derivative_array = flow_array/derivative_division_array
                        denominator = 0

                        id = 0
                        while id < len(sigma_L)-1:
                            denominator += derivative_array[sigma_L[id]][sigma_L[id+1]]
                            id += 1  
                        id = 0
                        while id < len(sigma_U)-1:
                            denominator += derivative_array[sigma_U[id]][sigma_U[id+1]]
                            id += 1

                        # set delta_h value
                        delta_h = numerator/denominator

                        # determine delta_h based on minimum constraints - only have to look at 
                        id = 0
                        while id < len(sigma_U)-1:
                            if delta_h <= flow_array[sigma_U[id]][sigma_U[id+1]]:
                                pass
                            else:
                                delta_h = flow_array[sigma_U[id]][sigma_U[id+1]]
                            id += 1

                        # Shift flows 
                        id = 0
                        while id < len(sigma_L)-1:
                            flow_array[sigma_L[id]][sigma_L[id+1]] += delta_h
                            id += 1
                        id = 0
                        while id < len(sigma_U)-1:
                            flow_array[sigma_U[id]][sigma_U[id+1]] -= delta_h
                            id += 1



                        # at end of loop, slice off last value and repeat
                        topo_order = topo_order[:-1]
                    
            
                    
                # update the weight array
                # FIXME : EQUATION WRITTEN: Change to BPR function
                weight_array_itter = (flow_array/division_array)**2+addition_array
                weight_array_itter = np.nan_to_num(weight_array_itter)





            # calculate termination criteria
            SPTT = 0
            for idx, x in enumerate(OD_matrix):
                origin = np.int(x[0])
                destination = np.int(x[2])
                flow = x[1]
                # Calculate shortest paths
                backnode, costlabel = shortestPath(
                    origin=origin, capacity_array=capacity_array, weight_array=weight_array_itter)
                path, cost = pathTo(
                    backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
                # For shortest path travel time calculation (termination criteria)
                SPTT += cost*flow

            # termination criteria
            TSTT = np.sum(np.multiply(flow_array, weight_array_itter))
            AEC = (TSTT-SPTT)/sum_d
            RG = TSTT/SPTT - 1
            print(RG)

            iter +=1
        #print(np.around(flow_array,1))

    test='output'
    return test
