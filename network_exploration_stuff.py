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
# TODO: In plot aoi - allow it to accept a list of res parcels
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



# TODO: Write a function that finds nearest points
# should put in one location a function that, based on a given input, that returns locations of nearest points on network
# this way can have in one location the multiple methodologies on how I think this will be done



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
        G.add_edge(x, 99999999, G_weight=0)
    # run the min_cost_flow function to retrieve FlowDict
    try:
        flow_dictionary = nx.min_cost_flow(
            G, demand=G_demand, weight=G_weight, capacity=G_capacity)
        # run the cost_of_flow function to retrieve total cost
        cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
        return flow_dictionary, cost_of_flow
    except nx.NetworkXUnfeasible:
        return None, None


def max_flow_parcels(G='', res_points='', dest_points='', G_capacity='capacity', G_weight='travel_time'):
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
    :param dest_points: Point locations of all of the resources the random residential parcel is to be routed to
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
    # add artifical origin node
    G.add_nodes_from([(0)])
    # find the number of unique origin nodes
    # sums will be equal to number of res parcels
    sums = 0
    unique_origin_nodes = np.unique(origins)
    for unique_node in unique_origin_nodes:
        count = np.count_nonzero(origins == unique_node)
        sums += count
        kwargs = {f"{G_weight}": 0, f"{G_capacity}": count}  # TODO: KWARGS
        # Create the aritifical source edges from the artificial origin, 0, to all origin nodes
        # The capacity of each edge is the maximum demand
        G.add_edge(0, unique_node, **kwargs)

    # determine the unique destination nodes
    unique_dest_nodes = np.unique(destinations)
    # Graph must be in digraph format: convert multidigraph to digraph
    G = nx.DiGraph(G)
    # create artificial destination node
    #    All sinks will go to this demand in order to balance equation
    G.add_nodes_from([(99999999)])
    # add edges from sinks to this artificial node
    for x in unique_dest_nodes:
        G.add_edge(x, 99999999, G_weight=0)

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
