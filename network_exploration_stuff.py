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
from rasterio.plot import show
import matplotlib.pylab as pl
from matplotlib import pyplot
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
from rasterstats import zonal_stats
from matplotlib_scalebar.scalebar import ScaleBar
import random
import copy
import heapq
from collections import defaultdict
import rasterio
import rasterio.mask
import fiona
import shapely
import time
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import Pool as Poolm
from multiprocessing import Process
import math
import sys
import matplotlib as mpl


from scipy import sparse

# set some panda display options
pd.set_option('display.max_columns', None)

# CAPACITY TODOS
# TODO: Change disruption impact on capacity - remove nodes that cannot support travel
# TODO: Potentially - add capacity limitation to min cost flow calculations (i.e., ideal travel functions )


# SPEED DISRUPTION TODOS
# TODO: input new equations for speed disruption, possibly add bining feature


# POTENTIAL LIMITATION TODOS
# TODO: Some of the closest origin nodes are the destination nodes - for now assuming that no paths need to be calculated (i.e., can just walk)


# IMPROVEMENT TODOS
# TODO: move digraph conversion into better spot and determine when it is needed
# TODO: Move CRS check into functions
# TODO: Min cost flow and max cost flow functions have some repeat code - I could condense a lot of the work
# TODO: Use kwargs methodology to create new edges with appropriate attributes
#   see # TODO: KWARGS below in code for example
# TODO: read_graph_from_disk - problem with preserving dtypes so I hardcoded a
#   temporary fix - right now it goes through the attributes and determines if
#    it could be converted to a float - if it can it does
#   this involves manually skipping oneway as well (i think b/c its bool)
#   this was needed b/c  inundation_G load not preserve datatypes
# TODO: update shortest path and path to functions to the new modified ones - much quicker
# TODO: add impact of rainfall intensity on reduction of speed and intensity
    # light rain reduces free flow speed 2-13% and capacity 4-11%
    # heavy rain reduces speed 3-16% and capacity 10-30%
    # see these sites for info: 
    # https://engrxiv.org/preprint/view/1322/2765
    # https://ops.fhwa.dot.gov/weather/best_practices/AMS2003_TrafficFlow.pdf 
    # https://ops.fhwa.dot.gov/weather/q1_roadimpact.htm#:~:text=Light%20rain%20can%20decrease%20freeway,12%20percent%20in%20low%20visibility.

# OTHER TODOS
# TODO: Need to use different method to determine random node variable name in the min cost flow function
# TODO: In plot aoi - allow it to accept a list of dest parcels
#   i.e., food and gas locations# TODO: convert parcel identificaiton (access vs. no access) to function




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


# THE FUNCTIONS #############################################################


def shape_2_graph(source):
    '''
    Extracts road network from shapefile.

    This function takes an input path to a .shp file (assuming the associated .shp, .shx, .dbf, etc. are all in the same folder)
    and extracts/returns the OpenStreetMap road network. This function automatically pulls the speed and travel time information,
    which are necessary components in further calculations. 

    **Important Notes**
    
    Ensure that the shapefile area has a sufficent buffer to capture enough of the road network.
    For example, if the AOI is a block group, residents might require driving outside bounadry to get to areas within
    and excluding roadnetworks, espcially those that border the AOI, can lead to errors.

    Only potential improvement is to include function to add a buffer to AOI instead of requiring that the user automatically do this, but this is not the most meaningful improvement


    :param source: Path to .shp file of the area of interest
    :type source: str
    :return: **G**, the resultant road network from OpenStreetMap
    :rtype: networkx.MultiDiGraph
    
    '''

    # extracting from OSMNx requires crs to be EPSG:4326(WGS84). This will preserve output crs if projected
    boundary = gpd.read_file(source)
    crs_proj = boundary.crs
    boundary=boundary.to_crs("EPSG:4326")

    G = ox.graph_from_polygon(
        boundary['geometry'][0],
        network_type='drive')

    # dictionary of highway speeds in kph
    hwy_speeds = {"motorway": 120,
                  "motorway_link": 120,
                  "trunk": 120,
                  "trunk_link": 120,
                  "primary": 65,
                  "primary_link": 65,
                  "secondary": 50,
                  "secondary_link": 50,
                  "tertiary": 50,
                  "tertiary_link": 50,
                  "residential": 40,
                  "minor": 40,
                  "unclassified": 40,
                  "living_street": 40}
    # add edge speeds where they exist
    ox.add_edge_speeds(G, hwy_speeds=hwy_speeds, fallback=30)
    # add travel times 
    ox.add_edge_travel_times(G)

    # some highways types have multiple inputs (e.g., ['residential','unclassified']), which messes with capacity and width assignments
    # determine which edges have multiple road types and select the first from the list to be that value
    edge_rtype = nx.get_edge_attributes(G, 'highway')
    edges = [k for k, v in edge_rtype.items() if not isinstance(v, str)]
    for edge in edges:
        list = G.get_edge_data(edge[0], edge[1], edge[2])['highway']
        G[edge[0]][edge[1]][edge[2]]['highway'] = list[0]
    
    # updated edge attributes variable to update lane capacities
    edge_rtype = nx.get_edge_attributes(G, 'highway')
    # set per lane capacities based on following scheme:
        # Motorway & its links                                  : 2300 pcu/hr/ln
        # Trunk & its links                                     : 2300 pcu/hr/ln
        # Primary & its links                                   : 1700 pcu/hr/ln
        # Secondary & its links                                 : 1500 pcu/hr/ln
        # Tertiary & its links                                  : 1000 pcu/hr/ln
        # Residential, minor, living street, and unclassified   :  600 pcu/hr/ln
    for k, v in edge_rtype.items():
        if v == 'motorway' or v == 'motorway_link' or v == 'trunk' or v == 'trunk_link':
            G[k[0]][k[1]][k[2]]['capacity'] = 2300
        elif v == 'primary' or v == 'primary_link':
            G[k[0]][k[1]][k[2]]['capacity'] = 1700
        elif v == 'secondary' or v == 'secondary_link':
            G[k[0]][k[1]][k[2]]['capacity'] = 1500
        elif v == 'tertiary' or v == 'tertiary_link':
            G[k[0]][k[1]][k[2]]['capacity'] = 1000
        elif v == 'residential' or v == 'minor' or v == 'living_street' or v == 'unclassified':
            G[k[0]][k[1]][k[2]]['capacity'] = 600


    # some roads gain/lose lanes during uninterupted segments, and therefore lanes are reported as a list (e.g., ['3', '4', '2'])
    # determine which edges have multiple number of lanes and replace with the smallest one available 
    edge_lane = nx.get_edge_attributes(G, 'lanes')
    edges = [k for k, v in edge_lane.items() if not isinstance(v, str)]
    for edge in edges:
        list = G.get_edge_data(edge[0], edge[1], edge[2])['lanes']
        list_int = [int(i) for i in list]
        minimum = min(list_int)       
        G[edge[0]][edge[1]][edge[2]]['lanes'] = minimum

    # convert all lanes to int, for some reason reported as strings
    edge_lane = nx.get_edge_attributes(G, 'lanes')
    edge_lane_int = dict([k, int(v)] for k, v in edge_lane.items())
    nx.set_edge_attributes(G, edge_lane_int, 'lanes')

    # Not all edges have lane information, need to fill in gaps and determine road widths
    # Factor of the type of road and if it is one way or not
    # Number of lanes ONE WAY will go as follows:
        # Trunk:                                                4
        # Motorway:                                             4
        # Primary:                                              3
        # Secondary:                                            2
        # Tertiary:                                             1
        # Residential, minor, unclassified, living street:      1
        # all link roads (typically on ramps or slip lanes):    1
    # AVERAGE LANE WIDTH IN METERS
    lane_width = 3.7
    edge_lane = nx.get_edge_attributes(G, 'lanes')
    edge_oneway = nx.get_edge_attributes(G, 'oneway')
    edge_rtype = nx.get_edge_attributes(G, 'highway')


    for k,v in edge_rtype.items():
        # for oneway roads...
        if G[k[0]][k[1]][k[2]]['oneway'] is True:
            # that have OSMNx lane information
            if 'lanes' in G[k[0]][k[1]][k[2]]:
                G[k[0]][k[1]][k[2]]['width'] = lane_width*G[k[0]][k[1]][k[2]]['lanes']
                G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']
            # need lane information
            else:
                if v == 'trunk' or v == 'motorway' :
                    G[k[0]][k[1]][k[2]]['lanes'] = 4
                    G[k[0]][k[1]][k[2]]['width'] = 4* lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']
                elif v == 'primary':
                    G[k[0]][k[1]][k[2]]['lanes'] = 3
                    G[k[0]][k[1]][k[2]]['width'] = 3 * lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']
                elif v == 'secondary':
                    G[k[0]][k[1]][k[2]]['lanes'] = 2
                    G[k[0]][k[1]][k[2]]['width'] = 2 * lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']
                elif v == 'tertiary' or v == 'residential' or v == 'minor' or v == 'unclassified' or v == 'living_street' or v == 'trunk_link' or v == 'motorway_link' or v == 'primary_link' or v == 'secondary_link' or v == 'tertiary_link':
                    G[k[0]][k[1]][k[2]]['lanes'] = 1
                    G[k[0]][k[1]][k[2]]['width'] = 1 * lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']
        # for roads not listed as one-way roads...
        else:
            # that have OSMNx lane information
            if 'lanes' in G[k[0]][k[1]][k[2]]:
                G[k[0]][k[1]][k[2]]['width'] = lane_width*G[k[0]][k[1]][k[2]]['lanes']
            # need lane information
            else:
                if v == 'trunk' or v == 'motorway' :
                    G[k[0]][k[1]][k[2]]['lanes'] = 8
                    G[k[0]][k[1]][k[2]]['width'] = 8* lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']/2
                elif v == 'primary':
                    G[k[0]][k[1]][k[2]]['lanes'] = 6
                    G[k[0]][k[1]][k[2]]['width'] = 6 * lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']/2
                elif v == 'secondary':
                    G[k[0]][k[1]][k[2]]['lanes'] = 4
                    G[k[0]][k[1]][k[2]]['width'] = 4 * lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']/2
                elif v == 'tertiary' or v == 'residential' or v == 'minor' or v == 'unclassified' or v == 'living_street':
                    G[k[0]][k[1]][k[2]]['lanes'] = 2
                    G[k[0]][k[1]][k[2]]['width'] = 2 * lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']/2
                elif v == 'trunk_link' or v == 'motorway_link' or v == 'primary_link' or v == 'secondary_link' or v == 'tertiary_link':
                    G[k[0]][k[1]][k[2]]['lanes'] = 1
                    G[k[0]][k[1]][k[2]]['width'] = 1 * lane_width
                    G[k[0]][k[1]][k[2]]['capacity'] *= G[k[0]][k[1]][k[2]]['lanes']/2
    

    # project back to original crs
    G=ox.project_graph(G, crs_proj)


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
    ox.io.save_graphml(G, f'{path}/{name}')


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
    
    G = ox.io.load_graphml(f'{path}/{name}')
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
 
    # append res_points with nearest network node
    longs=res_points.geometry.x
    lats=res_points.geometry.y
    origins = ox.distance.nearest_nodes(G, longs, lats)
    res_points['nearest_node'] = origins

    # will find the nearest node based solely on the centroid of resource parcel (i.e., dest_points)
    longs=dest_points.geometry.x
    lats=dest_points.geometry.y
    destinations = ox.distance.nearest_nodes(G, longs, lats)
    dest_points['nearest_node'] = destinations

    # Create the demand attribute and set all equal to 0
    nx.set_node_attributes(G, values=0, name=G_demand)

    # counting taking too long trying a new approach
    from collections import Counter
    #unique_origin_nodes_counts = dict(Counter(origins).items())
    unique_origin_nodes_counts = dict(zip(list(Counter(origins).keys()), [x*-1 for x in Counter(origins).values()]))
    
    # Create list of unique origins and destinations
    unique_origin_nodes = np.unique(origins)
    unique_dest_nodes = np.unique(destinations)
    
    # determine shared nodes and set demand equal to 0
    shared_nodes = list(set(unique_origin_nodes) & set(unique_dest_nodes))
    for x in shared_nodes:
        unique_origin_nodes_counts[x] = 0

    # add source information (the negative demand value prev. calculated)
    nx.set_node_attributes(G, unique_origin_nodes_counts, G_demand)

    # Calculate the positive demand: the sink demand
    positive_demand = sum(unique_origin_nodes_counts.values())*-1

    # TODO: I don't think returning res_points or dest_points is neccessary and can be removed - ACTUALLY MIGHT BE USED IN FLOW DECOMP
    return(G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points)


def nearest_nodes_vertices(G='', res_points='', dest_parcels='', G_demand='demand', dest_points=None):
    '''
    MODIFICATION OF NEAREST_NODES FUNCTION TO SNAP DEST PARCELS TO MULTIPLE NODES
    
    Return Graph with demand attribute based on snapping sources/sinks to nearest nodes (intersections).
    
    The purpose of this function is take a series of input locations and determine the nearest nodes (intersections) within a Graph network. It takes this information to create source and sink information. 

    **Important Note**: The output graph has an attribute labeled G_demand, which shows the number of closest residential parcels to each unique intersection. Because these are considered sources, they have a negative demand. Sink locations, or the intersections that are closest to the the dest_points, will not have a demand value (G_demand == 0) because we do not know how much flow is going to that node until after a min-cost flow algorithm is run. I.e., to route properly, we end up creating an artifical sink and then decompose the flow to determine how much flow is going to each destination.

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_parcels: Parcels of all of the resources the residential parcels are to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param G_demand: name of attribuet in G refering to node demand, *Default='demand*
    :type G_demand: string
    :param dest_points: Optional, if not None, will output with added columns on the nearest nodes
    :type dest_points: geopandas dataframe
    :returns: 
        - **G**, *networkx.DiGraph* with G_demand attribute showing source values (negative demand)
        - **unique_origin_nodes**, *lst* of unique origin nodes
        - **unique_dest_nodes**, *lst* of unqiue destination nodes
        - **positive_demand**, *int* of the total positive demand across the network. **Does not include the residents that share a closest intersection with a resource parcel.**
        - **Shared_nodes**, *lst* of nodes that res and destination parcels share same nearest 
        - **res_points**, *geopandas.GeoDataFrame*, residential points with appended attribute of 'nearest_node'
        - **dest_parcels**, *geopandas.GeoDataFrame*, destination/sink footprints with appended attribute of 'nearest_node'
        - **dest_points**, *geopandas.GeoDataFrame*, destination/sink points with appended attribute of 'nearest_node'
    :rtype: tuple

    '''

    # create empty lists for nodes of origins and destinations
    origins = []
    destinations = []

    # append res_points with nearest network node
    longs = res_points.geometry.x
    lats = res_points.geometry.y
    origins = ox.distance.nearest_nodes(G, longs, lats)
    res_points['nearest_node'] = origins

    # Will determine vertices of dest parcels and find multiple nearest nodes for destinations
    # simplify dest_parcels to reduce geometry
    dest_simp = dest_parcels.simplify(3, preserve_topology=False)
    # account for multipolygons by creating convex hulls of each parcel
    convex_hull = dest_simp.convex_hull
    # convert back to geodataframe
    dest_simp = gpd.GeoDataFrame({'geometry': convex_hull})
    # create nested list of vertices to examine for each parcel
    coords = [list(dest_simp.geometry.exterior[row_id].coords)
                for row_id in range(dest_simp.shape[0])]

    # initiate lists
    # nested list of destinations
    destinations_int = []
    # list of destinations in string form, each string is all detinations for one parcel
    destinations_str = []
    for idx, coord in enumerate(coords):
        # BUG maybe: in BPA, coords has 3 values, AN only has two. Not sure why, should examine, for now built in loop to bypass
        if len(coord[0]) == 2:
            longs = [i for i, j in coord]
            lats = [j for i, j in coord]
        elif len(coord[0]) == 3:
            longs = [i for i, j, z in coord]
            lats = [j for i, j, z in coord]
        dests = ox.distance.nearest_nodes(G, longs, lats)
        destinations_int.append(dests)
        destinations_str.append(' '.join(str(x) for x in dests))
    dest_parcels['nearest_nodes'] = destinations_str
    if dest_points is not None:
        dest_points=gpd.sjoin(dest_points, dest_parcels)
        dest_points.drop(columns=['index_right'], inplace=True)
    # create single list of destinations
    destinations = [element for innerList in destinations_int for element in innerList]

    # Create the demand attribute and set all equal to 0
    nx.set_node_attributes(G, values=0, name=G_demand)

    # counting taking too long trying a new approach
    from collections import Counter
    #unique_origin_nodes_counts = dict(Counter(origins).items())
    unique_origin_nodes_counts = dict(zip(list(Counter(origins).keys()), [
                                      x*-1 for x in Counter(origins).values()]))

    # Create list of unique origins and destinations
    unique_origin_nodes = np.unique(origins)
    unique_dest_nodes = np.unique(destinations)
    unique_dest_nodes_list = [np.unique(nodes) for nodes in destinations_int]

    # determine shared nodes and set demand equal to 0
    shared_nodes = list(set(unique_origin_nodes) & set(unique_dest_nodes))
    for x in shared_nodes:
        unique_origin_nodes_counts[x] = 0

    # add source information (the negative demand value prev. calculated)
    nx.set_node_attributes(G, unique_origin_nodes_counts, G_demand)

    # Calculate the positive demand: the sink demand
    positive_demand = sum(unique_origin_nodes_counts.values())*-1

    if dest_points is None:
        return(G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels)
    else:
        return(G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points)


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


def min_cost_flow_parcels(G='', res_points='', dest_points='', dest_parcels='', G_demand='demand', G_capacity='capacity', G_weight='travel_time', dest_method='single'):
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

    This function also does not factor in limits to a resource. In other words, there is no limit on the number of residential parcels
    that can go to any single resource location. **This feature also needs further exploration.**

    :param G: Graph network. Will be converted to DiGraph in function
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the random residential parcel is to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param dest_parcels: technically only required if dest_method == 'multiple', destination parcel shapefile
    :param dest_points: geopandas.GeoDataFrame
    :param G_demand: name of attribute in G refering to node demands, *Default='demand'* 
    :type G_demand: string
    :param G_capacity: name of attribute in G refering to edge capacities, *Default='capacity'* 
    :type G_capacity: string
    :param G_weight: name of attribute in G refering to edge weights, *Default='travel_time'* 
    :type G_weight: string  
    :param dest_method: either 'single', or 'multiple', determines which destination point methodology to use. Either connecting to nearest single node or multiple nearest nodes
    :type dest_method: string
    :returns: 
        - **flow_dictionary**, dictionary of dictionaries keyed by nodes for edge flows
        - **cost_of_flow**, integer, total cost of all flow. If weight/cost is travel time, this is the total time for everyone to reach the resource (seconds).
    :rtype: tuple
    :Raises: nx.NetworkXUnfeasible if all demand cannot be satisfied, i.e., all residential parcels cannot reach the resource.

    ''' 
    # for many functions to work, graph needs to be a digraph (NOT a multidigraph) i.e., no parallel edges
    # TODO factor in check of parallel edges. For now, just convert to digraph
    G = nx.DiGraph(G)
    # Travel times must be whole numbers -  round values if not whole numbers
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])
    
    # if snapping destination parcels to nearest singular node only
    if dest_method == 'single':
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

    elif dest_method == 'multiple':
        # Find the source and sink demands and append to graph G
        G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels = nearest_nodes_vertices(G=G, res_points=res_points, dest_parcels=dest_parcels, G_demand=G_demand)
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

        # add the super_sink that everything goes to 
        G.add_nodes_from([(99999999, {G_demand: positive_demand})])

        # identify the destination nodes, and create artifical sink edges
        # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to 
        for idx, dest_parcel in dest_parcels.iterrows():
            dest_node = dest_points[dest_points.geometry.within(dest_parcel.geometry)]
            #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
            for i, node in dest_node.iterrows():
                # add the dest node to the graph using OSMID as its ID
                # TODO: explore option about plotting with physical lines
                # x = node.geometry.x
                # y = node.geometry.y
                # G.add_nodes_from([(node['osmid'], {'demand': 0, 'x':x, 'y':y})])
                G.add_nodes_from([(node['osmid'], {G_demand: 0})])

                # add links from nearest intersections to parcel centroid
                for nearest_intersection in unique_dest_nodes_list[idx]:
                    kwargs = {G_weight: 0, G_capacity: 999999999}  # TODO: KWARGS
                    G.add_edge(nearest_intersection, node['osmid'], **kwargs)

                # add link from parcel centroid to super sink
                # TODO: This is where I can specifically add capacities for each specific grocery store
                # FOR NOW: just setting capacity to a high amount 
                kwargs = {G_weight: 0, G_capacity: 999999999}
                G.add_edge(node['osmid'], 99999999, **kwargs)


    # run the min_cost_flow function to retrieve FlowDict
    try:
        flow_dictionary = nx.min_cost_flow(
            G, demand=G_demand, weight=G_weight, capacity=G_capacity)
        # run the cost_of_flow function to retrieve total cost
        cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
        return flow_dictionary, cost_of_flow
    except nx.NetworkXUnfeasible:
        return None, None


def max_flow_parcels(G='', res_points='', dest_points='', G_capacity='capacity', G_weight='travel_time', G_demand='demand', dest_method='single', dest_parcels=None, ignore_capacity=False):
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
    :param dest_method: either 'multiple' or 'single', determines snapping mechanism for destination points, if it uses multiple vertices or just the nearest vertice
    :type dest_method: string
    :param dest_parcels: If dest_method =='multiple', need the destination resource parcels as an input
    :type dest_parcels: geopandas.GeoDataFrame
    :param ignore_capacity: if True, road link capacities will be set to infinity. If False, original road capacities are used
    :type ignore_capacity: bool
    :returns: 
        - **flow_dictionary**, dictionary of dictionaries keyed by nodes for edge flows
        - **cost_of_flow**, integer, total cost of all flow. If weight/cost is travel time, this is the total time for everyone to reach the resource (seconds).
        - **max_flow**, integer, the maximum flow (i.e., households) that can access the resource
        - **access**, string, either 'Complete' or'Partial', representing if every residential parcel can/cannot access the resource
    :rtype: tuple

    '''

    # for many functions to work, graph needs to be a digraph (NOT a multidigraph) i.e., no parallel edges
    # TODO factor in check of parallel edges. For now, just convert to digraph
    G = nx.DiGraph(G)

    # Based on input, change capacities 
    if ignore_capacity is True:
        nx.set_edge_attributes(G, 999999999999, name=G_capacity)

    # Travel times must be whole numbers - just round values
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])

    if dest_method == 'single':
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
    
    elif dest_method == 'multiple':
        # Find the source and sink demands and append to graph G
        G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels = nearest_nodes_vertices(G=G, res_points=res_points, dest_parcels=dest_parcels, G_demand=G_demand)
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

        # add the super_sink that everything goes to 
        G.add_nodes_from([(99999999, {G_demand: positive_demand})])

        # identify the destination nodes, and create artifical sink edges
        # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to 
        for idx, dest_parcel in dest_parcels.iterrows():
            dest_node = dest_points[dest_points.geometry.within(dest_parcel.geometry)]
            #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
            for i, node in dest_node.iterrows():
                # add the dest node to the graph using OSMID as its ID
                # TODO: explore option about plotting with physical lines
                # x = node.geometry.x
                # y = node.geometry.y
                # G.add_nodes_from([(node['osmid'], {'demand': 0, 'x':x, 'y':y})])
                G.add_nodes_from([(node['osmid'], {G_demand: 0})])

                # add links from nearest intersections to parcel centroid
                for nearest_intersection in unique_dest_nodes_list[idx]:
                    kwargs = {G_weight: 0, G_capacity: 999999999}  # TODO: KWARGS
                    G.add_edge(nearest_intersection, node['osmid'], **kwargs)

                # add link from parcel centroid to super sink
                # TODO: This is where I can specifically add capacities for each specific grocery store
                # FOR NOW: just setting capacity to a high amount 
                kwargs = {G_weight: 0, G_capacity: 999999999}
                G.add_edge(node['osmid'], 99999999, **kwargs)


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
                    background_edges=None,
                    edge_width=None, 
                    edge_color=None,
                    bbox=None, 
                    loss_access_parcels=None, 
                    scalebar=False,
                    inundation=None, 
                    insets=None, 
                    save_loc=None,
                    raster=None,
                    decomp_flow=False):
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
    :param edge_color: attribute in G to color edges. Designed to show max inundation on the road
    :type edge_color: string
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
    :param raster: raster of inundation extent
    :type raster: .tif file, use rasterio.open()
    :param decomp_flow: if True, plot residential parcels with color ramp symbolizing travel time
    :type decomp_flow: bool

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

    fig, ax = plt.subplots(facecolor='white', figsize=(12,12))
    ax.axis('off')
    
    # plot roads, residential parcels, and resource parcels
    # if decomp_flow is True, plot res parcels with scale bar to show travel times
    if decomp_flow is True:
        res_parcels.plot(ax=ax, color='tan') # alternative color is antiquewhite
        # # IQR METHOD to mask impact of outliers
        costs = sorted(res_parcels['cost_of_flow'].tolist())
        costs = [x for x in costs if math.isnan(x)==False]
        q1,q3,=np.percentile(costs, [25,75])
        iqr=q3-q1
        upper_bound=q3+(3*iqr)
        res_parcels.loc[res_parcels['cost_of_flow']>=upper_bound, ['cost_of_flow']]=upper_bound
        # res_parcels.plot(ax=ax, column='cost_of_flow',cmap='bwr', legend=True)
        res_parcels.plot(ax=ax, column='cost_of_flow',cmap='bwr', legend=False)
        # add colorbar
        # ax_cbar = fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=upper_bound), cmap='bwr'),ax=ax,
        #                         ticks=np.linspace(0,upper_bound,6))
        
        
        #res_parcels.plot(ax=ax, column='cost_of_flow',cmap='bwr', scheme='natural_breaks',k=10, legend=True)
    else:
        res_parcels.plot(ax=ax, color='tan')# alternative color is antiquewhite
    
    # plot the resource parcels
    resource_parcels.plot(ax=ax, color='mediumseagreen', edgecolor='darkgreen', linewidth=0.5)

    # #TODO: ADD FEATURE TO HIGHLIGHT NODES - USEFUL IN ORDER TO TRANSLATE RESULTS TO GRAPH QUICKLY WHILE TESTING 
    # G_gdf_nodes.loc[[6016], 'geometry'].plot(ax=ax, color='green')
    # G_gdf_nodes.loc[[10727], 'geometry'].plot(ax=ax, color='green')
    # G_gdf_nodes.loc[[7196], 'geometry'].plot(ax=ax, color='green')
    # #### END TEST

    # option to plot loss of access parcels
    if loss_access_parcels is None:
        pass
    else:
        loss_access_parcels.plot(ax=ax, color='saddlebrown')

    # plot background light gray roads
    # plot different background edges if value given
    if background_edges is None:
        ox.plot.plot_graph(G,
                            ax=ax,
                            node_size=0, node_color='black',
                            edge_color='lightgray',
                            show=False,
                            edge_linewidth=1,
                            bbox=total_bounds)
    else:
        ox.plot.plot_graph(background_edges,
                            ax=ax,
                            node_size=0, node_color='black',
                            edge_color='lightgray',
                            show=False,
                            edge_linewidth=1,
                            bbox=total_bounds)

    # plot edges based on color (typically max inundation on the road)
    if edge_color is None:
        pass
    else:
        # extract the values to plot colors to
        vals = pd.Series(nx.get_edge_attributes(G, edge_color))
        # set the bounds of the color
        bounds=[0,10,50,100,150,200,250,300,350,400,500,600,vals.dropna().max()]
        
        #define the colormap
        cmap = plt.cm.Blues  

        # extract all colors from the .jet map
        cmaplist = [cmap(i) for i in range(cmap.N)]
        # force the first color entry to be clear and only show underlying roads
        cmaplist[0] = (.5, .5, .5, 0.0)

        # create the new cmap
        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            'Custom cmap', cmaplist, cmap.N)

        bounds=[0,10,150,300,600,vals.dropna().max()]
        cmap = (mpl.colors.ListedColormap(['lightgray', 'paleturquoise','deepskyblue', 'royalblue', 'black']))

        # normalize the colors and create mapping color object
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        cm_map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        color_series=vals.map(cm_map.to_rgba)

        # TODO: incorporate into single plot, for now just create a second figure for the scale bar
        fig2, ax2 = plt.subplots(figsize=(6, 1))
        fig2.subplots_adjust(bottom=0.5)
        fig.colorbar( mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
                        cax=ax2,
                        ticks=bounds,
                        orientation='horizontal',
                        label='Max Inundation on Road (mm)')


        ox.plot_graph(G,
                      ax=ax,
                      node_size=0, node_color='black',
                      edge_color=color_series,
                      show=False,
                      edge_linewidth=1.,
                      bbox=total_bounds)

    # plot edges based on width
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
        new_max = 5

        # TODO: could add user inputs to assign values here
        # plotting each scaled width is too much, must bin into groups of 100
        bounds = np.linspace(1,max(unique_flows), num=new_max+1)
        widths = np.linspace(new_min, new_max, num=new_max)

        # previous width determination method => scale each to new min/max
        # for flow in unique_flows:
        #     new_value = (((flow - old_min) * (new_max - new_min)) /
        #                  (old_max - old_min)) + new_min
        
        for idx, width in enumerate(widths):

            # select edges of gdf based on flow value
            selected_edges = G_gdf_edges[(G_gdf_edges[edge_width] >= bounds[idx]) & (G_gdf_edges[edge_width] <= bounds[idx+1])]
            if len(selected_edges) > 0:
                sub = ox.graph_from_gdfs(gdf_edges=selected_edges,
                                        gdf_nodes=G_gdf_nodes)
                
                ox.plot.plot_graph(sub,
                                ax=ax,
                                node_size=0,
                                edge_color='black',
                                show=False,
                                edge_linewidth=width,
                                bbox=total_bounds)

    # optional plot inundation raster
    if inundation is None:
        pass
    else:
        # Choose colormap
        cmap = plt.cm.Blues
        # Get the colormap colors
        my_cmap = cmap(np.arange(cmap.N))
        # Set alpha
        my_cmap[:, -1] = np.linspace(0.5, 1, cmap.N)
        my_cmap[0][3] = 0
        # Create new colormap
        my_cmap = ListedColormap(my_cmap)
        # BUG: inundation plot kills figure sometimes - too big?
        print('plotting raster')
        show(inundation, ax=ax, cmap=my_cmap, zorder=100)
        print('raster plotted')
    
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


def inundate_network(G, path, inundation, G_capacity='capacity',G_width='width',G_speed='speed_kph',G_length='length',pp=4, name='AOI_Graph_Inundated'):
    '''
    Create a new graph network based on the imapct of an inundation layer.

    Currently - An inundation layer is overlayed on a graph network. The maximum depth that intersects each road segment (within a given distance based on the type of road). That maximum depth is then equated to a decrease in speed and capacity on that entire road segment.

    This function adds attributes to the graph network including including:

    - max_inundation_mm, the maximum inundation in milimeters
    - inundation_capacity_{mod}, the reduced capacity of the road segement
    - kph_{mod}, the reduced speed of the road segment
    - inundation_travel_time_{mod}, the new travel time of the road segment

    For inundation_capacity_, kph_, and inundation_travel_time_, there are three seperate occurances of each attribute with the suffix mod, agr, or con refering to the moderate, aggressive, or conservative depth-disruption equation is used in its calculation. All are quadratically decreasing depending on the maximum allowable road depth which are 300, 150, or 600 respectively. Any road segments with a greater depth of water on a road segement that these values has a travel time and capacity set equal to 0.

    This function also relies on parallel processing. Zonal statistics for each road segment is time consuming, but also independent of each other and therefore is easily parallelized. In one test, a single thread took 10 minutes to process 36,000 edges. The same study area took less than 3 minutes to process with 4 threads. 
    
    **Important Note** - There are a number of methodologies in which this can be accomplished, and others need to be explored further. The documenation listed must be updated regularly whenever changes are made to reflect current options of this function. Just a few things that need to be addressed include:

    - impact on partial segments compared to entire segment

    :param G: The graph network. Should have projected coordinates.
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param path: File path to save output graph network to. Will have the name 'AOI_Graph_Inundated'
    :type path: string
    :param inundation: File path to inundation raster .tif
    :type inundation: string
    :param eq: either 'conservative', 'moderate', or 'aggresive', refers to which inundation depth-disruption equation to implement *Default='moderate'*
    :type eq: string
    :param G_capacity: graph attribute refering to road capacities. Capacity has units of vehicles per lane per hour *Default='capacity'* 
    :type G_capacity: string
    :param G_width: graph attribute refering to road width. Width has units of meters *Default='width'* 
    :type G_width: string
    :param G_speed: graph attribute refering to free flow road speed limit. Speed has units of km per hour *Default='speed'* 
    :type G_speed: string
    :param G_length: graph attribute refering to length of road segments. Length has units of meters *Default='length'* 
    :type G_length: string
    :param pp: number of threads to use in parallel processing
    :type pp: int

    :return: **inundated_G**, the impacted graph network. Will have new edge attributes of 'max_inundation', 'inundation_capacity', 'inundation_speed_kph', and 'inundation_travel_time'
    :rtype: networkx.Graph [Multi, MultiDi, Di]
    
    '''
    # convert graph to nodes and edges
    nodes, edges = ox.graph_to_gdfs(G)
    # Set buffer geometry of edges gdf based on road width
    edges['buffer_geometry'] = edges.buffer(distance=edges[G_width], cap_style=2)
    edges.set_geometry(col='buffer_geometry', drop='geometry, inplace=True')
    # edges['buffer_geometry'].to_file(
    #     '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/test_edges.shp')

    # BEGIN PARALLEL PROCESSING
    # geometry column number
    geo_col_num = edges.columns.get_loc("geometry")
    # the individual inundation zonal statistic task
    def inundate_zs(n):
        return zonal_stats(edges.iat[n,geo_col_num], inundation, stats='max')
    # create the sequence of numbers for each edge
    n=range(0,len(edges))
    # create a pool of workers and run the function inundate_zs for each zone
    # map() function creates the background batches
    pool=Pool(pp)
    results=pool.map(inundate_zs, n)
    # close the pool
    pool.close()
    # END PARALLEL PROCESSING

    # convert results to list - if None, where raster doesn't intersect shapefile, replace with 0
    results_as_list = [0 if d[0]['max'] is None else round(d[0]['max']*1000/10)*10 for d in results]
    # relate back to edges geodataframe
    edges['max_inundation_mm'] = results_as_list
    
    # reduction equation options
    # conservative -> 0.0002415w**2 - 0.2898w + 86.94
    # moderate     -> 0.0009w**2 - 0.5529w + 86.9448
    # aggressive   -> 0.003864w**2 - 1.1592w + 86.94

    # CALCULATE THE REDUCTION IN SPEEDS AND SUBSEQUENT TRAVEL TIMES FOR CONSERVATIVE, MODERATE, AND AGGRESSIVE SCENARIOS
    # calcualte percentage of speed (PSR) remaining - replace with 0 if greater than the low points
    # Moderate
    edges.loc[edges['max_inundation_mm']>=300, 'PSR_mod'] = 0
    edges.loc[edges['max_inundation_mm'] < 300, 'PSR_mod'] = (0.0009*edges['max_inundation_mm']**2 - 0.5529*edges['max_inundation_mm'] + 86.9448)/86.9448
    # Conservative
    edges.loc[edges['max_inundation_mm']>=600, 'PSR_con'] = 0
    edges.loc[edges['max_inundation_mm'] < 600, 'PSR_con'] = (0.0002415*edges['max_inundation_mm']**2 - 0.2898*edges['max_inundation_mm'] + 86.94)/86.94
    # Aggressive
    edges.loc[edges['max_inundation_mm']>=150, 'PSR_agr'] = 0
    edges.loc[edges['max_inundation_mm'] < 150, 'PSR_agr'] = (0.003864*edges['max_inundation_mm']**2 - 1.1592*edges['max_inundation_mm'] + 86.94)/86.94

    opts=['mod','con','agr']
    for opt in opts:
        # calculate the inundation speed based on reduction equation
        edges[f'kph_{opt}']=edges[G_speed]*edges[f'PSR_{opt}']
        # calculate inundation travel time in seconds
        edges[f'inundation_travel_time_{opt}'] = edges[G_length]/(edges[f'kph_{opt}']/60/60*1000)
        # replace travel time infs with 0 (when speed reduced to 0)
        edges.loc[edges[f'inundation_travel_time_{opt}'] == np.inf,f'inundation_travel_time_{opt}'] = 0
        # Create inundation capacity. if speed is reduced to zero, setting capacity equal to 0 because it cannot support travel
        edges[f'inundation_capacity_{opt}'] = edges[G_capacity]
        edges.loc[edges[f'kph_{opt}'] == 0, f'inundation_capacity_{opt}'] = 0


    # edges.loc[edges['max_inundation_mm']>=307.167, 'PSR_mod'] = 0
    # edges.loc[edges['max_inundation_mm'] < 307.167, 'PSR_mod'] = (0.0009*edges['max_inundation_mm']**2 - 0.5529*edges['max_inundation_mm'] + 86.9448)/86.9448
    # # Conservative
    # edges.loc[edges['max_inundation_mm']>=600, 'PSR_con'] = 0
    # edges.loc[edges['max_inundation_mm'] < 600, 'PSR_con'] = (0.0002415*edges['max_inundation_mm']**2 - 0.2898*edges['max_inundation_mm'] + 86.94)/86.94
    # # Aggressive
    # edges.loc[edges['max_inundation_mm']>=150, 'PSR_agr'] = 0
    # edges.loc[edges['max_inundation_mm'] < 150, 'PSR_agr'] = (0.003864*edges['max_inundation_mm']**2 - 1.1592*edges['max_inundation_mm'] + 86.94)/86.94

    # # calculate the inundation speed based on reduction equation
    # edges['kph_mod'] = edges[G_speed]*edges['PSR_mod']
    # edges['kph_con'] = edges[G_speed]*edges['PSR_con']
    # edges['kph_agr'] = edges[G_speed]*edges['PSR_agr']
    
    # # calculate inundation travel time in seconds
    # edges['inundation_travel_time_mod'] = edges[G_length]/(edges['kph_mod']/60/60*1000)
    # edges['inundation_travel_time_con'] = edges[G_length]/(edges['kph_con']/60/60*1000)
    # edges['inundation_travel_time_agr'] = edges[G_length]/(edges['kph_agr']/60/60*1000)

    # # replace travel time infs with 0 (when speed reduced to 0)
    # edges.loc[edges['inundation_travel_time_mod'] == np.inf,'inundation_travel_time_mod'] = 0
    # edges.loc[edges['inundation_travel_time_con'] == np.inf,'inundation_travel_time_con'] = 0
    # edges.loc[edges['inundation_travel_time_agr'] == np.inf,'inundation_travel_time_agr'] = 0

    # # Create inundation capacity. if speed is reduced to zero, setting capacity equal to 0 because it cannot support travel
    # edges['inundation_capacity_mod'] = edges[G_capacity]
    # edges['inundation_capacity_con'] = edges[G_capacity]
    # edges['inundation_capacity_agr'] = edges[G_capacity]
    # edges.loc[edges['kph_mod'] == 0, 'inundation_capacity_mod']=0
    # edges.loc[edges['kph_con'] == 0, 'inundation_capacity_con']=0
    # edges.loc[edges['kph_agr'] == 0, 'inundation_capacity_agr'] = 0

    # save edited graph
    new_graph = ox.graph_from_gdfs(nodes, edges)
    save_2_disk(G=new_graph, path=path, name=name)
    return (new_graph)


def flow_decomposition(G='', res_points='', dest_points='',res_parcels='',dest_parcels='', G_demand='demand', G_capacity='capacity', G_weight='travel_time', dest_method='single'):
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
    :param res_parcels: Polygons of all residential parcels
    :type res_parcels: geopandas.GeoDataFrame
    :param dest_parcels: Polygons of all of a specific destination/resource
    :type dest_parcels: geopandas.GeodataFrame
    :param G_demand: name of attribuet in G refering to node demand, *Default='demand*
    :type G_demand: string
    :param G_capacity: name of attribute in G refering to edge capacities, *Default='capacity'* 
    :type G_capacity: string
    :param G_weight: name of attribute in G refering to edge weights, *Default='travel_time'* 
    :type G_weight: string      
    :param dest_method: either 'single' or 'multiple', determines nearest_node methodology to use whether it is the single closest node or the nearest nodes to corners of parcels
    :type dest_method: string

    :returns:
        - **decomposed_paths**, dictionary of dictionaries keyed by unique source nodes
        - **sink_insights**, dictionary of dictionaries containing information about each sink
        - **res_locs**, geopandas.GeoDataFrame, residential parcels with flow information
        - **dest_locs**, geopandas.GeoDataFrame, resource/destination parcels with flow information
    :rtype: tuple



    '''
    # Travel times must be whole numbers -  round values if not whole numbers
    for x in G.edges:
            G.edges[x][G_weight] = round(G.edges[x][G_weight])
    
    G=copy.deepcopy(G)

    if dest_method=='single':
        # Run max_flow_parcels to get flow dictionary
        flow_dict, flow_cost, max_flow, access = max_flow_parcels(G=G, 
                                                                    res_points=res_points, 
                                                                    dest_points=dest_points, 
                                                                    G_capacity=G_capacity, 
                                                                    G_weight=G_weight, 
                                                                    G_demand=G_demand)
        
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

    if dest_method == 'multiple':
        # Run max_flow_parcels to get flow dictionary
        flow_dict, flow_cost, max_flow, access = max_flow_parcels(G=G,
                                                                  res_points=res_points,
                                                                  dest_points=dest_points,
                                                                  G_capacity=G_capacity,
                                                                  G_weight=G_weight,
                                                                  G_demand=G_demand,
                                                                  dest_method='multiple',
                                                                  dest_parcels=dest_parcels)
        
        #Run nearest_nodes to get unique origins and destinations
        G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points = nearest_nodes_vertices(
            G=G, res_points=res_points, dest_parcels=dest_parcels, G_demand=G_demand, dest_points=dest_points)

        # add artifical source node
        G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
        # add edge from artifical source node to real source nodes with 0 weight and capacity equal to demand
        sums = 0
        for unique_node in unique_origin_nodes:
            sums -= G.nodes[unique_node][G_demand]
            # TODO: KWARGS
            kwargs = {f"{G_weight}": 0,
                      f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}
            G.add_edge(0, unique_node, **kwargs)
            # since we added an artificial source node, all original source nodes must have a zero demand
            G.nodes[unique_node][G_demand] = 0

        # add the super_sink that everything goes to
        G.add_nodes_from([(99999999, {G_demand: positive_demand})])

        # identify the destination nodes, and create artifical sink edges
        # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to
        for idx, dest_parcel in dest_parcels.iterrows():
            dest_node = dest_points[dest_points.geometry.within(
                dest_parcel.geometry)]
            #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
            for i, node in dest_node.iterrows():
                # add the dest node to the graph using OSMID as its ID
                # TODO: explore option about plotting with physical lines
                # x = node.geometry.x
                # y = node.geometry.y
                # G.add_nodes_from([(node['osmid'], {'demand': 0, 'x':x, 'y':y})])
                G.add_nodes_from([(node['osmid'], {G_demand: 0})])

                # add links from nearest intersections to parcel centroid
                for nearest_intersection in unique_dest_nodes_list[idx]:
                    # TODO: KWARGS
                    kwargs = {G_weight: 0, G_capacity: 999999999}
                    G.add_edge(nearest_intersection, node['osmid'], **kwargs)

                # add link from parcel centroid to super sink
                # TODO: This is where I can specifically add capacities for each specific grocery store
                # FOR NOW: just setting capacity to a high amount
                kwargs = {G_weight: 0, G_capacity: 999999999}
                G.add_edge(node['osmid'], 99999999, **kwargs)

        
    # Create intermediate graph from max_flow_parcels dictionary, consisting of only edges that have a flow going across them
    # Set the allowable flow of each edge in intermediate graph to the flow going across that edge in original max_flow_parcels solution
    G_inter = nx.DiGraph()
    G=nx.DiGraph(G)
    for i in flow_dict:
        for j in flow_dict[i]:
            if flow_dict[i][j] > 0:
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
            decomposed_paths[source]['Path'].append([' '.join(map(str, l)) for l in [path[1:-1]]])
            decomposed_paths[source]['Flow'].append(limiting_flow)
            decomposed_paths[source]['Cost Per Flow'].append(cost)
            decomposed_paths[source]['Number of Unique Paths'] += 1
            decomposed_paths[source]['Total Flow'] += limiting_flow  
        else:
            decomposed_paths[source] = {'Source': [source], 'Sink': [sink], 'Path': [' '.join(map(str, l)) for l in [path[1:-1]]], 'Flow': [limiting_flow], 'Cost Per Flow': [cost], 'Number of Unique Paths': 1, 'Total Flow': limiting_flow}

        # Add sink insights, keyed by the sink
        if dest_method == 'single': 
            if sink in sink_insights.keys():
                sink_insights[sink]['Number of unique paths'] +=1
                sink_insights[sink]['Total flow w/o walking']+= limiting_flow
                sink_insights[sink]['Total flow w/ walking']+= limiting_flow
            else:
                sink_insights[sink] = {'Number of unique paths': 1, 'Total flow w/o walking': limiting_flow, 'Total flow w/ walking': limiting_flow + len(res_points.loc[res_points['nearest_node'] == sink])}
           
        if dest_method == 'multiple':
            if sink in sink_insights.keys():
                sink_insights[sink]['Number of unique paths'] +=1
                sink_insights[sink]['Total flow w/o walking']+= limiting_flow
                sink_insights[sink]['Total flow w/ walking']+= limiting_flow 
            else:
                near_nodes = list(map(int, dest_points.loc[dest_points['osmid'] == sink, 'nearest_nodes'].iloc[0].split(' ')))
                sink_insights[sink] = {'Number of unique paths': 1, 'Total flow w/o walking': limiting_flow, 'Total flow w/ walking': limiting_flow + len(res_points.loc[res_points['nearest_node'].isin(near_nodes)])}

    # Begin relating decomposition results to outputs 
    # 1. first determine the unique origin nodes that don't have a path and mark accordingly
    # 2. Find all the unique origin nodes with a single path to a single destination, mark accordingly
    # 3. Find all unique origin nodes with mulitple paths to multiple/same desitnation, mark accordingly
    # 4. For shared nodes (nearest node is origin and destination), we consider these the walkable nodes

    # 1. Nodes with no path 
    nodes_to_pop = []
    for node in unique_origin_nodes:
        # Create empty lists of possible sinks, paths, and cost_per_flow
        # extract decomposition information -> if decomp_info doesn't exist, that means unique origin node has NO path to sink
        try:
            decomp_info = decomposed_paths[node]
        except: 
            missing_paths = True
        else:
            missing_paths=False
        # for each path from a unique origin to any sink, need to append lists of possible attributes
        if missing_paths is True:
            # If no decomp_info, source node is not serviceable -> unreachable
            res_points.loc[res_points['nearest_node'] == node,['service']] = 'no'
            res_points.loc[res_points['nearest_node'] == node, ['cost_of_flow']] = np.NaN
            nodes_to_pop.append(node)
            
    # 2. Nodes with a single path
    unique_origin_nodes = unique_origin_nodes.tolist()
    for node in nodes_to_pop:
        unique_origin_nodes.remove(node)

    # create subset of decomposed_paths whose number of unique paths is equal to 1
    decomp_subset = dict((k, decomposed_paths[k]) for k in unique_origin_nodes if decomposed_paths[k]['Number of Unique Paths'] == 1)

    # convert decomp_subset to dataframe
    decomp_subset_df = pd.DataFrame.from_dict(decomp_subset, orient='index',columns=['Source', 
                                                                                        'Sink', 
                                                                                        'Path',
                                                                                        'Flow', 
                                                                                        'Cost Per Flow',
                                                                                        'Number of Unique Paths', 
                                                                                        'Total Flow']).reset_index()
    
    # set the rows the appropriate dtypes before merge
    decomp_subset_df['source'] = decomp_subset_df['Source'].str[0]
    decomp_subset_df['sink'] = decomp_subset_df['Sink'].str[0]
    decomp_subset_df['cost_of_flow'] = decomp_subset_df['Cost Per Flow'].str[0]
    decomp_subset_df['walkable?'] = np.NaN
    decomp_subset_df['service'] = 'yes'
    decomp_subset_df['path'] = decomp_subset_df['Path'].str[0]

    # Merge the appropriate attributes
    res_points = res_points.merge(decomp_subset_df[['source','sink','cost_of_flow','walkable?','service','path']], how='left', left_on='nearest_node', right_on='source')

    # combine similar columns and remove unnecessary ones
    res_points['cost_of_flow'] = res_points['cost_of_flow_y'].combine_first(res_points['cost_of_flow_x'])
    res_points['service'] = res_points['service_y'].combine_first(res_points['service_x'])
    res_points.drop(columns=['cost_of_flow_y','cost_of_flow_x','service_y','service_x'], inplace=True)

    # 3. Nodes with multiple paths
    # create subset of decomposed_paths whose number of unique paths is greater than 1
    decomp_subset = dict((k, decomposed_paths[k]) for k in unique_origin_nodes if decomposed_paths[k]['Number of Unique Paths'] > 1)

    # convert decomp_subset to dataframe
    decomp_subset_df = pd.DataFrame.from_dict(decomp_subset, orient='index', columns=['Source',
                                                                                      'Sink',
                                                                                      'Path',
                                                                                      'Flow',
                                                                                      'Cost Per Flow',
                                                                                      'Number of Unique Paths',
                                                                                      'Total Flow']).reset_index()

    # for each unique origin that has multiple paths, extract each unique sink path and flow
    for ix, row in decomp_subset_df.iterrows():
        for idx, pathway in enumerate(row['Path']):
            flow = int(row['Flow'][idx])
            source = int(row['Source'][0])
            sink = row['Sink'][idx]
            path = pathway
            cost_per_flow = row['Cost Per Flow'][idx]

            # randomly select res_points with same source and no path information
            # the syntax for query is a bit odd: needed @ to reference variable and find nans using !=
            query = res_points.query("nearest_node == @source and sink != sink").sample(n=flow).index
            res_points.loc[query, ['sink']] = sink
            res_points.loc[query, ['path']] = path
            res_points.loc[query, ['cost_of_flow']] = cost_per_flow
            res_points.loc[query, ['walkable?']] = np.NaN
            res_points.loc[query, ['service']] = 'yes'

    #4. set shared nodes to walkable and set remaining features
    res_points.loc[res_points['nearest_node'].isin(shared_nodes), ['walkable?']] = 'yes'
    res_points.loc[res_points['nearest_node'].isin(shared_nodes), ['service']] = 'yes'
    res_points.loc[res_points['nearest_node'].isin(shared_nodes), ['cost_of_flow']] = 0
    res_points.loc[res_points['nearest_node'].isin(shared_nodes), ['sink']] = res_points['nearest_node']
    res_points.loc[res_points['nearest_node'].isin(shared_nodes), ['path']] = 0

    # convert sink_insights to dataframe
    sink_insights_df = pd.DataFrame.from_dict(sink_insights, orient='index', columns=['Number of unique paths',
                                                                                        'Total flow w/o walking',
                                                                                        'Total flow w/ walking']).reset_index()

    # append dest poitns with appropriate flow information
    if dest_method == 'single':
        dest_points = dest_points.merge(sink_insights_df, how='left', left_on='nearest_node', right_on='index')
    if dest_method == 'multiple':
        dest_points = dest_points.merge(sink_insights_df, how='left', left_on='osmid', right_on='index')
       
    # Spatial join the res_parcels/res_points and dest_parcels/dest_points data
    res_parcels = gpd.sjoin(res_parcels,res_points)
    dest_parcels = gpd.sjoin(dest_parcels, dest_points)
  
    # metrics
    # convert to list of flows
    cost_list = sorted(res_parcels['cost_of_flow'].tolist())
    #remove nans 
    cost_list = [x for x in cost_list if math.isnan(x)==False]
    # calculate interquartile range
    q1,q3,=np.percentile(cost_list, [25,75])
    iqr=q3-q1
    upper_bound=q3+(3*iqr)
    # calculate percentage of flows that can be considered a major outlier
    perc_outlier = sum(x>=upper_bound for x in cost_list)/positive_demand *100
    # create list with 'masked' outliers
    cost_list_wo = [upper_bound if x >= upper_bound else x for x in cost_list]
    # median flow cost(excluding including)
    med_cost_in = np.median(cost_list)
    # total flow cost(excluding including)
    total_cost_in = sum(cost_list)
    # average flow cost(excluding including)
    mean_cost_in = np.mean(cost_list)
    # median flow (excluding outliers)
    med_cost_ex = np.median(cost_list_wo)
    # total flow cost(excluding outliers)
    total_cost_ex = sum(cost_list_wo)
    # average flow cost(excluding outliers)
    mean_cost_ex = np.mean(cost_list_wo)

    print(f'percent of flow major outlier: {perc_outlier}')
    print(f'cost of flow with outliers: {total_cost_in}')
    print(f'cost of flow without outliers: {total_cost_ex}')
    print(f'average cost of flow with outliers: {mean_cost_in}')
    print(f'average cost of flow without outliers: {mean_cost_ex}')
    print(f'median cost of flow with outliers: {med_cost_in}')
    print(f'median cost of flow without outliers: {med_cost_ex}')


    # return the decomposed paths dict of dict, the sink_insights dict of dict, and the res_locs/dest_locs geodataframes
    return decomposed_paths, sink_insights, res_parcels, dest_parcels


def shortestPath(origin, capacity_array, weight_array):
    """
    This method finds the shortest path from origin to all other nodes
    in the network.  
    Algorithm utilized is Dijkstra's.s

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
        # Step 5: add i to f_label # TODO: I think this is indented correctly now but should test - was originally operating within last else statement
        f_label.append(node_i)

    backnode=np.asarray(backnode)
    costlabel=np.asarray(costlabel)
    # backnodes equal to list indicating previous node in shortest path
    # costlabel is list of costs associated with each node
    return (backnode, costlabel)


def shortestPath_heap(origin, capacity_array, weight_array, adjlist, destination=False):
    """
    This method finds the shortest path from origin to all other nodes in the network.
    Uses a binary heap priority que in an attempt to speed up the shortestPath function.
    Shortest path algorithm utilized is Dijkstra's.
    going to use built in function in python called heapq

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
    pq = [(0, origin)]
    # Step 2: Initialize SEL with origin, and turn into a heap
    # sel = [origin]

    # Create termination loop, if sel empty, terminate
    while len(pq) > 0:
        # Step 3: Choose a node i from SEL with minimum L
        # doing this by creating a list of cost labels for each node in SEL
        cost, node_i = heapq.heappop(pq)
        if node_i in f_label:
            pass
        else:

            # Step 4: for all arc (i,j) repeat following
            for node_j in adjlist[node_i]:
            # for node_j, matrix_value in enumerate(capacity_array[node_i]):
            #     # if matrix_value is equal to -1, no link exist
            #     if matrix_value == -1:
            #         pass
                #also need to skip nodes in f, because already finalized
                if node_j in f_label:
                    pass
                else:
                    # Step 4a: calculate costlabel temp
                    l_temp = cost + weight_array[node_i][node_j]
                    # Step 4b: if new path is shorter, update cost and backnode
                    if l_temp < costlabel[node_j]:
                        costlabel[node_j] = l_temp
                        backnode[node_j] = node_i
                    # TODO: Step 4c: add j tos heap ->  no way of checking if already in or not, so added if statement at top to filter
                    heapq.heappush(pq, (costlabel[node_j], node_j))
        # Step 5: add i to f_label # TODO: I think this is indented correctly now but should test - was originally operating within last else statement
        f_label.append(node_i)

        # STEP 6: if destination found AND finalized 
        if destination is not False:
            if node_j == destination and node_j in f_label:
                break

    backnode = np.asarray(backnode)
    costlabel = np.asarray(costlabel)
    # backnodes equal to list indicating previous node in shortest path
    # costlabel is list of costs associated with each node
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


def pathTo_mod(backnode, costlabel, origin, destination):
    """test pathTo function to output node-node pairs instead of list of nodes to avoid having to do loops elsewhere
    
    """
    cost = costlabel[destination]
    path = []
    nxt1=[]
    nxt2=[]
    idx = destination
    counter=0
    while idx > 0:
        idx = backnode[idx]
        if counter==0:
            path.append([idx,destination])
            nxt1.append(idx)
        elif counter %2 !=0:
            nxt2.append(idx)
            nxt1.insert(0,idx)
            path.append(nxt1)
            nxt1=[]
        else:
            nxt2.insert(0,idx)
            path.append(nxt2)
            nxt2=[]
            nxt1.append(idx)
        counter+=1
        
    path.append(origin)
    path.reverse()
    path.pop(0)
    path.pop(0)
    #path=np.asarray(path)
    cost=np.asarray(cost)
    return(path, cost)


###### I DON'T THINK I NEED maxFlow OR mod_search ######
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
###### I DON'T THINK I NEED maxFlow OR mod_search ######



# PARALLELIZATION FUNCTIONS USED IN TRAFFIC ASSIGNMENT
def initial_feasible_flow_parallel(n):

    backnode, costlabel = shortestPath_heap(origin=n[0], 
                                                capacity_array=capacity_array, 
                                                weight_array=weight_array, 
                                                adjlist=adjlist, 
                                                destination=n[2])

    path, cost = pathTo_mod(backnode=backnode, 
                            costlabel=costlabel, 
                            origin=n[0], 
                            destination=n[2])
    return path
    

def SPTT_parallel(n):
    backnode, costlabel = shortestPath_heap(origin=n[0], 
                                            capacity_array=capacity_array, 
                                            weight_array=weight_array_iter,
                                            adjlist=adjlist, 
                                            destination=n[2])
    path, cost = pathTo_mod(backnode=backnode, 
                                costlabel=costlabel, 
                                origin=n[0], 
                                destination=n[2])
    return (path, cost)


def traffic_assignment(G,
                        res_points, 
                        dest_points,
                        dest_parcels,  
                        G_capacity='capacity', 
                        G_weight='travel_time', 
                        algorithm='path_based', 
                        method='CFW',
                        link_performance='BPR',
                        termination_criteria=['iter',0],
                        dest_method='single'):
    """
    Solves the static traffic assignment problem based using algorithms and termination criteria from user inputs.

    Taking in a network of roads, shapefile of residential points, and shapefile of destination points, all residential parcels are snapped to the nearest intesection and routed to the nearest destination (also snapped to nearest intersection). Each link has a theoretical capacity (G_capacity argument) and free flow travel time (G_weight argument) associated with it that are used as inputs into the link performance function (link_performance argument). Currently, only one link performance function is available, the user equilibrium BPR function (Bureau of Public Roads function). Other functions will be added in the future.

    Furthermore, the user can define the type of algorithm (algorithm argument) as either link_based, path_based, or bush_based. Currently for link_based arguments only, there are multiple methods that could be used including: MSA (method of successive averages), Bisection, Newton's Method, and CFW (Conjugate Franke-Wolfe). 

    The termination criteria (termination_criteria argument) is a list of length two, with the first argument identifying the type of criteria (e.g., the number of iterations, the relative gap value, etc.) and the second argument being the value to compare the criteria to.
    
    To summarize the algorithms and methods currently avaialble:
            - Link Based
                - Method of Successive Averages (MSA)
                - Bisection
                - Newton's Method
                - Conjugate Frank-Wolfe (CFW)
            
            - Path Based
                - Gradient Projection
            
            -Bush Based 
                - Algorithm B 


    :param G: The graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Location of the residential parcel points
    :type res_points: geopandas Geodataframe
    :param dest_points: Location of the destination parcel points
    :type dest_points: geopanda. Geodataframe
    :param dest_parcels: Technically only needed if dest_method == 'multiple', geodataframe of destination parcels
    :type dest_parcels: gepandas.Geodataframe
    :param G_capacity: Attribute of G that is theoretical road link capacity
    :type G_capacity: string
    :param G_weight: Attribute of G that is free flow travel time across link
    :type G_weight: string
    :param algorithm: Algorithm to be used to solve traffic assignment. Can be link_based, path_based, or bush_based
    :type algorithm: string
    :param method: Method of algorithm to be used. Currenlty only applies to link_based, with options of MSA, Bisection, Newton's Method, or CFW
    :type method: string
    :param link_performance: The link performance function to solve. Currenlty only BPR is available (user equilrium BPR to be specific)
    :type link_performance: string
    :param termination_criteria: The criteria type and comparision number to stop iterations. Available criteria include AEC, iters, and RG
    :type termination_criteria: list, first element string and second element number
    :param dest_method: either 'single' or 'multiple' determines snapping methodology for destination parcels, either single nearest node or multiple nearest nodes
    :type dest_method: string
    :returns:
        - **G_output**, nx.Multigraph of road network with traffic assignment flow attribute, 'TA_Flow'
        - **AEC_list**, list of average access cost, AEC, values for each iteration
        - **TSTT_list**, list of total system travel times, TSTT, values for each iteration
        - **SPTT_list**, list of shortest path travel times, SPTT, values for each iteration
        - **RG_list**, list of relative gap, RG, values for each iteration
        - **iter**, number of iterations completed
    :rtype: tuple

    """

    # Right now, only focused on BPR but could expand to include the following:
    #         - Davidson's delay model
    #         - Akcelik function
    #         - Conical delay model

    # Convert graph to digraph format
    # TODO: Need to see if this is having an impact on the output
    G = nx.DiGraph(G)

    # BUG fix? Remove all edges with 0 capacity
    G=copy.deepcopy(G)
    no_cap_edges = list(filter(lambda e: e[2] == 0, (e for e in G.edges.data(G_capacity))))
    no_cap_ids = list(e[:2] for e in no_cap_edges)
    G.remove_edges_from(no_cap_ids)

    # Travel times must be whole numbers -  round values if not whole numbers
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])

    # Set variables that will be used as constants throughout algorithm 
    super_sink = 99999999
    super_origin = 0
    # Meets int32 byte criteria
    artificial_weight = 999999999
    artificial_capacity = 999999999

    # Set BPR alpha and beta constants
    nx.set_edge_attributes(G, values=0.15, name='alpha')
    nx.set_edge_attributes(G, values=4, name='beta')

    # original node list - used to relate flow back to graph representation
    og_node_list = [num for num in range(0,len(G.nodes))]

    # PREPROCESSING STEPS: 
    # a.   Determine sources/sinks, 
    # b.   set aritificial components of network, 
    # c.   create OD Matrix
    # d.   create list of ordered nodes
    # e.   Convert capacity, weight, alpha, and beta arrays

    if dest_method == 'single':
        # a. Determine sources and sinks
        G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points = nearest_nodes(G=G, 
                                                                                                                            res_points=res_points, 
                                                                                                                            dest_points=dest_points, 
                                                                                                                            G_demand='demand')

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
        #   also rename the nodes in graph
        G = nx.convert_node_labels_to_integers(G, 0, ordering="sorted", label_attribute='old_label')
        nodes_in_order = sorted(G.nodes())

    elif dest_method == 'multiple':
        # a. Determine sources and sinks
        G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points = nearest_nodes_vertices(
            G=G, res_points=res_points, dest_parcels=dest_parcels, dest_points=dest_points, G_demand='demand')

        # b_1. Add an artifical 0 node to maintain other nodes positions
        #   Doesn't actually need attributes since not connected, but filling for continuity sake
        G.add_node(super_origin, **{G_weight: artificial_weight, G_capacity: 0})

        # b_2. Add artifical edges from vertices around destinations to aggregate points than to artifical sink with 0 cost and max capacity
        # This is to allow minimum cost routing to whatever resource - destination of OD pair which can change based on cost
        # identify the destination nodes, and create artifical sink edges
        # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to 
        
        # artificial node id tracker - useful in maintianing compatability with dtype
        dest_node_ids=99999998

        for idx, dest_parcel in dest_parcels.iterrows():
            dest_node = dest_points[dest_points.geometry.within(dest_parcel.geometry)]
            #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
            for i, node in dest_node.iterrows():
                dest_node_ids-=1
                # add the dest node to the graph using OSMID as its ID
                G.add_nodes_from([(dest_node_ids, {'demand': 0})])

                # add links from nearest intersections to parcel centroid
                for nearest_intersection in unique_dest_nodes_list[idx]:
                    kwargs = {G_weight: 0, G_capacity: artificial_capacity}
                    G.add_edge(nearest_intersection, dest_node_ids, **kwargs)

                # add link from parcel centroid to super sink
                # TODO: This is where I can specifically add capacities for each specific grocery store
                # FOR NOW: just setting capacity to a high amount 
                kwargs = {G_weight: 0, G_capacity: artificial_capacity}
                G.add_edge(dest_node_ids, super_sink, **kwargs)

        
        # c. Create OD_matrix and add artifical edge from each origin to super_sink with extreme weight and capacity to capture loss of service 
        #    2 reasons to add origin->super_sink edge: 
        #   (1): Can build in elastic demand cost if people are priced out of accessign resource, or can set to arbitrarily high value to
        #   (2): By having artifical route, can always send max flow w/o considering if accessible - if path goes through artifical link, than no access when cost is arbitrarily high
        # initialize OD_matrix        
        OD_matrix = np.zeros([len(unique_origin_nodes),3])

        for idx,x in enumerate(unique_origin_nodes):
            # BUG: before adding to OD matrix, determine if a path even exists
            if nx.has_path(G, x, super_sink) is False:
                pass
            else:
                OD_matrix[idx][0] = x                            # origin nodes
                OD_matrix[idx][1] = G.nodes[x]['demand']*-1      # flow associated with the origin node
                OD_matrix[idx][2] = len(G.nodes())-1             # destination node (when nodes reordered for numpy array, it has a value of length-1 (b/c of artifical origin 0 ))
                G.add_edge(x, super_sink, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

        # delete OD_matrix rows that we no longer need
        indexList = np.where(~OD_matrix.any(axis=1))[0]
        OD_matrix = np.delete(OD_matrix, indexList, axis=0)


        # d. Sort nodes in proper order
        #   This should in theory then preserve saying node '2' in array is node '2' in network g, especially with adding the artifical node, 0
        #   also rename the nodes in graph
        G = nx.convert_node_labels_to_integers(G, 0, ordering="sorted", label_attribute='old_label')
        nodes_in_order = sorted(G.nodes())

    # set global variables for parallel processing
    global capacity_array
    global weight_array
    global adjlist
    global weight_array_iter
    print('making arrays')
    # convert to capacity and weight arrays
    capacity_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight=G_capacity, nonedge=-1, dtype=np.int32)    # array of theoretical link capacities
    weight_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight=G_weight, nonedge=-1, dtype=np.int32)        # array of free flow weights (costs)
    #alpha_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='alpha', nonedge=-1) 
    #beta_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='beta', nonedge=-1)
    print('arrays made')

    # #TODO - needs to be implmeneted ASAP    
    # # convert arrays to sparse arrays to see if that reduces memory usage
    # weight_array_s = sparse.coo_array(weight_array)
    # capacity_array_s = sparse.csr_array(capacity_array)
    # print(f'weight_array bytes: {weight_array.nbytes}, {weight_array_s.data.nbytes}')
    # print(f'capacity_array bytes: {capacity_array.nbytes}, {capacity_array_s.data.nbytes}')


    # TODO: Check the timing on this -> I believe if I actually used a full array (above) instead of constant (below) the compute time would be significantly longer
    alpha_array = 0.15
    beta_array = 4
    # create adjacency list - faster for some functions to use this instead of matrix
    # adjlist = defaultdict(list)
    # for i in range(capacity_array.shape[0]):
    #     for j in range(capacity_array.shape[0]):
    #         if capacity_array[i][j] != -1:
    #             adjlist[i].append(j)

    # alternative method for creating adjacency list
    adjlist = defaultdict(list)
    row, col = np.where(capacity_array!=-1)
    for i in range(len(row)):
        adjlist[row[i]].append(col[i])
    
    # set variables
    sum_d = positive_demand                             # sum of the demand
    num_nodes = weight_array.shape[0]                   # number of noads
    node_list = [num for num in range(0, num_nodes)]    # list of all nodes
    ##################################################################################################################################

    # WRITE FUNCTIONS FOR CALCULATING BPR FUNCTIONS FOR WEIGHT ARRAY AND WEIGHT ARRAY DERIVATIVE
    # weight_array variable should remain the free flow time array
        # weight_array_itter (or weight_array_iter, need to fix spelling for different algorithms) will be the changing weight array with each iteration

    # LINK PERFORMANCE HELPER FUNCTIONS 
    # So far, only link_performance_function and link_performance derivative require equation inputs 

    def link_performance_function(capacity_array,flow_array,weight_array,eq=link_performance,alpha_array=0.15,beta_array=4):
        """ Function returning the weight_array_iter array based on user selected equation
        
        Universival arguments:
            capacity_array == node-node array of link capacities
            flow_array     == node-node array of link flows
            weight_array   == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
            alpha_array == either node-node array of alpha values for each link OR single value 
            beta_array  == either node-node array of beta values for each link OR single value     

        returns: weight_array_iter
        """
        
        if eq=="BPR":
            weight_array_iter = weight_array*(1+alpha_array*(flow_array/capacity_array)**beta_array)
        # elif eq == "Davidsons":
            # pass
        
        return weight_array_iter


    def link_performance_derivative(capacity_array, flow_array, weight_array, eq=link_performance, alpha_array=0.15, beta_array=4):
        """ Function returning the derivative array based on user selected equation
        
        Universival arguments:
            capacity_array == node-node array of link capacities
            flow_array     == node-node array of link flows
            weight_array   == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
            alpha_array == either node-node array of alpha values for each link OR single value 
            beta_array  == either node-node array of beta values for each link OR single value     

        returns: link_performance_derivative
        """

        if eq == "BPR":
            link_performance_derivative = (beta_array*weight_array*alpha_array/capacity_array**beta_array)*flow_array**(beta_array-1)
        # elif eq == "Davidsons":
            # pass

        return link_performance_derivative


    def bisection_zeta(lambda_val, flow_array_star, capacity_array, flow_array, weight_array, eq=link_performance, alpha_array=0.15, beta_array=4):
        """ Function returning the zeta value when using the bisection link-based method
        
        Universival arguments:
            lambda_val      == lambda value to shift flow
            flow_array_star == updated flow array
            capacity_array  == node-node array of link capacities
            flow_array      == node-node array of link flows
            weight_array    == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
            alpha_array == either node-node array of alpha values for each link OR single value 
            beta_array  == either node-node array of beta values for each link OR single value     

        returns: bisection zeta value
        """
        print('calc x hat')
        x_hat = lambda_val*flow_array_star+(1-lambda_val)*flow_array
        
        print('calc x hat link performance')
        x_hat_link_performance = link_performance_function(capacity_array=capacity_array, 
                                                           flow_array=x_hat,
                                                           weight_array=weight_array,
                                                            eq=eq, 
                                                            alpha_array=alpha_array, 
                                                            beta_array=beta_array)
        print('calc zeta')

        zeta = np.sum(x_hat_link_performance*(flow_array_star-flow_array))

        return(zeta)
        

    def newton_f(lambda_val, flow_array_star, capacity_array, flow_array, weight_array, eq=link_performance, alpha_array=0.15, beta_array=4):
        """ Function returning the f and f_prime values when using the newton's and CFW link-based methods
        
        Universival arguments:
            lambda_val      == lambda value to shift flow
            flow_array_star == updated flow array
            capacity_array  == node-node array of link capacities
            flow_array      == node-node array of link flows
            weight_array    == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
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


    def cfw_alpha(flow_array_star_old, flow_array_star, capacity_array, flow_array, weight_array, eq=link_performance, alpha_array=0.15, beta_array=4):
        """ Function returning the numerator and denominator of alpha value for CFW link-based method

        Universival arguments:
            flow_array_star_old == previous updated flow array
            flow_array_star     == updated flow array (the AON, or All-Or-Nothing assignment)
            capacity_array      == node-node array of link capacities
            flow_array          == node-node array of link flows
            weight_array        == node-node array of free-flow weights (i.e., zero flow travel time)

        BPR arguments:
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


    def topo_order_function(bush, origin, weight_array, num_nodes = num_nodes, node_list=node_list):
        """Function to calculate topographic order on a bush."""

        # Find shortest path from origin to all locations - need to determine what nodes are unreachable 
        # BUG: check if shortestpath_heap actually works here
        #backnode, costlabel = shortestPath(origin, bush, weight_array)
        backnode, costlabel = shortestPath_heap(origin, capacity_array, weight_array, adjlist)
        # where backnode == -1, unreachable, and won't be considered in topological ordering 
        unreachable = [i for i,j in enumerate(backnode) if j == -1]
        unreachable.pop(unreachable.index(origin))
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
                adj_table[node] = [i for i, j in enumerate(bush[node]) if j > 0]
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

        #topo_order.append(origin)
        topo_order = topo_order[::-1]
        return topo_order


    def termination_function(termination_criteria, iters, AEC, RG):
        """Function to determine if algorithm loop should continue or not.
        The input terminatio_criteria is a list with a string and numerical value. Once numerical value is reached, the loop should stop.
        """
        val=True
        if termination_criteria[0] == 'iter':
            if iters >= termination_criteria[1]:
                val=False

        if termination_criteria[0] == 'AEC':
            if AEC <= termination_criteria[1]:
                val=False

        if termination_criteria[0] == 'RG':
            if RG <= termination_criteria[1]:
                val=False

        return val


    # ALGORITHMS
    if algorithm == "link_based":

        # INITIALIZE LINK BASED ALGORITHM
        # a.    Create empty flow array - same size as the capacity_array (b/c node-node relationship)
        # b.    Create initial feasible flow array with all-or-nothing minimum cost flow assignment (free flow travel cost)

        # a. Create empty flow array 
        flow_array = np.zeros_like(capacity_array)

        # BEGIN ORIGINAL NON PARALLEL VERSION 
        # b. Create initial feasible flow array
        # for x in OD_matrix:
        #     origin = np.int(x[0])
        #     flow = x[1]
        #     destination = np.int(x[2])
        #     # calculate shortest path from an origin to all other locations
        #     #backnode, costlabel = shortestPath(origin=origin, capacity_array=capacity_array, weight_array=weight_array)
        #     backnode, costlabel = shortestPath_heap(origin=origin, capacity_array=capacity_array, weight_array=weight_array, adjlist=adjlist, destination=destination)
        #     # determine path and cost from origin to super sink
        #     #path, cost = pathTo(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
        #     path, cost = pathTo_mod(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
        #     # update the flow array by iterating through the shortest paths just determined
        #     for i in path:
        #         flow_array[i[0]][i[1]] += flow
        # END ORIGINAL NON PARALLEL VERSION
  
        # BEGIN PARALLEL PROCESSING
        # set parallel processing variables
        print('starting first parallel')
        pp=8
        n = [[np.int(OD_matrix[n][0]), OD_matrix[n][1], np.int(OD_matrix[n][2])] for n in range(0, len(OD_matrix))]
        # initialize the pool
        pool=Poolm(pp)
        # map function
        results = pool.map(initial_feasible_flow_parallel, n)
        # close the pool
        pool.close()
        pool.join()
        # edit the flow_array matrix
        for idx, result in enumerate(results):
            for i in result:
                flow_array[i[0]][i[1]] += OD_matrix[idx][1]     
        # END PARALLEL PROCESSING 
        print('finished first parallel')
        
        
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
        # set iter values
        iter=0
        iter_val = True
        # begin loop
        while iter_val is True:
            print(f'Begin iteration: {iter}')
            # 1. recalculate weight_array
            weight_array_iter = link_performance_function(capacity_array=capacity_array,
                                                            flow_array=flow_array,
                                                          weight_array=weight_array,
                                                            eq=link_performance,
                                                            alpha_array=alpha_array,
                                                          beta_array=beta_array)
            # 2. Create empty kappa array and flow_array_star
            # TODO: I DON'T THINK KAPPA IS ACTUALLY USED-> I DON'T THINK KAPPA IS NEEDED??????? 
            kappa = np.zeros([1, np.shape(OD_matrix)[0]])
            flow_array_star = np.zeros_like(capacity_array)
            # For shortest path travel time calculation (termination criteria)
            SPTT=0      
            # #BEGIN ORIGINAL CALC
            # for idx, x in enumerate(OD_matrix):
            #     origin = np.int(x[0])
            #     destination = np.int(x[2])
            #     flow = x[1]
            #     # 3. Calculate shortest paths
            #     backnode, costlabel = shortestPath_heap(origin=origin, 
            #                                             capacity_array=capacity_array, 
            #                                             weight_array=weight_array, 
            #                                             adjlist=adjlist, 
            #                                             destination=destination)
            #     path, cost = pathTo_mod(backnode=backnode, 
            #                                 costlabel=costlabel, 
            #                                 origin=origin, 
            #                                 destination=destination)
            #     # 4. Fill shortest path cost array, kappa
            #     kappa[0][idx]=cost      
            #     # For shortest path travel time calculation (termination criteria)
            #     SPTT += cost*flow
                # 5. Update the flow array
                # for i in path:
                #     flow_array_star[i[0]][i[1]] += flow
            # END ORIGINAL CALC
            print(f'start iteration {iter} parallel process')
            # BEGIN PARALLEL CALC
            n = [[np.int(OD_matrix[n][0]), OD_matrix[n][1], np.int(OD_matrix[n][2])] for n in range(0, len(OD_matrix))]
            pool=Poolm(pp)
            results = pool.map(SPTT_parallel, n)
            for idx, result in enumerate(results):
                path, cost = result
                flow = OD_matrix[idx][1]
                SPTT += cost*flow
                for i in path:
                    flow_array_star[i[0]][i[1]] += flow
            pool.close()
            pool.join()
            # END PARALLEL CALC
            print(f'finish iteration {iter} parallel process')
            # Calculate termination variables
            TSTT = np.sum(np.multiply(flow_array, weight_array_iter))
            AEC = (TSTT-SPTT)/sum_d
            RG = TSTT/SPTT - 1
            
            print(f'calc lambda')
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
                    zeta = bisection_zeta(lambda_val = lambda_val, 
                                            flow_array_star=flow_array_star, 
                                             capacity_array=capacity_array, 
                                            flow_array=flow_array,
                                            weight_array=weight_array,
                                            eq=link_performance, 
                                            alpha_array=alpha_array, 
                                            beta_array=beta_array)

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
                # # calculate x_hat
                # x_hat = lambda_val*flow_array_star + (1-lambda_val)*flow_array
                # calculate f and f_prime
                f, f_prime = newton_f(lambda_val=lambda_val, 
                                        flow_array_star=flow_array_star, 
                                        capacity_array=capacity_array, 
                                        flow_array=flow_array, 
                                        weight_array=weight_array, 
                                        eq=link_performance, 
                                        alpha_array=alpha_array, 
                                      beta_array=beta_array)


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
                epsilon=0.01
                if CFW_marker > 0:
                    # solve for alpha, eq (6.3) in boyles textbook, page 176
                    num, den = cfw_alpha(flow_array_star_old=flow_array_star_old, 
                                            flow_array_star=flow_array_star, 
                                            capacity_array=capacity_array, 
                                            flow_array=flow_array, 
                                            weight_array=weight_array, 
                                            eq=link_performance, 
                                            alpha_array=alpha_array, 
                                            beta_array=beta_array)

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
                # calculate f and f_prime
                f, f_prime = newton_f(lambda_val, 
                                        flow_array_star=flow_array_star, 
                                        capacity_array=capacity_array, 
                                        flow_array=flow_array,
                                        weight_array=weight_array,
                                        eq=link_performance, 
                                        alpha_array=alpha_array, 
                                        beta_array=beta_array)

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

            print(f'shift flows')
            # 7. Shift flows
            flow_array = lambda_val*flow_array_star + (1-lambda_val)*flow_array

            # Fill termination variable lists
            AEC_list.append(AEC)
            TSTT_list.append(TSTT)
            SPTT_list.append(SPTT)
            RG_list.append(RG)
            iter += 1

            print(f'termination criteria')
            # determine if iterations should continue
            iter_val = termination_function(termination_criteria = termination_criteria, 
                                            iters = iter, 
                                            AEC = AEC, 
                                            RG = RG)
    
    elif algorithm == 'path_based':
        # Specifically using Newton's Method of Gradient Projection (not to be confused with Projected Gradient)

        # Termination criteria variables I am interested in keeping track off:
        AEC_list = []
        TSTT_list = []
        SPTT_list = []
        RG_list = []

        # Initialize the flow and weight arrays
        flow_array = np.zeros_like(capacity_array)

        # Initialize nested list to store OD path information 
        paths_array = []

        # populate paths_array matrix
        for OD_pair in OD_matrix:
            # The paths_array matrix is set up as follows:
                # paths_array[x][0]: [O, D]
                # paths_array[x][1]: [list of arrays that are paths for OD-pair]
                # paths_array[x][2]: [list of flows (associated with paths)]
            
            origin = OD_pair[0].astype(int)
            flow = OD_pair[1]
            destination = OD_pair[2].astype(int)
            # Calculate shortest path for OD pair
            backnode, costlabel = shortestPath_heap(origin=origin, 
                                                    capacity_array=capacity_array, 
                                                    weight_array=weight_array, 
                                                    adjlist=adjlist)
            path, cost = pathTo_mod(backnode=backnode, 
                                        costlabel=costlabel, 
                                        origin=origin, 
                                        destination=destination)
            
            # update the flow array in order to update the weight array
            for i in path:
                flow_array[i[0]][i[1]]+=flow
            # append the paths_array
            paths_array.append([[origin,destination],[path],[flow]])

        # create and update the weight_array for this itteration
        weight_array_iter = link_performance_function(capacity_array=capacity_array,
                                                        flow_array=flow_array,
                                                        weight_array=weight_array,
                                                        eq=link_performance,
                                                        alpha_array=alpha_array,
                                                        beta_array=beta_array)

        # calculate initial termination criteria variables
        SPTT = 0
        for idx, x in enumerate(OD_matrix):
            origin = np.int(x[0])
            destination = np.int(x[2])
            flow = x[1]
            # Calculate shortest paths
            backnode, costlabel = shortestPath_heap(origin=origin, 
                                                    capacity_array=capacity_array, 
                                                    weight_array=weight_array_iter, 
                                                    adjlist=adjlist)
            # don't actually need shortest path, just need the cost of the path
            cost = costlabel[destination]
            # For shortest path travel time (SPTT) calculation 
            SPTT += cost*flow

        # termination criteria
        TSTT = np.sum(np.multiply(flow_array, weight_array_iter))
        AEC = (TSTT-SPTT)/sum_d
        RG = TSTT/SPTT - 1

        # append the termination criteria lists
        AEC_list.append(AEC)
        TSTT_list.append(TSTT)
        SPTT_list.append(SPTT)
        RG_list.append(RG)

        # set iter values
        iter = 0
        iter_val = True
        # begin loop
        while iter_val is True:

            for OD_pair in paths_array:
                origin = OD_pair[0][0]
                destination = OD_pair[0][1]
                paths = OD_pair[1]
                flows = OD_pair[2]
                # # Find the new shortest path 
                backnode, costlabel = shortestPath_heap(origin=origin, 
                                                        capacity_array=capacity_array, 
                                                        weight_array=weight_array_iter, 
                                                        adjlist=adjlist)
                path_hat, cost_hat = pathTo_mod(backnode=backnode, 
                                                costlabel=costlabel, 
                                                origin=origin, 
                                                destination=destination)
                # if path already in, don't append, just move it, and associated information to the end
                if path_hat in paths:
                    index = paths.index(path_hat)
                    paths.append(paths.pop(index))
                    flows.append(flows.pop(index))

                # if not, append the new path 
                else:
                    paths.append(path_hat)
                    # set to 0 to act as place holder for upcoming calculations
                    flows.append(0)

                # perform functions on each path except the shortest path (in the last index of list)
                for idx, path in enumerate(paths[:-1]):
                    # Determine unique_links between shortest path (path_hat) and current path under investigation
                    all_links = [item for sublst in zip(path, path_hat) for item in sublst]
                    unique_links = unique_links = [list(x) for x in set(tuple(x) for x in all_links)]

                    # calculate numerator of delta_h
                    cost = sum([weight_array_iter[x[0]][x[1]] for x in path])
                    num = cost-cost_hat

                    # determine the sum of derivatives (denominator of delta_h equation)
                    den = 0
                    derivative = link_performance_derivative(capacity_array=capacity_array, 
                                                                flow_array=flow_array, 
                                                                weight_array=weight_array,
                                                                eq=link_performance, 
                                                                alpha_array=alpha_array, 
                                                                beta_array=beta_array)

                    den=sum([derivative[x[0]][x[1]] for x in unique_links])

                    # calculate delta_h
                    delta_h = num/den

                    # adjusted path flow value
                    adj = min(flows[idx], delta_h)

                    # update paths_array flow values
                    flows[idx] -= adj
                    flows[-1] += adj

                    # for the path under examination and path_hat, need to update flow array to calc new weight array
                    for i in path:
                        flow_array[i[0]][i[1]] -= adj
                    for i in path_hat:
                        flow_array[i[0]][i[1]] += adj

                    # Recalculate weight_array for next itteration
                    weight_array_iter = link_performance_function(capacity_array=capacity_array,
                                                                    flow_array=flow_array,
                                                                    weight_array=weight_array,
                                                                    eq=link_performance,
                                                                    alpha_array=alpha_array,
                                                                    beta_array=beta_array)
 
                # TODO: eliminate paths where flow is 0 -> could incorporate this into loop above 
                marker=0
                for flow in flows:
                    if flow == 0:
                        flows.pop(marker)
                        paths.pop(marker)
                    else:
                        marker+=1    

            # calculate  termination criteria variables
            SPTT = 0
            for idx, x in enumerate(OD_matrix):
                origin = np.int(x[0])
                destination = np.int(x[2])
                flow = x[1]
                # # Calculate shortest paths
                backnode, costlabel = shortestPath_heap(origin=origin, 
                                                        capacity_array=capacity_array, 
                                                        weight_array=weight_array_iter, 
                                                        adjlist=adjlist)
                # don't actually need shortest path, just need the cost of the path
                cost = costlabel[destination]
                # For shortest path travel time calculation (termination criteria)
                SPTT += cost*flow

            TSTT = np.sum(np.multiply(flow_array, weight_array_iter))
            AEC = (TSTT-SPTT)/sum_d
            RG = TSTT/SPTT - 1

            # Fill termination variable lists
            AEC_list.append(AEC)
            TSTT_list.append(TSTT)
            SPTT_list.append(SPTT)
            RG_list.append(RG)
            iter += 1

            # determine if iterations should continue
            iter_val = termination_function(termination_criteria=termination_criteria, 
                                            iters = iter, 
                                            AEC=AEC,
                                            RG=RG)

    elif algorithm == 'bush_based':
        
 
        # Convert network into a set of bushes for each OD pair
        #topo_orders = []
        bushes = []
        bush_flows = [] 

        # Termination criteria variables I am interested in keeping track off:
        AEC_list = []
        TSTT_list = []
        SPTT_list = []
        RG_list = []

        for pair in OD_matrix:
            origin = pair[0].astype(int) 
            flow = pair[1]
            destination = pair[2].astype(int)

            # Create an empty  bush that is the same shape as weight_array
            bush = np.full_like(weight_array, -1)
            # Find shortest path from origin to all locations - need to determine what nodes are unreachable 
            backnode, costlabel = shortestPath(origin, capacity_array, weight_array)
       
            # create alternative to prim's just using the backnode array? not minimum spanning tree but shortest path graph?
            #   TODO: Did i fix this??
            for i, j in enumerate(backnode):
                if j ==-1:
                    pass
                else:
                    # I think J and I are placed properly?
                    bush[j][i]=1

            # Set initial flow for each bush - flow along shortest path
            bush_flow = np.zeros_like(weight_array)
            backnode, costlabel = shortestPath(origin, bush, weight_array)
            path, cost = pathTo(backnode=backnode, costlabel=costlabel, origin=origin, destination=destination)
            # path_mod, cost_mod = pathTo_mod(backnode, costlabel, origin, destination)

            for i in range(len(path)-1):
                bush_flow[path[i]][path[i+1]] += flow

            # append the trackers with topo orders and bushes
            bushes.append(bush)
            bush_flows.append(bush_flow)

        # sum flow on bushes for initial network flow
        flow_array = np.zeros_like(weight_array)
        for bush in bush_flows:
            flow_array += bush

        # have to replace the negative flows in flow_array with 0s
        flow_array[flow_array<0]=0

        # # calculate initial network travel times
        weight_array_itter = link_performance_function(capacity_array=capacity_array,
                                                        flow_array=flow_array,
                                                        weight_array=weight_array,
                                                        eq=link_performance,
                                                        alpha_array=alpha_array,
                                                        beta_array=beta_array)

        def label_function(topo_order, bush, weight_array_itter):
            '''function to calculate all of the L and U labels for algorithm B'''
            # Set L and U variables
            L_link = np.full_like(weight_array_itter, np.inf)
            U_link = np.full_like(weight_array_itter, -np.inf)
            L_node = np.full(len(topo_order), np.inf)
            U_node = np.full(len(topo_order), -np.inf)
            # L and U at r (origin) are equal to 0
            L_node[0] = 0
            U_node[0] = 0
            
            # determine L and U labels in forward topological ordering
            id = 0
            while id <= len(topo_order)-1:

                # for first in topological order
                # i == topo_order value & id == index
                i = topo_order[id]
                if id == 0:
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
                            u_val = max(u_val, U_link[h][i])
                    L_node[id] = l_val
                    U_node[id] = u_val

                # for all others
                else:
                    l_val = L_node[id]
                    u_val = U_node[id]
                    for h in topo_order[:id+1]:
                        if bush[h][i] == 1:
                            l_val = min(l_val, L_link[h][i])
                            u_val = max(u_val, U_link[h][i])

                    if l_val == np.inf:
                        for j in topo_order[id:]: 
                            if bush[i][j] == 1:
                                L_link[i][j] = weight_array_itter[i][j]
                                U_link[i][j] = weight_array_itter[i][j]
                        L_node[id] = 0
                        U_node[id] = 0

                    else:
                        L_node[id] = l_val
                        U_node[id] = u_val
                        for j in topo_order[id:]:
                            if bush[i][j] == 1:
                                L_link[i][j] = min(L_link[i][j], L_node[id]+weight_array_itter[i][j])
                                U_link[i][j] = max(U_link[i][j], U_node[id]+weight_array_itter[i][j])

                id += 1

            # remove infs for clarity
            U_link[np.isinf(U_link)] = 0
            L_link[np.isinf(L_link)] = 0

            return L_node, U_node, L_link, U_link


        # Begin the loops of, for each bush,:
            # 1. find if shortcuts exist (i.e., find new shortest path from O-D, if it uses links not on bush, add them)
                # 1a. find shortest path
                # 1b. add shortcuts to bush
            # 2. Calculate L and U labels in forward topological order
            # 3. conduct divergence node determination loop to shift flows using Newton's method
            # 4. After all Newton's adjustments for this bush, calculate new network flow and adjust travel times
            # repeat for each bush

        # set iter values
        iter = 0
        iter_val = True
        # begin loop
        while iter_val is True:

            #print('#######################################################################################')
            for idx, bush in enumerate(bushes):
                print(idx)
                # set variables
                origin = OD_matrix[idx][0].astype(int)
                flow = OD_matrix[idx][1]
                destination = OD_matrix[idx][2].astype(int)
                bush_flow = bush_flows[idx]

                # Calculate initial L and U Labels
                # conduct topological ordering from the origin node
                topo_order = topo_order_function(bush, origin=origin, weight_array=weight_array_itter, num_nodes=num_nodes, node_list=node_list)

                # calculate initial L and U labels
                L_node, U_node, L_link, U_link = label_function(topo_order, 
                                                                    bush, 
                                                                    weight_array_itter)

                # 1. FIND SHORTCUTS AND ADD TO BUSH
                # I don't like adding every possible shortcut - this doesn't seem efficent to me
                # BUG : furthermore, I am not positive but I believe this can inadvertently add a loop, thus invalidating the topo_order function later 
                #       THIS NEEDS TO BE ADDRESSED
                # from textbook and slides, proper way is to add any link where Ui + tij < Uj # BUG L or U??? POWERPOINT AND BOOK SAY OPPOSITE
                for i in topo_order:
                    for j in topo_order:
                        if (bush[i][j] == -1) and (capacity_array[i][j] > 0):      # doesen't exist in bush but does exist in weight arary (i.e., whole network)
                            if L_node[topo_order.index(i)] + weight_array_itter[i][j] <= L_node[topo_order.index(j)]:
                                bush[i][j] = 1
                                
   
                # 2. CALCULATE L AND U LABELS IN FORWARD TOPOLOGICAL ORDER WITH NEW SHORTCUTS
                # conduct topological ordering from the origin node
                topo_order = topo_order_function(bush, origin=origin, weight_array=weight_array_itter, num_nodes=num_nodes, node_list=node_list)
                # calculate new L and U labels
                L_node, U_node, L_link, U_link = label_function(topo_order, bush, weight_array_itter)

                # instead of being interested in last node topoligcally, maybe we are interested in destination node?
                # BUG: this is helping calculate delta_hs but it really is supposed to be last node topologically from notes
                topo_order=topo_order[:topo_order.index(num_nodes-1)+1]

                # DIVERGENCE NODE LOOP
                loop = True
                while loop is True:
                    if len(topo_order) == 1:
                        loop=False
                    else:
                        # determine topoligcally last node in the bush
                        i = topo_order[-1]
                        
                        # determine the longest and shortest paths from U and L Labels
                        # max loop calculation 
                        # because we plan to subtract flows from max path, if any flow on max path is equal to 0, we can skip this iteration 
                        max_val_flow = np.inf
                        max_val = -np.inf
                        max_path = [i]
                        for h in topo_order[:topo_order.index(i)+1]:
                            if bush[h][i] == 1:
                                if max_val < U_link[h][i]:
                                    max_val = U_link[h][i]
                                    next_node_max = h
                        max_path.append(next_node_max)
            
                        # loop for max path
                        while max_path[-1] != origin:
                            max_val = -np.inf
                            for h in topo_order[:topo_order.index(next_node_max)+1]:
                                if bush[h][next_node_max].astype(int) == 1: 
                                    if max_val < U_link[h][next_node_max]:
                                        max_val = U_link[h][next_node_max]
                                        max_val_flow = min(max_val_flow, bush_flow[h][next_node_max])
                                        next_node_holder = h
                            max_path.append(next_node_holder)
                            next_node_max = next_node_holder

                        # check if max path is feasible, otherwise restart search
                        if max_val_flow == 0:
                            topo_order = topo_order[:-1]
                        else:
                            
                            # min loop calculation         
                            min_val = np.inf
                            min_path = [i]
                            for h in topo_order:
                                if bush[h][i] == 1:
                                    if min_val > L_link[h][i]:
                                        min_val = L_link[h][i]
                                        next_node_min = h
                            min_path.append(next_node_min)
                            # loop for min path
                            while min_path[-1] != origin:
                                min_val = np.inf
                                for h in topo_order[:topo_order.index(next_node_min)+1]:
                                    if bush[h][next_node_min].astype(int) == 1:
                                        if min_val > L_link[h][next_node_min]:
                                            min_val = L_link[h][next_node_min]
                                            next_node_holder = h
                                min_path.append(next_node_holder)
                                next_node_min = next_node_holder
                        
                            backnode, costlabel = shortestPath(origin, bush, weight_array_itter)
                            path_mod, cost_mod = pathTo_mod(backnode, costlabel, origin, destination)

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
                                num1 = U_node[topo_order.index(i)] - U_node[topo_order.index(a)]
                                num2 = L_node[topo_order.index(i)] - L_node[topo_order.index(a)]
                                numerator = num1-num2
                                #numerator = (U_node[topo_order.index(i)] - U_node[topo_order.index(a)]) - (L_node[topo_order.index(i)] - L_node[topo_order.index(a)])

                                # BPR function derivative
                                derivative_array = link_performance_derivative(capacity_array=capacity_array, 
                                                                                flow_array=flow_array, 
                                                                                weight_array=weight_array, 
                                                                                eq=link_performance, 
                                                                                alpha_array=alpha_array, 
                                                                                beta_array=beta_array)

                                denominator = 0
                                for id in range(len(sigma_L)-1):
                                    denominator += derivative_array[sigma_L[id]][sigma_L[id+1]]
                                for id in range(len(sigma_U)-1):
                                    denominator += derivative_array[sigma_U[id]][sigma_U[id+1]]

                                # set delta_h value
                                delta_h = numerator/denominator

                                # determine delta_h based on maximum flow that can be shifted 
                                for id in range(len(sigma_U)-1):
                                    if delta_h <= bush_flow[sigma_U[id]][sigma_U[id+1]]:
                                        pass
                                    else:
                                        delta_h = bush_flow[sigma_U[id]][sigma_U[id+1]]
                                
                                # Shift flows
                                for id in range(len(sigma_L)-1):
                                    flow_array[sigma_L[id]][sigma_L[id+1]] += delta_h
                                    bush_flow[sigma_L[id]][sigma_L[id+1]] += delta_h
                                for id in range(len(sigma_U)-1):
                                    flow_array[sigma_U[id]][sigma_U[id+1]] -= delta_h
                                    bush_flow[sigma_U[id]][sigma_U[id+1]] -= delta_h

                                # at end of loop, slice off last value and repeat
                                topo_order = topo_order[:-1]
                    
                # update the weight array
                weight_array_itter = link_performance_function(capacity_array,flow_array,weight_array,eq=link_performance,alpha_array=alpha_array,beta_array=beta_array)
               
                # TODO Won't work right now but need to figure out how to remove unused eges
                # # Remove edges from bush that have 0 flow while preserving conectivity 
                # backnode, costlabel = shortestPath(origin, bush, weight_array)
                # reachable = np.where(backnode != -1)
                # for i in topo_order:
                #     for j in topo_order[topo_order.index(i):]:
                #         if (bush[i][j] == 1) and (bush_flow[i][j] == 0):
                #             # remove the edge
                #             bush[i][j] = -1
                #             backnode_n, costlabel_n = shortestPath(origin, bush, weight_array)
                #             reachable_n = np.where(backnode_n != -1)
                #             print(f'reachable reachable_n {reachable, reachable_n, [i], [j]}')
                #             if np.array_equal(reachable_n, reachable):
                #                 print(f'removed {[i], [j], idx}')
                #                 pass
                #             else:
                #                 bush[i][j] = 1
        
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

            # Fill termination variable lists
            AEC_list.append(AEC)
            TSTT_list.append(TSTT)
            SPTT_list.append(SPTT)
            RG_list.append(RG)
            iter += 1

            # determine if iterations should continue
            iter_val = termination_function(termination_criteria=termination_criteria,
                                            iters=iter,
                                            AEC=AEC,
                                            RG=RG)


    # EDIT THE FINAL OUTPUTS
    nx.set_edge_attributes(G, nx.get_edge_attributes(G, 'travel_time'), name='Weight_Array_Iter')
    # relate flow array back to network
    new_data= {}
    for i in og_node_list:
        for j in og_node_list:
            if flow_array[i][j] > 0:
                new_data.update({(i, j): {'TA_Flow': flow_array[i][j], 'Weight_Array_Iter': weight_array_iter[i][j]}})

    # set edge attributes with new data
    nx.set_edge_attributes(G, new_data)

    # remove artificial super source and sink
    G.remove_node(len(G.nodes())-1)
    G.remove_node(super_origin)

    # if multiple dest method, remove the intermediate artifical links
    if dest_method == 'multiple':
        while len(list(G.nodes)) > len(og_node_list):
            G.remove_node(len(G.nodes()))
    #TODO:  has to be multidigraph to save properly should adjust save function to represent this requirement
    G_output = nx.MultiDiGraph(G)

    # #TODO: Fix outputs - either add save option or remove and keep as seperate function - in scratch. This function should return the graph
    save_2_disk(G=G_output, path='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs',
                name='AN_Graph_Traffic_Assignment_No_Inundation')

    # delete global variables 
    del capacity_array
    del weight_array
    del adjlist
    del weight_array_iter

    return G_output, AEC_list, TSTT_list, SPTT_list, RG_list, iter
