# network exploration stuff
# using a number of websites for background information
# https://automating-gis-processes.github.io/2017/lessons/L7/network-analysis.html
# https://networkx.org/documentation/stable/tutorial.html
# https://geoffboeing.com/2016/11/osmnx-python-street-networks/

# Packages actively using
import sys
import copy
import math
import heapq
from heapq import heappop, heappush
from itertools import count
from itertools import islice
from collections import defaultdict
import rasterio
import numpy as np
import osmnx as ox
import pandas as pd
import rasterio.mask
import networkx as nx
import geopandas as gpd
import matplotlib as mpl
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDirectionArrows
import matplotlib.patches as mpatches
from typing import Union
from scipy import sparse
from scipy.stats import iqr
from rasterio.plot import show
import matplotlib.pyplot as plt
from collections import Counter
from rasterstats import zonal_stats
from collections import defaultdict
from multiprocessing import Pool as Poolm
from matplotlib.colors import ListedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from multiprocessing.pool import ThreadPool as Pool

from networkx.algorithms.connectivity import build_auxiliary_edge_connectivity
from networkx.algorithms.flow import build_residual_network



# testing sknetwork for faster shortest path calcualtions
# import sknetwork as skn

from concurrent.futures import ThreadPoolExecutor

# packages I use for testing 
import time
# use following for testing time:
    # start=time.time()
    # ...some code...
    # end=time.time()
    # print(end-start)

# packages I may be able to remove
import fiona
import random
import shapely
from scipy import sparse
import matplotlib.pylab as pl
from matplotlib import pyplot
import matplotlib.colors as colors
from multiprocessing import Process

# set some temporary panda display options
pd.set_option('display.max_columns', None)

############# NOTES SECTION: contains organized table of todos and bugs that I am actively and planning on working on
    # Bugs should be seen as major issues that need to be addressed ASAP, TODOS are smaller improvements to improve functionality

# BUGS
# BUG1: In traffic assignment function, convert numpy arrays to sparse arrays to save memory and actually be able to calculate 


# MAJOR TODOS: I think the most important items that need to be addressed ASAP
# TODO13: Multiple resource assignment versus individual resources at once
# TODO15: Convert entire code base into OOP - convert all the individual functions to a class

# CAPACITY TODOS
# TODO1: Change disruption impact on capacity: fully disrupted nodes/edges need to just be completly removed (i.e., BPR can't divide by zero)
# TODO2: Add resource capacity. e.g., a grocery store can only support so many people. 
    # can make this a hard capacity (no one else gets served)
    # or a soft capacity, where an additional cost is incured, which grows exponentially with each new person

# SPEED DISRUPTION TODOS
# 

# POTENTIAL LIMITATION TODOS
# TODO3: Some of the closest origin nodes are the destination nodes - for now assuming that no paths need to be calculated (i.e., can just walk)
    # Could probably set a buffer, e.g., origins with XX meters are within walking distance and can be ignored from travelling network
        # this will be related to resource capacity, see #TODO2


# IMPROVEMENT TODOS
# TODO4: move digraph conversion into better spot and determine when it is needed
    # should be within each function NOT in the save or load function, because we should keep the original multidigraph saved
# TODO5: Move CRS check into functions
# TODO6: Min cost flow and max cost flow functions have some repeat code - I could condense a lot of the work
# TODO7: Use kwargs methodology to create new edges with appropriate attributes
#   see # TODO7: KWARGS below in code for example
# TODO8: read_graph_from_disk - problem with preserving dtypes so I hardcoded a
    # temporary fix - right now it goes through the attributes and determines if it could be converted to a float - if it can it does. 
    # this involves manually skipping oneway as well (i think b/c its bool)
    # this was needed b/c  inundation_G load not preserve datatypes
    # I think this can just stay how it is and isn't really an improvement that needs to be made
# TODO9: update shortest path and path to functions to the new modified ones - much quicker
# TODO10: add impact of rainfall intensity on reduction of speed and intensity
    # light rain reduces free flow speed 2-13% and capacity 4-11%
    # heavy rain reduces speed 3-16% and capacity 10-30%
    # see these sites for info: 
    # https://engrxiv.org/preprint/view/1322/2765
    # https://ops.fhwa.dot.gov/weather/best_practices/AMS2003_TrafficFlow.pdf 
    # https://ops.fhwa.dot.gov/weather/q1_roadimpact.htm#:~:text=Light%20rain%20can%20decrease%20freeway,12%20percent%20in%20low%20visibility.
# TODO11: modallity component: i.e., factoring in public transportation and access to vehicle - stochastic component?
# TODO12: add a node snapping filter to correct when parcels are snapped to the wrong nearest intersection
    # a few different options exist:
        # choose to ignore roads that are a certain classification (i.e., freeway)
        # see if parcel has information that aligns with OSM data (i.e., match street names)
        # examine neighbors nearest nodes (i.e., nearest parcels within radius, where do they snap? if pattern exists change snapping for incorrect one)
# TODO14: Identify and double check bridges that are/are not flooding
# TODO15: In some functions, we add a demand attribute to the network. The network isn't an output so demand is not saved. May need to consider changing this
# TODO16: In many calculations, artifical node is set as node 99999999 - may want to reconsider this to avoid the small chance of duplicating a node
# TODO17: Could explore option of graphing artifical nodes/edges
        # x = node.geometry.x
        # y = node.geometry.y
# TODO18: Add option to change color of specific nodes in plot_aoi
    # example: G_gdf_nodes.loc[[6016], 'geometry'].plot(ax=ax, color='green')


# RESOURCES TO EXAMINE
    #https://osmnx.readthedocs.io/en/stable/osmnx.html?highlight=interpolate#osmnx.utils_geo.interpolate_points 
    #https://stackoverflow.com/questions/64104884/osmnx-project-point-to-street-segments 
    #https://github.com/gboeing/osmnx/issues/269 


# ASSUMPTIONS THAT NEED TO BE DOCUMENTED: Tried to capture all assumptions, regardless of how small, for transparency
# 1. Origin node snapping - access to transportation network means nearest intersection isn't flooded
# 2. Destination node snapping - if you can access the nearest corner of a destination, you can access that resource
# 3. OpenStreetMap edges with multiples type attributes (e.g., ['residential','unclassified']) are set to the first in the list
# 4. OpenStreetMap edges with multiple lanes listed (e.g., ['3', '4', '2']) are set to the smallest number in the list
# 5. All of the predetermined road attributes (e.g., number of lanes, speeds, width of lane, capacities, etc.)
# 6. When multidigraph is converted to digraph, if parallel edges exist, only one set of attributes remains the same. Little to no impact on results so far.




# THE FUNCTIONS #############################################################

def crs_check(source: str, shapefiles: list[str]=[], crs: int=None):
    '''
    Resave shapefiles with a matching coordinate reference system

    This function rewrites a list of of shapefiles to have the same coordinate reference system as the shapefile given in the source argument. Alternatively, if a *crs*, is passed, the source and the list of shapefiles will all be written with the given coordinate reference system. To pass a crs, use the EPSG integer convention. For example, 'crs=4326' is WGS84 coordinate reference system. For some OSMNx operations to function properly, WGS84 is required, and this is taken care of in the functions themselves.

    :param source: Path to .shp file of 
    :type source: str
    :param source: Paths to .shp files that need their coordinate reference system checked
    :type source: list of strings
    :param source: if not None, rewrite all inputs with given crs
    :type source: int or None
    
    '''

    if crs is not None:
        crs_proj = f"EPSG:{crs}"
        gpd.read_file(source).to_crs(crs_proj).to_file(source)
    else:
        crs_proj = gpd.read_file(source).crs
    for shapefile in shapefiles:
        gpd.read_file(shapefile).to_crs(crs_proj).to_file(shapefile)


def shape_2_graph(source: str):
    '''
    Extracts OSM road network from shapefile.

    This function takes an input string path to a .shp file of the area of interest (assuming the associated .shp, .shx, .dbf, etc. are all in the same folder) and returns the OpenStreetMap road network for that area.  It is important to consider an area of interest with an appropriate buffer to capture a large enough portion of the road network. For example, if the area of interest is a Census block group, residents will likely require driving outside of that boundary to get to another area within the bounadry. This function automatically cleans and sets the lane capacities (passenger car units per hour), number of lanes, lane width (standard 3.7 meters), speed (kilometers per hour), and travel time (seconds) information, all necessary components in further calculations. It should be noted that output lane capacities have already been converted to passenger car units per hour per road segment, meaning that the number of lanes is accounted for in the capacity attribute. For road segments with missing attributes in OpenStreetMap, values are set to standard values from appropriate highway manuals based on road type. 
    
    It is recomended that the input shapefile have a projected coordinate reference system. Regardless, the coordiante reference system of the input shapefile will be preserved in the output road network. If this function is being run multiple times and the user needs a new and updated network, delete the cache or the same network will be used every time.

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
    

    # Need to check/add bridge data. When inundating network, bridge nodes are ignored and assumed to be above flood waters
    for k, v in G.edges.items():
        if 'bridge' not in G[k[0]][k[1]][k[2]]:
            G[k[0]][k[1]][k[2]]['bridge'] = 'No'

    # project back to original crs
    G=ox.project_graph(G, crs_proj)


    return(G)


def save_2_disk(G: nx.MultiDiGraph, path: str, name='AOI_Graph'):
    '''
    Save network as graphxml file.

    Maintaining an original copy of the road network used during an analysis is beneficial because edits to the network (e,g., deletion of inundated edges ir addition of artificial edges) can inadvertnely impact or prevent other analyses or plotting. Users can continue to reload and access the original network under consideration. Typically the input graph, G, will be the graph obtained from :func:`shape_2_graph`. This function can only save a network if it is a MultiDiGraph.
    
    :param G: Input graph to be saved
    :type G: networkx.MultiDiGraph 
    :param path: Path to folder where graph will be saved
    :type path: str
    :param name: Name the graph will be saved as. *Default='AOI_Graph*
    :type name: str

    '''
    ox.io.save_graphml(G, f'{path}/{name}')


def rename(G: nx.Graph,
           start=1):
    '''
    Rename nodes from 1 to N, the total number of nodes.

    Renaming nodes in this order makes it easier to access nodes/edges of interest. It is also necessary when comparing nodes/edges before/after :func:`traffic_assignment` which converts networks to/from arrays.

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param start: the number to start the numbering
    :type start: int
    :return: **G**, renamed graph network
    :rtype: networkx.Graph [Multi, MultiDi, Di] 
    '''

    # Number of nodes
    num_nodes = len(G)
    # list from start to [num_nodes]
    lst = list(range(start, num_nodes+1))
    # Create a dictionary (value for each key is None) from sorted G nodes
    mapping = dict.fromkeys(sorted(G))
    # iterate through this dictionary (with enumerate to have count var.)
    for idx, key in enumerate(mapping):
        # set each key value (previously None) to new value
        mapping[key] = lst[idx]
    # nx.relabel_nodes relabels the nodes in G based on mapping dict
    G = nx.relabel_nodes(G, mapping)
    return (G)


def read_graph_from_disk(path: str, name='AOI_Graph', rename_nodes=True):
    '''
    Open a saved graphxml file.

    This works in conjunction with :func:`save_2_disk`, opening a saved graphxml file.

    :param path: Path to folder where graphxml file is saved.
    :type path: str
    :param name: Name of the graphxml file. *Default='AOI_Graph'*
    :type name: str
    :param rename_nodes: Whether or not nodes should be renamed from 1 to N, where N is the total number of nodes.
    :type rename_nodes: bool. *Default=True*
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
    
    if rename_nodes is True:
        # rename nodes from 1 to N, the total number of nodes
        G = rename(G)
    else:
        pass
    return(G)


def parallel_edges(G: nx.Graph):
    '''
    Find parallel edges, and self loop edges, in a graph network.

    A lot of network algorithms in NetworkX require that networks be DiGraphs, or directional graphs. When extracting road networks with OSMnx, they are often returned as MultiDiGraphs, meaning they could potentially have parallel edges (multiple edges with the same start and finish node) or self loop edges (an edge with the same start and end node). This function determines if these instances exists. There could potentially be an issue in end results depending on the attributes of of the parallel/loop edge and each should be examined in a study area if they exist. Typically, the do not have an impact.

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :returns: 
        - **parallel_edges**, list of parallel edges, *[[u,v],...]*
        - **self_loop_edges**, list of self loop edges, *[[u,u],...]*
    :rtype: tuple

    '''
   
    parallel_edges=[]
    self_loop_edges=[]
    for u,v,k in G.edges(keys=True):
        if k>0:
            parallel_edges.append([u,v])
        if u==v:
            self_loop_edges.append([u, v])

    return parallel_edges, self_loop_edges


def nearest_nodes(G: nx.Graph, res_points: gpd.GeoDataFrame, dest_points: gpd.GeoDataFrame, G_demand='demand'):
    '''
    Returns graph with demand attribute based on snapping sources/sinks to each's nearest node (intersection).
    
    The purpose of this function is take a series of input locations (residential points and destination points) and determine each's nearest node (i.e., roadway intersections) within a Graph network. It takes this information to append demand (i.e., source and sink) information. 

    **Important Note**: The output graph has an attribute labeled based on the *G_demand argument*, which shows the number of closest residential parcels to each unique intersection. Because these are considered sources, they have a negative demand. Sink locations, or the intersections that are closest to the the dest_points, will not have a demand value (G_demand == 0) because we do not know how much flow is going to that node until after a minimum-cost flow algorithm is run. In the subsequent flow routing algorithms, we create an artificial sink and then decompose the flow to determine how much flow is going to each unique destination point. The graph that is returned is maintains the type of graph that was used as an input.
 
    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the residential parcels are to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param G_demand: name of attribuet in G refering to node demand, *Default='demand*
    :type G_demand: string
    :returns: 
        - **G**, *networkx.Graph [Multi, MultiDi, Di]* with G_demand attribute showing source values (negative demand)
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

    # Create the demand attribute and set all equal to 0
    nx.set_node_attributes(G, values=0, name=G_demand)
    
    # Create list of unique origins 
    unique_origin_nodes = np.unique(origins)

    # seperate instances of when a single geodataframe are passed or a list are
    if isinstance(dest_points, gpd.GeoDataFrame):
        # Count the number of unique origin nodes
        # create a dictionary of origin nodes and demand (i.e., number of origins at that node for *-1)
        unique_origin_nodes_counts = dict(zip(list(Counter(origins).keys()), [x*-1 for x in Counter(origins).values()]))
        # will find the nearest node based solely on the centroid of resource parcel (i.e., dest_points)
        longs=dest_points.geometry.x
        lats=dest_points.geometry.y
        
        destinations = ox.distance.nearest_nodes(G, longs, lats)
        dest_points['nearest_node'] = destinations
        # create list of unique destinations
        unique_dest_nodes = np.unique(destinations)
        
        # determine shared nodes and set demand equal to 0
        shared_nodes = list(set(unique_origin_nodes) & set(unique_dest_nodes))
        for x in shared_nodes:
            unique_origin_nodes_counts[x] = 0

        # add source information (the negative demand value prev. calculated)
        nx.set_node_attributes(G, unique_origin_nodes_counts, G_demand)

        # Calculate the positive demand: the sink demand
        positive_demand = sum(unique_origin_nodes_counts.values())*-1
        return(G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points)
    
    elif isinstance(dest_points, list):
        # Count the number of unique origin nodes
        # create a dictionary of origin nodes and demand (i.e., number of origins at that node for *-1)
        unique_origin_nodes_counts = dict(zip(list(Counter(origins).keys()), [x*-1*len(dest_points) for x in Counter(origins).values()]))
        unique_origin_nodes_counts_copy = dict(zip(list(Counter(origins).keys()), [x*-1*len(dest_points) for x in Counter(origins).values()]))
        # create variables to fill
        unique_dest_nodes_final = []
        shared_nodes_final = []
        # iterate through shapefiles, appending the appropriate outputs
        for idx, shp in enumerate(dest_points):
            unique_origin_nodes_counts_iter = dict(zip(list(Counter(origins).keys()), [x*-1 for x in Counter(origins).values()]))
            # Create the singular demand attribute and set all equal to 0
            nx.set_node_attributes(G, values=0, name=f'{G_demand}{idx}')
            # will find the nearest node based solely on the centroid of resource parcel (i.e., dest_points)
            longs = shp.geometry.x
            lats = shp.geometry.y
            destinations = ox.distance.nearest_nodes(G, longs, lats)
            shp['nearest_node'] = destinations

            # create list of unique dest nodes
            unique_dest_nodes = np.unique(destinations)
            unique_dest_nodes_final.append(list(unique_dest_nodes))
            
            # determine shared nodes and reduce demand for that origin node by the number of origins that share with that specific resource
            shared_nodes = list(set(unique_origin_nodes) & set(unique_dest_nodes))
            shared_nodes_final.append(shared_nodes)
            for x in shared_nodes:
                unique_origin_nodes_counts[x] += unique_origin_nodes_counts_copy[x]*-1/len(dest_points)
                # convert to make sure int
                unique_origin_nodes_counts[x] = int(unique_origin_nodes_counts[x])
                # adjust and append the iter counts for individual demand values
                unique_origin_nodes_counts_iter[x]=0
            nx.set_node_attributes(G, unique_origin_nodes_counts_iter, f'{G_demand}{idx}')

            
        # Calculate the positive demand: the sink demand
        positive_demand = sum(unique_origin_nodes_counts.values())*-1
        # add source information (the negative demand value prev. calculated)
        nx.set_node_attributes(G, unique_origin_nodes_counts, G_demand)

        return (G, unique_origin_nodes, unique_dest_nodes_final, positive_demand, shared_nodes_final, res_points, dest_points)
    
    
def nearest_nodes_vertices(G: nx.Graph, 
                           res_points: gpd.GeoDataFrame, 
                           dest_parcels: gpd.GeoDataFrame, 
                           dest_points: gpd.GeoDataFrame, 
                           G_demand='demand'):
    '''
    Returns graph with demand attribute based on snapping sources/sinks to nearest nodes (intersections).
    
    This is a modification of the :func:`nearest_nodes` function, where instead of snapping destination parcels to their nearest singular intersection, it finds the nearest intersection to each corner of the parcel boundary. This more accuretly reflects multiple entrences to destinations/resources.
    
    **Important Note**: The output graph has an attribute labeled based on the *G_demand argument*, which shows the number of closest residential parcels to each unique intersection. Because these are considered sources, they have a negative demand. Sink locations, or the intersections that are closest to the the dest_points, will not have a demand value (G_demand == 0) because we do not know how much flow is going to each location until after a minimum-cost flow algorithm is run. In the subsequent flow routing algorithms, we create an artificial sink and then decompose the flow to determine how much flow is going to each unique destination point. The graph that is returned is maintains the type of graph that was used as an input.

    Destination points (i.e., resource centroids) is passed with the *dest_points* argument to output information regarding the nearest intersections being used for each destination point. This parameter is also used in determining if multiple destination files are being pased

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_parcels: Parcels of all of the resources the residential parcels are to be routed to
    :type dest_parcels: geopandas.GeoDataFrame
    :param G_demand: name of attribute in G refering to node demand, *Default='demand*
    :type G_demand: string
    :param dest_points: Will output with added columns on the nearest nodes
    :type dest_points: geopandas.GeoDataFrame
    :returns: 
        - **G**, *networkx.DiGraph* with G_demand attribute showing source values (negative demand)
        - **unique_origin_nodes**, *lst* of unique origin nodes
        - **unique_dest_nodes**, *lst* of unqiue destination nodes
        - **positive_demand**, *int* of the total positive demand across the network. **Does not include the residents that share a closest intersection with a resource parcel.**
        - **Shared_nodes**, *lst* of nodes that res and destination parcels share same nearest 
        - **res_points**, *geopandas.GeoDataFrame*, residential points with appended attribute of 'nearest_node'
        - **dest_parcels**, *geopandas.GeoDataFrame*, destination/sink footprints with appended attribute of 'nearest_node'
        - **dest_points**, *geopandas.GeoDataFrame*, *Optional,* destination/sink points with appended attribute of 'nearest_node'
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
    # Create the demand attribute and set all equal to 0
    nx.set_node_attributes(G, values=0, name=G_demand)

    # Create list of unique origins
    unique_origin_nodes = np.unique(origins)

    # seperate instance of when a single geodataframe are passed or a list are
    if isinstance(dest_points, gpd.GeoDataFrame):
        # create a dictionary of origin nodes and demand (i.e., number of origins at that node for *-1 demand)
        unique_origin_nodes_counts = dict(zip(list(Counter(origins).keys()), [x*-1 for x in Counter(origins).values()]))

        # Will determine vertices of dest parcels and find multiple nearest nodes for destinations
        # simplify dest_parcels to reduce geometry
        # Preserve topology of the parcels
        dest_simp = dest_parcels.simplify(3, preserve_topology=True)
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

        # Create list of unique destinations
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
        
        return(G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points)

    elif isinstance(dest_points, list):
        # create a dictionary of origin nodes and demand (i.e., number of origins at that node for *-1 demand)
        unique_origin_nodes_counts = dict(zip(list(Counter(origins).keys()), [x*-1*len(dest_points) for x in Counter(origins).values()]))
        unique_origin_nodes_counts_copy = dict(zip(list(Counter(origins).keys()), [x*-1*len(dest_points) for x in Counter(origins).values()]))
        # create variables to fill
        unique_dest_nodes_list_final=[]
        shared_nodes_final=[]
        
        # iterate through shapefiles, appending the appropriate outputs
        for step, shp in enumerate(dest_parcels):
            unique_origin_nodes_counts_iter = dict(zip(list(Counter(origins).keys()), [x*-1 for x in Counter(origins).values()]))
            # Will determine vertices of dest parcels and find multiple nearest nodes for destinations
            # simplify dest_parcels to reduce geometry
            # Preserve topology of the parcels
            dest_simp = shp.simplify(3, preserve_topology=True)
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
            shp['nearest_nodes'] = destinations_str
            
            dest_points[step] = gpd.sjoin(dest_points[step], shp)
            dest_points[step].drop(columns=['index_right'], inplace=True)
            # create single list of destinations
            destinations = [element for innerList in destinations_int for element in innerList]

            # Create list of unique destinations
            unique_dest_nodes = np.unique(destinations)
            unique_dest_nodes_list = [np.unique(nodes) for nodes in destinations_int]
            unique_dest_nodes_list_final.append(unique_dest_nodes_list)

            # determine shared nodes and set demand equal to 0
            shared_nodes = list(set(unique_origin_nodes) & set(unique_dest_nodes))
            shared_nodes_final.append(shared_nodes)
            for x in shared_nodes:
                unique_origin_nodes_counts[x] += unique_origin_nodes_counts_copy[x]*-1/len(dest_parcels)
                # convert to make sure int
                unique_origin_nodes_counts[x] = int(unique_origin_nodes_counts[x])
                # adjust and append the iter c ounts for individual demand values
                unique_origin_nodes_counts_iter[x]=0
            nx.set_node_attributes(G, unique_origin_nodes_counts_iter, f'{G_demand}{step}')

        # add source information (the negative demand value prev. calculated)
        nx.set_node_attributes(G, unique_origin_nodes_counts, G_demand)

        # Calculate the positive demand: the sink demand
        positive_demand = sum(unique_origin_nodes_counts.values())*-1

        return(G, unique_origin_nodes, unique_dest_nodes_list_final, positive_demand, shared_nodes_final, res_points, dest_parcels, dest_points)
    

def random_shortest_path(G: nx.Graph, res_points: gpd.GeoDataFrame, dest_points: gpd.GeoDataFrame, plot=False):
    '''
    Shortest path between a random residential parcel and all of a given resource.

    This function takes a set of residential points, chooses one at random, and then routes it to every destination point (i.e., some given resource). This function is capable of producing a figure that also shows how that residential parcel would be routed to each of the destination points. 

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


def min_cost_flow_parcels(G: nx.Graph, 
                          res_points: gpd.GeoDataFrame, 
                          dest_points: gpd.GeoDataFrame, 
                          dest_parcels: gpd.GeoDataFrame=None, 
                          G_demand='demand', 
                          G_capacity='capacity', 
                          G_weight='travel_time', 
                          dest_method='single') -> tuple[dict, int]:
    '''
    Find the minimum cost flow of all residential parcels to a given resource.

    Given a set of sources (residential parcels) and sinks (destination parcels) this function determines the minimum cost flow while accounting for edge capacities. Sources send demand and therefore have a negative demand value. Sinks recieve demand and therefore have a positive demand value. This notation comes from from the NetworkX package. All sinks are connected to an artificial sink in order to route sources to their nearest available destination.

    Residential and resource parcels are snapped to the road network using :func:`nearest_nodes` and :func:`nearest_nodes_vertices`. If the *dest_method* argument is set to 'single', destination parcels are snapped to their nearest nodes. If the *dest_method* argument is set to 'multiple', than destination parcels are snapped to the nearest nodes for each vertex of the parcel boundary. 

    Similar to other functions, source demand is the sum of the nearest residential parcels that share that node as their nearest and if a residential parcel and shares a nearest node with a destination parcel this demand is ignored and they are assumed to be within walking distance.

    If this function returns an error, than every source is not able to reach the sink. We recommend using :func:`max_flow_parcels` in order to examine the network.

    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the residential parcels will be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param dest_parcels: technically only required if *dest_method* == 'multiple', destination parcel shapefile
    :param dest_points: geopandas.GeoDataFrame
    :param G_demand: name of attribute in G refering to node demands, *Default='demand'* 
    :type G_demand: string
    :param G_capacity: name of attribute in G refering to edge capacities, *Default='capacity'* 
    :type G_capacity: string
    :param G_weight: name of attribute in G refering to edge weights, *Default='travel_time'* 
    :type G_weight: string  
    :param dest_method: either 'single', or 'multiple', determines which destination point methodology to use. Either connecting destinations/sinks to nearest single node or multiple nearest nodes
    :type dest_method: string
    :returns: 
        - **flow_dictionary**, dictionary of dictionaries keyed by nodes for edge flows
        - **cost_of_flow**, integer, total cost of all flow. If weight/cost is travel time, this is the total time for everyone to reach the resource (seconds).
    :rtype: tuple
    :Raises: nx.NetworkXUnfeasible if all demand cannot be satisfied, i.e., all residential parcels cannot reach the resource(s).

    ''' 
    # For algorithm to function properly, need to remove parallel and self loops by converting to digraph
    
    G = nx.DiGraph(G)
    # Travel times must be whole numbers -  round values if not whole numbers
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])
    
    # if only going to single destination file, proceed normally
    if isinstance(dest_points, gpd.GeoDataFrame):
        
    # if snapping destination parcels to nearest singular node only
        if dest_method == 'single':
            # Find the source and sink demands and append to graph G
            G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points = nearest_nodes(
                G=G, res_points=res_points, dest_points=dest_points, G_demand=G_demand)
            
            # create artificial sink node with calculated total demand
            #    All sinks will go to this demand in order to balance equation
            G.add_nodes_from([(99999999, {G_demand: positive_demand})])
            # add edges from sinks to this artificial node
            #TODO2: this is where I could add destination capacities if I wanted to
            for x in unique_dest_nodes:
                kwargs = {f"{G_weight}": 0}
                G.add_edge(x, 99999999, **kwargs)

        # if snapping destination parcels to nearest multiple nodes, based on parcel vertices
        elif dest_method == 'multiple':
            
            # Find the source and sink demands and append to graph G
            G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points = nearest_nodes_vertices(G=G, res_points=res_points, dest_parcels=dest_parcels, dest_points=dest_points, G_demand=G_demand)
            # add artifical source node
            G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
            # add edge from artifical source node to real source nodes with 0 weight and capacity equal to demand
            sums=0
            for unique_node in unique_origin_nodes:
                sums -= G.nodes[unique_node][G_demand]
                kwargs = {f"{G_weight}": 0, f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}  
                G.add_edge(0, unique_node, **kwargs)
                # since we added an artificial source node, all original source nodes must have a zero demand
                G.nodes[unique_node][G_demand]=0

            # add the super_sink that everything goes to 
            G.add_nodes_from([(99999999, {G_demand: positive_demand})])

            # artificial node id tracker - useful in maintianing compatability with dtype
            dest_node_ids = 99999998
            # identify the destination nodes, and create artifical sink edges
            # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to 
            for idx, dest_parcel in dest_parcels.iterrows():
                dest_node = dest_points[dest_points.geometry.within(dest_parcel.geometry)]
                #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
                for i, node in dest_node.iterrows():
                    dest_node_ids -= 1
                    # add the dest node to the graph using OSMID as its ID
                    G.add_nodes_from([(dest_node_ids, {'demand': 0})])

                    # add links from nearest intersections to parcel centroid
                    for nearest_intersection in unique_dest_nodes_list[idx]:
                        kwargs = {G_weight: 0, G_capacity: 999999999}
                        G.add_edge(nearest_intersection, dest_node_ids, **kwargs)

                    # add link from parcel centroid to super sink
                    # TODO2: This is where I can specifically add capacities for each specific grocery store
                    # FOR NOW: just setting capacity to a high amount 
                    kwargs = {G_weight: 0, G_capacity: 999999999}
                    G.add_edge(dest_node_ids, 99999999, **kwargs)

        # run the min_cost_flow function to retrieve FlowDict
        try:
            flow_dictionary = nx.min_cost_flow(G, demand=G_demand, weight=G_weight, capacity=G_capacity)
            # run the cost_of_flow function to retrieve total cost
            cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
            return flow_dictionary, cost_of_flow
        except nx.NetworkXUnfeasible:
            return None, None

    # if multiple destinations, simply nx.run min_cost_flow and nx.cost_of_flow multiple times and combine outputs
    # congestion isn't considered, so don't need to "run them all at the same time"
    elif isinstance(dest_points, list) or isinstance(dest_parcels, list):

        dictionaries=[]
        final_cost_of_flow = 0

        for idx, dest_point in enumerate(dest_points):

            if dest_method == 'single':
                # Find the source and sink demands and append to graph G
                G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_point_out = nearest_nodes(
                    G=G, res_points=res_points, dest_points=dest_point, G_demand=G_demand)
                # Create artificial sink node with calculated total demand
                # All sinks will go to this demand in order to balance equation
                G.add_nodes_from([(99999999, {G_demand: positive_demand})])
                # add edges from sinks to this artificial node
                for x in unique_dest_nodes:
                    kwargs = {f"{G_weight}": 0}
                    G.add_edge(x, 99999999, **kwargs)

            elif dest_method == 'multiple':
                # keep track of edges added to make it easier to remove them each iteration
                edges_added=[]
                # Find the source and sink demands and append to graph G
                G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcel, dest_point_out = nearest_nodes_vertices(
                    G=G, res_points=res_points, dest_parcels=dest_parcels[idx],dest_points=dest_point, G_demand=G_demand)
                # add artifical source node
                G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
                # add edge from artifical source node to real source nodes with 0 weight and capacity equal to demand
                sums = 0
                for unique_node in unique_origin_nodes:
                    sums -= G.nodes[unique_node][G_demand]
                    kwargs = {f"{G_weight}": 0,
                            f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}
                    G.add_edge(0, unique_node, **kwargs)
                    edges_added.append((0,unique_node))
                    # since we added an artificial source node, all original source nodes must have a zero demand
                    G.nodes[unique_node][G_demand] = 0

                # add the super_sink that everything goes to
                G.add_nodes_from([(99999999, {G_demand: positive_demand})])

                # artificial node id tracker - useful in maintianing compatability with dtype
                dest_node_ids = 99999998
                # identify the destination nodes, and create artifical sink edges
                # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to
                for idx, parcel in dest_parcel.iterrows():
                    dest_node = dest_point[dest_point.geometry.within(parcel.geometry)]
                    # since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
                    for i, node in dest_node.iterrows():
                        dest_node_ids -= 1
                        # add the dest node to the graph using OSMID as its ID
                        G.add_nodes_from([(dest_node_ids, {'demand': 0})])

                        # add links from nearest intersections to parcel centroid
                        for nearest_intersection in unique_dest_nodes_list[idx]:
                            kwargs = {G_weight: 0, G_capacity: 999999999}
                            G.add_edge(nearest_intersection, dest_node_ids, **kwargs)
                            edges_added.append((nearest_intersection, dest_node_ids))

                        # add link from parcel centroid to super sink
                        # TODO2: This is where I can specifically add capacities for each specific grocery store
                        # FOR NOW: just setting capacity to a high amount
                        kwargs = {G_weight: 0, G_capacity: 999999999}
                        G.add_edge(dest_node_ids, 99999999, **kwargs)
                        edges_added.append((dest_node_ids, 99999999))

            # run the min_cost_flow function to retrieve FlowDict
            try:
                flow_dictionary = nx.min_cost_flow(G, demand=G_demand, weight=G_weight, capacity=G_capacity)
                dictionaries.append(flow_dictionary)
                # run the cost_of_flow function to retrieve total cost
                cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
                final_cost_of_flow += cost_of_flow

                # before next iteration, remove the edges that were added 
                # TODO: don't know if i should ignore capacity (almost always feasible) or include capacity constraint (order of dest_points matters)
                # TODO: FOR NOW, ignoring compounding capacity on each iteration
                if dest_method == 'single':
                    for x in unique_dest_nodes:
                        G.remove_edge(x, 99999999)
                    G.remove_node(99999999)
                if dest_method == 'multiple':
                    G.remove_edges_from(edges_added)
                    G.remove_node(99999999)

            except nx.NetworkXUnfeasible:
                print('Error in calculating one of the min cost flow paths')

        # sum to get final flow dictionary
        final_flow_dictionary = {}
        for data in dictionaries:
            for k in data.keys():
                    for x, y in data[k].items():
                        if not k in final_flow_dictionary.keys():
                            final_flow_dictionary[k] = {x: y}
                        else:
                            if not x in final_flow_dictionary[k].keys():
                                final_flow_dictionary[k].update({x: y})
                            else:
                                final_flow_dictionary[k].update({x: final_flow_dictionary[k][x] + y})

        return final_flow_dictionary, final_cost_of_flow


def max_flow_parcels(G: nx.Graph, 
                     res_points: gpd.GeoDataFrame, 
                     dest_points: gpd.GeoDataFrame, 
                     dest_parcels: gpd.GeoDataFrame=None, 
                     G_capacity='capacity', 
                     G_weight='travel_time', 
                     G_demand='demand', 
                     dest_method='single', 
                     ignore_capacity=False) -> tuple[dict, int, int, str]:
    '''
    Find maximum flow with minimum cost of a network between residential parcels and a given resource.
    
    This function calculates the maximum feasible flow from all sources (residential parcels) and sinks (destination parcels). This function is independent of whether or not every residential parcel can actually reach a resource. Residential and resource parcels are snapped to the road network using :func:`nearest_nodes` and :func:`nearest_nodes_vertices`. If the *dest_method* argument is set to 'single', destination parcels are snapped to their nearest nodes. If the *dest_method* argument is set to 'multiple', than destination parcels are snapped to the nearest nodes for each vertex of the parcel boundary. Similar to other functions, source demand is the sum of the nearest residential parcels that share that node as their nearest and if a residential parcel and shares a nearest node with a destination parcel this demand is ignored and they are assumed to be within walking distance.


    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_points: Point locations of all of the residential parcels
    :type res_points: geopandas.GeoDataFrame
    :param dest_points: Point locations of all of the resources the residential parcels are to be routed to
    :type dest_points: geopandas.GeoDataFrame
    :param dest_parcels: echnically only required if *dest_method* == 'multiple', destination parcel shapefile
    :type dest_parcels: geopandas.GeoDataFrame
    :param G_capacity: name of attribute in G refering to edge capacities, *Default='capacity'* 
    :type G_capacity: string
    :param G_weight: name of attribute in G refering to edge weights, *Default='travel_time'* 
    :type G_weight: string  
    :param dest_method: either 'multiple' or 'single', determines snapping mechanism for destination points, if it uses multiple vertices or just the nearest vertice
    :type dest_method: string
    :param ignore_capacity: if True, road link capacities will be set to infinity. If False, original road capacities are used
    :type ignore_capacity: bool
    :returns: 
        - **flow_dictionary**, dictionary of dictionaries keyed by nodes for edge flows
        - **cost_of_flow**, integer, total cost of all flow. If weight/cost is travel time, this is the total time for everyone to reach the resource (seconds).
        - **max_flow**, integer, the maximum flow (i.e., households) that can access the resource
        - **access**, string, either 'Complete' or 'Partial', representing if every residential parcel can/cannot access the resource
    :rtype: tuple

    '''

    # for many functions to work, graph needs to be a digraph (NOT a multidigraph) i.e., no parallel edges
    G = nx.DiGraph(G)

    # Based on input, change capacities 
    if ignore_capacity is True:
        nx.set_edge_attributes(G, 999999999999, name=G_capacity)

    # Travel times must be whole numbers - just round values
    for x in G.edges:
        G.edges[x][G_weight] = round(G.edges[x][G_weight])

    # if only going to single destination file, proceed normally
    if isinstance(dest_points, gpd.GeoDataFrame):
        # if snapping destination parcels to nearest singular node only
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
                kwargs = {f"{G_weight}": 0, f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}  
                G.add_edge(0, unique_node, **kwargs)
                # since we added an artificial source node, all original source nodes must have a zero demand
                G.nodes[unique_node][G_demand]=0
                
            # create artificial sink node 
            G.add_nodes_from([(99999999, {G_demand: positive_demand})])
            # add edges from sinks to this artificial node
            for x in unique_dest_nodes:
                kwargs={f"{G_weight}":0}
                G.add_edge(x, 99999999, **kwargs)
        # if snapping destionation parcels to nearest mulitple nodes, based on parcel vertices
        elif dest_method == 'multiple':
            # Find the source and sink demands and append to graph G
            G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points = nearest_nodes_vertices(G=G, res_points=res_points, dest_parcels=dest_parcels, dest_points=dest_points, G_demand=G_demand)
            # add artifical source node
            G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
            # add edge from artifical source node to real source nodes with 0 weight and capacity equal to demand
            sums=0
            for unique_node in unique_origin_nodes:
                sums -= G.nodes[unique_node][G_demand]
                kwargs = {f"{G_weight}": 0, f"{G_capacity}": G.nodes[unique_node][G_demand]*-1} 
                G.add_edge(0, unique_node, **kwargs)
                # since we added an artificial source node, all original source nodes must have a zero demand
                G.nodes[unique_node][G_demand]=0

            # add the super_sink that everything goes to 
            G.add_nodes_from([(99999999, {G_demand: positive_demand})])


            # artificial node id tracker - useful in maintianing compatability with dtype
            dest_node_ids = 99999998

            # identify the destination nodes, and create artifical sink edges
            # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to 
            for idx, dest_parcel in dest_parcels.iterrows():
                dest_node = dest_points[dest_points.geometry.within(dest_parcel.geometry)]
                #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
                for i, node in dest_node.iterrows():
                    dest_node_ids -= 1
                    # add the dest node to the graph using OSMID as its ID
                    G.add_nodes_from([(dest_node_ids, {G_demand: 0})])

                    # add links from nearest intersections to parcel centroid
                    for nearest_intersection in unique_dest_nodes_list[idx]:
                        kwargs = {G_weight: 0, G_capacity: 999999999} 
                        G.add_edge(nearest_intersection, dest_node_ids, **kwargs)

                    # add link from parcel centroid to super sink
                    # TODO2: This is where I can specifically add capacities for each specific grocery store
                    # FOR NOW: just setting capacity to a high amount 
                    kwargs = {G_weight: 0, G_capacity: 999999999}
                    G.add_edge(dest_node_ids, 99999999, **kwargs)

        # run the max_flow_min_cost function to retrieve FlowDict
        flow_dictionary = nx.max_flow_min_cost(G=G, s=0, t=99999999, weight=G_weight, capacity=G_capacity)
        # run the cost_of_flow function to retrieve total cost
        cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
        max_flow, max_dict = nx.maximum_flow(G,_s=0,_t=99999999,capacity=G_capacity)
                                    
        if max_flow == sums:
            access = 'Complete'
        else:
            access = 'Partial'
        return flow_dictionary, cost_of_flow, max_flow, access
    
    # if multiple dsetinations, simply run nx functions multiple times and combine outputs
    # congestion isn't considered, so no need to "run them all at the same time"
    if isinstance(dest_points, list) or isinstance(dest_parcels, list):

        dictionaries=[]
        final_cost_of_flow=0
        final_max_flow=0
        accesses=[]

        for idx, dest_point in enumerate(dest_points):

            if dest_method == 'single':
                # Find the source and sink demands and append to graph G
                G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_point_out = nearest_nodes(
                    G=G, res_points=res_points, dest_points=dest_point, G_demand=G_demand)

                # add artifical source node
                G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
                # add edge from artifical source node to real source nodes with 0 weight and capacity equal to demand
                sums = 0
                for unique_node in unique_origin_nodes:
                    sums -= G.nodes[unique_node][G_demand]
                    kwargs = {f"{G_weight}": 0,
                            f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}
                    G.add_edge(0, unique_node, **kwargs)
                    # since we added an artificial source node, all original source nodes must have a zero demand
                    G.nodes[unique_node][G_demand] = 0

                # create artificial sink node
                G.add_nodes_from([(99999999, {G_demand: positive_demand})])
                # add edges from sinks to this artificial node
                for x in unique_dest_nodes:
                    kwargs = {f"{G_weight}": 0}
                    G.add_edge(x, 99999999, **kwargs)

            elif dest_method == 'multiple':
                # keep track of edges added to make it easier to remove them each iteration
                edges_added = []
                # Find the source and sink demands and append to graph G
                G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcel, dest_points_out = nearest_nodes_vertices(
                    G=G, res_points=res_points, dest_parcels=dest_parcels[idx], dest_points=dest_point, G_demand=G_demand)
                # add artifical source node
                G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
                # add edge from artifical source node to real source nodes with 0 weight and capacity equal to demand
                sums = 0
                for unique_node in unique_origin_nodes:
                    sums -= G.nodes[unique_node][G_demand]
                    kwargs = {f"{G_weight}": 0,
                            f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}
                    G.add_edge(0, unique_node, **kwargs)
                    edges_added.append((0, unique_node))
                    # since we added an artificial source node, all original source nodes must have a zero demand
                    G.nodes[unique_node][G_demand] = 0

                # add the super_sink that everything goes to
                G.add_nodes_from([(99999999, {G_demand: positive_demand})])
                # artificial node id tracker - useful in maintianing compatability with dtype
                dest_node_ids = 99999998
                # identify the destination nodes, and create artifical sink edges
                # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to
                for idx, dest_parcel in dest_parcel.iterrows():
                    dest_node = dest_point[dest_point.geometry.within(dest_parcel.geometry)]
                    # since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
                    for i, node in dest_node.iterrows():
                        dest_node_ids -= 1
                        # add the dest node to the graph using OSMID as its ID
                        G.add_nodes_from([(dest_node_ids, {G_demand: 0})])

                        # add links from nearest intersections to parcel centroid
                        for nearest_intersection in unique_dest_nodes_list[idx]:
                            kwargs = {G_weight: 0, G_capacity: 999999999}
                            G.add_edge(nearest_intersection,dest_node_ids, **kwargs)
                            edges_added.append((nearest_intersection,dest_node_ids))

                        # add link from parcel centroid to super sink
                        # TODO2: This is where I can specifically add capacities for each specific grocery store
                        # FOR NOW: just setting capacity to a high amount
                        kwargs = {G_weight: 0, G_capacity: 999999999}
                        G.add_edge(dest_node_ids, 99999999, **kwargs)
                        edges_added.append((dest_node_ids, 99999999))

            # run the max_flow_min_cost function to retrieve FlowDict
            flow_dictionary = nx.max_flow_min_cost(G=G, s=0, t=99999999, weight=G_weight, capacity=G_capacity)
            # run the cost_of_flow function to retrieve total cost
            cost_of_flow = nx.cost_of_flow(G, flow_dictionary, weight=G_weight)
            max_flow, max_dict = nx.maximum_flow(G,_s=0,_t=99999999,capacity=G_capacity)
                                        
            if max_flow == sums:
                access = 'Complete'
            else:
                access = 'Partial'
            # append the outputs
            dictionaries.append(flow_dictionary)
            final_cost_of_flow += cost_of_flow
            final_max_flow += max_flow
            accesses.append(access)

            # before next iteration, remove the edges that were added 
            # TODO: don't know if i should ignore capacity (almost always feasible) or include capacity constraint (order of dest_points matters)
            # TODO: FOR NOW, ignoring compounding capacity on each iteration
            if dest_method == 'single':
                for x in unique_dest_nodes:
                    G.remove_edge(x, 99999999)
                G.remove_node(99999999)
            if dest_method == 'multiple':
                G.remove_edges_from(edges_added)
                G.remove_node(99999999)


        # sum to et final flow dictionary
        final_flow_dictionary = {}
        for data in dictionaries:
            for k in data.keys():
                    for x, y in data[k].items():
                        if not k in final_flow_dictionary.keys():
                            final_flow_dictionary[k] = {x: y}
                        else:
                            if not x in final_flow_dictionary[k].keys():
                                final_flow_dictionary[k].update({x: y})
                            else:
                                final_flow_dictionary[k].update({x: final_flow_dictionary[k][x] + y})

        return final_flow_dictionary, final_cost_of_flow, final_max_flow, accesses


def plot_aoi(G: nx.Graph,
            res_parcels: gpd.GeoDataFrame, 
            resource_parcels: gpd.GeoDataFrame, 
            background_edges: gpd.GeoDataFrame=None,
            edge_width: str = None, 
            edge_color: str = None,
            bbox: gpd.GeoDataFrame = None, 
            loss_access_parcels: gpd.GeoDataFrame = None, 
            scalebar: bool = False,
            inundation: rasterio.io.DatasetReader = None, 
            insets: list[str]=None, 
            save_loc: str = None,
            decomp_flow: bool = False,
            rotation=False,
            bg_water=None,
            default=True):
    '''
    Create a plot with commonly used features.

    This function has a variety of parameters to create quick figures that I commonly look at. This can be actively changed to add more/less features.
    **Returns final figure, use plt.show() to see plot.**


    :param G: Graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param res_parcels: The residential parcels to be plotted. Default color is tan.
    :type res_parcels: geopandas.GeoDataframe
    :param resource_parcels: The resource parcels to be plotted. Default color is blue.
    :type resource_parcels: geopandas.GeoDataFrame
    :param background_edges: *Optional*, background edges of plot if missing edges in input graph (*Default=None*)
    :type background_edges: geopandas.GeoDataFrame
    :param edge_width: *Optional*, attribute of *G* to scale road edges by (*Default=None*)
    :type edge_width: string
    :param edge_color: attribute in G to color edges. Originally designed to show max inundation on the road, other attributes can be used though (*Default=None*)
    :type edge_color: string
    :param bbox: set boundary of figures. If *None*, bbox set to max boundary of *G* (*Default=None*) 
    :type bbox: geopandas.GeoDataframe
    :param lose_access_parcels: The residential parcels that can no longer access a resource (i.e., due to flooding) (*Default=None*)
    :type lose_access_parcels: geopandas.GeoDataframe
    :param scalebar: If *True*, plots a scalebar (*Default=False*)
    :type scalebar: bool
    :param inundation: Raster of inundation to plot over figure. Use rasterio.open
    :type inundadtion: rasterio Object
    :param insets: File locations of areas to create inset box outlines (*Default=None*) 
    :type insets: list of strings
    :param save_loc: Location to save figure (*Default=None*)
    :type save_loc: string
    :param decomp_flow: if True, plot residential parcels with color ramp symbolizing travel time (*Default=False*)
    :type decomp_flow: bool
    :param default: if true, creates fig and ax object itself. if not, input [fig, ax object]
    type default: bool or list of [fig, ax] object

    :return: **fig**, the produced figure
    :rtype: matplotlib figure

    '''
    parcel_linewidth=0.5
    linewidth=0.5
    linewidth_max=2


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
        center_point=((total_bounds[2]-total_bounds[3])/2 + total_bounds[3],
                    (total_bounds[0]-total_bounds[1])/2 + total_bounds[1])
    else:
        # set bbox to localized area
        total_bounds = bbox.total_bounds
        # convert bounds to N, S, E, W
        total_bounds = (total_bounds[3],
                        total_bounds[1],
                        total_bounds[2],
                        total_bounds[0])
        center_point=((total_bounds[2]-total_bounds[3])/2 + total_bounds[3],
                      (total_bounds[0]-total_bounds[1])/2 + total_bounds[1])

    # calculate new bounds if rotated
    if rotation is not None:
        if bbox is None:
            newbbox=gpd.GeoSeries(g for g in G_gdf_edges['geometry'])
        else:
            newbbox = gpd.GeoSeries(g for g in bbox['geometry'])
        G_gdf_edges_rot = newbbox.rotate(rotation, origin=center_point)
        total_bounds = G_gdf_edges_rot.total_bounds
        total_bounds = (total_bounds[3],
                        total_bounds[1],
                        total_bounds[2],
                        total_bounds[0])

    if default is True:
        # create subplot with background color and remove axis
        fig, ax = plt.subplots(facecolor='white', figsize=(12,12))
        # fig.set_facecolor('none')
        # Hide X and Y axes label marks
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        # Hide X and Y axes tick marks
        ax.set_xticks([])
        ax.set_yticks([])
        plt.axis('off')
    
    else:
        fig=default[0]
        ax=default[1]

    # set facecolor
    ax.set_facecolor('whitesmoke')

    # set bounding limits 
    ax.set_xlim((total_bounds[3], total_bounds[2]))
    ax.set_ylim((total_bounds[1], total_bounds[0]))
    
    # plot background water
    if bg_water is None:
        pass
    else:
        if rotation is False:
            bg_water.plot(ax=ax, color='midnightblue')
        else:
            bg_water_rot = gpd.GeoSeries(g for g in bg_water['geometry'])
            bg_water_rot = bg_water_rot.rotate(rotation, origin=center_point)
            bg_water_rot.plot(ax=ax, color='midnightblue')

    # plot roads, residential parcels, and resource parcels
    # if decomp_flow is True, plot res parcels with scale bar to show travel times
    if decomp_flow is True:
        if rotation is False:        
            res_parcels.plot(ax=ax, color='tan')
        else:
            res_parcels_rot = gpd.GeoSeries(g for g in res_parcels['geometry'])
            res_parcels_rot = res_parcels_rot.rotate(rotation, origin=center_point)
            res_parcels_rot.plot(ax=ax, color='tan') 
        # # IQR METHOD to mask impact of outliers
        costs = sorted(res_parcels['cost_of_flow'].tolist())
        costs = [x for x in costs if math.isnan(x)==False]
        q1,q3,=np.percentile(costs, [25,75])
        iqr=q3-q1
        upper_bound=q3+(3*iqr)
        res_parcels.loc[res_parcels['cost_of_flow']>=upper_bound, ['cost_of_flow']]=upper_bound
        if rotation is False:        
            res_parcels.plot(ax=ax, color='tan')
        else:
            res_parcels_rot = gpd.GeoSeries(g for g in res_parcels['geometry'])
            res_parcels_rot = res_parcels_rot.rotate(rotation, origin=center_point)
            res_parcels_rot.plot(ax=ax, column='cost_of_flow',cmap='bwr', legend=False)      
    else:
        if rotation is False:        
            res_parcels.plot(ax=ax, color='tan')
        else:
            res_parcels_rot = gpd.GeoSeries(g for g in res_parcels['geometry'])
            res_parcels_rot = res_parcels_rot.rotate(rotation, origin=center_point)
            res_parcels_rot.plot(ax=ax, color='tan')
    
    # plot the resource parcels
    if resource_parcels is None:
        pass
    elif isinstance(resource_parcels, gpd.GeoDataFrame):
        if rotation is False:        
            resource_parcels.plot(ax=ax, color='mediumseagreen', edgecolor='darkgreen', linewidth=parcel_linewidth)
        else:
            resource_parcels_rot = gpd.GeoSeries(g for g in resource_parcels['geometry'])
            resource_parcels_rot = resource_parcels_rot.rotate(rotation, origin=center_point)
            resource_parcels_rot.plot(ax=ax, color='mediumseagreen', edgecolor='darkgreen', linewidth=parcel_linewidth)
    elif isinstance(resource_parcels, list):
        colors = ['#ff7f00', '#1f78b4',  '#33a02c','#b2df8a',
                  '#fb9a99', '#e31a1c', '#fdbf6f', '#a6cee3']
        edgecolors=['#663300','#0c3048','#144012','#46711f',
                    '#9b0806','#5b0a0b','#905202','#265b78']
        for i, parcel in enumerate(resource_parcels):
            if rotation is False:        
                parcel.plot(
                    ax=ax, color=colors[i], edgecolor=edgecolors[i], linewidth=parcel_linewidth)
            else:
                parcel_rot = gpd.GeoSeries(g for g in parcel['geometry'])
                parcel_rot = parcel_rot.rotate(rotation, origin=center_point)
                parcel_rot.plot(ax=ax, color=colors[i], edgecolor=edgecolors[i], linewidth=parcel_linewidth)

    # option to plot loss of access parcels
    if loss_access_parcels is None:
        pass
    else:
        if rotation is False:        
            loss_access_parcels.plot(ax=ax, color='saddlebrown')
        else:
            loss_access_parcels_rot = gpd.GeoSeries(g for g in loss_access_parcels['geometry'])
            loss_access_parcels_rot = loss_access_parcels_rot.rotate(rotation, origin=center_point)
            loss_access_parcels_rot.plot(ax=ax, color='saddlebrown')

    # plot background light gray roads
    # plot different background edges if value given
    if background_edges is None:
        G_gdf_edges = ox.graph_to_gdfs(G=G, nodes=False)
        if rotation is False:
            G_gdf_edges.plot(ax=ax, edgecolor='lightgray',
                             linewidth=linewidth, zorder=-1)
        else:
            G_gdf_edges_rot = gpd.GeoSeries(g for g in G_gdf_edges['geometry'])
            G_gdf_edges_rot = G_gdf_edges_rot.rotate(rotation, origin=center_point)
            G_gdf_edges_rot.plot(ax=ax, edgecolor='lightgray', linewidth=linewidth,zorder=-1)

    else:
        G_gdf_edges_c = ox.graph_to_gdfs(G=background_edges, nodes=False)
        if rotation is False:
            G_gdf_edges_c.plot(ax=ax, edgecolor='lightgray', linewidth=linewidth,zorder=-1)
        else:
            G_gdf_edges_rot = gpd.GeoSeries(
                g for g in G_gdf_edges_c['geometry'])
            G_gdf_edges_rot = G_gdf_edges_rot.rotate(rotation, origin=center_point)
            G_gdf_edges_rot.plot(ax=ax, edgecolor='lightgray', linewidth=linewidth,zorder=-1)

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
        cmap = (mpl.colors.ListedColormap(['lightgray', 'paleturquoise','deepskyblue', 'royalblue', 'mediumblue']))

        # normalize the colors and create mapping color object
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        cm_map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        color_series=vals.map(cm_map.to_rgba)

        # TODO: incorporate into single plot, for now just create a second figure for the scale bar
        fig2, ax2 = plt.subplots(figsize=(6, 1))
        fig2.subplots_adjust(bottom=0.5)
        fig2.colorbar( mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
                        cax=ax2,
                        ticks=bounds,
                        orientation='horizontal',
                        label='Max Inundation on Road (mm)')


        G_gdf_edges = ox.graph_to_gdfs(G=G, nodes=False)
        if rotation is False:
            G_gdf_edges.plot(ax=ax, edgecolor=color_series,
                             linewidth=linewidth)
        else:
            G_gdf_edges_rot = gpd.GeoSeries(g for g in G_gdf_edges['geometry'])
            G_gdf_edges_rot = G_gdf_edges_rot.rotate(rotation, origin=center_point)
            G_gdf_edges_rot.plot(ax=ax, edgecolor=color_series, linewidth=linewidth)

    # plot edges based on width
    if edge_width is None:
        pass
    else:
        flow_dict = nx.get_edge_attributes(G, edge_width)
        flows = []
        for x in flow_dict:
            flows.append(flow_dict[x])
        unique_flows = np.unique(flows).tolist()
        # unique_flows.pop(0)
        old_min = min(unique_flows)
        old_max = max(unique_flows)
        new_min = linewidth
        new_max = linewidth_max

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

                
                if rotation is False:
                    selected_edges.plot(ax=ax, edgecolor='black', linewidth=width)
                else:
                    G_gdf_edges_rot = gpd.GeoSeries(g for g in selected_edges['geometry'])
                    G_gdf_edges_rot = G_gdf_edges_rot.rotate(rotation, origin=center_point)
                    G_gdf_edges_rot.plot(ax=ax, edgecolor='black', linewidth=width)

    # optional plot inundation raster
    if inundation is None:
        pass
    else:
        print('PROBABLY BUG: NEED TO ADD ROTATION')
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
            
            if rotation is False:
                inset.plot(ax=ax,edgecolor='forestgreen',linewidth=3,zorder=1000)
            else:
                inset_rot = gpd.GeoSeries(g for g in inset['geometry'])
                inset_rot = inset_rot.rotate(rotation, origin=center_point)
                inset_rot.plot(ax=ax,edgecolor='forestgreen',linewidth=3,zorder=1000)




            inset = gpd.read_file(inset)
            inset.boundary.plot(ax=ax,
                                edgecolor='forestgreen',
                                linewidth=3,
                                zorder=1000)

    # plot tight layout
    fig.tight_layout()

    # optional add scale bar
    if scalebar is False:
        pass
    else:
        scale=ax.add_artist(ScaleBar(1,
                                frameon=True,
                                box_color='lightgray',
                                location='lower left'))
        scale.zorder=10000
        # to determine width and length of head, need to convert display coordinates to data coordinates
        # display coordinates is in inches
        # want an arrow that has a width of half an inch and a length of 1 inch
        inv=ax.transData.inverted()
        trans_size = inv.transform([(0,0),(0.5,1.5)])
        head_width = trans_size[1][0]-trans_size[0][0]
        head_length = trans_size[1][1]-trans_size[0][1]

        x_tail = 0.90
        y_tail = 0.04
        x_head = 0.90
        y_head = 0.105
        old_center = ((x_tail+x_head)/2,(y_tail+y_head)/2)

        tails=(x_tail,y_tail)
        heads=(x_head,y_head)

        if rotation is not False:
            # rotate points around center - rotation is postive when counterclockwise
            angle=rotation*np.pi/180
            x_tail_r = (x_tail-old_center[0])*np.cos(angle)-(y_tail-old_center[1])*np.sin(angle)+old_center[0]
            y_tail_r = (x_tail-old_center[0])*np.sin(angle)+(y_tail-old_center[1])*np.cos(angle)+old_center[1]        
            x_head_r = (x_head-old_center[0])*np.cos(angle)-(y_head-old_center[1])*np.sin(angle)+old_center[0]
            y_head_r = (x_head-old_center[0])*np.sin(angle)+(y_head-old_center[1])*np.cos(angle)+old_center[1]
            tails=(x_tail_r,y_tail_r)
            heads=(x_head_r,y_head_r)
        arrow = mpatches.FancyArrowPatch(tails, 
                                         heads,
                                         mutation_scale=25,
                                         transform=ax.transAxes,
                                         facecolor='black',arrowstyle='fancy',
                                         zorder=1000)
        ax.add_patch(arrow)

        # variables for adding N
        xy = (0.95, 0.1)
        xy_text=(0.95, 0.01)

        if rotation is not False:
            if abs(rotation) >= 45:
                xy_text = (0.95, 0.02)

        # add N
        ax.annotate('N', 
                    xy=xy, 
                    xycoords='axes fraction',
                    xytext=xy_text,
                    textcoords='axes fraction',
                    va='bottom', 
                    ha='center',
                    fontsize=14)

    if save_loc is None:
        pass
    else:
        plt.savefig(save_loc)
    
    return(fig, ax)


def summary_function(G: nx.Graph):
    '''
    Return quick summary information regarding a network.

    This function can be useful for displaying quick information describing the attributes of a network. It should be expanded upon to show more relevent information compared to its current state. 

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


def inundate_network(G: nx.Graph, 
                     path: str, 
                     inundation: str, 
                     G_capacity: str = 'capacity',
                     G_width: str = 'width',
                     G_speed: str = 'speed_kph', 
                     G_length: str = 'length',
                     pp: int = 8, 
                     name: str = 'AOI_Graph_Inundated'):
    '''
    Create a new graph network based on the imapct of an inundation layer.

    An inundation layer is overlayed on a graph network to determine the maximum depths that intersect each road segment. The maximum depth is then used to calculate the reduction in speed and capacity on that entire road segment using a conservative, moderate, and aggressive equation. This function adds attributes to the graph network including:

    - max_inundation_mm, the maximum inundation in milimeters
    - inundation_capacity_{mod}, the reduced capacity of the road segement
    - kph_{mod}, the reduced speed of the road segment
    - inundation_travel_time_{mod}, the new travel time of the road segment

    For inundation_capacity_, kph_, and inundation_travel_time_, there are three seperate occurances of each attribute with the suffix 'mod', 'agr', or 'con' refering to the moderate, aggressive, or conservative depth-disruption equation that is used in its calculation. All are quadratically decreasing depending on the maximum allowable road depth which are 300, 150, or 600 respectively. Any road segments with a greater depth of water on a road segement that these values has a travel time and capacity set equal to 0.

    This function also relies on parallel processing. Zonal statistics for each road segment is time consuming, but also independent of each other and therefore is easily parallelized. In one test, a single thread took 10 minutes to process 36,000 edges. The same study area took less than 3 minutes to process with 4 threads. 

    The inundation raster and graph network should have the same coordinate reference system. 

    :param G: The graph network
    :type G: networkx.Graph [Multi, MultiDi, Di]
    :param path: File path to save output graph network to. Will have name based on name argument
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


    with rasterio.open(inundation) as src:
        array=src.read(1)
        transform =src.transform

        # BEGIN PARALLEL PROCESSING
        # geometry column number
        geo_col_num = edges.columns.get_loc("geometry")
        # the individual inundation zonal statistic task
        def inundate_zs(n):
            return zonal_stats(edges.iat[n,geo_col_num], 
                               raster=array,
                               affine=transform,
                               stats='max')
        # create the sequence of numbers for each edge
        n=range(0,len(edges))
        # create a pool of workers and run the function inundate_zs for each zone
        # map() function creates the background batches
        pool=Pool(pp)
        results=pool.map(inundate_zs, n)
        # close the pool
        pool.close()
        pool.join()
        # END PARALLEL PROCESSING

    # convert results to list - if None, where raster doesn't intersect shapefile, replace with 0
    results_as_list = [0 if d[0]['max'] is None else round(d[0]['max']*1000/10)*10 for d in results]
    # relate back to edges geodataframe
    edges['max_inundation_mm'] = results_as_list
    
    # any edge that is a bridge, set max_inundation_mm to 0
    edges.loc[edges['bridge'] == 'yes', 'max_inundation_mm'] = 0

    
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

    # save edited graph
    new_graph = ox.graph_from_gdfs(nodes, edges)
    save_2_disk(G=new_graph, path=path, name=name)
    # eliminate variables to save memory
    n=None
    nodes=None
    edges=None
    results=None
    geo_col_num=None
    results_as_list=None
    

    return (new_graph)


def flow_decomposition(G: nx.Graph, 
                       res_points: gpd.GeoDataFrame, 
                       dest_points: gpd.GeoDataFrame,
                       res_parcels: gpd.GeoDataFrame,
                       dest_parcels: gpd.GeoDataFrame, 
                       G_demand: str ='demand', 
                       G_capacity: str ='capacity', 
                       G_weight: str ='travel_time', 
                       dest_method: str = 'single'):
    '''
    Given a network, G, decompose the maximum flow (w/ minimum cost) into individual flow paths.

    This function determines the maximum feasible flow (with a minimum cost) and then decomposes the resultant network solution into a set of individual feasible flow paths. The decomposed information including path, sink, and cost of each path is related back to and returned in the residential and destination parcel boundaries outputs. The returned decomposed flow is just one possible solution. Flows are decomposed by looking at the shortest path from a super sink to a super source, resulting in a "greedy" algorithm that routes closer parcels first, and furthest parcels last. Therefore, the decomposed network may have a higher variance, as compared to other possible decomposition solutions. 

    If multiple diverging flow paths exist (i.e., a residential parcel has the same sourec but a different path to a potentially different resource parcel), the function becomes stochastic. This means that each residential parcel sharing a common source node is assigned a random feasible sink-path-cost from the feasible set.


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
    :param G_demand: name of attribute in G refering to node demand, *Default='demand*
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
    G = copy.deepcopy(G)
    for x in G.edges:
            G.edges[x][G_weight] = round(G.edges[x][G_weight])
    
    # if only using a single destination file, proceed normally
    if isinstance(dest_points, gpd.GeoDataFrame):

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
                kwargs = {f"{G_weight}": 0,
                        f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}
                G.add_edge(0, unique_node, **kwargs)
                # since we added an artificial source node, all original source nodes must have a zero demand
                G.nodes[unique_node][G_demand] = 0

            # add the super_sink that everything goes to
            G.add_nodes_from([(99999999, {G_demand: positive_demand})])

            # artificial node id tracker - useful in maintianing compatability with dtype
            dest_node_ids = 99999998
            # identify the destination nodes, and create artifical sink edges
            # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to
            for idx, dest_parcel in dest_parcels.iterrows():
                dest_node = dest_points[dest_points.geometry.within(dest_parcel.geometry)]
                #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
                for i, node in dest_node.iterrows():
                    dest_node_ids -= 1
                    # add the dest node to the graph using OSMID as its ID
                    G.add_nodes_from([(dest_node_ids, {G_demand: 0})])
                    # for sink insights, need to set dest_node_ids to attribute within dest_points
                    dest_points.loc[dest_points.geometry == node['geometry'], 'dest_node_ids'] = dest_node_ids


                    # add links from nearest intersections to parcel centroid
                    for nearest_intersection in unique_dest_nodes_list[idx]:
                        kwargs = {G_weight: 0, G_capacity: 999999999}
                        G.add_edge(nearest_intersection, dest_node_ids, **kwargs)

                    # add link from parcel centroid to super sink
                    # TODO2: This is where I can specifically add capacities for each specific grocery store
                    # FOR NOW: just setting capacity to a high amount
                    kwargs = {G_weight: 0, G_capacity: 999999999}
                    G.add_edge(dest_node_ids, 99999999, **kwargs)

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
        # CURRENTLY-> set up as greedy algorithm, will produce highest variance in final flows if flow paths from common sources split
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
                    near_nodes = list(map(int, dest_points.loc[dest_points['dest_node_ids'] == sink, 'nearest_nodes'].iloc[0].split(' ')))
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
        
        # combine similar columns and remove unnecessary ones if necessary
        try:
            res_points['cost_of_flow'] = res_points['cost_of_flow_y'].combine_first(res_points['cost_of_flow_x'])
            res_points['service'] = res_points['service_y'].combine_first(res_points['service_x'])
            res_points.drop(columns=['cost_of_flow_y','cost_of_flow_x','service_y','service_x'], inplace=True)
        except:
            pass

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

        # append dest points with appropriate flow information
        if dest_method == 'single':
            dest_points = dest_points.merge(sink_insights_df, how='left', left_on='nearest_node', right_on='index')
        if dest_method == 'multiple':
            dest_points = dest_points.merge(sink_insights_df, how='left', left_on='dest_node_ids', right_on='index')
        
        # Spatial join the res_parcels/res_points and dest_parcels/dest_points data
        res_parcels = gpd.sjoin(res_parcels, res_points)
        dest_parcels = gpd.sjoin(dest_parcels, dest_points)

        # sjoin is messing with dataframe attributes, editing so returns only attributes we are interested in from this function
        # to avoid data issues, for now only returning attributes created during this funciton
        res_attributes = ['geometry', 'nearest_node','source', 'sink', 'walkable?', 'path', 'cost_of_flow', 'service']
        dest_attributes = ['geometry', 'dest_node_ids', 'Number of unique paths','Total flow w/o walking','Total flow w/ walking','nearest_nodes']
        # dest_parcels['nearest_nodes'] = dest_parcels['nearest_nodes_right']
        res_parcels.drop(columns=[col for col in res_parcels if col not in res_attributes], inplace=True)
        dest_parcels.drop(columns=[col for col in dest_parcels if col not in dest_attributes], inplace=True)

        # METRICS: Could output in some way in the future if necessary 
        # # Right now, no need to print, but could export in some manner
        # # convert to list of flows
        # cost_list = sorted(res_parcels['cost_of_flow'].tolist())
        # #remove nans 
        # cost_list = [x for x in cost_list if math.isnan(x)==False]
        # # calculate interquartile range
        # q1,q3,=np.percentile(cost_list, [25,75])
        # iqr=q3-q1
        # upper_bound=q3+(3*iqr)
        # # calculate percentage of flows that can be considered a major outlier
        # perc_outlier = sum(x>=upper_bound for x in cost_list)/positive_demand *100
        # # create list with 'masked' outliers
        # cost_list_wo = [upper_bound if x >= upper_bound else x for x in cost_list]
        # # median flow cost(excluding including)
        # med_cost_in = np.median(cost_list)
        # # total flow cost(excluding including)
        # total_cost_in = sum(cost_list)
        # # average flow cost(excluding including)
        # mean_cost_in = np.mean(cost_list)
        # # median flow (excluding outliers)
        # med_cost_ex = np.median(cost_list_wo)
        # # total flow cost(excluding outliers)
        # total_cost_ex = sum(cost_list_wo)
        # # average flow cost(excluding outliers)
        # mean_cost_ex = np.mean(cost_list_wo)

        # # print(f'percent of flow major outlier: {perc_outlier}')
        # # print(f'cost of flow with outliers: {total_cost_in}')
        # # print(f'cost of flow without outliers: {total_cost_ex}')
        # # print(f'average cost of flow with outliers: {mean_cost_in}')
        # # print(f'average cost of flow without outliers: {mean_cost_ex}')
        # # print(f'median cost of flow with outliers: {med_cost_in}')
        # # print(f'median cost of flow without outliers: {med_cost_ex}')


        # return the decomposed paths dict of dict, the sink_insights dict of dict, and the res_locs/dest_locs geodataframes
        return decomposed_paths, sink_insights, res_parcels, dest_parcels

    # if multiple destinations, simply run function mulitple times and combine outputs
    # congestion isn't considered, so don't need to run them all at the same time
    elif isinstance(dest_points, list) or isinstance(dest_parcels, list):
        
        decomposed_paths_final=[]
        sink_insights_final=[]

        res_attributes = ['geometry', 'nearest_node']

        for idx, dest_point in enumerate(dest_points):
            
            if dest_method=='single':
                # Run max_flow_parcels to get flow dictionary
                flow_dict, flow_cost, max_flow, access = max_flow_parcels(G=G, 
                                                                            res_points=res_points, 
                                                                            dest_points=dest_point, 
                                                                            G_capacity=G_capacity, 
                                                                            G_weight=G_weight, 
                                                                            G_demand=G_demand)
                
                #Run nearest_nodes to get unique origins and destinations
                G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points_out = nearest_nodes(
                    G=G, res_points=res_points, dest_points=dest_point, G_demand=G_demand)
                
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
                # keep track of edges added to make it easier to remove them after each iteration
                edges_added=[]
                # Run max_flow_parcels to get flow dictionary
                flow_dict, flow_cost, max_flow, access = max_flow_parcels(G=G,
                                                                        res_points=res_points,
                                                                        dest_points=dest_point,
                                                                        G_capacity=G_capacity,
                                                                        G_weight=G_weight,
                                                                        G_demand=G_demand,
                                                                        dest_method='multiple',
                                                                        dest_parcels=dest_parcels[idx])
                
                #Run nearest_nodes to get unique origins and destinations
                G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcel, dest_point = nearest_nodes_vertices(
                    G=G, res_points=res_points, dest_parcels=dest_parcels[idx], G_demand=G_demand, dest_points=dest_point)

                # add artifical source node
                G.add_nodes_from([(0, {G_demand: positive_demand*-1})])
                # add edge from artifical source node to real source nodes with 0 weight and capacity equal to demand
                sums = 0
                for unique_node in unique_origin_nodes:
                    sums -= G.nodes[unique_node][G_demand]
                    kwargs = {f"{G_weight}": 0,
                            f"{G_capacity}": G.nodes[unique_node][G_demand]*-1}
                    G.add_edge(0, unique_node, **kwargs)
                    edges_added.append((0,unique_node))
                    # since we added an artificial source node, all original source nodes must have a zero demand
                    G.nodes[unique_node][G_demand] = 0

                # add the super_sink that everything goes to
                G.add_nodes_from([(99999999, {G_demand: positive_demand})])

                # artificial node id tracker - useful in maintianing compatability with dtype
                dest_node_ids = 99999998
                # identify the destination nodes, and create artifical sink edges
                # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to
                for i, dest_parcel in dest_parcel.iterrows():
                    dest_node = dest_point[dest_point.geometry.within(dest_parcel.geometry)]
                    #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
                    for x, node in dest_node.iterrows():
                        dest_node_ids -= 1
                        # add the dest node to the graph using OSMID as its ID
                        G.add_nodes_from([(dest_node_ids, {G_demand: 0})])
                        # for sink insights, need to set dest_node_ids to attribute within dest_points
                        dest_point.loc[dest_point.geometry == node['geometry'], 'dest_node_ids'] = dest_node_ids

                        # add links from nearest intersections to parcel centroid
                        for nearest_intersection in unique_dest_nodes_list[i]:
                            kwargs = {G_weight: 0, G_capacity: 999999999}
                            G.add_edge(nearest_intersection, dest_node_ids, **kwargs)
                            edges_added.append((nearest_intersection, dest_node_ids))

                        # add link from parcel centroid to super sink
                        # TODO2: This is where I can specifically add capacities for each specific grocery store
                        # FOR NOW: just setting capacity to a high amount
                        kwargs = {G_weight: 0, G_capacity: 999999999}
                        G.add_edge(dest_node_ids, 99999999, **kwargs)
                        edges_added.append((nearest_intersection, dest_node_ids))
            # create a new copy of res_points
            res_points_c = res_points.copy()

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
            # CURRENTLY-> set up as greedy algorithm, will produce highest variance in final flows if flow paths from common sources split
            while len(G_inter.edges)>0:
                # Find the shortest path from the artifical source to artificial sink
                path = ox.distance.shortest_path(G_inter,0,99999999,weight=G_weight)
                # BUG: there are times where edges still exist but no viable path, likely due to rounding errors?
                if path is None:
                    break
                else:
                    # create empty list of flows -> allowable flow through each edge along shortest path
                    flows = []
                    # length of the shortest path
                    len_path=len(path)
                    # calculate the cost per unit flow along this path
                    cost = nx.path_weight(G_inter, path, weight=G_weight)
                    # for each edge in the path, append flows with the allowable flow
                    for i, x in enumerate(path):
                        if i==len_path-1:
                            pass
                        else:
                            flows.append(G_inter[path[i]][path[i+1]]['allowable_flow'])
                    # limiting flow is the minimum allowable flow along the edges in the path
                    limiting_flow=min(flows)
                    # for each edge along the path, subtract the limiting flow from allowable flow
                    # if the allowable flow is reduced to zero, remove the edge
                    for i, x in enumerate(path):
                        if i == len_path-1:
                            pass
                        else:
                            G_inter[path[i]][path[i+1]]['allowable_flow'] -= limiting_flow
                            if G_inter[path[i]][path[i+1]]['allowable_flow'] == 0:
                                G_inter.remove_edge(path[i],path[i+1])
                    
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
                            sink_insights[sink] = {'Number of unique paths': 1, 'Total flow w/o walking': limiting_flow, 'Total flow w/ walking': limiting_flow + len(res_points_c.loc[res_points_c['nearest_node'] == sink])}
                    
                    if dest_method == 'multiple':
                        if sink in sink_insights.keys():
                            sink_insights[sink]['Number of unique paths'] +=1
                            sink_insights[sink]['Total flow w/o walking']+= limiting_flow
                            sink_insights[sink]['Total flow w/ walking']+= limiting_flow 
                        else:
                            near_nodes = list(map(int, dest_point.loc[dest_point['dest_node_ids'] == sink, 'nearest_nodes'].iloc[0].split(' ')))
                            sink_insights[sink] = {'Number of unique paths': 1, 'Total flow w/o walking': limiting_flow, 'Total flow w/ walking': limiting_flow + len(res_points_c.loc[res_points_c['nearest_node'].isin(near_nodes)])}

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
                    res_points_c.loc[res_points_c['nearest_node'] == node,['service']] = 'no'
                    res_points_c.loc[res_points_c['nearest_node'] == node, ['cost_of_flow']] = np.NaN
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
            res_points_c = res_points_c.merge(decomp_subset_df[['source','sink','cost_of_flow','walkable?','service','path']], how='left', left_on='nearest_node', right_on='source')
            # combine similar columns and remove unnecessary ones if necessary
            try:
                res_points_c['cost_of_flow'] = res_points_c['cost_of_flow_y'].combine_first(res_points_c['cost_of_flow_x'])
                res_points_c['service'] = res_points_c['service_y'].combine_first(res_points_c['service_x'])
                res_points_c.drop(columns=['cost_of_flow_y','cost_of_flow_x','service_y','service_x'], inplace=True)
            except:
                pass

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
                for i, pathway in enumerate(row['Path']):
                    flow = int(row['Flow'][i])
                    source = int(row['Source'][0])
                    sink = row['Sink'][i]
                    path = pathway
                    cost_per_flow = row['Cost Per Flow'][i]

                    # randomly select res_points with same source and no path information
                    # the syntax for query is a bit odd: needed @ to reference variable and find nans using !=
                    query = res_points_c.query("nearest_node == @source and sink != sink").sample(n=flow).index
                    res_points_c.loc[query, ['sink']] = sink
                    res_points_c.loc[query, ['path']] = path
                    res_points_c.loc[query, ['cost_of_flow']] = cost_per_flow
                    res_points_c.loc[query, ['walkable?']] = np.NaN
                    res_points_c.loc[query, ['service']] = 'yes'

            #4. set shared nodes to walkable and set remaining features
            res_points_c.loc[res_points_c['nearest_node'].isin(shared_nodes), ['walkable?']] = 'yes'
            res_points_c.loc[res_points_c['nearest_node'].isin(shared_nodes), ['service']] = 'yes'
            res_points_c.loc[res_points_c['nearest_node'].isin(shared_nodes), ['cost_of_flow']] = 0
            res_points_c.loc[res_points_c['nearest_node'].isin(shared_nodes), ['sink']] = res_points['nearest_node']
            res_points_c.loc[res_points_c['nearest_node'].isin(shared_nodes), ['path']] = 0

            # convert sink_insights to dataframe
            sink_insights_df = pd.DataFrame.from_dict(sink_insights, orient='index', columns=['Number of unique paths',
                                                                                                'Total flow w/o walking',
                                                                                                'Total flow w/ walking']).reset_index()

            # append dest points with appropriate flow information
            if dest_method == 'single':
                dest_point = dest_point.merge(sink_insights_df, how='left', left_on='nearest_node', right_on='index')
            if dest_method == 'multiple':
                dest_point = dest_point.merge(sink_insights_df, how='left', left_on='dest_node_ids', right_on='index')
            
            # Spatial join the res_parcels/res_points and dest_parcels/dest_points data
            res_parcels = gpd.sjoin(res_parcels, res_points_c)
            dest_parcels[idx] = gpd.sjoin(dest_parcels[idx], dest_point)

            # sjoin is messing with dataframe attributes, editing so returns only attributes we are interested in from this function
            # to avoid data issues, for now only returning attributes created during this funciton
            
            # since iteratively adding to res_parcels, need to rename columns and add to res_attributes not to remove
            res_attributes.append(f'source_r{idx}')
            res_attributes.append(f'sink_r{idx}')
            res_attributes.append(f'walkable?_r{idx}')
            res_attributes.append(f'path_r{idx}')
            res_attributes.append(f'cost_of_flow_r{idx}')
            res_attributes.append(f'service_r{idx}')
            res_parcels.rename(columns={'source': f'source_r{idx}',
                                        'sink': f'sink_r{idx}',
                                        'walkable?' : f'walkable?_r{idx}',
                                        'path' : f'path_r{idx}',
                                        'cost_of_flow' : f'cost_of_flow_r{idx}',
                                        'service' : f'service_r{idx}'}, inplace=True)

            dest_attributes = ['geometry', 'dest_node_ids', 'Number of unique paths','Total flow w/o walking','Total flow w/ walking','nearest_nodes']
            # dest_parcels[idx]['nearest_nodes'] = dest_parcels[idx]['nearest_nodes_right']
            res_parcels.drop(columns=[col for col in res_parcels if col not in res_attributes], inplace=True)
            dest_parcels[idx].drop(columns=[col for col in dest_parcels[idx] if col not in dest_attributes], inplace=True)

            decomposed_paths_final.append(decomposed_paths)
            sink_insights_final.append(sink_insights)

            # before next iteration, remove the edges that were added 
            # TODO: don't know if i should ignore capacity (almost always feasible) or include capacity constraint (order of dest_points matters)
            # TODO: FOR NOW, ignoring compounding capacity on each iteration
            if dest_method == 'single':
                for x in unique_dest_nodes:
                    G.remove_edge(x, 99999999)
                G.remove_node(99999999)
            if dest_method == 'multiple':
                G.remove_edges_from(edges_added)
                G.remove_node(99999999)

        # return the decomposed paths dict of dict, the sink_insights dict of dict, and the res_locs/dest_locs geodataframes
        return decomposed_paths_final, sink_insights_final, res_parcels, dest_parcels



# HELPER FUNCTIONS:
#  utilized in traffic assignment

def shortestPath(origin: int, 
                 capacity_array: np.array, 
                 weight_array: np.array):
    """
    This method finds the shortest path from origin to all other nodes in the network. Uses Dijkstra's algorithm.
    
    :param origin: origin node
    :type origin: int
    :param capacity_array: node-node array of edge capacities
    :type capacity_array: np.array
    :param weight_array: node-node array of edge weights (i.e., travel times)
    :type weight_array: np.array

    :returns:
        - **backnodes**, np.array indiciating previous node in shorest path
        - **costlabel**, list of costs associated with each node
    :rtype: tuple

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

    # Step 2: Initialize SEL (scan eligable list) with origin
    sel = [origin]
    # TODO: Can i delete this x=0?
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

def shortestPath_heap(origin: int, 
                      capacity_array: Union[np.array, sparse.lil_array], 
                      weight_array: Union[np.array, sparse.lil_array], 
                      adjlist: list[list], 
                      destination: int = False, 
                      sparse_array: bool = True,
                      limit=np.inf):
    """
    This method finds the shortest path from origin to all other nodes in the network. 
    Uses a binary heap priority que in an attempt to speed up the shortestPath function. Uses Dijkstra's algorithm and built in python function heapq.
    
    This function replaces shorestPath because it is a faster alternative due to the binary heap priority que.

    :param origin: origin node
    :type origin: int
    :param capacity_array: node-node array of edge capacities
    :type capacity_array: np.array or sparse.csr_array
    :param weight_array: node-node array of edge weights (i.e., travel times)
    :type weight_array: np.array or sparse.csr_array
    :param adjlist: adjacency list of network edges
    :type adjlist: nested list
    :param destination: if not False, set to destination node, terminates algorithm when shortest path to destination node is found
    :type destination: int
    :param sparse_array: If True, input arrays are sparse.csr_arrays, if False, input arrays are np.arrays
    :type sparse_array: bool


    :returns:
        - **backnodes**, np.array indiciating previous node in shorest path
        - **costlabel**, list of costs associated with each node
    :rtype: tuple

    # """

    # UNTESTED: chatgpt vesrion start
    # Initialize the data structures
    # lil_array=weight_array
    # source=origin
    # target=destination
    # num_nodes = lil_array.shape[0]
    # visited = np.zeros(num_nodes, dtype=bool)
    # costlabel = np.full(num_nodes, np.inf)
    # backnode = np.full(num_nodes, -1, dtype=int)

    # # Set the cost of the source node to zero
    # costlabel[source] = 0

    # # Initialize the priority queue with the source node
    # pq = [(0, source)]

    # # Main loop
    # while pq:
    #     # Pop the node with the smallest tentative distance from the priority queue
    #     current_dist, current_node = heapq.heappop(pq)

    #     # If the node has already been visited or is the target node, continue to the next iteration
    #     if visited[current_node] or (target is not None and current_node == target):
    #         continue

    #     # Mark the current node as visited
    #     visited[current_node] = True

    #     # Update the costlabel and backnode of the neighbors of the current node
    #     for neighbor, weight in zip(lil_array.rows[current_node], lil_array.data[current_node]):
    #         new_cost = costlabel[current_node] + weight
    #         if new_cost < costlabel[neighbor]:
    #             costlabel[neighbor] = new_cost
    #             backnode[neighbor] = current_node
    #             # Add the neighbor to the priority queue with its tentative distance as the priority
    #             heapq.heappush(pq, (new_cost, neighbor))
    # chat gpt version end


    # # ORIGINAL algorithm
    # # Need the number of nodes for calculations
    # numnodes = np.shape(capacity_array)[0]

    # # # Set up backnode and cost lists to return
    # backnode = [-1] * numnodes
    # costlabel = [np.inf] * numnodes

    # # I am going to be implementing Dijkstra's Algorithm

    # # Step 1: Initialize all labels Li=infinity origin Ls=0
    # # costlabel already initialized, need to change cost label for origin
    # costlabel[origin] = 0
    # # Also creating f_label for final labels
    # f_label = []
    # pq = [(0, origin)]
    # # Step 2: Initialize SEL with origin, and turn into a heap
    # # sel = [origin]

    # # Create termination loop, if sel empty, terminate
    # while len(pq) > 0:
    #     # Step 3: Choose a node i from SEL with minimum L
    #     # doing this by creating a list of cost labels for each node in SEL
    #     cost, node_i = heapq.heappop(pq)
    #     if node_i in f_label:
    #         pass
    #     else:

    #         # Step 4: for all arc (i,j) repeat following
    #         for node_j in adjlist[node_i]:
    #             #also need to skip nodes in f, because already finalized
    #             if node_j in f_label:
    #                 pass
    #             else:
    #                 # Step 4a: calculate costlabel temp
    #                 if sparse_array is True:
    #                     l_temp = cost + weight_array[node_i, node_j]            
    #                 else:
    #                     l_temp = cost + weight_array[node_i][node_j]

    #                 # Step 4b: if new path is shorter, update cost and backnode
    #                 if l_temp < costlabel[node_j]:
    #                     costlabel[node_j] = l_temp
    #                     backnode[node_j] = node_i
    #                 # TODO: Step 4c: add j tos heap ->  no way of checking if already in or not, so added if statement at top to filter
    #                 heapq.heappush(pq, (costlabel[node_j], node_j))
    #         # Step 5: add i to f_label # TODO: I think this is indented correctly  but should test - was originally operating within outside of last else statement
    #         f_label.append(node_i)

    #     # STEP 6: if destination found AND finalized 
    #     if destination is not False:
    #         if node_j == destination and node_j in f_label:
    #             break

    # backnode = np.asarray(backnode)
    # costlabel = np.asarray(costlabel)
    # backnodes equal to list indicating previous node in shortest path
    # costlabel is list of costs associated with each node
    # # ORIGINAL algorithm end


    #  scipy sparse function
    # BUG: Sometimes dijkstras gets caught in infinite loop -> not sure why 
    # Bellman-Ford works but its incredibly slow
    # costlabel,backnode = sparse.csgraph.shortest_path(weight_array,
    #                                                   indices=origin,
    #                                                   return_predecessors=True,
    #                                                   method='D')
    costlabel,backnode = sparse.csgraph.dijkstra(weight_array,
                                                 indices=origin,
                                                 return_predecessors=True,
                                                 limit=limit)
    
    return (backnode, costlabel)

def pathTo(backnode, costlabel, origin, destination):
    """
    Given backnode costlabel, origin, and destination, return a list of nodes that is the shortest path.    
    """
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
    """
    Modified pathTo function to output node-node pairs instead of list of nodes to avoid having to do loops elsewhere
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

# PARALLELIZATION HELPER FUNCTIONS USED IN TRAFFIC ASSIGNMENT
def SPTT_parallel(n):
    """
    Calculate initial feasible flow in parallel. 

    # TODO: Needs to be much better documented, at a later date. 
    n is array version of OD_matrix where:
        n[0]: origin
        n[1]: flow
        n[2]: destination
        n[3]: what we should return
        n[4]: which weight array to use

    n[3] can be either 'cost', 'backnode, 'path', or 'path_cost', or 'path_cost_backnode'
    n[4] can be either 'weight_array' or 'weight_array_iter'

    # maybe following as well:
    n[5] = capacity_array
    n[6] = weight_array

    """
    if n[4] == 'weight_array':
        if len(n) == 5:
            backnode, costlabel = shortestPath_heap(origin=n[0],
                                                    capacity_array=capacity_array,
                                                    weight_array=weight_array,
                                                    adjlist=adjlist,
                                                    destination=n[2])
        else:
            backnode, costlabel = shortestPath_heap(origin=n[0],
                                                    capacity_array=n[5],
                                                    weight_array=n[6],
                                                    adjlist=adjlist,
                                                    destination=n[2])                                              
    elif n[4] == 'weight_array_iter':
        if len(n) == 5:
            backnode, costlabel = shortestPath_heap(origin=n[0],
                                                    capacity_array=capacity_array,
                                                    weight_array=weight_array_iter,
                                                    adjlist=adjlist,
                                                    destination=n[2])
        else:
            backnode, costlabel = shortestPath_heap(origin=n[0],
                                                    capacity_array=n[5],
                                                    weight_array=n[6],
                                                    adjlist=adjlist,
                                                    destination=n[2])

    if n[3] == 'cost':
        cost = costlabel[n[2]]
        return cost
    
    elif n[3] == 'backnode':
        return backnode

    else:
        path, cost = pathTo_mod(backnode=backnode, 
                                costlabel=costlabel, 
                                origin=n[0], 
                                destination=n[2])
        if n[3] == 'path':
            return path
        if n[3] == 'path_cost':
            return path, cost 
        if n[3] == 'path_cost_backnode':
            return path, cost, backnode
    

# PARALLELIZATION HELPER FUNCTION FOR REDUNDANCY METRICS
def redundancy_parallel(n):
    """
    Calculate individual redundancy score for
    n is OD_matrix where:
        n[0]: origin
        n[1]: nested list of destinations (by resource)

    need global variables as well:
        G_copy:             copy of network Graph
        G_weight_copy       copy of weight array 
        inc_copy            increment of isochrone
        max_length_copy     maximum acceptable travel time
    """

    source=n[0]
    targets=n[1]

    # calculate shortest path matrix for source
    distance = nx.single_source_dijkstra_path_length(G=G_copy, source=source, weight=G_weight_copy)

    # initialize output [source, target, val, resource number]
    output=[]
    # find distances
    for i, resource in enumerate(targets):
        for target in resource:
            length = distance.get(target)
            # if no path exists, give a score of 0
            if length is None:
                output.append([source,target,0, int(i)])
            # if path does exist,
            else:
                # determine score based on determined length
                # if length is greater than maximum
                if length > max_length_copy:
                    output.append([source, target, 0, int(i)])
                # if length is equal to the maximum length
                elif length == max_length_copy:
                    output.append([source, target, 1, int(i)])
                else:
                    # round to next highest multiple of increment
                    val = inc_copy * math.ceil(length/inc_copy)
                    # subtract from max length, divide by increment to determine score
                    val = (max_length_copy-val)/inc_copy
                    output.append([source,target,val,int(i)])
    return output


def edge_disjoint_paths(n):
    """
    Parallel processing of edge_disjoint paths for network redundancy metric
    n: origin
    output takes form of:
    [origin, score for resource1, score for resource2, etc.]
    """

    # first find the neighbors of neighbors of origins
    origin_nn=[]
    origin_n1=[node for node in G_copy.neighbors(n)]
    origin_nn = origin_nn + origin_n1
    for node in origin_n1:
        origin_n2=[neighbor for neighbor in G_copy.neighbors(node)]
        # if origin is a neighbor remove it
        try: 
            origin_n2.remove(n)
        except ValueError:
            pass
        origin_nn += origin_n2
        

        # # attempt at neighbors-neighbors-neighbors
        # for node2 in origin_n2:
        #     origin_n3 = [neighbor for neighbor in G_copy.neighbors(node2)]
        #     # if origin is a neighbor remove it
        #     try: 
        #         origin_n3.remove(n)
        #     except ValueError:
        #         pass
        #     origin_nn += origin_n3

    # remove duplicates 
    origin_nn = list(set(origin_nn))

    # add artificial edges from origin to the neighbors-of-neighbors values
    for node in origin_nn:
        kwargs = {G_weight_copy: 1,
                  G_capacity_copy: 999999999}
        G_copy.add_edge(n, node, **kwargs)

    H = build_auxiliary_edge_connectivity(G_copy)
    R = build_residual_network(H, 'capacity')

    output=[]
    # for each resource type,
    for idx, node in enumerate(sink_nodes_copy):
        output.append(0)

        # if edge exists between origin and resource than they are connected and should be given maximum score
        # if G_copy.has_edge(n,node):
        # if edge resource node is in nn, then they are connected
        if node in origin_nn:
            # will filter out values that are -1 and set to the maximum
            output[idx]=-1
        else:
            try:
                k = len(list(nx.edge_disjoint_paths(G_copy,n,node, auxiliary=H,residual=R)))
            # if no path exists,
            except nx.NetworkXNoPath:
                k=0
            output[idx]=k
    output.insert(0, n)

    return output


# Can probably remove below functions
def get_path_length(G, path, weight):
    length = 0
    if len(path) > 1:
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]

            length += G[u][v][weight]
    return (length)

def k_shortest_paths(G, source, target, k, weight=None):
    return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))

def k_shortest_paths(n):
    """
    #TODO:
    Has major bugs needs to be addressed

    Calculate the number of viable paths

    # TODO: Needs to be much better documented, at a later date. 
    n is array version of OD_matrix where:
        n[0]: source
        n[1]: target
        n[2]: shared nodes for that resource
    """
    source=n[0]
    target=n[1]
    shared_nodes=n[2]

    # # k-shortest paths
    # TODO: make a seperate function
    # # determine if path exists
    # if nx.has_path(G_copy,source,target) is False:
    #     pass
    # else:
    #     # calculate shortest path
    #     cost = nx.shortest_path_length(G_copy,
    #                                     source=source,
    #                                     target=target,
    #                                     weight=G_weight_copy)

    #     # if the shortest path cost is greater than acceptable cost:
    #     if cost > max_tt_copy:
    #         output = [source, target, 0]

    #     # if source shares a node with that resources location, set to max paths
    #     if cost <= max_tt_copy and source in shared_nodes:
    #          output = [source, target, max_paths_copy]

    #     # all other paths, which should have acceptable length and not share node 
    #     if cost <= max_tt_copy and source not in shared_nodes:
    #         # calculate k shortest paths (k == max_paths_copy)
    #         s_paths=[]
    #         for s_path in k_shortest_paths(G_copy, source, target, max_paths_copy, weight=None):
    #             s_paths.append(s_path)
    #         #total number of acceptable paths var
    #         num_paths = max_paths_copy
    #         # in reverse order, find how many of the calculated paths have an acceptable time
    #         for s_path in reversed(s_paths):
    #             cost = nx.path_weight(G_copy, s_path, weight=G_weight_copy)
    #             if cost <= max_tt_copy:
    #                 output = [source, target, num_paths]
    #                 break
    #             else:
    #                 num_paths-=1
    #         output=[source, target,num_paths]
    
    # old k-shortest path method, native python function
    # # a source must have at least one acceptable path to be passed to the k-shortest path algorithm
    # # a path is acceptable if it is below the maximum allowable travel time 
    # if cost <= max_tt_copy and source in shared_nodes:
    #     output=[source, target, max_paths_copy]
    # elif cost <= max_tt_copy and source not in shared_nodes:
    #     length, path = nx.single_source_dijkstra(G_copy, source, target, weight=G_weight_copy)
    #     lengths = [length]
    #     paths = [path]
    #     c = count()
    #     B = []
    #     G_original = G_copy.copy()
    #     i = 1
    #     iter = True

    #     while iter is True:
    #         for j in range(len(paths[-1]) - 1):

    #             spur_node = paths[-1][j]
    #             root_path = paths[-1][:j + 1]

    #             edges_removed = []
    #             for c_path in paths:
    #                 if len(c_path) > j and root_path == c_path[:j + 1]:
    #                     u = c_path[j]
    #                     v = c_path[j + 1]
    #                     if G_copy.has_edge(u, v):
    #                         # edge_attr = G.edge[u][v]
    #                         edge_attr = G_copy.get_edge_data(u, v)
    #                         G_copy.remove_edge(u, v)
    #                         edges_removed.append((u, v, edge_attr))

    #             for n in range(len(root_path) - 1):
    #                 node = root_path[n]
    #                 # out-edges
    #                 # for u, v, edge_attr in G.edges_iter(node, data=True):
    #                 edgelist = list(G_copy.edges(node, data=True))
    #                 for u, v, edge_attr in edgelist:
    #                     G_copy.remove_edge(u, v)
    #                     edges_removed.append((u, v, edge_attr))

    #                 edgelist_in = list(G_copy.in_edges(node, data=True))
    #                 if G_copy.is_directed():
    #                     # in-edges
    #                     for u, v, edge_attr in edgelist_in:
    #                         # for u, v, edge_attr in G.edges(data=True):
    #                         G_copy.remove_edge(u, v)
    #                         edges_removed.append((u, v, edge_attr))

    #             try:
    #                 spur_path_length, spur_path = nx.single_source_dijkstra(
    #                     G_copy, spur_node, target, weight=G_weight_copy)
    #                 if target in spur_path:
    #                     total_path = root_path[:-1] + spur_path
    #                     total_path_length = get_path_length(
    #                         G_original, root_path, G_weight_copy) + spur_path_length
    #                     heappush(B, (total_path_length, next(c), total_path))
    #             except:
    #                 pass

    #             for e in edges_removed:
    #                 u, v, edge_attr = e
    #                 G_copy.add_edge(u, v, **edge_attr)

    #         if B:
    #             (l, _, p) = heappop(B)
    #             lengths.append(l)
    #             paths.append(p)
    #         else:
    #             break
    #         if l <= max_tt_copy:
    #             i += 1
    #         else:
    #             iter = False
    #         if len(lengths) >= max_paths_copy:
    #             iter = False
        
    #     output = [source, target, len(lengths)]
    # else:
    #     # still want to keep track of 0s
    #     output = [source, target, 0]
    # return(output)


###### HELPER FUNCTIONS CURRENTLY NOT BEING USED ######
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
    reachable, backnode = mod_search(
        source=source, flow=flow, numnodes=numnodes, capacity_array=capacity_array)
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
###### HELPER FUNCTIONS CURRENTLY NOT BEING USED ######


def traffic_assignment(G: nx.Graph,
                        res_points: gpd.GeoDataFrame, 
                        dest_points: gpd.GeoDataFrame,
                        dest_parcels: gpd.GeoDataFrame,  
                        G_capacity: str = 'capacity', 
                        G_weight: str = 'travel_time', 
                        algorithm: str = 'path_based', 
                        method: str = 'CFW',
                        link_performance: str = 'BPR',
                        termination_criteria: list[str, int] = ['iter',0],
                        dest_method: str = 'single',
                        sparse_array: bool = True):
    """
    Solves the static traffic assignment problem using algorithms and termination criteria from user inputs.

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
    :param dest_method: either 'single' or 'multiple' determines snapping methodology for destination parcels, either nearest node or multiple nearest nodes
    :type dest_method: string
    :param sparse_array: If True, uses sparse.csr_array in calculations, if False, use np.array in calculations
    :type sparse_array: bool
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


    # commonly used variables and their meaning:
        # weight_array:                     free flow time 
        # weight_array_iter (or _iter):     weight array for that specific itteration
        # OD_matrix:                        



    # Convert graph to digraph format
    G = nx.DiGraph(G)

    # TODO: remove this deepcopy?
    # G=copy.deepcopy(G)

    # Remove all edges with 0 capacity
    no_cap_edges = list(filter(lambda e: e[2] == 0, (e for e in G.edges.data(G_capacity))))
    no_cap_ids = list(e[:2] for e in no_cap_edges)
    G.remove_edges_from(no_cap_ids)

    # Travel times must be whole numbers -  round values up to nearest whole number
    for x in G.edges:
        G.edges[x][G_weight] = math.ceil(G.edges[x][G_weight])

    # Set variables that will be used as constants throughout algorithm - set keeping int32 byte criteria in mind
    super_sink = 99999999
    super_origin = 0
    artificial_weight = 999999
    artificial_weight_min=1
    artificial_capacity = 999999999

    # original node list - used to relate flow back to graph representation
    #TODO: Potential Bug: this is only og_node_list if nodes are renamed, which is done when network is loaded
    # og_node_list is used to relate 
    og_node_list = [num for num in range(0,len(G.nodes))]
    # PREPROCESSING STEPS: 
    # a.   Determine sources/sinks, 
    # b.   set aritificial components of network, 
    # c.   create OD Matrix
    # d.   create list of ordered nodes
    # e.   Convert capacity, weight, alpha, and beta arrays

    if dest_method == 'single':
        # determine if single destination file or mulitple
        if isinstance(dest_points, gpd.GeoDataFrame):
            # a. Determine sources and sinks
            G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points = nearest_nodes(G=G, 
                                                                                                                                res_points=res_points, 
                                                                                                                                dest_points=dest_points, 
                                                                                                                                G_demand='demand')

            # b_1. Add an artifical 0 node to maintain other nodes positions
            #   Doesn't actually need attributes since not connected, but filling for continuity sake
            G.add_node(super_origin, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

            # b_2. Add artifical edge from each destination node to the artifical sink with zero cost and maximum capacity
            # This is to allow minimum cost routing to whatever resource - destination of OD pair can change based on cost
            for idx, sink in enumerate(unique_dest_nodes):
                G.add_edge(sink, super_sink, **{G_weight: artificial_weight_min, G_capacity: artificial_capacity})
        
            # c. Create OD_matrix and add artifical edge from each origin to super_sink with extreme weight and capacity to capture loss of service 
            #    2 reasons to add origin->super_sink edge: 
            #   (1): Can build in elastic demand cost if people are priced out of accessign resource, or can set to arbitrarily high value to
            #   (2): By having artifical route, can always send max flow w/o considering if accessible - if path goes through artifical link, than no access when cost is arbitrarily high
            OD_matrix = np.empty([len(unique_origin_nodes),3])
            for idx, x in enumerate(unique_origin_nodes):
                OD_matrix[idx][0] = x                        # origin nodes
                OD_matrix[idx][1] = G.nodes[x]['demand']*-1  # flow associated with the origin node
                OD_matrix[idx][2] = len(G.nodes())-1         # destination node (when nodes reordered, it has a value of length-1 (b/c of artifical origin 0 ))

                # need to add edges from each origin to super_sink with extreme weight and capacity to capture loss of service
                # HOWEVER, issue arrises if shared origin-destination, therefore not adding this artifical edge of x is also in unique_dest_nodes
                if x not in unique_dest_nodes:
                    G.add_edge(x, super_sink, **{G_weight: artificial_weight, G_capacity: artificial_capacity})
                    
            # d. Sort nodes in proper order
            #   This should in theory then preserve saying node '2' in array is node '2' in network g, especially with adding the artifical node, 0
            #   also rename the nodes in graph
            G = nx.convert_node_labels_to_integers(G, 0, ordering="sorted", label_attribute='old_label')
            nodes_in_order = sorted(G.nodes())

        if isinstance(dest_points, list):
            # a. Determine sources and sinks
            G, unique_origin_nodes, unique_dest_nodes, positive_demand, shared_nodes, res_points, dest_points = nearest_nodes(G=G,
                                                                                                                              res_points=res_points,
                                                                                                                              dest_points=dest_points,
                                                                                                                              G_demand='demand')

            # b_1. Add an artifical 0 node to maintain other nodes positions
            #   Doesn't actually need attributes since not connected, but filling for continuity sake
            G.add_node(super_origin, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

            # b_2. Add artifical edge from each destination node to the artifical sink with zero cost and maximum capacity
            # This is to allow minimum cost routing to whatever resource - destination of OD pair can change based on cost
            # need a len variable before adding anymore nodes
            temp_len = len(G.nodes())
            for idx, sinks in enumerate(unique_dest_nodes):
                for sink in sinks:
                    G.add_edge(sink, super_sink+idx, **{G_weight: artificial_weight_min, G_capacity: artificial_capacity})

            # c. Create OD_matrix and add artifical edge from each origin to super_sink with extreme weight and capacity to capture loss of service 
            #    2 reasons to add origin->super_sink edge: 
            #   (1): Can build in elastic demand cost if people are priced out of accessign resource, or can set to arbitrarily high value to
            #   (2): By having artifical route, can always send max flow w/o considering if accessible - if path goes through artifical link, than no access when cost is arbitrarily high
            OD_matrix = np.empty([len(unique_origin_nodes)*len(unique_dest_nodes),3])

            for i, resource in enumerate(unique_dest_nodes):
                for idx, x in enumerate(unique_origin_nodes):
                    dest_node_val = temp_len+(i)                  # similar to single dest file, except accounting for multiple super sinks
                    OD_matrix[idx+i*len(unique_origin_nodes)][0] = x                           # origin nodes
                    OD_matrix[idx+i*len(unique_origin_nodes)][1] = G.nodes[x][f'demand{i}']*-1 # flow associated with the origin node
                    OD_matrix[idx+i*len(unique_origin_nodes)][2] = dest_node_val         

                    # need to add edges from each origin to super_sink with extreme weight and capacity to capture loss of service
                    # HOWEVER, issue arrises if shared origin-destination, therefore not adding this artifical edge of x is also in unique_dest_nodes
                    if x not in unique_dest_nodes[i]:
                        G.add_edge(x, super_sink+i, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

            # d. Sort nodes in proper order
            #   This should in theory then preserve saying node '2' in array is node '2' in network g, especially with adding the artifical node, 0
            #   also rename the nodes in graph
            G = nx.convert_node_labels_to_integers(G, 0, ordering="sorted", label_attribute='old_label')
            nodes_in_order = sorted(G.nodes())


    elif dest_method == 'multiple':
        # determine if single destination file or multiple
        if isinstance(dest_points,gpd.GeoDataFrame):
            # a. Determine sources and sinks
            G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points = nearest_nodes_vertices(G=G, res_points=res_points, dest_parcels=dest_parcels, dest_points=dest_points, G_demand='demand')

            # b_1. Add an artifical 0 node to maintain other nodes positions
            #   Doesn't actually need attributes since not connected, but filling for continuity sake
            G.add_node(super_origin, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

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
                        kwargs = {G_weight: artificial_weight_min, G_capacity: artificial_capacity}
                        G.add_edge(nearest_intersection, dest_node_ids, **kwargs)

                    # add link from parcel centroid to super sink
                    # TODO2: This is where I can specifically add capacities for each specific grocery store
                    # FOR NOW: just setting capacity to a high amount 
                    kwargs = {G_weight: artificial_weight_min, G_capacity: artificial_capacity}
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
                    OD_matrix[idx][0] = x                        # origin nodes
                    OD_matrix[idx][1] = G.nodes[x]['demand']*-1  # flow associated with the origin node
                    OD_matrix[idx][2] = len(G.nodes())-1         # destination node, when nodes reordered, it has a value of length-1 (b/c of artifical origin 0 ))
                    
                    # need to add edges from each origin to super_sink with extreme weight and capacity to acpture loss of service
                    # HOWEVER, issue arrises if shared origin-destination, therefore not adding this artifical edge of x is also in inque_dest_nodes
                    if x not in shared_nodes:
                        G.add_edge(x, super_sink, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

            # delete OD_matrix rows that we no longer need
            indexList = np.where(~OD_matrix.any(axis=1))[0]
            OD_matrix = np.delete(OD_matrix, indexList, axis=0)

            # d. Sort nodes in proper order
            #   This should in theory then preserve saying node '2' in array is node '2' in network g, especially with adding the artifical node, 0
            #   also rename the nodes in graph
            G = nx.convert_node_labels_to_integers(G, 0, ordering="sorted", label_attribute='old_label')
            nodes_in_order = sorted(G.nodes())

        if isinstance(dest_points, list):
            # a. Determine sources and sinks
            G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points = nearest_nodes_vertices(G=G, res_points=res_points, dest_parcels=dest_parcels, dest_points=dest_points, G_demand='demand')

            # b_1. Add an artifical 0 node to maintain other nodes positions
            #   Doesn't actually need attributes since not connected, but filling for continuity sake
            G.add_node(super_origin, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

            # b_2. Add artifical edges from vertices around destinations to aggregate points than to artifical sink with 0 cost and max capacity
            # This is to allow minimum cost routing to whatever resource - destination of OD pair which can change based on cost
            # identify the destination nodes, and create artifical sink edges
            # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to 
            # artificial node id tracker - useful in maintianing compatability with dtype
            dest_node_ids=99999998
            # need a len variable before adding anymore nodes - this is for offsetting dest node in OD-matrix
            temp_len = len(G.nodes())
            for i, parcels in enumerate(dest_parcels):
                for idx, dest_parcel in parcels.iterrows():
                    dest_node = dest_points[i][dest_points[i].geometry.within(dest_parcel.geometry)]
                    #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
                    for a, node in dest_node.iterrows():
                        dest_node_ids-=1
                        # add the dest node to the graph 
                        G.add_nodes_from([(dest_node_ids, {'demand': 0})])
                        temp_len+=1
                        # add links from nearest intersections to parcel centroid
                        for nearest_intersection in unique_dest_nodes_list[i][idx]:
                            kwargs = {G_weight: artificial_weight_min, G_capacity: artificial_capacity}
                            G.add_edge(nearest_intersection, dest_node_ids, **kwargs)

                        # add link from parcel centroid to super sink
                        # TODO2: This is where I can specifically add capacities for each specific resource
                        # FOR NOW: just setting capacity to a high amount 
                        kwargs = {G_weight: artificial_weight_min, G_capacity: artificial_capacity}
                        G.add_edge(dest_node_ids, super_sink+i, **kwargs)

            # c. Create OD_matrix and add artifical edge from each origin to super_sink with extreme weight and capacity to capture loss of service 
            #    2 reasons to add origin->super_sink edge: 
            #   (1): Can build in elastic demand cost if people are priced out of accessign resource, or can set to arbitrarily high value to
            #   (2): By having artifical route, can always send max flow w/o considering if accessible - if path goes through artifical link, than no access when cost is arbitrarily high
            # initialize OD_matrix        
            OD_matrix = np.empty([len(unique_origin_nodes)*len(unique_dest_nodes_list),3])
            
            indexList=[]
            for i, resource in enumerate(unique_dest_nodes_list):
                for idx,x in enumerate(unique_origin_nodes):
                    # BUG: before adding to OD matrix, determine if a path even exists
                    if nx.has_path(G, x, super_sink+1) is False:
                        # save index to be removed
                        indexList.append(idx+i*len(unique_origin_nodes))
                        pass
                    else:
                        dest_node_val = temp_len+(i)    # similar to single dest file, except accounting for multiple super sinks
                        OD_matrix[idx+i*len(unique_origin_nodes)][0] = x                        # origin nodes
                        OD_matrix[idx+i*len(unique_origin_nodes)][1] = G.nodes[x][f'demand{i}']*-1  # flow associated with the origin node
                        OD_matrix[idx+i*len(unique_origin_nodes)][2] = dest_node_val         
                        
                        # need to add edges from each origin to super_sink with extreme weight and capacity to capture loss of service
                        # HOWEVER, issue arrises if shared origin-destination, therefore not adding this artifical edge of x is also in inque_dest_nodes
                        if x not in shared_nodes[i]:
                            G.add_edge(x, super_sink+i, **{G_weight: artificial_weight, G_capacity: artificial_capacity})

            # delete OD_matrix rows that we no longer need - where values are 0 or nans
            # indexList = np.where(~OD_matrix.any(axis=1))[0]
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

    # convert capacities and weights to numpy arrays
    if sparse_array is True:
        capacity_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight=G_capacity, nonedge=0, dtype=np.float64)    
        weight_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight=G_weight, nonedge=0, dtype=np.float64)
        # create adjacency list - faster for some functions to use this instead of matrix
        adjlist = defaultdict(list)
        row, col = np.where(capacity_array != 0)
        for i in range(len(row)):
            adjlist[row[i]].append(col[i])
    else:
        capacity_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight=G_capacity, nonedge=-1, dtype=np.float64)    
        weight_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight=G_weight, nonedge=-1, dtype=np.float64)
        # create adjacency list - faster for some functions to use this instead of matrix
        adjlist = defaultdict(list)
        row, col = np.where(capacity_array!=-1)
        for i in range(len(row)):
            adjlist[row[i]].append(col[i])        
    
    # Set BPR alpha and beta constants -> Currently not using array values for this - too time consuming
    nx.set_edge_attributes(G, values=0.15, name='alpha')
    # nx.set_edge_attributes(G, values=4, name='beta')
    # alpha_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='alpha', nonedge=-1) 
    # beta_array = nx.to_numpy_array(G, nodelist=nodes_in_order, weight='beta', nonedge=-1)
    alpha_array = 0.15
    beta_array = 4

    # array to create empty copies later on
    copy_array = np.zeros_like(capacity_array, dtype=np.float32)

    if sparse_array is True:    
        # convert weight and capacity arrays to sparse arrays in list-of-list format
        weight_array = sparse.lil_array(weight_array, dtype=np.float32)
        capacity_array = sparse.lil_array(capacity_array, dtype=np.float32)
        # need an array of same shape of capacity of just ones for link performance function function
        temp_cap = capacity_array.tocsr()
        addition = sparse.csr_array(([1 for i in temp_cap.data], temp_cap.indices, temp_cap.indptr))
        temp_cap = None
        # TODO: COULD also precompute capacity**beta_array, which is used in link_performance_function_derivative
    
    # set variables
    sum_d = positive_demand                             # sum of the demand
    num_nodes = capacity_array.shape[0]                  # number of nodes
    node_list = [num for num in range(0, num_nodes)]    # list of all nodes
    pp=8                                                # set parallel processing variables
    
    # increase the factor of the weight array
    weight_array=weight_array.multiply(100)
    ####################################################################################################################

    # HELPER FUNCTIONS 
    # link_performance_function and link_performance_derivative require link performance equations
        # just a note if future link performance functions need to be added

    def link_performance_function(capacity_array,
                                flow_array,
                                weight_array,
                                eq,
                                sparse_array,
                                alpha_array=0.15,
                                beta_array=4,
                                addition=addition,
                                update_arr = None,
                                update_links = None):
        """ 
        Function returning the weight_array_iter array based on user selected equation.

        If update_arr is not None, this is a weight_array_iter that only has a few edges that need to be updated and the entire array doesn't need to be calculated
        
        Universival arguments:
            capacity_array == array of link capacities (either node-node np.array or sparse.csr_array)
            flow_array     == array of link flows (either node-node np.array or sparse.csr_array)
            weight_array   == array of free-flow weights (i.e., zero flow travel time, either node-node np.array or sparse.csr_array)
            sparse_array   == Bool, if True, inputs are sparse.csr_array, if False, inputs are np.array

        BPR arguments:
            alpha_array == BPR alpha value (currently only single value, potential for node-node array in future)
            beta_array  == BPR beta value (currently only single value, potential for node-node array in future)
            addition    == an array that hasa value of 1 where there is a link (i.e., capacity > 0), used in BPR calcualtion

        other arguments:
            update_arr   == If not None, array that needs to be updated based on just a few links 
            update_links == The links that need to be updated 

        returns: weight_array_iter
        """

        if eq == "BPR":
            if sparse_array is True:
                # Using scipy.sparse functions much faster than using just *, /, and **
                # I believe this is because these operators convert to numpy arrays and the conversions are slow
                # limited documentation on how scipy.sparse methods work
                # convert arrays to csr arrays for faster math calculations

                if update_arr is not None:
                    for unique_links in update_links:
                        for link in unique_links:
                            # update_arr[link[0], link[1]] = np.round(weight_array[link[0], link[1]] * (
                            #     1+alpha_array*(flow_array[link[0], link[1]]/capacity_array[link[0], link[1]])**beta_array),0)
                            update_arr[link[0], link[1]] = weight_array[link[0], link[1]] * (
                                1+alpha_array*(flow_array[link[0], link[1]]/capacity_array[link[0], link[1]])**beta_array)
                    return update_arr

                else:
                    weight_array = sparse.csr_array(weight_array)
                    flow_array = sparse.csr_array(flow_array)
                    capacity_array = sparse.csr_array(capacity_array)
                    
                    # calculate BPR function
                    weight_array_iter = weight_array.multiply(flow_array.multiply(capacity_array.power(-1)).power(beta_array).multiply(alpha_array)+addition)
                    weight_array_iter = weight_array_iter.tolil()
                    # BUG: need to add rounding to be able to converge???
                    # weight_array_iter = sparse.lil_matrix(np.round(weight_array_iter.todense(), decimals=0))
                
                # old method of calculating BPR function - easier to read but much slower
                # weight_array_iter = weight_array * (1+alpha_array*(flow_array/capacity_array)**beta_array)
                # weight_array_iter = sparse.lil_array(weight_array_iter)
                
            # If not using sparse_arrays:
            else:
                weight_array_iter = weight_array * \
                    (1+alpha_array*(flow_array/capacity_array)**beta_array)
        # elif eq == "Davidsons":
            # pass

        return weight_array_iter

    def link_performance_derivative(capacity_array,
                                    flow_array,
                                    weight_array,
                                    eq,
                                    sparse_array,
                                    alpha_array=0.15,
                                    beta_array=4,
                                    update_arr=None,
                                    update_links=None):
        """ 
        Function returning the derivative array based on user selected equation
        
        Universival arguments:
            capacity_array == array of link capacities (either node-node np.array or sparse.csr_array)
            flow_array     == array of link flows (either node-node np.array or sparse.csr_array)
            weight_array   == array of free-flow weights (i.e., zero flow travel time, either node-node np.array or sparse.csr_array)
            sparse_array   == Bool, if True, inputs are sparse.csr_array, if False, inputs are np.array

        BPR arguments:
            alpha_array == BPR alpha value (currently only single value, potential for node-node array in future)
            beta_array  == BPR beta value (currently only single value, potential for node-node array in future) 

        other arguments:
            update_arr   == If not None, array that needs to be updated based on just a few links 
            update_links == The links that need to be updated 
            
        returns: link_performance_derivative
        """

        if eq == "BPR":
            if sparse_array is True:

                if update_arr is not None:
                    for link in update_links:
                        # update_arr[link[0], link[1]] = np.round((
                        #     beta_array*weight_array[link[0], link[1]]*alpha_array/capacity_array[link[0], link[1]]**beta_array)*flow_array[link[0], link[1]]**(beta_array-1),0)
                        update_arr[link[0], link[1]] = np.round((
                            beta_array*weight_array[link[0], link[1]]*alpha_array/capacity_array[link[0], link[1]]**beta_array)*flow_array[link[0], link[1]]**(beta_array-1),0)
                    return update_arr
                
                else:
                    # convert to csr_arrays for faster calculations
                    weight_array = sparse.csr_array(weight_array)
                    flow_array = sparse.csr_array(flow_array)
                    capacity_array = sparse.csr_array(capacity_array)
                    
                    # calculate BPR derivative
                    link_performance_derivative = flow_array.power(
                        beta_array-1).multiply(capacity_array.power(-1).power(beta_array).multiply(alpha_array).multiply(weight_array).multiply(beta_array))

                    # old method of calculating BPR function - easier to read but much slower
                    # link_performance_derivative = (beta_array*weight_array*alpha_array/capacity_array**beta_array)*flow_array**(beta_array-1)
                    
                    link_performance_derivative = link_performance_derivative.tolil()
                    # BUG: need to add rounding to be able to converge???
                    # link_performance_derivative = sparse.lil_matrix(np.round(link_performance_derivative.todense(), decimals=0))
            else:
                link_performance_derivative = (
                    beta_array*weight_array*alpha_array/capacity_array**beta_array)*flow_array**(beta_array-1)
        # elif eq == "Davidsons":
            # pass



        return link_performance_derivative

    def bisection_zeta(lambda_val, 
                       flow_array_star, 
                       capacity_array, 
                       flow_array, 
                       weight_array, 
                       eq=link_performance, 
                       alpha_array=0.15, 
                       beta_array=4,
                       sparse_array=sparse_array):
        """ 
        Function returning the zeta value when using the bisection link-based method
        
        Universival arguments:
            lambda_val      == lambda value to shift flow
            flow_array_star == array of updated link flows (i.e., zero flow travel time, either node-node np.array or sparse.csr_array)
            capacity_array  == array of link capacities (either node-node np.array or sparse.csr_array)
            flow_array      == array of link flows (either node-node np.array or sparse.csr_array)
            weight_array    == array of free-flow weights (i.e., zero flow travel time, either node-node np.array or sparse.csr_array)
            sparse_array    == Bool, if True, inputs are sparse.csr_array, if False, inputs are np.array

        BPR arguments:
            alpha_array == BPR alpha value (currently only single value, potential for node-node array in future)
            beta_array  == BPR beta value (currently only single value, potential for node-node array in future)     

        returns: bisection zeta value
        """
        x_hat = lambda_val*flow_array_star+(1-lambda_val)*flow_array
        
        x_hat_link_performance = link_performance_function(capacity_array=capacity_array, 
                                                           flow_array=x_hat,
                                                           weight_array=weight_array,
                                                            eq=eq, 
                                                            alpha_array=alpha_array, 
                                                            beta_array=beta_array,
                                                            sparse_array=sparse_array)

        zeta = np.sum(x_hat_link_performance*(flow_array_star-flow_array))

        return(zeta)
        
    def newton_f(lambda_val, 
                 flow_array_star, 
                 capacity_array, 
                 flow_array, 
                 weight_array, 
                 eq=link_performance, 
                 alpha_array=0.15, 
                 beta_array=4,
                 sparse_array=sparse_array):
        """ Function returning the f and f_prime values when using the newton's and CFW link-based methods
        
        Universival arguments:
            lambda_val      == lambda value to shift flow
            flow_array_star == array of updated link flows (i.e., zero flow travel time, either node-node np.array or sparse.csr_array)
            capacity_array  == array of link capacities (either node-node np.array or sparse.csr_array)
            flow_array      == array of link flows (either node-node np.array or sparse.csr_array)
            weight_array    == array of free-flow weights (i.e., zero flow travel time, either node-node np.array or sparse.csr_array)
            sparse_array    == Bool, if True, inputs are sparse.csr_array, if False, inputs are np.array

        BPR arguments:
            alpha_array == BPR alpha value (currently only single value, potential for node-node array in future)
            beta_array  == BPR beta value (currently only single value, potential for node-node array in future)     

        returns: f and f_prime
        """
        x_hat = lambda_val*flow_array_star+(1-lambda_val)*flow_array

        x_hat_link_performance = link_performance_function(capacity_array=capacity_array,
                                                           flow_array=x_hat,
                                                           weight_array=weight_array,
                                                           eq=eq,
                                                           alpha_array=alpha_array,
                                                           beta_array=beta_array,
                                                           sparse_array=sparse_array)

        x_hat_link_performance_derivative = link_performance_derivative(capacity_array=capacity_array, 
                                                                        flow_array=x_hat,
                                                                        weight_array=weight_array, 
                                                                        eq=eq, 
                                                                        alpha_array=alpha_array, 
                                                                        beta_array=beta_array,
                                                                        sparse_array=sparse_array)

        f = np.sum(x_hat_link_performance*(flow_array_star-flow_array))
        f_prime = np.sum(x_hat_link_performance_derivative*(flow_array_star-flow_array)**2)

        return f, f_prime

    def cfw_alpha(flow_array_star_old, 
                  flow_array_star, 
                  capacity_array, 
                  flow_array, 
                  weight_array, 
                  eq=link_performance, 
                  alpha_array=0.15, 
                  beta_array=4,
                  sparse_array=sparse_array):
        """ Function returning the numerator and denominator of alpha value for CFW link-based method

        Universival arguments:
            flow_array_star_old == ARRAY OF previously updated flow array (either node-node np.array or sparse.csr_array)
            flow_array_star     == array of updated link flows (i.e., zero flow travel time, either node-node np.array or sparse.csr_array)
            capacity_array      == array of link capacities (either node-node np.array or sparse.csr_array)
            flow_array          == array of link flows (either node-node np.array or sparse.csr_array)
            weight_array        == array of free-flow weights (i.e., zero flow travel time, either node-node np.array or sparse.csr_array)
            sparse_array        == Bool, if True, inputs are sparse.csr_array, if False, inputs are np.array

        BPR arguments:
            alpha_array == BPR alpha value (currently only single value, potential for node-node array in future)
            beta_array  == BPR beta value (currently only single value, potential for node-node array in future)       

        returns: f and f_prime
        """

        deriv = link_performance_derivative(capacity_array=capacity_array, 
                                            flow_array=flow_array, 
                                            weight_array=weight_array, 
                                            eq=eq, 
                                            alpha_array=alpha_array, 
                                            beta_array=beta_array,
                                            sparse_array=sparse_array)

        #eq (6.3) in boyles textbook page 176
        item1 = flow_array_star_old - flow_array
        item2 = flow_array_star - flow_array
        item3 = flow_array_star - flow_array_star_old
        

        num = np.sum(deriv*item1*item2)
        den = np.sum(deriv*item1*item3)
        return num, den

    def topo_order_function(bush, 
                            origin, 
                            weight_array, 
                            num_nodes = num_nodes, 
                            node_list=node_list,
                            sparse_array=sparse_array):
        """
        Function that calculates topographic order on a bush.
        """

        # Find shortest path from origin to all locations - need to determine what nodes are unreachable 
        backnode, costlabel = shortestPath_heap(origin=origin, 
                                                capacity_array=capacity_array, 
                                                weight_array=weight_array, 
                                                adjlist=adjlist,
                                                destination=False,
                                                sparse_array=sparse_array)

        # where backnode == -1, unreachable, and won't be considered in topological ordering 
        unreachable = [i for i,j in enumerate(backnode) if j == -1]
        unreachable.pop(unreachable.index(origin))
        # creat visited labels 
        visited = [False]* (num_nodes)
        # set the unreachable nodes and origin nodes to True in order to skip them
        for val in unreachable:
            visited[val]=True
        
        # create adjacency table of bush array
        adj_table = {}
        for node in node_list:
            if visited[node] is True:
                pass
            else:
                if sparse_array is True:
                    adj_table[node]=[i for i, j in enumerate(bush.A[[node], :][0]) if j > 0]
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

        topo_order = topo_order[::-1]
        return topo_order

    def termination_function(termination_criteria, 
                             iters, 
                             AEC, 
                             RG):
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

    ########################################################################################################################
    # ALGORITHMS
    if algorithm == "link_based":

        # INITIALIZE LINK BASED ALGORITHM
        # a.    Create empty flow array - same size as the capacity_array (b/c node-node relationship)
        # b.    Create initial feasible flow array with all-or-nothing minimum cost flow assignment (free flow travel cost)

        # a. Create empty flow array 
        flow_array = np.zeros_like(copy_array)

        # BEGIN ORIGINAL NON PARALLEL VERSION 
        # b. Create initial feasible flow array
        # for x in OD_matrix:
        #     origin = int(x[0])
        #     flow = x[1]
        #     destination = int(x[2])
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
        # n = [origin, flow, destination, what to return, which weight array to use]
        # I think OD_matrix can be reformated to start in this format, and we can therefore save time 
        n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]), 'path', 'weight_array'] for n in range(0, len(OD_matrix))]
        # initialize the pool
        pool=Poolm(pp)
        # map function
        results = pool.map(SPTT_parallel, n)
        # close the pool
        pool.close()
        pool.join()
        # edit the flow_array matrix
        for idx, result in enumerate(results):
            for i in result:
                flow_array[i[0]][i[1]] += OD_matrix[idx][1]     
        # END PARALLEL PROCESSING 
        
        # convert flow_array to sparse_array if necessary
        if sparse_array is True:
            flow_array = sparse.lil_array(flow_array)

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
            # 1. recalculate weight_array
            weight_array_iter = link_performance_function(capacity_array=capacity_array,
                                                          flow_array=flow_array,
                                                          weight_array=weight_array,
                                                          eq=link_performance,
                                                          alpha_array=alpha_array,
                                                          beta_array=beta_array,
                                                          sparse_array=sparse_array)

            # 2. Create empty kappa array and flow_array_star
            # TODO: I DON'T THINK KAPPA IS ACTUALLY USED
            kappa = np.zeros([1, np.shape(OD_matrix)[0]])
            if sparse_array is True:
                flow_array_star = sparse.lil_array(capacity_array.shape)
            else:
                flow_array_star = np.zeros_like(capacity_array)
            # For shortest path travel time calculation (termination criteria)
            SPTT=0    

            # #BEGIN ORIGINAL CALC
            # for idx, x in enumerate(OD_matrix):
            #     origin = int(x[0])
            #     destination = int(x[2])
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
            #     # 5. Update the flow array
            #     for i in path:
            #         flow_array_star[i[0]][i[1]] += flow
            # END ORIGINAL CALC
           
            # BEGIN PARALLEL CALC
            # n = [origin, flow, destination, what to return, which weight array to use]
            n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]), 'path_cost', 'weight_array_iter'] for n in range(0, len(OD_matrix))]
            pool=Poolm(pp)
            results = pool.map(SPTT_parallel, n)
            pool.close()
            pool.join()
            for idx, result in enumerate(results):
                path, cost = result
                flow = OD_matrix[idx][1]
                SPTT += cost*flow
                if sparse_array is True:
                    for i in path:
                        flow_array_star[i[0], i[1]] += flow
                else:
                    for i in path:
                        flow_array_star[i[0]][i[1]] += flow
            
            # END PARALLEL CALC

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
                                      beta_array=beta_array,
                                      sparse_array=sparse_array)


                # Update Lambda
                #TODO: see PPT, add criteria about hitting the same endpoint twice in a row to terminate
                lambda_val -= f/f_prime
                # maximum lambda_value -> shift all flow
                if lambda_val > 1:
                    lambda_val = 1
                # minimum lambda_value -> shift none of the flow
                elif lambda_val < 0:
                    lambda_val = 0
                # if f_prime is 0 -> shift none of the flow
                else:
                    lambda_val=0

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
                                            beta_array=beta_array,
                                            sparse_array=sparse_array)
                    if den == 0:
                        alpha=0
                    else:
                        alpha = num/den
                        if alpha >= 1:
                            alpha = 1 - epsilon
                        if alpha < 0:
                            alpha = 0
                    if sparse_array is True:
                        flow_array_star = sparse.lil_array(flow_array_star_old.multiply(alpha) + flow_array_star.multiply(1-alpha))
                    else:
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
                                        beta_array=beta_array,
                                        sparse_array=sparse_array)
                
                # Update Lambda
                lambda_val -= f/f_prime
                # maximum lambda_value -> shift all flow
                if lambda_val > 1:
                    lambda_val = 1
                # minimum lambda_value -> shift none of the flow
                elif lambda_val < 0:
                    lambda_val = 0
                # if f_prime is 0 -> shift none of the flow
                else:
                    lambda_val=0

                # flow_array_star_old value for the next iteration
                if sparse_array is True:
                    flow_array_star_old = sparse.lil_array.copy(flow_array_star)
                else:
                    flow_array_star_old = np.copy(flow_array_star)
                # set marker to pass for subsequent interations
                CFW_marker = 1

            # 7. Shift flows
            if sparse_array is True:
                flow_array = sparse.lil_array(flow_array_star.multiply(lambda_val) + flow_array.multiply(1-lambda_val))
            else:
                flow_array = lambda_val*flow_array_star + (1-lambda_val)*flow_array 

            # Fill termination variable lists
            print(SPTT, TSTT, AEC, RG)
            AEC_list.append(AEC)
            TSTT_list.append(TSTT)
            SPTT_list.append(SPTT)
            RG_list.append(RG)
            iter += 1

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
        if sparse_array is True:
            flow_array = sparse.lil_array(np.shape(copy_array))
        else:
            flow_array = np.zeros_like(copy_array)
        # Initialize nested list to store OD path information 
            # The paths_array matrix is set up as follows:
            # paths_array[x][0]: [O, D]
            # paths_array[x][1]: [list of arrays that are paths for OD-pair]
            # paths_array[x][2]: [list of flows (associated with paths)]
        paths_array = []

        # populate paths_array matrix
        # BEGIN PARALLEL PROCESS
        # n = [origin, flow, destination, what to return, which weight array to use]
        n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]), 'path', 'weight_array',capacity_array, weight_array] for n in range(0, len(OD_matrix))]
        # n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]), 'path', 'weight_array'] for n in range(0, len(OD_matrix))]
        pool=Poolm(pp)
        results = pool.map(SPTT_parallel, n)
        pool.close()
        pool.join()
        for idx, result in enumerate(results):
            path = result
            flow = OD_matrix[idx][1]
            origin = OD_matrix[idx][0].astype(int)
            destination = OD_matrix[idx][2].astype(int)
            if sparse_array is True:
                for i in path:
                    flow_array[i[0], i[1]] += flow
            else:
                for i in path:
                    flow_array[i[0]][i[1]] += flow
            # append the paths_array
            paths_array.append([[origin,destination],[path],[flow]])
        # END PARALLEL PROCESS

        # Initiate the weight_array_iter
        weight_array_iter = link_performance_function(capacity_array=capacity_array,
                                                        flow_array=flow_array,
                                                        weight_array=weight_array,
                                                        eq=link_performance,
                                                        alpha_array=alpha_array,
                                                        beta_array=beta_array,
                                                        sparse_array=sparse_array)
        # weight_array_iter=weight_array_iter.astype('float32')

        # calculate initial termination criteria variables
        SPTT=0
        # BEGIN PARALLEL PROCESS
        # n = [origin, flow, destination, what to return, which weight array to use]
        n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]), 'cost', 'weight_array_iter',capacity_array, weight_array_iter] for n in range(0, len(OD_matrix))]

        pool=Poolm(pp)
        results = pool.map(SPTT_parallel, n)
        pool.close()
        pool.join()
        for idx, result in enumerate(results):
            cost = result
            flow = OD_matrix[idx][1]
            SPTT += cost*flow
        # END PARALLEL PROCESS

        # termination criteria
        TSTT = np.sum(np.multiply(flow_array, weight_array_iter))
        AEC = (TSTT-SPTT)/sum_d
        RG = TSTT/SPTT - 1

        # append the termination criteria lists
        AEC_list.append(AEC)
        TSTT_list.append(TSTT)
        SPTT_list.append(SPTT)
        RG_list.append(RG)
 
        # initialize a derivate flow array
        derivative = link_performance_derivative(capacity_array=capacity_array, 
                                                                flow_array=flow_array, 
                                                                weight_array=weight_array,
                                                                eq=link_performance, 
                                                                alpha_array=alpha_array, 
                                                                beta_array=beta_array,
                                                                sparse_array=sparse_array)

        # set iter values
        iter = 0
        iter_val = True
        # begin loop
        # have to initative adj_links_d to maintain track of adjacency links used in derivative
        adj_links_d=[] 
        
        while iter_val is True:
            print(f'iter: {iter}')
            timea=time.time()
            # Because flow shifts with each OD_pair, we cannot precalculate in parallel all of the shortest paths
            for OD_pair in paths_array:
                origin = OD_pair[0][0]
                destination = OD_pair[0][1]
                paths = OD_pair[1]
                flows = OD_pair[2]

                # calculate the current cost of known paths to set as limit of new shortest paths
                new_limits=[np.inf]
                for known_path in paths:
                    value=0
                    for edges in known_path:
                        value += weight_array_iter[edges[0],edges[1]]
                    new_limits.append(value)
                limit = np.min(new_limits)
                # limit=np.inf

                # # Find the new shortest path 
                backnode, costlabel = shortestPath_heap(origin=origin, 
                                                        capacity_array=capacity_array, 
                                                        weight_array=weight_array_iter, 
                                                        adjlist=adjlist,
                                                        sparse_array=sparse_array,
                                                        destination=destination,
                                                        limit=limit)
                
                
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

                # perform functions on each path except the shortest path (in the last position of the list)
                # adj_marker keeps track of if flows were shifted for that OD pair
                # adj_links are the links that were changed and need to be updated
                adj_marker = False
                adj_links=[]  # specifically for weight_array_iter
                # adj_links_d=[] # specifically for derivative

                for idx, path in enumerate(paths[:-1]):
                    # Determine unique_links between shortest path (path_hat) and current path under investigation
                    # Previous incorrect method
                    # all_links = [item for sublst in zip(path, path_hat) for item in sublst]
                    # unique_links = [list(x) for x in set(tuple(x) for x in all_links)]
                    set1 = set(tuple(x) for x in path)
                    set2 = set(tuple(x) for x in path_hat)
                    unique_links=[list(x) for x in (set1 - set2)]+[list(x) for x in (set2 - set1)]
                    
                    # calculate numerator of delta_h
                    if sparse_array is True:
                        cost = sum([weight_array_iter[x[0],x[1]] for x in path])
                    else:
                        cost = sum([weight_array_iter[x[0]][x[1]] for x in path])
                    num = cost-cost_hat
                    
                    # determine the sum of derivatives (denominator of delta_h equation)
                    den = 0

                    # Update link_performance_derivative
                    # BUG : where does adj_links_d need to be initialized?
                    # Only need to recalculate edges if there have been adjustments made to edges
                    if len(adj_links_d) >= 1:
                        derivative = link_performance_derivative(capacity_array=capacity_array, 
                                                                    flow_array=flow_array, 
                                                                    weight_array=weight_array,
                                                                    eq=link_performance, 
                                                                    alpha_array=alpha_array, 
                                                                    beta_array=beta_array,
                                                                    sparse_array=sparse_array,
                                                                    update_arr=derivative,
                                                                    update_links=adj_links_d)
                        adj_links_d=[]
                    
                    if sparse_array is True:
                        den=sum([derivative[x[0],x[1]] for x in unique_links])
                    else:
                        den=sum([derivative[x[0]][x[1]] for x in unique_links])
                    
                    # calculate delta_h
                    delta_h = num/den
                    if delta_h < 0:
                        # There should never be a negative delta_h value, if this is printed there is some sort of error that needs to be ID'd
                        print(f'NEGATIVE delta_h: {num}/{den}={delta_h}')
                     
                    # adjusted path flow value
                    adj = min(flows[idx], delta_h)
                    # There should never be a negative adj (i.e., projected delta_h). if this is printed there is some sort of error that needs to be ID'd
                    if adj < 0:
                        print(f'NEGATIVE adj! Delta_h:{num}/{den}={delta_h}, flow:{flows[idx]}')
                        print(unique_links)
                        print(path)
                        print(path_hat)
                    if adj > 0:
                        # only need to update weight_array_iter if values actually changed
                        adj_marker=True
                        adj_links.append(unique_links)
                        adj_links_d=copy.deepcopy(unique_links)

                    # update paths_array flow values
                    flows[idx] -= adj
                    flows[-1] += adj
                    # for the path under examination and path_hat, need to update flow array to calc new weight array
                    if sparse_array is True:
                        for i in path:
                            flow_array[i[0],i[1]] -= adj
                        for i in path_hat:
                            flow_array[i[0],i[1]] += adj
                        if flow_array[i[0],i[1]] < 0:
                            # there should never be a negative flow, if this prints there is an error that needs to be ID'd
                            print('NEGATIVE FLOW')
                    else:
                        for i in path:
                            flow_array[i[0]][i[1]] -= adj
                        for i in path_hat:
                            flow_array[i[0]][i[1]] += adj
                    
                # Recalculate weight_array for next itteration
                # only need to calculate a new weight array 2 conditions are met:
                # 1. we have multiple paths (captured with the adj_marker)
                # 2. there was a positive change to a flow (captured with the adj_marker)
                # And we only need to update the weight_array_iter links that experienced a change
                if adj_marker is True:
                    weight_array_iter = link_performance_function(capacity_array=capacity_array,
                                                                    flow_array=flow_array,
                                                                    weight_array=weight_array,
                                                                    eq=link_performance,
                                                                    alpha_array=alpha_array,
                                                                    beta_array=beta_array,
                                                                    sparse_array=sparse_array,
                                                                    update_arr=weight_array_iter,
                                                                    update_links=adj_links)      
                
                # TODO: eliminate paths where flow is 0 -> could incorporate this into loop above 
                marker=0
                for flow in flows:
                    if flow == 0:
                        flows.pop(marker)
                        paths.pop(marker)
                    else:
                        marker+=1   


            time1=time.time()
            # TODO: To avoid having to calculate this every iteration, if we know we are going to do more, add if statement to skip 
            #   for example, only calculate these vars every other iteration, or every five iterations
            # BEGIN NEW SPTT PARALLEL CALC
            # For shortest path travel time calculation (termination criteria)
            SPTT=0  
            # BEGIN PARALLEL CALC
            n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]), 'cost', 'weight_array_iter',capacity_array, weight_array_iter] for n in range(0, len(OD_matrix))]
            # n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]), 'cost', 'weight_array_iter'] for n in range(0, len(OD_matrix))]
            pool=Poolm(pp)
            results = pool.map(SPTT_parallel, n)
            pool.close()
            pool.join()
            for idx, cost in enumerate(results):
                flow = OD_matrix[idx][1]
                SPTT += cost*flow
            # END NEW SPTT PARALLEL CALC
            time2=time.time()
            print(time2-time1)

            # # calculate  termination criteria variables
            # SPTT = 0
            # for idx, x in enumerate(OD_matrix):
            #     origin = int(x[0])
            #     destination = int(x[2])
            #     flow = x[1]
            #     # # Calculate shortest paths
            #     backnode, costlabel = shortestPath_heap(origin=origin, 
            #                                             capacity_array=capacity_array, 
            #                                             weight_array=weight_array_iter, 
            #                                             adjlist=adjlist,
            #                                             sparse_array=sparse_array)
            #     # don't actually need shortest path, just need the cost of the path
            #     cost = costlabel[destination]
            #     # For shortest path travel time calculation (termination criteria)
            #     SPTT += cost*flow

            TSTT = np.sum(np.multiply(flow_array, weight_array_iter))
            AEC = (TSTT-SPTT)/sum_d
            RG = TSTT/SPTT - 1
            print(SPTT, TSTT, AEC, RG)

            # Fill termination variable lists
            AEC_list.append(AEC)
            TSTT_list.append(TSTT)
            SPTT_list.append(SPTT)
            RG_list.append(RG)
            timeb=time.time()
            print(f'iter {iter} time: {timeb-timea}')
            iter += 1

            # determine if iterations should continue
            iter_val = termination_function(termination_criteria=termination_criteria, 
                                            iters = iter, 
                                            AEC=AEC,
                                            RG=RG)
            
    elif algorithm == 'bush_based':
        
        # Convert network into a set of bushes for each OD pair
        bushes = []
        bush_flows = [] 

        # Termination criteria variables I am interested in keeping track off:
        AEC_list = []
        TSTT_list = []
        SPTT_list = []
        RG_list = []

        # BUG: Don't know why, but this caused a bug for it to get caught when calclating alternative paths
        # ended up impacting speed very little, can try again later or just ignore for now
        # Initizalize bushes for each OD pair
        # BEGIN PARALLEL
        # n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]),'path_cost_backnode', 'weight_array'] for n in range(0, len(OD_matrix))]
        # pool=Poolm(pp)
        # results = pool.map(SPTT_parallel, n)
        # for idx, result in enumerate(results):
        #     path = result[0]
        #     cost = result[1]
        #     backnode = result[2]
        #     flow = OD_matrix[idx][1]
        #     origin = OD_matrix[idx][0].astype(int)
        #     destination = OD_matrix[idx][2].astype(int)

        #     # Create an empty bush that is the same shape as weight_array
        #     if sparse_array is True:
        #         bush = sparse.lil_array(np.shape(copy_array))
        #     else:
        #         bush = np.full_like(copy_array, -1)

        #     # add appropriate edges to bush, using a shortest path tree
        #     for i, j in enumerate(backnode):
        #         if j ==-1:
        #             pass
        #         else:
        #             # I think J and I are placed properly?
        #             if sparse_array is True:
        #                 bush[[j],[i]]=1
        #             else:
        #                 bush[j][i]=1

        #     # Set initial flow for each bush equal to flow along shortest path
        #     # create empty bush flow array
        #     if sparse_array is True:
        #         bush_flow = sparse.lil_array(np.shape(copy_array))
        #     else:
        #         bush_flow = np.zeros_like(copy_array)

        #     # update the bush flow array
        #     for i in path:
        #         if sparse_array is True:
        #             bush_flow[i[0], i[1]] += flow
        #         else:
        #             bush_flow[i[0]][i[1]] += flow

        #     # append the trackers with topo orders and bushes
        #     bushes.append(bush)
        #     bush_flows.append(bush_flow)

        # pool.close()
        # pool.join()
        # #END PARALLEL


        # UNPARALLELIZED BUSHES INITIALIZING
        for pair in OD_matrix:
            origin = pair[0].astype(int) 
            flow = pair[1]
            destination = pair[2].astype(int)
            # Create an empty bush that is the same shape as weight_array
            if sparse_array is True:
                bush = sparse.lil_array(np.shape(copy_array))
            else:
                bush = np.full_like(copy_array, -1)
            # Find shortest path from origin to all locations - need to determine what nodes are unreachable 
            backnode, costlabel = shortestPath_heap(origin=origin, 
                                                    capacity_array=capacity_array, 
                                                    weight_array=weight_array, 
                                                    adjlist=adjlist,
                                                    sparse_array=sparse_array)
            # create alternative to prim's just using the backnode array? not minimum spanning tree but shortest path graph?
            #   TODO: Did i fix this??
            for i, j in enumerate(backnode):
                if j ==-1:
                    pass
                else:
                    # I think J and I are placed properly?
                    if sparse_array is True:
                        bush[[j],[i]]=1
                    else:
                        bush[j][i]=1
            # Set initial flow for each bush equal to flow along shortest path
            # create empty bush flow array
            if sparse_array is True:
                bush_flow = sparse.lil_array(np.shape(copy_array))
            else:
                bush_flow = np.zeros_like(copy_array)
            # determine shortest paths
            backnode, costlabel = shortestPath_heap(origin=origin,
                                                    capacity_array=capacity_array,
                                                    weight_array=weight_array,
                                                    adjlist=adjlist,
                                                    sparse_array=sparse_array)
            # extract path data
            path, cost = pathTo_mod(backnode=backnode, 
                                    costlabel=costlabel, 
                                    origin=origin, 
                                    destination=destination)
            # update the bush flow array 
            for i in path:
                if sparse_array is True:
                    bush_flow[i[0], i[1]] += flow
                else:
                    bush_flow[i[0]][i[1]] += flow

            # append the trackers with topo orders and bushes
            bushes.append(bush)
            bush_flows.append(bush_flow)
            # UNPARALLELIZED BUSH INITIALIZATION


        # create empty flow array
        if sparse_array is True:
            flow_array = sparse.lil_array(np.shape(copy_array))
        else:
            flow_array = np.zeros_like(copy_array, dtype=np.float64)
        
        # sum flow on bushes for initial network flow
        for bush in bush_flows:
            flow_array += bush

        if sparse_array is False:
            # have to replace the negative flows in flow_array with 0s, only applicable in numpy calcautions
            flow_array[flow_array<0]=0

        # # calculate initial network travel times
        weight_array_iter = link_performance_function(capacity_array=capacity_array,
                                                        flow_array=flow_array,
                                                        weight_array=weight_array,
                                                        eq=link_performance,
                                                        alpha_array=alpha_array,
                                                        beta_array=beta_array,
                                                        sparse_array=sparse_array)


        def label_function(topo_order, bush, weight_array_iter):
            '''
            function to calculate all of the L and U labels for algorithm B
            '''
            # Set L and U variables
            L_link = np.full_like(weight_array_iter, np.inf)
            U_link = np.full_like(weight_array_iter, -np.inf)
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
                            L_link[i][j] = weight_array_iter[i][j]
                            U_link[i][j] = weight_array_iter[i][j]

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
                                L_link[i][j] = weight_array_iter[i][j]
                                U_link[i][j] = weight_array_iter[i][j]
                        L_node[id] = 0
                        U_node[id] = 0

                    else:
                        L_node[id] = l_val
                        U_node[id] = u_val
                        for j in topo_order[id:]:
                            if bush[i][j] == 1:
                                L_link[i][j] = min(L_link[i][j], L_node[id]+weight_array_iter[i][j])
                                U_link[i][j] = max(U_link[i][j], U_node[id]+weight_array_iter[i][j])

                id += 1

            # remove infs for clarity
            U_link[np.isinf(U_link)] = 0
            L_link[np.isinf(L_link)] = 0

            return L_node, U_node, L_link, U_link

        def label_function_sparse(topo_order, bush, weight_array_iter):
            '''
            function to calculate all of the L and U labels for algorithm B using sparse arrays
            '''

            # Set L and U variables
            infs = [np.inf]*sparse.lil_array(weight_array_iter)
            neginfs=[-np.inf]*sparse.lil_array(weight_array_iter)
            L_link = sparse.lil_array(infs)
            U_link = sparse.lil_array(neginfs)
            
            L_node = np.full(len(topo_order), np.inf)
            U_node = np.full(len(topo_order), -np.inf)
            # L and U at r (origin) are equal to 0
            L_node[0] = 0
            U_node[0] = 0

            # determine L and U labels in forward topological ordering
            id = 0
            
            for id, i in enumerate(topo_order):
                # i == topo_order node & id == index
                # for first in topological order
                if id == 0:
                    for j in topo_order:
                        if bush[i,j] == 1:
                            L_link[i,j] = weight_array_iter[i,j]
                            U_link[i,j] = weight_array_iter[i,j]

                # for last in topological order
                elif id == len(topo_order)-1:
                    l_val = L_node[id]
                    u_val = U_node[id]
                    for h in topo_order:
                        if bush[h,i] == 1:
                            l_val = min(l_val, L_link[h,i])
                            u_val = max(u_val, U_link[h,i])
                    L_node[id] = l_val
                    U_node[id] = u_val

                # for all others
                else:
                    l_val = L_node[id]
                    u_val = U_node[id]
                    for h in topo_order[:id+1]:
                        if bush[h,i] == 1:
                            l_val = min(l_val, L_link[h,i])
                            u_val = max(u_val, U_link[h,i])

                    if l_val == np.inf:
                        for j in topo_order[id:]:
                            if bush[i,j] == 1:
                                L_link[i,j] = weight_array_iter[i,j]
                                U_link[i,j] = weight_array_iter[i,j]
                        L_node[id] = 0
                        U_node[id] = 0
                    
                    else:
                        L_node[id] = l_val
                        U_node[id] = u_val
                        for j in topo_order[id:]:
                            if bush[i,j] == 1:
                                L_link[i,j] = min(L_link[i,j], L_node[id]+weight_array_iter[i,j])
                                U_link[i,j] = max(U_link[i,j], U_node[id]+weight_array_iter[i,j])
            
            # removes infs from U_link and L_link for clarity
            U_link = sparse.csr_array(U_link)
            U_link.data = np.nan_to_num(U_link.data, posinf=0,neginf=0, copy=False)
            U_link.eliminate_zeros()
            U_link = sparse.lil_array(U_link)

            L_link = sparse.csr_array(L_link)
            L_link.data = np.nan_to_num(
                L_link.data, posinf=0, neginf=0, copy=False)
            L_link.eliminate_zeros()
            L_link = sparse.lil_array(L_link)



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
            
            for idx, bush in enumerate(bushes):
                print(idx)
                start=time.time()
                # set variables
                origin = OD_matrix[idx][0].astype(int)
                flow = OD_matrix[idx][1]
                destination = OD_matrix[idx][2].astype(int)
                bush_flow = bush_flows[idx]
                # Calculate initial L and U Labels
                # conduct topological ordering from the origin node
                topo_order = topo_order_function(bush, 
                                                 origin=origin, 
                                                 weight_array=weight_array_iter, 
                                                 num_nodes=num_nodes, 
                                                 node_list=node_list,
                                                 sparse_array=sparse_array)
                
                # calculate initial L and U labels
                # BUG: I think this can be rewritten ignoring topo-order - just use min and max paths, which should be much quicker to calc
                if sparse_array is True:
                    L_node, U_node, L_link, U_link = label_function_sparse(topo_order=topo_order,
                                                                           bush=bush,
                                                                           weight_array_iter=weight_array_iter)
                else:
                    L_node, U_node, L_link, U_link = label_function(topo_order=topo_order, 
                                                                    bush=bush, 
                                                                    weight_array_iter=weight_array_iter)

                # 1. FIND SHORTCUTS AND ADD TO BUSH
                # TODO: I don't like adding every possible shortcut - this doesn't seem efficent to me, might be adding in loops to network: Need to explore more
                # need to find where link doesn't exist on bush but does exist on weight array
                # #  Method below speeds up sparse_array shortcut addition by approx. 0.5 seconds per bush iteration 
                if sparse_array is True:
                    for i, rows in enumerate(capacity_array.rows):
                        for j in rows:
                            if (bush[i, j] == 0): 
                                try:
                                    if (L_node[topo_order.index(i)] + weight_array_iter[i, j] <= L_node[topo_order.index(j)]):
                                        bush[i,j] = 1
                                except:
                                    pass
                # if sparse_array is True:
                #     for i in topo_order:
                #         for j in topo_order:
                #             if (bush[i,j] == 0) and (capacity_array[i,j] > 0):      # doesen't exist in bush but does exist in weight array (i.e., whole network)
                #                 if L_node[topo_order.index(i)] + weight_array_iter[i,j] <= L_node[topo_order.index(j)]:
                #                     bush[i,j] = 1             
                else:
                    for i in topo_order:
                        for j in topo_order:
                            if (bush[i][j] == -1) and (capacity_array[i][j] > 0):      # doesen't exist in bush but does exist in weight arary (i.e., whole network)
                                if L_node[topo_order.index(i)] + weight_array_iter[i][j] <= L_node[topo_order.index(j)]:
                                    bush[i][j] = 1



                # 2. CALCULATE L AND U LABELS IN FORWARD TOPOLOGICAL ORDER WITH NEW SHORTCUTS
                # conduct topological ordering from the origin node
                topo_order = topo_order_function(bush, 
                                                 origin=origin, 
                                                 weight_array=weight_array_iter, 
                                                 num_nodes=num_nodes, 
                                                 node_list=node_list, 
                                                 sparse_array=sparse_array)
                # calculate new L and U labels
                if sparse_array is True:
                    L_node, U_node, L_link, U_link = label_function_sparse(topo_order=topo_order,
                                                                           bush=bush,
                                                                           weight_array_iter=weight_array_iter)
                else:
                    L_node, U_node, L_link, U_link = label_function(topo_order=topo_order, 
                                                                    bush=bush, 
                                                                    weight_array_iter=weight_array_iter)

                
                # instead of being interested in last node topoligcally, maybe we are interested in destination node?
                # BUG: this is helping calculate delta_hs but it really is supposed to be last node topologically from notes
                topo_order=topo_order[:topo_order.index(num_nodes-1)+1]

                # DIVERGENCE NODE LOOP
                loop = True
                if sparse_array is True:
                    while loop is True:
                        if len(topo_order) == 1:
                            loop = False
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
                                if bush[h,i] == 1:
                                    if max_val < U_link[h,i]:
                                        max_val = U_link[h,i]
                                        next_node_max = h
                            max_path.append(next_node_max)

                            # loop for max path
                            while max_path[-1] != origin:
                                max_val = -np.inf
                                for h in topo_order[:topo_order.index(next_node_max)+1]:
                                    if bush[h,next_node_max].astype(int) == 1:
                                        if max_val < U_link[h,next_node_max]:
                                            max_val = U_link[h,next_node_max]
                                            max_val_flow = min(max_val_flow, bush_flow[h,next_node_max])
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
                                    if bush[h,i] == 1:
                                        if min_val > L_link[h,i]:
                                            min_val = L_link[h,i]
                                            next_node_min = h
                                min_path.append(next_node_min)
                                # loop for min path
                                while min_path[-1] != origin:
                                    min_val = np.inf
                                    for h in topo_order[:topo_order.index(next_node_min)+1]:
                                        if bush[h,next_node_min].astype(int) == 1:
                                            if min_val > L_link[h,next_node_min]:
                                                min_val = L_link[h,next_node_min]
                                                next_node_holder = h
                                    min_path.append(next_node_holder)
                                    next_node_min = next_node_holder

                                # reverse order of lists for better comprehension
                                min_path = min_path[::-1]
                                max_path = max_path[::-1]

                                # if max and min path are the same, there is only one path and therefore we do not need to shift flow
                                if min_path == max_path:
                                    loop = False

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
                                                a = val
                                        separation = min_path.index(i) - min_path.index(a)

                                    # extract sigma_l and sigma_u, the list of links from a to i
                                    sigma_U = max_path[max_path.index(a):(max_path.index(i)+1)]
                                    sigma_L = min_path[min_path.index(a):(min_path.index(i)+1)]

                                    # calculate delta h, equation 6.61 in the textbook
                                    num1 = U_node[topo_order.index(i)] - U_node[topo_order.index(a)]
                                    num2 = L_node[topo_order.index(i)] - L_node[topo_order.index(a)]
                                    numerator = num1-num2

                                    # BPR function derivative
                                    derivative_array = link_performance_derivative(capacity_array=capacity_array,
                                                                                   flow_array=flow_array,
                                                                                   weight_array=weight_array,
                                                                                   eq=link_performance,
                                                                                   alpha_array=alpha_array,
                                                                                   beta_array=beta_array,
                                                                                   sparse_array=sparse_array)

                                    denominator = 0
                                    for id in range(len(sigma_L)-1):
                                        denominator += derivative_array[[sigma_L[id]],[sigma_L[id+1]]]
                                    for id in range(len(sigma_U)-1):
                                        denominator += derivative_array[[sigma_U[id]],[sigma_U[id+1]]]

                                    # set delta_h value
                                    delta_h = numerator/denominator

                                    # determine delta_h based on maximum flow that can be shifted
                                    for id in range(len(sigma_U)-1):
                                        if delta_h <= bush_flow[sigma_U[id],sigma_U[id+1]]:
                                            pass
                                        else:
                                            delta_h = bush_flow[sigma_U[id],sigma_U[id+1]]

                                    # Shift flows
                                    for id in range(len(sigma_L)-1):
                                        flow_array[sigma_L[id],sigma_L[id+1]] += delta_h
                                        bush_flow[sigma_L[id],sigma_L[id+1]] += delta_h
                                    for id in range(len(sigma_U)-1):
                                        flow_array[sigma_U[id],sigma_U[id+1]] -= delta_h
                                        bush_flow[sigma_U[id],sigma_U[id+1]] -= delta_h

                                    # at end of loop, slice off last value and repeat
                                    topo_order = topo_order[:-1]

                else:
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

                                    # BPR function derivative
                                    derivative_array = link_performance_derivative(capacity_array=capacity_array, 
                                                                                    flow_array=flow_array, 
                                                                                    weight_array=weight_array, 
                                                                                    eq=link_performance, 
                                                                                    alpha_array=alpha_array, 
                                                                                    beta_array=beta_array,
                                                                                    sparse_array=sparse_array)

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
                weight_array_iter = link_performance_function(capacity_array,
                                                              flow_array,
                                                              weight_array,
                                                              eq=link_performance,
                                                              alpha_array=alpha_array,
                                                              beta_array=beta_array,
                                                              sparse_array=sparse_array)


                # TODO Won't work right now but need to figure out how to remove unused edges
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
                end=time.time()
                print(end-start)

            # calculate termination criteria
            SPTT = 0
            for idx, x in enumerate(OD_matrix):
                origin = int(x[0])
                destination = int(x[2])
                flow = x[1]
                # Calculate shortest paths
                backnode, costlabel = shortestPath_heap(origin=origin,
                                                       capacity_array=capacity_array,
                                                       weight_array=weight_array_iter,
                                                       adjlist=adjlist,
                                                       destination=False,
                                                       sparse_array=sparse_array)

                # don't actually need shortest path, just need the cost of the path
                cost = costlabel[destination]
                
                # For shortest path travel time calculation (termination criteria)
                SPTT += cost*flow

            print(SPTT)
            # termination criteria
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
                                            iters=iter,
                                            AEC=AEC,
                                            RG=RG)

    ########################################################################################################################
    
    # EDIT THE FINAL OUTPUTS
    # Graph edge attributes that are added:
        # 'Weight_Array_Iter':  The new travel time given final flow solution
        # TA_flow':             Final flow on that edge

    # Set the background Weight_Array_Iter value to travel_time and update accordingly
    nx.set_edge_attributes(G, nx.get_edge_attributes(G, 'travel_time'), name='Weight_Array_Iter')

    # relate flow array back to network
    if sparse_array is True:
        # set background TA flow value and update accordingly
        nx.set_edge_attributes(G, 0, name='TA_Flow')
        # for i in og_node_list:
        #     for j in og_node_list:
        for i, row in enumerate(flow_array.rows):
            for j in row:
                if flow_array[i,j] > 0:
                    nx.set_edge_attributes(G, {(i, j): {'TA_Flow': flow_array[i,j], 'Weight_Array_Iter': weight_array_iter[i,j]}})
                    # new_data.update({(i, j): {'TA_Flow': flow_array[i,j], 'Weight_Array_Iter': weight_array_iter[i,j]}})
                # else:
                #     nx.set_edge_attributes(G, {(i, j): {'TA_Flow': 0}})
                    # new_data.update({(i, j): {'TA_Flow': 0}})
    
    else:
        new_data= {}
        for i in og_node_list:
            for j in og_node_list:
                if flow_array[i][j] > 0:
                    new_data.update({(i, j): {'TA_Flow': flow_array[i][j], 'Weight_Array_Iter': weight_array_iter[i][j]}})
                else:new_data.update({(i, j): {'TA_Flow': 0}})

    # set edge attributes with new data
    # nx.set_edge_attributes(G, new_data)

    # remove artificial super source and super sink
    if isinstance(dest_points, gpd.GeoDataFrame):
        G.remove_node(len(G.nodes())-1)
        G.remove_node(super_origin)
    elif isinstance(dest_points, list):
        temp_len=len(G.nodes())-1
        for x in range(0, len(dest_points)):
            G.remove_node(temp_len-x)
        G.remove_node(super_origin)

    # if multiple dest method, remove the intermediate artifical links
    if dest_method == 'multiple':
        while len(list(G.nodes)) > len(og_node_list):
            G.remove_node(len(G.nodes()))
    
    #TODO:  has to be multidigraph to save properly should adjust save function to represent this requirement
    G_output = nx.MultiDiGraph(G)

    # delete global variables created for parallel processing
    del capacity_array
    del weight_array
    del adjlist
    del weight_array_iter



    return G_output, AEC_list, TSTT_list, SPTT_list, RG_list, iter


def redundancy_metric(G,
                      res_points: gpd.GeoDataFrame,
                      res_parcels: gpd.GeoDataFrame,
                      dest_points: list[gpd.GeoDataFrame],
                      dest_parcels: list[gpd.GeoDataFrame],
                      inc=10,
                      max_length=180,
                      G_capacity: str = 'capacity', 
                      G_weight: str = 'travel_time',
                      norm=True):
    """
    Calculates the individual network redundancy metric.

    This metric is based on the k-shortest path algorithm. The k-shortest paths are calculated for each residential parcel, up until the time to go from the origin to destination is within a specific "time elongation factor" of the shortest routing time. The time-elongation factor default is 1.6 (i.e., the longest acceptable time to reach a resource can be, at most, X1.6 the shortest available time). See research for more (https://doi.org/10.1016/j.trd.2022.103275).

    The final metric quantified as follows:
        -determine number of paths that go to each specific resource location (e.g., grocery store 1, grocery store 2, etc.)
        -multiply number of paths to each specific resource location by the number of resource locations that can be reached (e.g., if a person can go to 2 different locations, multiply the number of paths to each store by 2)
        -sum the value across all the resources for a final score
        -repeat for each resource type

    FOR NOW: Assuming multiple list of destinations and nearest_nodes_vertices, can add other functionality later
    
    #TODO: units check for inc and max_length
    # inc: incremental length for each iterative isoline for scoring purposes
    # max_length: maximimum allowable travel time
    # norm: normalize each resources travel times before summation, to balance effects of resource with numerous locations vs. those with few 
    """



    # Convert graph to digraph format
    G = nx.DiGraph(G)

    # Travel times must be whole numbers -  round values up to nearest whole number
    for x in G.edges:
        G.edges[x][G_weight] = math.ceil(G.edges[x][G_weight])

    # Begin OD-matrix construction
        # takes form of [origin, destination, minimum cost]
        # TODO: This is heavily copied from traffic_assignment function - should consider putting all in a seperate helper 
        # steps are as follows:
            #a. determine sources and sinks
            #b. add artifical edges
            #c. create OD-matrix
    
    # The OD-matrix needs include paths from all origins to those destinations that are within the time-elongation factor of the shortest path time
    # What this means, is that each resource option that is within time_elongation factor of the absolute shortest path, needs an entry in OD-matrix

    # Set variables that will be used as constants throughout algorithm - set keeping int32 byte criteria in mind
    artificial_weight_min=1
    artificial_capacity = 999999999

    # a. Determine sources and sinks
    G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points = nearest_nodes_vertices(G=G, res_points=res_points, dest_parcels=dest_parcels, dest_points=dest_points, G_demand='demand')

    # b_2. Add artifical edges from vertices around destinations to aggregate points than to artifical sink with 0 cost and max capacity
    # This is to allow minimum cost routing to whatever resource - destination of OD pair which can change based on cost
    # identify the destination nodes, and create artifical sink edges
    # need to relate the nodes that are nearest to corners of parcels with the dest_points to associate the appropriate capacity to 
    # artificial node id tracker - useful in maintianing compatability with dtype
    dest_node_ids=99999998

    # need a list of the sink nodes for each resource
    sink_nodes=[]

    for i, parcels in enumerate(dest_parcels):
        #add new entry to sink_nodes
        sink_nodes.append([])
        for idx, dest_parcel in parcels.iterrows():
            dest_node = dest_points[i][dest_points[i].geometry.within(dest_parcel.geometry)]
            #since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
            for a, node in dest_node.iterrows():
                dest_node_ids-=1
                # add dest_node_ids to the sink_nodes list
                sink_nodes[i].append(dest_node_ids)
                # add the dest node to the graph 
                G.add_nodes_from([(dest_node_ids, {'demand': 0})])

                # add links from nearest intersections to parcel centroid
                for nearest_intersection in unique_dest_nodes_list[i][idx]:
                    kwargs = {G_weight: artificial_weight_min, G_capacity: artificial_capacity}
                    G.add_edge(nearest_intersection, dest_node_ids, **kwargs)
    
    # c. Create OD_matrix 
    # create redundnacy_dictionary output
    redundancy_dict = {}
    output_dict={}

    # OD_matrix needs to contain path from origin to specific resource location (not just super_sink)
    # this way we can assess redundancy to specific resource locations, not just the closest

    # create copies of G, G_weight, and other intput variables to make global for parallelization
    global G_copy
    global G_weight_copy
    global inc_copy
    global max_length_copy
    G_copy = copy.deepcopy(G)
    G_weight_copy = copy.deepcopy(G_weight)
    inc_copy=copy.deepcopy(inc)
    max_length_copy=copy.deepcopy(max_length)


    # need to stop calculating shortest path for each and every Origin-Destionation
    # Instead, just find shortest path from origin to all destinations, find cost of specific paths, and that's it

    # initialize OD_matrix: origin, sink_nodes
    OD_matrix = []

    # for each origin
    for origin in unique_origin_nodes:
        OD_matrix.append([origin,sink_nodes])

    # number of parallel processes 
    pp=8
    pool=Poolm(pp)
    results = pool.map(redundancy_parallel, OD_matrix)
    pool.close()
    pool.join()

    # create list of max_values to keep track for normalizing purposes
    max_values=[0]* len(unique_dest_nodes_list)
    # append the appropriate key-value pairs into output_dict
    for items in results:
        for result in items:
            #format of result: result= [source, target, value, resource number]
            i=result[3]
            
            # determine if the key alreadu exists in output_dict. If it does, append results. If not, create new key entry
            var = output_dict.get(result[0])
            # if key (origin) doesn't exist
            if var is None:
                output_dict[result[0]] = {f'Res {i} list of redundancy scores': [result[2]]}     # list of scores

                redundancy_dict[result[0]] = {'redundancy scores': [result[2]]}
                max_values[i] = max(max_values[i], result[2])
            # if key does exist
            else:
                # if key exists but the value doesn't
                if f'Res {i} list of redundancy scores' not in output_dict[result[0]]:
                    output_dict[result[0]][f'Res {i} list of redundancy scores'] = [result[2]]

                    redundancy_dict[result[0]]['redundancy scores'].append(result[2])
                    max_values[i] = max(max_values[i], result[2])

                else:
                    output_dict[result[0]][f'Res {i} list of redundancy scores'].append(result[2])

                    # update the redundancy score
                    score = 0
                    for val in output_dict[result[0]][f'Res {i} list of redundancy scores']:
                        score += val* np.count_nonzero(output_dict[result[0]][f'Res {i} list of redundancy scores'])

                    redundancy_dict[result[0]]['redundancy scores'][-1] = score
                    max_values[i] = max(max_values[i], score)


    if norm is True:
        # normalize redundancy scores, and find redundancy score total
        for source in redundancy_dict:
            total=0
            for i, val in enumerate(redundancy_dict[source]['redundancy scores']):
                norm_val = val / max_values[i]
                redundancy_dict[source]['redundancy scores'][i] = norm_val
                total += norm_val
            redundancy_dict[source]['redundancy score total'] = total
    else:
        for source in redundancy_dict:
            total = 0
            for val in redundancy_dict[source]['redundancy scores']:
                total += val
            redundancy_dict[source]['redundancy score total'] = total
    
    # convert dictionary to dataframe, and merge with res_points geodataframe
    data_df = pd.DataFrame(redundancy_dict).T
    data_df = data_df.astype({'redundancy score total':'float'})
    merged_gdf = res_points.merge(data_df, left_on='nearest_node', right_index=True)

    # Spatial join the res_parcels/res_points 
    res_parcels = gpd.sjoin(res_parcels, merged_gdf)

    # need to convert redundancy scores to individual attribute columns in order to export as a shapefile
    new_columns = []
    for i in range(len(dest_points)):
        new_columns.append(f'r{i}_redun')
    res_parcels[new_columns] = res_parcels['redundancy scores'].apply(pd.Series)
    res_parcels=res_parcels.drop('redundancy scores', axis=1)
    # edit the attributes of shapefile so only attributes we are intereseted in from this function are returned
    res_attributes = ['geometry', 
                      'nearest_node',
                      'redundancy score total']
    for col in new_columns:
        res_attributes.append(col)
    res_parcels.drop(columns=[col for col in res_parcels if col not in res_attributes], inplace=True)

    # delete global variables
    del G_copy
    del G_weight_copy
    del inc_copy
    del max_length_copy

    # TODO: determine value in returning output dict, or just remove it
    # FOR NOW: Just return res_parcels
    return res_parcels


def network_redundancy_metric(G,
                                res_points: gpd.GeoDataFrame,
                                res_parcels: gpd.GeoDataFrame,
                                dest_points: list[gpd.GeoDataFrame],
                                dest_parcels: list[gpd.GeoDataFrame],
                                G_capacity: str = 'capacity', 
                                G_weight: str = 'travel_time'):
    """
    Modified edge disjoint paths algorithm to determine alternative paths from origins to destinations
    """
    
    # Convert graph to digraph format
    G = nx.DiGraph(G)

    # Travel times must be whole numbers -  round values up to nearest whole number
    for x in G.edges:
        G.edges[x][G_weight] = math.ceil(G.edges[x][G_weight])

    # remove all edges with 0 capacity
    remove = [(u, v, d) for u, v, d in G.edges(data=True) if d[G_capacity] == 0]
    G.remove_edges_from(remove)

    # Determine sources and sinks
    G, unique_origin_nodes, unique_dest_nodes_list, positive_demand, shared_nodes, res_points, dest_parcels, dest_points = nearest_nodes_vertices(G=G, res_points=res_points, dest_parcels=dest_parcels, dest_points=dest_points, G_demand='demand')

    # find the neighbors-of-neighbors of unique_dest_nodes_list to crate artifical edges
    final_dest_nns=[]
    # for each resource type,
    for dest_nodes_list in unique_dest_nodes_list:
        int_dest_nns=[]
        # for each option of that resource type,
        for dest in dest_nodes_list:
            # for each nearest_node for that option of that type,
            for node in dest:
                dest_nn=[]
                # first set of neighbors
                dest_n1 = [a for a in G.neighbors(node)]
                for neighbor in dest_n1:
                    dest_n2 = [a for a in G.neighbors(neighbor)]
                    # if neighbor is a new neighbor, remove it
                    try: 
                        dest_n2.remove(neighbor)
                    except ValueError:
                        pass
                    dest_nn = dest_nn + dest_n2
            
                # remove duplicates 
                dest_nn = list(set(dest_nn))

            int_dest_nns.append(dest_nn)
        final_dest_nns.append(int_dest_nns)

    
    # Set variables that will be used as constants throughout algorithm - set keeping int32 byte criteria in mind
    artificial_weight_min = 1
    artificial_capacity = 999999999

    # Add artifical edges from vertices around destinations (their neighbors-of-neighbors nodes) to aggregate points 
    # identify the destination nodes, and create artifical sink edges
    # artificial node id tracker - useful in maintianing compatability with dtype
    dest_node_ids = 99999998
    # need a list of the sink nodes for each resource
    sink_nodes = []
    for i, parcels in enumerate(dest_parcels):
        # add new entry to sink_nodes
        # THIS IS DIFFERENT THEN OTHER FUNCTIONS
        # sink_nodes.append([])
        dest_node_ids -= 1
        sink_nodes.append(dest_node_ids)
        
        for idx, dest_parcel in parcels.iterrows():
            dest_node = dest_points[i][dest_points[i].geometry.within(dest_parcel.geometry)]
            
            # add new sink node
            G.add_nodes_from([(dest_node_ids, {'demand': 0})])

            # since we could have multiple dest nodes within a single boundary (multiple resources located at same parcel) need to iterate through dest_node
            for a, node in dest_node.iterrows():
                # dest_node_ids -= 1
                # # add dest_node_ids to the sink_nodes list
                # sink_nodes[i].append(dest_node_ids)
                # # add the dest node to the graph
                # G.add_nodes_from([(dest_node_ids, {'demand': 0})])

                # add links from NEIGHBORS-OF-NEIGHBORS of nearest nodes of parcel to the parcel centroid
                for nearest_intersection in final_dest_nns[i][idx]:
                    kwargs = {G_weight: artificial_weight_min,
                              G_capacity: artificial_capacity}
                    G.add_edge(nearest_intersection, dest_node_ids, **kwargs)


    # make other necessary variables global for parallel processing
    global G_copy
    global G_weight_copy
    global G_capacity_copy
    global sink_nodes_copy
    G_copy = copy.deepcopy(G)
    G_weight_copy = copy.deepcopy(G_weight)
    G_capacity_copy = copy.deepcopy(G_capacity)
    sink_nodes_copy = copy.deepcopy(sink_nodes)

    # number of parallel processes 
    pp=8
    pool=Poolm(pp)
    results = pool.map(edge_disjoint_paths, unique_origin_nodes)
    pool.close()
    pool.join()

    # convert results to dataframe
    # create column names
    columns=['origin']
    for idx, dest in enumerate(dest_parcels):
        columns.append(f'r{idx}_redun')
    results_df = pd.DataFrame(results, columns=columns)
    results_df.set_index('origin', inplace=True)
    # for each column, find values equal to -1 and set to maximum for that column
    for column in results_df.columns:
        max_value = results_df[column].max()  # Get the maximum value of the column
        results_df[column] = results_df[column].replace(-1, max_value)  # Replace -1 with the maximum value
        results_df[column] = results_df[column]/max_value
    #normalizing dataframe
    # results_df = results_df.apply(lambda iterator: ((iterator/iterator.max())))
    results_df['net_redun'] = results_df.sum(axis=1)

    # merge with res_points geodataframe
    merged_gdf = res_points.merge(results_df, left_on='nearest_node', right_index=True)

    # # Spatial join the res_parcels/res_points 
    res_parcels = gpd.sjoin(res_parcels, merged_gdf)

    # edit the attributes of shapefile so only attributes we are intereseted in from this function are returned
    res_attributes = ['geometry', 
                      'nearest_node',
                      'net_redun']
    for col in columns:
        res_attributes.append(col)
    res_parcels.drop(columns=[col for col in res_parcels if col not in res_attributes], inplace=True)



    del G_copy
    del G_weight_copy
    del G_capacity_copy
    del sink_nodes_copy

    return res_parcels