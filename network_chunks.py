# removed bits from exploration stuff file to clean up codes


import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
import geopandas as gpd

# Importing OSM Data Examples ###############################################
# can import OSM street level data by using a:
#   bounding box
#   lat-long point and radius
#   address
#   input polygon area
#   place name, i.e., 'Los Angeles, California'

# Address example
# G = ox.graph_from_address('2208 Pennsylvania Ave, Austin, Texas', network_type='drive')

# Place example
# G = ox.graph_from_place('Austin, Texas', network_type='drive')

# # Can easily plot graphs
# ox.plot_graph(G)


# Extract Edge attributes #####################################################
# retrieve the dictionaries of speed and travel time
speed_dict = nx.get_edge_attributes(G, 'speed_kph')
travel_time_dict = nx.get_edge_attributes(G, 'travel_time')

# retrieve speed and travel time values
speeds = []
for x in speed_dict:
    speeds.append(speed_dict[x])
travel_times = []
for x in travel_time_dict:
    travel_times.append(travel_time_dict[x])


#########
# more complicated plot function
# # More Difficult, but can also do manually
# print(G.edges)
# print(routes[0])
# edge_color = ['w']*len(G.edges)
# edge_line_width = [1]*len(G.edges)
# for idx, node in enumerate(routes[0]):
#     if idx == len(routes[0])-1:
#         pass
#     else:
#         num = 0
#         for u, v, k in G.edges:
#             if (u == routes[0][idx] and v == routes[0][idx+1]) or (v == routes[0][idx] and u == routes[0][idx+1]):
#                 edge_color[num] = 'r'
#                 edge_line_width[num] = 3
#             num += 1
#
#
# # node_sizes = [75 if u == routes[0][-1] else 0 for u in G.nodes]
#
# fig, ax = ox.plot_graph(G,
#                         node_size=node_sizes,
#                         edge_color=edge_color,
#                         edge_linewidth=edge_line_width,
#                         bgcolor='k')

# plotting
# node_sizes = [75 if u == routes[0][-1] else 0 for u in G.nodes]

# can plot using built function, but limited customizability of routes
# fig, ax = ox.plot_graph_routes(G, routes, route_colors=['r', 'b', 'g'],
#                               route_linewidth=6,
#                               orig_dest_size=0,
#                               node_size=node_sizes,
#                               bgcolor='k')
