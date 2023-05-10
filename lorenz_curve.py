# This is where I am testing out a lorenz curve function
import numpy as np
import networkx as nx
import osmnx as ox
import pandas as pd
import geopandas as gpd
import math
from matplotlib import pyplot as plt
import matplotlib
from matplotlib import patches
from matplotlib.lines import Line2D
from scipy.interpolate import make_interp_spline
import network_exploration_stuff as mynet
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)



plt.rcParams['font.family'] = "ubuntu"
# plt.rcParams['font.serif'] = "Times"
# plt.rcParams['font.sans-serif'] = 'Helvetica'
# NOTES:
# 177 BGs - double check

# FUNCTIONS
# define lorenz curve and gini coefficent functions
def gini(arr):
    # arr should already be sorted how we want it 
    # first sort
    sorted_arr = arr.copy()
    # sorted_arr.sort()
    n = arr.size
    coef_ = 2. / n
    const_ = (n + 1.) / n
    weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
    return coef_*weighted_sum/(sorted_arr.sum()) - const_

def lorenz_curve(X):
    """
    Creates Lorenz curve plot. X is a numpy array in the order of lowest to highest SVI
    """
    lorenz = X.cumsum() / X.sum()
    lorenz = np.insert(lorenz, 0, 0)
    vals = np.arange(lorenz.size)/(lorenz.size-1)

    fig, ax = plt.subplots(figsize=[6, 6])
    # Equality line
    eq_line, = ax.plot([0, 1], [0, 1], color='k',
                       linestyle='--', label='Equality Line')

    # # create interpoliated spline line
    # X_Y_Spline = make_interp_spline(lorenz, Y)
    # X_ = np.linspace(lorenz.min(), lorenz.max(), 1000)
    # Y_ = X_Y_Spline(X_)
    # lorenz_plot=ax.plot(X_,Y_, color='r')
    lorenz_plot, = ax.plot(vals, lorenz, label='tsamp', color='red')

    plt.title("Lorenz Curve of Impacted Roads")
    ax.set_xlabel("Block Group SVI")
    ax.set_ylabel("Normalized Cumulative Sum of Impacted Roads by Block Group")

    # add legend
    ax.legend(handles=[eq_line, lorenz_plot], loc='lower right')

    # add annotation box
    # annotation value
    val1 = gini(X)
    # annotation box text
    textstr = '\n'.join((f'Gini: {round(val1,4)}',))
    # patch properties
    props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
    # text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

def lorenz_curve_multiplot(X, row, column, axes, idx, names, method, annotations,tstamp=False):
    """
    Creates Lorenz curve plot. X is a numpy array in the order of lowest to highest SVI
    inputs:
    row: int, row number of plot
    col: int, col number of plot
    axes: axes of function
    idx: enumerate number of plot
    names: list of names of the resources
    method: x axis, should be 'SVI_scaled', 'Factor 1', or 'Factor 2', used only to determine line color
    tstamp: if False, just use normal red,blue, green, if not False, select appropriate value
    """
    # set plot color based on method
    if method == 'SVI_scaled':
        if tstamp is False:
            c_line='red'
        else:
            colors = ['#ffffcc', '#ffeda0', '#fed976', '#feb24c',
                      '#fd8d3c', '#fc4e2a', '#e31a1c', '#bd0026', '#800026']
            c_line = colors[tstamp]
    elif method == 'Factor 1':
        if tstamp is False:
            c_line='blue'
        else:
            colors = ['#ffffd9', '#edf8b1', '#c7e9b4', '#7fcdbb',
                      '#41b6c4', '#1d91c0', '#225ea8', '#253494', '#081d58']
            c_line = colors[tstamp]
    elif method == 'Factor 2':
        if tstamp is False:
            c_line='green'
        else:
            colors = ['#ffffe5', '#f7fcb9', '#d9f0a3', '#addd8e',
                      '#78c679', '#41ab5d', '#238443', '#006837', '#004529']
            c_line = colors[tstamp]

    lorenz = X.cumsum() / X.sum()
    lorenz = np.insert(lorenz, 0, 0)
    val = np.arange(lorenz.size)/(lorenz.size-1)

    # Equality line
    eq_line, = axes[row,column].plot([0, 1], [0, 1], color='k',linestyle='--', label='Equality Line')
    # lorenze curve
    lorenz_plot, = axes[row,column].plot(val, lorenz, label='tsamp', color=c_line, linewidth=2)
    # set title to resource name
    axes[row,column].set_title(names[idx], style='italic')

    # Set x and y tick markers
    axes[row, column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%','25%','50%','75%','100%'])
    axes[row, column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'])

    # set x and y limits
    axes[row, column].set_xlim(0,1)
    axes[row, column].set_ylim(0,1)


    if annotations is True:
        # add annotation box
        # annotation value
        val1 = gini(X)
        # annotation box text
        textstr = '\n'.join((f'Gini: {round(val1,4)}',))
        # patch properties
        props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
        # text box in upper left in axes coords
        axes[row,column].text(0.025, 0.975, 
                            textstr,
                            transform=axes[row,column].transAxes, 
                            fontsize=10,
                            verticalalignment='top', 
                            bbox=props)
        
        # add proportion of weight annotations
        # find the 25th and 75th percentile
        perc_25=np.percentile(lorenz,25)
        perc_75=np.percentile(lorenz,75)

        # plot vertical and horizontal lines for these locations
        perc_25_line_v, = axes[row, column].plot([0.25, 0.25], [0, perc_25], color='dimgray', linestyle='dotted')
        perc_25_line_h, = axes[row, column].plot([0, 0.25], [perc_25, perc_25], color='dimgray',linestyle='dotted')
        perc_75_line_v, = axes[row, column].plot([0.75, 0.75], [0, perc_75], color='dimgray',linestyle='dotted')
        perc_75_line_h, = axes[row, column].plot([0, 0.75], [perc_75, perc_75], color='dimgray',linestyle='dotted')

        # Plot double arrow annotations
        arrow1 = patches.FancyArrowPatch((0.025, 0), (0.025, perc_25), 
                                        arrowstyle='<->', 
                                        mutation_scale=10, 
                                        color='dimgray')
        arrow2 = patches.FancyArrowPatch((0.2, perc_25), (0.2, perc_75),
                                        arrowstyle='<->',
                                        mutation_scale=10,
                                        color='dimgray')
        arrow3 = patches.FancyArrowPatch((0.45, perc_75), (0.45, 1), 
                                        arrowstyle='<->', 
                                        mutation_scale=10,
                                        color='dimgray')
        axes[row,column].add_patch(arrow1)
        axes[row,column].add_patch(arrow2)
        axes[row,column].add_patch(arrow3)

        # Plot percentage values
        axes[row,column].annotate(text=f'{round(perc_25*100)}%',xy=(0.025,perc_25+.01))
        axes[row,column].annotate(text=f'{round((perc_75-perc_25)*100)}%',xy=(0.07,(perc_75-perc_25)/2+perc_25))
        axes[row,column].annotate(text=f'{round((1-perc_75)*100)}%', xy=(0.3, (1-perc_75)/2+perc_75))

#############################################################################################################
# read in svi shapefile
# OLD SVI
# svi = gpd.read_file(
#     '/home/mdp0023/Desktop/external/Data/Austin_North_BGs_SVI/Austin_North_BGs_SVI_projected.shp')
# STUDY AREA SPECIFIC SVI
# svi = gpd.read_file(
#     '/home/mdp0023/Documents/Codes_Projects/SVI_Code/Travis_County/SVI_Shapefiles/Travis_county_Austin_AOI_svi_2020_proj.shp')
# 2015 SPECIFIC -> relevant to the flood
svi = gpd.read_file(
    '/home/mdp0023/Documents/Codes_Projects/SVI_Code/Travis_County/SVI_Shapefiles/Travis_county_svi_2015_selected.shp')

# read in network and convert edges to shapefile
# BUG: has to be non TA graph for inundation lorenz curves
network = mynet.read_graph_from_disk(path='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs',
                                     name='AN_Graph_2015052522_inundation')
gdf_edges = ox.graph_to_gdfs(G=network, nodes=False)

# read in TA res parcels
res_parcels = gpd.read_file(
    '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs/Flow_decomp/2015052522_inundation_res_parcel_flow_decomp.shp')

# read in 2015 vars, grab ppunit data, merge with svi
vars2015 = pd.read_csv(
    '/home/mdp0023/Documents/Codes_Projects/SVI_Code/Travis_County/Calculated_Variables/Travis_county_var_2015.csv')
vars2015=vars2015[['PPUNIT','QBLACK','QSPANISH','QPOVTY','GEOID']].copy()
svi=svi.merge(vars2015,left_on='GEOID', right_on='GEOID')


#############################################################################################
# IMPACTED ROADS LORENZ PLOTS
def lorenz_curve_road_impacts_timeseries(fig, 
                                         axes,
                                        method='SVI_scaled',
                                        variable='agr_no_cap',
                                        times=[], 
                                        annotations=True, 
                                        columns=3,
                                        ftypes=['inundation']):
    """
    Function to plot lorenz curve of impact on roads, with options to show compound, fluvial, and/or pluvial impacts in a time series or not

    fig: matplotlib figure
    axes: matplotlib axes object
    method: string, either 'SVI_scaled', 'Factor 1', or 'Factor 2', what the x-axis is of final plot 
    variable: string, can be 'agr_no_cap', 'perc_no_cap', 'agr_increased_travel_time', 'perc_increased_travel_time', 'agr_impact', or 'perc_perc_impact'
    times: list of strings, the tstamps to include
    ftypes: list of strings, can contain 'inundation' for compound, 'inundation_fluvial' for fluvial, and 'inundation_pluvial' for pluvial
    """
    column = 0
    row = 0

    # read in appropriate network and convert to geopandas dataframe
    for a,time in enumerate(times):
        table_vals=[[0,0,0],[0,0,0],[0,0,0]]
        for b,ftype in enumerate(ftypes):

            network = mynet.read_graph_from_disk(path='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs',
                                            name=f'AN_Graph_{time}_{ftype}')
            gdf_edges = ox.graph_to_gdfs(G=network, nodes=False)

            # spatial join
            sjoined_data = gpd.sjoin(left_df=svi,
                                    right_df=gdf_edges,
                                    how='left')

            # count the number of roads within each block group and relate back to svi geodataframe
            count_dict = sjoined_data['GEOID'].value_counts().to_dict()
            svi["count"] = svi["GEOID"].apply(lambda x: count_dict.get(x))


            # count the number of roads within each block group with 0 capacity under agressive flood relationship
            subset_df = sjoined_data.loc[sjoined_data['inundation_capacity_agr'] == 0]
            count_dict = subset_df['GEOID'].value_counts().to_dict()
            svi["agr_no_cap"] = svi["GEOID"].apply(lambda x: count_dict.get(x))
            svi['agr_no_cap'] = svi['agr_no_cap'].fillna(0)

            # count the number of roads within each block group with an increased travel time under agressive flood relationship
            subset_df = sjoined_data.loc[sjoined_data['inundation_travel_time_agr']
                                        > sjoined_data['travel_time']]
            count_dict = subset_df['GEOID'].value_counts().to_dict()
            svi["agr_increased_travel_time"] = svi["GEOID"].apply(
                lambda x: count_dict.get(x))
            svi['agr_increased_travel_time'] = svi['agr_increased_travel_time'].fillna(0)

            # count the number of roads impacted (increased travel time & 0 capacity)
            svi['agr_impact'] = svi['agr_increased_travel_time'] + svi['agr_no_cap']

            # calc percentage of roads within a BG with no capacity
            svi['perc_no_cap'] = svi['agr_no_cap']/svi['count']*100

            # calc percentage of roads within a BG with an increased travel time
            svi['perc_increased_travel_time'] = svi['agr_increased_travel_time']/svi['count']*100

            # calc percentage of roads within a BG impacted
            svi['perc_impact'] = svi['agr_impact']/svi['count']*100

            # sort svi dataframe by appropriate column
            if method == 'Factor 2':
                svi.sort_values('Factor 2', axis=0, inplace=True, ascending=False)
            else:
                svi.sort_values(method, axis=0, inplace=True, ascending=True)

            # TODO: weight value by PPUNIT
            svi[variable] = svi[variable]*svi["PPUNIT"]

            # LORENZ CURVE
            # extract values to use in Lorenz curve
            if method == 'Factor 2':
                array = svi[variable].values*-1
            else:
                array = svi[variable].values

            # calculate lorenz values
            lorenz = array.cumsum() / array.sum()
            lorenz = np.insert(lorenz, 0, 0)
            val = np.arange(lorenz.size)/(lorenz.size-1)

            if ftype == 'inundation':
                c_line = '#cc57a4'
            elif ftype == 'inundation_pluvial':
                c_line = '#ff772e'
            elif ftype == 'inundation_fluvial':
                c_line = '#216aa6'

            # lorenze curve
            lorenz_plot, = axes[row, column].plot(val, lorenz, label='tsamp', color=c_line, linewidth=2)
            
            perc_25=np.percentile(lorenz,25)
            perc_75=np.percentile(lorenz,75)

            table_vals[b][0]=f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_75*100)-round(perc_25*100)}%'
            table_vals[b][2]=f'{100-round(perc_75*100)}%'

          
        col_labels = ['$\mathregular{B25^{th}}$', '$\mathregular{M50^{th}}$', '$\mathregular{T25^{th}}$']
        row_labels = ['Compound', 'Pluvial', 'Fluvial']
        # the rectangle is where I want to place the table
        if annotations is True:
            the_table = axes[row, column].table(cellText=table_vals,
                                colWidths=[0.12]*3,
                rowColours=['#cc57a4', '#ff772e', '#4d88b8'],
                                rowLabels=row_labels,
                                colLabels=col_labels,
                                loc='lower right')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(9)

        # Equality line
        eq_line, = axes[row,column].plot([0, 1], [0, 1], color='k',linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes[row, column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%','25%','50%','75%','100%'])
        axes[row, column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'])

        # set x and y limits
        axes[row, column].set_xlim(0,1)
        axes[row, column].set_ylim(0,1)
        
        # set axes title
        axes[row,column].set_title(f'{time[-2:]}:00 {time[6:8]}/{time[4:6]}/{time[:4]}', style='italic')

        # update row and column names
        column += 1
        if column == columns:
            column = 0
            row += 1


def lorenz_curve_road_impacts_individual(fig,
                                         axes,
                                         method='SVI_scaled',
                                         variable='agr_no_cap',
                                         times=[],
                                         annotations=True,
                                         columns=3,
                                         ftypes=['inundation']):
        

# read in appropriate network and convert to geopandas dataframe
    for a,time in enumerate(times):
        ginis = [0, 0, 0]
        table_vals=[[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        for b,ftype in enumerate(ftypes):

            network = mynet.read_graph_from_disk(path='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs',
                                            name=f'AN_Graph_{time}_{ftype}')
            gdf_edges = ox.graph_to_gdfs(G=network, nodes=False)

            # spatial join
            sjoined_data = gpd.sjoin(left_df=svi,
                                    right_df=gdf_edges,
                                    how='left')

            # count the number of roads within each block group and relate back to svi geodataframe
            count_dict = sjoined_data['GEOID'].value_counts().to_dict()
            svi["count"] = svi["GEOID"].apply(lambda x: count_dict.get(x))


            # count the number of roads within each block group with 0 capacity under agressive flood relationship
            subset_df = sjoined_data.loc[sjoined_data['inundation_capacity_agr'] == 0]
            count_dict = subset_df['GEOID'].value_counts().to_dict()
            svi["agr_no_cap"] = svi["GEOID"].apply(lambda x: count_dict.get(x))
            svi['agr_no_cap'] = svi['agr_no_cap'].fillna(0)

            # count the number of roads within each block group with an increased travel time under agressive flood relationship
            subset_df = sjoined_data.loc[sjoined_data['inundation_travel_time_agr']
                                        > sjoined_data['travel_time']]
            count_dict = subset_df['GEOID'].value_counts().to_dict()
            svi["agr_increased_travel_time"] = svi["GEOID"].apply(
                lambda x: count_dict.get(x))
            svi['agr_increased_travel_time'] = svi['agr_increased_travel_time'].fillna(0)

            # count the number of roads impacted (increased travel time & 0 capacity)
            svi['agr_impact'] = svi['agr_increased_travel_time'] + svi['agr_no_cap']

            # calc percentage of roads within a BG with no capacity
            svi['perc_no_cap'] = svi['agr_no_cap']/svi['count']*100

            # calc percentage of roads within a BG with an increased travel time
            svi['perc_increased_travel_time'] = svi['agr_increased_travel_time']/svi['count']*100

            # calc percentage of roads within a BG impacted
            svi['perc_impact'] = svi['agr_impact']/svi['count']*100

            # sort svi dataframe by appropriate column
            if method == 'Factor 2':
                svi.sort_values('Factor 2', axis=0, inplace=True, ascending=False)
            else:
                svi.sort_values(method, axis=0, inplace=True, ascending=True)

            # TODO: weight value by PPUNIT
            svi[variable] = svi[variable]*svi["PPUNIT"]

            # LORENZ CURVE
            # extract values to use in Lorenz curve
            if method == 'Factor 2':
                array = svi[variable].values*-1
            else:
                array = svi[variable].values

            # calculate lorenz values
            lorenz = array.cumsum() / array.sum()
            lorenz = np.insert(lorenz, 0, 0)
            val = np.arange(lorenz.size)/(lorenz.size-1)

            if ftype == 'inundation':
                c_line = '#cc57a4'
            elif ftype == 'inundation_pluvial':
                c_line = '#ff772e'
            elif ftype == 'inundation_fluvial':
                c_line = '#216aa6'

            # lorenze curve
            lorenz_plot, = axes.plot(val, lorenz, label='tsamp', color=c_line, linewidth=2)
            
            perc_25 = np.percentile(lorenz,25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz,75)

            table_vals[b][0]=f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3]=f'{100-round(perc_75*100)}%'

            # calc gini and save
            ginis[b] = gini(array)

        # plot ginis
        textstr = '\n'.join((f'Compound Gini: {round(ginis[0],2)}',
                            f'Pluvial Gini: {round(ginis[1],2)}',
                            f'Fluvial Gini: {round(ginis[2],2)}'))
        # patch properties
        props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
        # text box in upper left in axes coords
        axes.text(0.05, 0.95, textstr, transform=axes.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
          
        col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
        row_labels = ['Compound', 'Pluvial', 'Fluvial']
        # the rectangle is where I want to place the table
        if annotations is True:
            the_table = axes.table(cellText=table_vals,
                                colWidths=[0.10]*4,
                rowColours=['#cc57a4', '#ff772e', '#4d88b8'],
                                rowLabels=row_labels,
                                colLabels=col_labels,
                colColours=['lightgrey', 'lightgrey',
                                    'lightgrey', 'lightgrey'],
                                loc='lower right')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(14)
            the_table.scale(1, 1.5)

        # Equality line
        eq_line, = axes.plot([0, 1], [0, 1], color='k',linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes.set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%','25%','50%','75%','100%'],fontsize=14)
        axes.set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'],fontsize=14)

        # set x and y limits
        axes.set_xlim(0,1)
        axes.set_ylim(0,1)
        
        # set axes title
        axes.set_title(f'{time[-2:]}:00 {time[6:8]}/{time[4:6]}/{time[:4]}', style='italic', fontsize=14)


#############################################################################################
# TIME SERIES OF IMPACTED ROADS
# fig, axes = plt.subplots(2, 3, figsize=(15, 10))
# fig.set_facecolor('none')
# # fig.suptitle('Lorenz Curves of Impacted Roads',
# #              fontsize=14, weight='bold')
# fig.supxlabel('Social Vulnerability Index Percentile',
#               weight='bold',
#               x=0.5,
#               fontsize=20)
# fig.supylabel('Normalized Cumulative Sum of PPUNIT Weighted Closed Roads',
#               weight='bold',
#               fontsize=20)

# times = ['2015052521', '2015052522',
#          '2015052523', '2015052600',
#          '2015052601', '2015052602']

# lorenz_curve_road_impacts_timeseries(fig=fig,
#                                      axes=axes,
#                                      method='SVI_scaled',
#                                      variable='agr_no_cap',
#                                      times=times,
#                                      annotations=True,
#                                      columns=3,
#                                      ftypes=['inundation', 'inundation_pluvial','inundation_fluvial'])

# plt.tight_layout()
# plt.subplots_adjust(top=0.95,
#                     bottom=0.088,
#                     left=0.076,
#                     right=0.973,
#                     hspace=0.206,
#                     wspace=0.255)

# plt.show()


#############################################################################################
# # # INDIVIDUAL IMPACTED ROADS
# fig, axes = plt.subplots(1, 1, figsize=(6, 6))
# fig.set_facecolor('none')
# # fig.suptitle('Lorenz Curves of Impacted Roads',
# #              fontsize=14, weight='bold')
# axes.set_xlabel('Social Vulnerability Index Percentile',
#               weight='bold',
#               x=0.5,
#               fontsize=16)
# axes.set_ylabel('Normalized Cumulative Sum of Closed Roads',
#               weight='bold',
#               fontsize=16)

# times = ['2015052602']

# lorenz_curve_road_impacts_individual(fig=fig,
#                                      axes=axes,
#                                      method='SVI_scaled',
#                                      variable='agr_no_cap',
#                                      times=times,
#                                      annotations=True,
#                                      columns=3,
#                                      ftypes=['inundation', 'inundation_pluvial','inundation_fluvial'])

# plt.tight_layout()
# plt.subplots_adjust(top=0.95,
#                     bottom=0.093,
#                     left=0.151,
#                     right=0.953,
#                     hspace=0.206,
#                     wspace=0.255)

# plt.show()


#############################################################################################
# Residential access lorenz curve


def lorenz_curve_resource_accessibility(fig,
                                        axes, 
                                        method='SVI_scaled', 
                                        annotations=True, 
                                        weight=10,
                                        aggregate=True,
                                        columns=3):
    """
    fig: matplotlib figure
    axes: axes associated with that matplotlib figure
    method: string, either 'SVI_scaled', 'Factor 1', or 'Factor 2', what the x-axis is of final plot 
    annotations: bool, if True, plot annotations 
    weight: int, how much to multiply max value for no access parcels
    aggregate, bool, if True, include an aggregate plot
    columns, int, number of columns in final graph
    """
    column=0
    row=0

    # list of cost attributes to  consider:
    cost_atrs = ['cost_of_fl',
                'cost_of__1',
                'cost_of__2',
                'cost_of__3',
                'cost_of__4',
                'cost_of__5',
                'cost_of__6',
                'cost_of__7']

    # list of resource names
    res_names = ['Supermarket',
                'Emergency Room',
                'Pharmacy',
                'Police Station',
                'Convenience Store',
                'Fire Station',
                'Ambulance Station',
                'Gas Station',
                'Aggregate']

    # create blank array of all_arrays for the aggregate calcuation
    all_arrays = np.zeros_like(svi.shape[0], dtype='float')

    # Reproject cost of acces values into a "feasible set"
    # removing the impact of outliers by setting them all equal to the 3rd quantile plus 3*the IQR
    for cost_atr in cost_atrs:
        # # IQR METHOD to mask impact of outliers
        costs = sorted(res_parcels[cost_atr].tolist())
        costs = [x for x in costs if math.isnan(x)==False]
        q1,q3,=np.percentile(costs, [25,75])
        iqr=q3-q1
        upper_bound=q3+(3*iqr)
        res_parcels.loc[res_parcels[cost_atr] >= upper_bound, [cost_atr]] = upper_bound


    # BEGIN MAIN LOOP
    for idx, cost_atr in enumerate(cost_atrs):
        # Spatial join
        sjoined_data = gpd.sjoin(left_df=svi,
                                right_df=res_parcels,
                                how='left')

        # Fill nans with maximum from that column, and return back to GeoDataFrame
        filled_data = sjoined_data[cost_atr].replace(np.nan, sjoined_data[cost_atr].max()*weight)
        sjoined_data[cost_atr] = filled_data

        # count the number of residential parcels within each block group and relate back to  geodataframe
        count_dict = sjoined_data['GEOID'].value_counts().to_dict()
        sjoined_data["count"] = sjoined_data["GEOID"].apply(lambda x: count_dict.get(x))

        # Aggregate results
        summed_gdf = sjoined_data.dissolve(by='GEOID', aggfunc={cost_atr: 'sum', 
                                                                'SVI_scaled': 'first', 
                                                                'Factor 1': 'first',
                                                                'Factor 2': 'first',
                                                                'PPUNIT':'first',
                                                                'QBLACK':'first',
                                                                'QSPANISH': 'first',
                                                                'QPOVTY': 'first',
                                                                'count':'first'})

        # TODO: create ppunit AND/OR count weighted column
        summed_gdf[cost_atr] = summed_gdf[cost_atr]*summed_gdf["PPUNIT"]

        # SORT BY appropriate column
        if method == 'Factor 2':
            summed_gdf.sort_values(method, axis=0, inplace=True,ascending=False)
        else:
            summed_gdf.sort_values(method, axis=0, inplace=True)
        # summed_gdf.sort_values('SVI_scaled', axis=0, inplace=True)
        # summed_gdf.sort_values('Factor 1', axis=0, inplace=True)
        # summed_gdf.sort_values('Factor 2', axis=0, inplace=True,ascending=False)

        # TODO: HAVE TO MULTIPLY ARRAY by -1 IF FACTOR 2 DUE TO INVERSE WEALTH
        # single resource specific flood conditions
        if method == 'Factor 2':
            array = summed_gdf[cost_atr].values*-1
        else:
            array = summed_gdf[cost_atr].values

        # aggregate all arrays for combined lorenz
        all_arrays=np.add(array,all_arrays)
        
        lorenz_curve_multiplot(X=array,
                                row=row,
                                column=column,
                                axes=axes, 
                                idx=idx,
                                names=res_names,
                                method=method,
                               annotations=annotations)
        
        # update row and column names
        column += 1
        if column == columns:
            column = 0
            row += 1

    idx+=1
    if aggregate is True:
        # add the final aggregate plot, summing across all of the resources
        lorenz_curve_multiplot(X=all_arrays,
                                row=row,
                                column=column,
                                axes=axes, 
                                idx=idx,
                                names=res_names,
                                method=method,
                                annotations=annotations)


def lorenz_curve_resource_accessibility_time_series(fig, 
                                                    axes, 
                                                    method='SVI_scaled', 
                                                    annotations=True,
                                                    times=[],
                                                    weight=10):
    """
    fig: matplotlib figure
    axes: axes associated with that matplotlib figure
    method: string, either 'SVI_scaled', 'Factor 1', or 'Factor 2', what the x-axis is of final plot 
    annotations: bool, if True, plot annotations 
    times: list of times to refering to shapefiles to open
    weight, int, amount to multiply max value to fill no access parcels
    """

    for tstamp_index, tstamp in enumerate(times):
        res_parcels = gpd.read_file(
            f'/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs/Flow_decomp/{tstamp}_inundation_res_parcel_flow_decomp.shp')
        column = 0
        row = 0

        # list of cost attributes to  consider:
        cost_atrs = ['cost_of_fl',
                    'cost_of__1',
                    'cost_of__2',
                    'cost_of__3',
                        'cost_of__4',
                        'cost_of__5',
                        'cost_of__6',
                        'cost_of__7']

        # list of resource names
        res_names = ['Supermarket',
                    'Emergency Room',
                    'Pharmacy',
                        'Police Station',
                        'Convenience Store',
                        'Fire Station',
                        'Ambulance Station',
                        'Gas Station',
                        'Aggregate']

        # create blank array of all_arrays for the aggregate calcuation
        all_arrays = np.zeros_like(svi.shape[0], dtype='float')

        # Reproject cost of acces values into a "feasible set"
        # removing the impact of outliers by setting them all equal to the 3rd quantile plus 3*the IQR
        for cost_atr in cost_atrs:
                # # IQR METHOD to mask impact of outliers
                costs = sorted(res_parcels[cost_atr].tolist())
                costs = [x for x in costs if math.isnan(x) == False]
                q1, q3, = np.percentile(costs, [25, 75])
                iqr = q3-q1
                upper_bound = q3+(3*iqr)
                res_parcels.loc[res_parcels[cost_atr] >=
                                upper_bound, [cost_atr]] = upper_bound

        # BEGIN MAIN LOOP
        for idx, cost_atr in enumerate(cost_atrs):
            # Spatial join
            sjoined_data = gpd.sjoin(left_df=svi,
                                    right_df=res_parcels,
                                    how='left')

            # Fill nans with maximum from that column, and return back to GeoDataFrame
            filled_data = sjoined_data[cost_atr].replace(np.nan, sjoined_data[cost_atr].max()*weight)
            sjoined_data[cost_atr] = filled_data

            # count the number of residential parcels within each block group and relate back to  geodataframe
            count_dict = sjoined_data['GEOID'].value_counts().to_dict()
            sjoined_data["count"] = sjoined_data["GEOID"].apply(lambda x: count_dict.get(x))

            # Aggregate results
            summed_gdf = sjoined_data.dissolve(by='GEOID', aggfunc={cost_atr: 'sum',
                                                                    'SVI_scaled': 'first',
                                                                    'Factor 1': 'first',
                                                                    'Factor 2': 'first',
                                                                    'PPUNIT': 'first',
                                                                    'QBLACK': 'first',
                                                                    'QSPANISH': 'first',
                                                                    'QPOVTY': 'first',
                                                                    'count': 'first'})

            # TODO: create ppunit  weighted column
            summed_gdf[cost_atr] = summed_gdf[cost_atr]*summed_gdf["PPUNIT"]

            # SORT BY appropriate column
            if method == 'Factor 2':
                summed_gdf.sort_values(method, axis=0, inplace=True, ascending=False)
            else:
                summed_gdf.sort_values(method, axis=0, inplace=True)
            # summed_gdf.sort_values('SVI_scaled', axis=0, inplace=True)
            # summed_gdf.sort_values('Factor 1', axis=0, inplace=True)
            # summed_gdf.sort_values('Factor 2', axis=0, inplace=True,ascending=False)

            # TODO: HAVE TO MULTIPLY ARRAY by -1 IF FACTOR 2 DUE TO INVERSE WEALTH
            # single resource specific flood conditions
            if method == 'Factor 2':
                array = summed_gdf[cost_atr].values*-1
            else:
                array = summed_gdf[cost_atr].values

            # aggregate all arrays for combined lorenz
            all_arrays = np.add(array, all_arrays)

            lorenz_curve_multiplot(X=array,
                                row=row,
                                column=column,
                                axes=axes,
                                idx=idx,
                                names=res_names,
                                method=method,
                                annotations=annotations,
                                tstamp=tstamp_index)

            # update row and column names
            column += 1
            if column == 3:
                column = 0
                row += 1

        idx += 1
        # add the final aggregate plot, summing across all of the resources
        lorenz_curve_multiplot(X=all_arrays,
                            row=row,
                            column=column,
                            axes=axes,
                            idx=idx,
                            names=res_names,
                            method=method,
                            annotations=annotations,
                               tstamp=tstamp_index)


def lorenz_curve_resource_accessibility_aggregate_only(fig,
                                        axes,  
                                        weight=10):
    # plot the eq line
    # Equality line
    eq_line, = axes.plot([0, 1], [0, 1], color='k',
                         linestyle='--', label='Equality Line')
    
    # Set x and y tick markers
    axes.set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%','25%','50%','75%','100%'],fontsize=14)
    axes.set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'],fontsize=14)

    # set x and y limits
    axes.set_xlim(0,1)
    axes.set_ylim(0,1)


    # list of cost attributes to  consider:
    cost_atrs = ['cost_of_fl',
                 'cost_of__1',
                 'cost_of__2',
                 'cost_of__3',
                 'cost_of__4',
                 'cost_of__5',
                 'cost_of__6',
                 'cost_of__7']

    # list of resource names
    res_names = ['Supermarket',
                 'Emergency Room',
                 'Pharmacy',
                 'Police Station',
                 'Convenience Store',
                 'Fire Station',
                 'Ambulance Station',
                 'Gas Station',
                 'Aggregate']

    # Reproject cost of acces values into a "feasible set"
    # removing the impact of outliers by setting them all equal to the 3rd quantile plus 3*the IQR
    for cost_atr in cost_atrs:
        # # IQR METHOD to mask impact of outliers
        costs = sorted(res_parcels[cost_atr].tolist())
        costs = [x for x in costs if math.isnan(x) == False]
        q1, q3, = np.percentile(costs, [25, 75])
        iqr = q3-q1
        upper_bound = q3+(3*iqr)
        res_parcels.loc[res_parcels[cost_atr] >=
                        upper_bound, [cost_atr]] = upper_bound

    # BEGIN MAIN LOOP
    methods=['SVI_scaled','Factor 1', 'Factor 2']
    table_vals=[[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    ginis=[0,0,0]
    for b, method in enumerate(methods):
        # create blank array of all_arrays for the aggregate calcuation
        all_arrays = np.zeros_like(svi.shape[0], dtype='float')     
        for idx, cost_atr in enumerate(cost_atrs):
            # Spatial join
            sjoined_data = gpd.sjoin(left_df=svi,
                                    right_df=res_parcels,
                                    how='left')

            # Fill nans with maximum from that column, and return back to GeoDataFrame
            filled_data = sjoined_data[cost_atr].replace(
                np.nan, sjoined_data[cost_atr].max()*weight)
            sjoined_data[cost_atr] = filled_data

            # count the number of residential parcels within each block group and relate back to  geodataframe
            count_dict = sjoined_data['GEOID'].value_counts().to_dict()
            sjoined_data["count"] = sjoined_data["GEOID"].apply(
                lambda x: count_dict.get(x))

            # Aggregate results
            summed_gdf = sjoined_data.dissolve(by='GEOID', aggfunc={cost_atr: 'sum',
                                                                    'SVI_scaled': 'first',
                                                                    'Factor 1': 'first',
                                                                    'Factor 2': 'first',
                                                                    'PPUNIT': 'first',
                                                                    'QBLACK': 'first',
                                                                    'QSPANISH': 'first',
                                                                    'QPOVTY': 'first',
                                                                    'count': 'first'})

            # TODO: create ppunit AND/OR count weighted column
            summed_gdf[cost_atr] = summed_gdf[cost_atr]*summed_gdf["PPUNIT"]

            # SORT BY appropriate column
            if method == 'Factor 2':
                summed_gdf.sort_values(
                    method, axis=0, inplace=True, ascending=False)
            else:
                summed_gdf.sort_values(method, axis=0, inplace=True)

            # TODO: HAVE TO MULTIPLY ARRAY by -1 IF FACTOR 2 DUE TO INVERSE WEALTH
            # single resource specific flood conditions
            if method == 'Factor 2':
                array = summed_gdf[cost_atr].values*-1
            else:
                array = summed_gdf[cost_atr].values

            # aggregate all arrays for combined lorenz
            all_arrays = np.add(array, all_arrays)

        # Plot the aggregate
        lorenz = all_arrays.cumsum() / all_arrays.sum()
        lorenz = np.insert(lorenz, 0, 0)
        vals = np.arange(lorenz.size)/(lorenz.size-1)

        if method == 'SVI_scaled':
            c='red'
        elif method == 'Factor 1':
            c='blue'
        elif method == 'Factor 2':
            c='green'

        # calc gini and save
        ginis[b] = gini(all_arrays)

        # plot lorenz
        lorenz_plot, = axes.plot(vals, lorenz, label='tsamp', color=c,linewidth=2)

        perc_25 = np.percentile(lorenz, 25)
        perc_50 = np.percentile(lorenz, 50)
        perc_75=np.percentile(lorenz,75)

        table_vals[b][0]=f'{round(perc_25*100)}%'
        table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
        table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
        table_vals[b][3]=f'{100-round(perc_75*100)}%'

          
    col_labels = ['Q1', 'Q2', 'Q3' ,'Q4']
    row_labels = ['SVI', 'SS', 'ES']
    # the rectangle is where I want to place the table
    the_table = axes.table(cellText=table_vals,
                        colWidths=[0.10]*4,
        rowColours=['red', 'blue', 'green'],
                        rowLabels=row_labels,
                        colLabels=col_labels,
        colColours=['lightgrey', 'lightgrey', 'lightgrey','lightgrey'],
                        loc='lower right')
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(14)
    the_table.scale(1,1.5)
    the_table[1,-1].get_text().set_color('white')
    the_table[2,-1].get_text().set_color('white')
    the_table[3,-1].get_text().set_color('white')

    # plot ginis
    # annotation box text
    textstr = '\n'.join((f'SVI Gini: {round(ginis[0],2)}',
                         f'SS Gini: {round(ginis[1],2)}',
                         f'ES Gini: {round(ginis[2],2)}'))
    # patch properties
    props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
    # text box in upper left in axes coords
    axes.text(0.05, 0.95, textstr, transform=axes.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    


#
##########################################################################################################################################
# SINGLE PLOT OF ALL INDEPENDENT VARIABLES OF AGGREGATE TRAVEL TIMES
# fig, axes = plt.subplots(1,1,figsize=(7,7))
# fig.set_facecolor('none')
# axes.set_title('Lorenz Curves of Residential Accessability',
#              fontsize=18, weight='bold',style='italic')
# axes.set_xlabel('Indicator Percentile',
#               weight='bold',
#               x=0.5,
#               fontsize=16)
# axes.set_ylabel('Normalized Cumulative Sum of Travel Times',
#               weight='bold',
#               fontsize=16)
# lorenz_curve_resource_accessibility_aggregate_only(fig=fig, axes=axes, weight=10)
# plt.tight_layout()
# plt.show()

########################################################################################################################################
# 3X3 PLOT OF ALL INDEPENDENT VARIABLES
# fig, axes = plt.subplots(3, 3, figsize=(10, 10))
# fig.suptitle('Lorenz Curves of Residential Access to Critical Resources',
#              fontsize=14, weight='bold')
# fig.supxlabel('Social Vulnerability Index Percentile',
#               weight='bold',
#               x=0.5)
# fig.supylabel('Normalized Cumulative Sum of PPUNIT Weighted Travel Times',
#               weight='bold')
# lorenz_curve_resource_accessibility(fig=fig, axes=axes, method='SVI_scaled', annotations=True)
# lorenz_curve_resource_accessibility(fig=fig, axes=axes, method='Factor 1', annotations=True)
# lorenz_curve_resource_accessibility(fig=fig, axes=axes, method='Factor 2', annotations=True)
# add custom legend
# legend_elements = [Line2D([0], [0], color='red', lw=4, label='SVI'),
#                    Line2D([0], [0], color='blue', lw=4, label='Social Status Factor'),
#                    Line2D([0], [0], color='green', lw=4, label='Economic  Factor')]
# plt.figlegend(handles=legend_elements, loc='lower right',ncol=3)

#########################################################################################################################################
# 4X2 SVI PLOT
# fig, axes = plt.subplots(2,4, figsize=(16, 8))
# fig.set_facecolor('none')
# # fig.suptitle('Lorenz Curves of Residential Access to Critical Resources',
# #              fontsize=14, weight='bold')
# fig.supxlabel('Social Vulnerability Index Percentile',
#               weight='bold',
#               x=0.5,
#               fontsize=18)
# fig.supylabel('Normalized Cumulative Sum of PPUNIT Weighted Travel Times',
#               weight='bold',
#               fontsize=18)
# lorenz_curve_resource_accessibility(fig=fig,
#                                     axes=axes,
#                                     method='SVI_scaled',
#                                     annotations=True,
#                                     weight=10,
#                                     aggregate=False,
#                                     columns=4)

# plt.tight_layout()
# plt.subplots_adjust(top=0.925,
#                     bottom=0.075,
#                     left=0.09,
#                     right=0.975,
#                     hspace=0.25,
#                     wspace=0.2)
# plt.show()





##############################################
# PLOT EXAMPLE OF LORENZ CURVE

fig, ax = plt.subplots(figsize=[5, 5])
fig.set_facecolor('none')
# fig.suptitle('Example Lorenz Curve',
#              fontsize=16, weight='bold')
fig.supxlabel('Social Vulnerability Index Percentile',
              weight='bold',
              x=0.5,
              fontsize=16)
fig.supylabel('Normalized Cumulative Sum of Variable',
              weight='bold',
              fontsize=16)
# Equality line
eq_line, = ax.plot([0, 1], [0, 1], color='k',
                    linestyle='--', label='Equality Line')


def hanging_line(point1, point2):
    import numpy as np

    a = (point2[1] - point1[1])/(np.cosh(point2[0]) - np.cosh(point1[0]))
    b = point1[1] - a*np.cosh(point1[0])
    x = np.linspace(point1[0], point2[0], 100)
    y = a*np.cosh(x) + b

    return (x, y)

point1 = [0, 0]
point2 = [1, 1]
x, y = hanging_line(point1, point2)
ax.plot(x, y, color='red',linewidth=2, label='Lorenz Curve')

# Set x and y tick markers
ax.set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=['0%','25%','50%','75%','100%'], fontsize=14)
ax.set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
              '0%', '25%', '50%', '75%', '100%'], fontsize=14)

# set x and y limits
ax.set_xlim(0,1)
ax.set_ylim(0,1)

# add annotation box
# annotation value
val1 = gini(x)
# annotation box text
textstr = '\n'.join((f'Gini: {round(val1,4)}',))
# patch properties
props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
# text box in upper left in axes coords
ax.text(0.025, 0.975, 
                    textstr,
                    transform=ax.transAxes, 
                    fontsize=14,
                    verticalalignment='top', 
                    bbox=props)

# add proportion of weight annotations
# find the 25th and 75th percentile
perc_25=np.percentile(y,25)
perc_50=np.percentile(y,50)
perc_75=np.percentile(y,75)

# plot vertical and horizontal lines for these locations
perc_25_line_v, = ax.plot([0.25, 0.25], [0, perc_25], color='dimgray', linestyle='dotted')
perc_25_line_h, = ax.plot([0, 0.25], [perc_25, perc_25], color='dimgray',linestyle='dotted')
perc_50_line_v, = ax.plot([0.50, 0.50], [0, perc_50], color='dimgray',linestyle='dotted')
perc_50_line_h, = ax.plot([0, 0.50], [perc_50, perc_50],color='dimgray', linestyle='dotted')
perc_75_line_v, = ax.plot([0.75, 0.75], [0, perc_75], color='dimgray',linestyle='dotted')
perc_75_line_h, = ax.plot([0, 0.75], [perc_75, perc_75],color='dimgray', linestyle='dotted')

# Plot double arrow annotations
arrow1 = patches.FancyArrowPatch((0.025, 0), (0.025, perc_25), 
                                arrowstyle='<->', 
                                mutation_scale=10, 
                                color='dimgray')
arrow2 = patches.FancyArrowPatch((0.2, perc_25), (0.2, perc_50),
                                arrowstyle='<->',
                                mutation_scale=10,
                                color='dimgray')
arrow3 = patches.FancyArrowPatch((0.45, perc_50), (0.45, perc_75),
                                 arrowstyle='<->',
                                 mutation_scale=10,
                                 color='dimgray')
arrow4 = patches.FancyArrowPatch((0.65, perc_75), (0.65, 1), 
                                arrowstyle='<->', 
                                mutation_scale=10,
                                color='dimgray')
ax.add_patch(arrow1)
ax.add_patch(arrow2)
ax.add_patch(arrow3)
ax.add_patch(arrow4)

# Plot percentage values
ax.annotate(text=f'{round(perc_25*100)}%',xy=(0.06,perc_25/2-0.01),fontsize=14)
ax.annotate(text=f'{round((perc_50-perc_25)*100)}%',xy=(0.2,(perc_50-perc_25)/2+perc_25),fontsize=14)
ax.annotate(text=f'{round((perc_75-perc_50)*100)}%',xy=(0.46, (perc_75-perc_50)/2+perc_50),fontsize=14)
ax.annotate(text=f'{round((1-perc_75)*100)}%', xy=(0.67, (1-perc_75)/2+perc_75),fontsize=14)

# add table
col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
row_labels = ['Lorenz']
table_vals = [[f'{round(perc_25*100)}%', 
               f'{round((perc_50-perc_25)*100)}%',
               f'{round((perc_75-perc_50)*100)}%',
               f'{round((1-perc_75)*100)}%']]

the_table = ax.table(cellText=table_vals,colWidths=[0.13]*4,
                     rowColours=['r'],
                     rowLabels=row_labels,
                     colLabels=col_labels, 
                     colColours=['lightgray', 'lightgray',
                                 'lightgray', 'lightgray'],
                     loc='lower right')
the_table.auto_set_font_size(False)
the_table.set_fontsize(14)
the_table.set_zorder(100)
the_table.scale(1, 1.5)

plt.tight_layout()
plt.show()
