# This contains lorenz curve functions to measure equity of different resiliency functions
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

# set font family
plt.rcParams['font.family'] = "ubuntu"


# filepath variables
network_data_fp = '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North'
svi_data_fp = '/home/mdp0023/Documents/Codes_Projects/SVI_Code/Travis_County'
graphs_fp = '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs'

# Load variables
# SVI Shapefile: specific to 2015, relevant to the flood
svi = gpd.read_file(f'{svi_data_fp}/SVI_Shapefiles/Travis_county_svi_2015_selected.shp')

# read in background network, convert eges to shapefile 
# BUG: has to be no TA network
# network = mynet.read_graph_from_disk(path=f'{network_data_fp}/AN_Graphs',
#                                      name='AN_Graph_2015052522_inundation')
network = mynet.read_graph_from_disk(path=f'{network_data_fp}/AN_Graphs',
                                     name='AN_Graph_2015052518_inundation')
gdf_edges = ox.graph_to_gdfs(G=network, nodes=False)

# Read in the TA parcels
# res_parcels = gpd.read_file(f'{network_data_fp}/AN_Graphs/Flow_decomp/2015052522_inundation_res_parcel_flow_decomp.shp')
res_parcels = gpd.read_file(f'{network_data_fp}/AN_Graphs/Flow_decomp/2015052518_inundation_res_parcel_flow_decomp.shp')

# Read in 2015 specific socio-demographic variables, grab relevent variables data and merge with SVI
vars2015 = pd.read_csv(f'{svi_data_fp}/Calculated_Variables/Travis_county_var_2015.csv')
vars2015=vars2015[['PPUNIT','GEOID']].copy()
svi=svi.merge(vars2015,left_on='GEOID', right_on='GEOID')

# figure size converstion
mm=1/25.4
##############################################################################################################################################
# FUNCTIONS
# Gini coefficent
def gini(arr):
    """
    Calculates the Gini coefficient.
    
    A Gini coefficeint is the measure deviation the lorenz curve has from a perfectly equitable (1:1) relationship

    :param arr: list of values, sorted from lowest to highest SVI (or whatever independent variable is being analyzed)
    :type arr: np.array
    :returns:
        - **gini**, Value of deviation between lorenz curve and line of perfect equity 
    :rtype: float
    """
    # arr should already be sorted how we want it 
    sorted_arr = arr.copy()
    n = arr.size
    coef_ = 2. / n
    const_ = (n + 1.) / n
    weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
    return coef_*weighted_sum/(sorted_arr.sum()) - const_

# example lorenz curve for presentations
def example_lorenz():
    """
    Plots a descriptive example of a lorenz curve with appropraite annotations. This isn't meant to plot any real data, just a presentation quality example.
    
    returns: matplotlib fig
    """
    # set figure titles
    fig, ax = plt.subplots(figsize=[5, 5])
    fig.set_facecolor('none')
    fig.supxlabel('Social Vulnerability Index Percentile',
                  weight='bold',
                  x=0.5,
                  fontsize=16)
    fig.supylabel('Normalized Cumulative Sum of Variable',
                  weight='bold',
                  fontsize=16)

    # Plot perfect equity line
    eq_line, = ax.plot([0, 1], [0, 1], color='k',
                       linestyle='--', label='Equality Line')

    # Set x and y tick markers
    ax.set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                  '0%', '25%', '50%', '75%', '100%'], fontsize=14)
    ax.set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
        '0%', '25%', '50%', '75%', '100%'], fontsize=14)

    # set x and y limits
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # function to plot curve
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
    ax.plot(x, y, color='red', linewidth=2, label='Lorenz Curve')

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
    perc_25 = np.percentile(y, 25)
    perc_50 = np.percentile(y, 50)
    perc_75 = np.percentile(y, 75)

    # plot vertical and horizontal lines for these locations
    perc_25_line_v, = ax.plot(
        [0.25, 0.25], [0, perc_25], color='dimgray', linestyle='dotted')
    perc_25_line_h, = ax.plot(
        [0, 0.25], [perc_25, perc_25], color='dimgray', linestyle='dotted')
    perc_50_line_v, = ax.plot(
        [0.50, 0.50], [0, perc_50], color='dimgray', linestyle='dotted')
    perc_50_line_h, = ax.plot(
        [0, 0.50], [perc_50, perc_50], color='dimgray', linestyle='dotted')
    perc_75_line_v, = ax.plot(
        [0.75, 0.75], [0, perc_75], color='dimgray', linestyle='dotted')
    perc_75_line_h, = ax.plot(
        [0, 0.75], [perc_75, perc_75], color='dimgray', linestyle='dotted')

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
    ax.annotate(text=f'{round(perc_25*100)}%',
                xy=(0.06, perc_25/2-0.01), fontsize=14)
    ax.annotate(text=f'{round((perc_50-perc_25)*100)}%',
                xy=(0.2, (perc_50-perc_25)/2+perc_25), fontsize=14)
    ax.annotate(text=f'{round((perc_75-perc_50)*100)}%',
                xy=(0.46, (perc_75-perc_50)/2+perc_50), fontsize=14)
    ax.annotate(text=f'{round((1-perc_75)*100)}%',
                xy=(0.67, (1-perc_75)/2+perc_75), fontsize=14)

    # add table
    col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
    row_labels = ['Lorenz']
    table_vals = [[f'{round(perc_25*100)}%',
                   f'{round((perc_50-perc_25)*100)}%',
                   f'{round((perc_75-perc_50)*100)}%',
                   f'{round((1-perc_75)*100)}%']]

    the_table = ax.table(cellText=table_vals, colWidths=[0.13]*4,
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
    return fig

# Single basic Lorenz curve: isn't actually called anywhere, either remove, or modify other functions to call this 
def lorenz_curve(X):
    """
    Creates a Lorenz curve plot. 

    Creates Lorenz curve plot based on given datasets of ordered variables.


    :param X: list of values, sorted from lowest to highest SVI (or whatever independent variable is being analyzed)
    :type X: np.array
    :returns:
        - **lorenz plot**, figure showing equitable/inequitable distribution of a variable across a population
    :rtype: matplotlib fig

    """
    lorenz = X.cumsum() / X.sum()
    lorenz = np.insert(lorenz, 0, 0)
    vals = np.arange(lorenz.size)/(lorenz.size-1)

    fig, ax = plt.subplots(figsize=[6, 6])
    # Equality line
    eq_line, = ax.plot([0, 1], [0, 1], color='k',
                       linestyle='--', label='Equality Line')

    lorenz_plot, = ax.plot(vals, lorenz, label='tsamp', color='red')

    plt.title("Lorenz Curve")
    ax.set_xlabel("Independent Variable (i.e., ordered SVI")
    ax.set_ylabel("Impact Variable")

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

# Lorenz curve multiplot: also not called??? 
def lorenz_curve_multiplot(X, row, column, axes, idx, names, method, annotations, tstamp=False):
    """
    Manipulates existing instances of a matplotlib figure and axes to plot multiple lorenz curves on a single figure.

    Each time this function is called, it manipulates a different subplot within an existing instance of a figure created from matplotlib.subplots
    
    :param X: list of values, sorted from lowest to highest SVI (or whatever independent variable is being analyzed)
    :type X: np.array
    :param row: row number within existing figure
    :type row: int
    :param col: column number within existing figure
    :type col: int
    :param axes: axes within the figure being manipulated
    :type axes: array of matplotlib.axes
    :param idx: number used to call elements from names, the "plot number" within the figure
    :type idx: int
    :param names: the name that gets plotted with each subplot (e.g., the resource names)
    :type names": list of strings
    :param method: the x-axis variable, one of the following: 'SVI_scaled', 'Factor 1', 'Factor 2', determines plotting color
    :type method: string
    :param tstamp: if True, plots time series data for each subplot, and adjusts colors accordingly
    :type tstamp: bool
    
    inputs:
    X: np.array, in the order of lowest to highest SVI
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

# Network reliability, impacted roads Lorenz curve
def network_reliability(axes,
                        fpath,
                        prefix,
                        svi,
                        method='SVI_scaled',
                        variable='agr_no_cap',
                        times=[],
                        annotations=True,
                        columns=3,
                        ftypes=['inundation']):
    """
    Function to plot lorenz curve of impact on roads with options to show compound, fluvial, and/or pluvial impacts in time series or not.

    This should automate the creation of network reliability lorenz curves, with a highly customizable plot tool. Some important usage notes:

    fpath and prefix:
        - used to find the graphs that need to be analyzed
        - graphs should have the following naming convention: {fpath}/{prefix}_{time}_{ftype}, matching the 'times' and 'ftypes' inputs
    
    method:
        -the x-axis variable the Lorenz curve is plotted against
        -acceptable values: 'SVI_scaled', 'Factor 1', or 'Factor 2' (refers to attributes of svi geodataframe)
        -Currently, Factor 2 is inversed, because it is the economic factor (inversely related to vulnerability)
        -TODO: should be a seperate function method to allow input on whether variable is inversed or not

    variable:
        - The y-axis varaible being plotted. The acceptable values and their meaning are as follows:
            -'agr_no_cap': aggressive inundation, number of roads with no capacity per block group
            -'perc_no_cap': aggressive inundation, percentage of roads with no capacity per block group
            -'agr_increased_travel_time': aggresive inundation, number of roads with an increased travel time per block group
            -'perc_increased_travel_time': aggressive inundation, percentage of roads with an increased travel time per block group
            -'agr_impact': aggressive inundation, number of roads impacted (increased travel time & no capacity) per block group
            -'perc_perc_impact':aggressive inundation, percentage of roads impacted (increased travel time & no capacity per block group)
        -default: 'agr_no_cap'
        - regardless of the variable, it is also weighted by PPUNIT (people per unit) Census data, to account for population density differences

    times:
        - should be in the format YYYYMMDDHH
        - used to plot subplot titles as well as point to which network needs to be loaded

    columns:
        - the number of columns in the figure, used to iterate through plotting

    ftypes:
        - Which floods should be considered, list of strings
        - Acceptable values are 'inundation' for compound, 'inundation_fluvial' for fluvial, and 'inundation_pluvial' for pluvial


    :param axes: axes of figure being manipulated. If only plotting a single instance, put axes in [[]]
    :type axes: array of matplotlib.axes
    :param fpath: the folder path where the graphs are stored
    :type fpath: string
    :param prefix: the graph network prefixes
    :type prefix: string
    :param method: the x-axis variable
    :type method: string
    :param variable: y-axis variable
    :type variable: string
    :param times: the times associated with each network to be loaded and analyzed 
    :type times: list of strings
    :param annotations: If True, plots annotation box containing quartile burden information and Gini coefficients
    :type annotations: bool
    :param columns: number of columns in figure being manipulated
    :type columns: int
    :param ftypes: which flood types should be plotted 
    :type ftypes: list of strings
    """
    column = 0
    row = 0

    # read in appropriate network and convert to geopandas dataframe
    for a, time in enumerate(times):
        # create variable to hold information for percentile table
        table_vals = [[0, 0, 0, 0] for _ in range(len(ftypes))]
        # create a variable to hold information regarding gini coefficeints
        ginis = [0 for _ in range(len(ftypes))]
        # for each flood type of interest
        for b, ftype in enumerate(ftypes):

            network = mynet.read_graph_from_disk(path=fpath,
                                                 name=f'{prefix}_{time}_{ftype}')
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
            svi['agr_increased_travel_time'] = svi['agr_increased_travel_time'].fillna(
                0)

            # count the number of roads impacted (increased travel time & 0 capacity)
            svi['agr_impact'] = svi['agr_increased_travel_time'] + svi['agr_no_cap']

            # calc percentage of roads within a BG with no capacity
            svi['perc_no_cap'] = svi['agr_no_cap']/svi['count']*100

            # calc percentage of roads within a BG with an increased travel time
            svi['perc_increased_travel_time'] = svi['agr_increased_travel_time'] / \
                svi['count']*100

            # calc percentage of roads within a BG impacted
            svi['perc_impact'] = svi['agr_impact']/svi['count']*100

            # sort svi dataframe by appropriate column
            if method == 'Factor 2':
                svi.sort_values('Factor 2', axis=0,
                                inplace=True, ascending=False)
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
            lorenz_plot, = axes[row][column].plot(val, lorenz, label='tsamp', color=c_line, linewidth=2)

            perc_25 = np.percentile(lorenz, 25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz, 75)

            table_vals[b][0] = f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3] = f'{100-round(perc_75*100)}%'

            # calc gini and save
            ginis[b] = gini(array)

        col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
        # determine appropriate row labels
        row_labels = []
        for ftype in ftypes:
            if ftype == 'inundation':
                row_labels.append('Cmpd.')
            elif ftype == 'inundation_fluvial':
                row_labels.append('Fluvial')
            elif ftype == 'inundation_pluvial':
                row_labels.append('Pluvial')

        # determine appropriate row colors
        row_colors=[]
        for ftype in ftypes:
            if ftype == 'inundation':
                row_colors.append('#cc57a4')
            elif ftype == 'inundation_fluvial':
                row_colors.append('#216aa6')
            elif ftype == 'inundation_pluvial':
                row_colors.append('#ff772e')

        # the rectangle is where I want to place the table
        if annotations is True:

            # plot ginis
            textstr = []
            for i, row_label in enumerate(row_labels):
                textstr.append(f'{row_label} G: {round(ginis[i],2)}')
            textstr = '\n'.join(textstr)

            # Gini patch properties
            props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
            # text box in upper left in axes coords
            axes[row][column].text(0.05, 0.95, textstr, transform=axes[row][column].transAxes, fontsize=8,
                                   verticalalignment='top', bbox=props)

            # Plot quartile information
            the_table = axes[row][column].table(cellText=table_vals,
                                                # colWidths=[0.10]*4,
                                                rowColours=row_colors,
                                                colColours=['lightgrey', 'lightgrey',
                                                            'lightgrey', 'lightgrey'],
                                                rowLabels=row_labels,
                                                colLabels=col_labels,
                                                loc='bottom')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(8)
            the_table.scale(1, 0.8)

        # Equality line
        eq_line, = axes[row][column].plot(
            [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[' ', ' ', ' ', ' ', ' '])
                                    #labels=[' ', '25%', '50%', '75%', '100%'])
        axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0],
                                      labels=['0%', '25%', '50%', '75%', '100%'])

        # adjust ticks
        axes[row][column].tick_params(axis="x", direction="in")
        
        # set x and y limits
        axes[row][column].set_xlim(0, 1)
        axes[row][column].set_ylim(0, 1)

        # set axes title
        axes[row][column].set_title(
            f'{time[-2:]}:00 {time[6:8]}/{time[4:6]}/{time[:4]}', style='italic')

        # update row and column names
        column += 1
        if column == columns:
            column = 0
            row += 1

# Individual reliability, cost to access resources Lorenz curve
def individual_reliability(axes,
                           cost_atrs,
                           res_names,
                           svi,
                           res_parcels,
                           methods=['SVI_scaled'],
                           weight=4,
                           aggregate=True,
                           annotations=True,
                           columns=3):
    """
    Function to plot Lorenz curves of individual reliability for resource accessibility, analyzing the cost to access resources. Some important usage notes:

    cost_atrs:
        - this is a list of strings associated with the res_parcels geodataframe attributes that are the variables of interest. They should be the attributes that represent the cost to access each resource.
        - each variable is weighted by the PPUNIT Census variable to account for population density differences

    res_names:
        - These are the real resource names associated with each cost attribute variable in the same order. 

    methods:
        -the x-axis variable the Lorenz curve is plotted against
        -acceptable values: 'SVI_scaled', 'Factor 1', and/or 'Factor 2' (refers to attributes of svi geodataframe)
        -Currently, Factor 2 is inversed, because it is the economic factor (inversely related to vulnerability)
        -TODO: should be a seperate function method to allow input on whether variable is inversed or not

    outlier masking:
        - Due to the occurance of severe outliers occuring in travel times, (e.g., a residential parcel located at the end of a long roadway that has an extremely high travel time), the cost to access resources are "projected" onto a feasible set. This is done by calculating the interquartile range (IQR). A value is considered an outlier if it is equal to or greater than 3 times the IQR. These values are set to be equal to 3*IQR.


    :param axes: axes of figure being manipulated
    :type axes: array of matplotlib.axes
    :param cost_atrs: the attributes within the res_parcels geodataframe that need to be considered
    :type cost_atrs: list of strings
    :param res_names: the names associated with each element in cost_atrs (i.e., the resources actual name)
    :type res_names: list of strings
    :param svi: the social vulnerabilty dataset to aggregate results to
    :type svi: geodataframe 
    :param res_parcels: the residential parcels with flow decomposition information (i.e., the cost to get to each resource type)
    :type res_parcels: geodataframe
    :param methods: the x-axis variable
    :type methods: list of strings
    :param weight: weight to apply to parcels that cannot access the resource. The value will be set to the maximum travel time for that resource times weight
    :type weight: int
    :param aggregate: if True, also plots an aggregate Lorenz curve. 
    :type aggregate: bool
    :param annotations: If True, plots annotation box containing quartile burden information and Gini coefficients
    :type annotations: bool
    :param columns: number of columns in figure being manipulated
    :type columns: int
    
    """
    column = 0
    row = 0

    # create blank array of all_arrays for the aggregate calcuation
    if aggregate is True:
        all_arrays = [np.zeros_like(svi.shape[0], dtype='float') for _ in range(len(methods))]

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
    # for each resource,
    for idx, cost_atr in enumerate(cost_atrs):
        # create variable to hold information for percentile table
        table_vals = [[0, 0, 0, 0] for _ in range(len(methods))]
        # create a variable to hold information regarding gini coefficeints
        ginis = [0 for _ in range(len(methods))]
        
        # for each method,
        for b, method in enumerate(methods):
            
            # Spatial join
            sjoined_data = gpd.sjoin(left_df=svi,
                                    right_df=res_parcels,
                                    how='left')

            # Fill nans with maximum from that column, and return back to GeoDataFrame
            filled_data = sjoined_data[cost_atr].replace(np.nan, sjoined_data[cost_atr].max()*weight)
            sjoined_data[cost_atr] = filled_data

            # count the number of residential parcels within each block group and relate back to geodataframe
            count_dict = sjoined_data['GEOID'].value_counts().to_dict()
            sjoined_data["count"] = sjoined_data["GEOID"].apply(
                lambda x: count_dict.get(x))

            # Aggregate results
            summed_gdf = sjoined_data.dissolve(by='GEOID', aggfunc={cost_atr: 'sum',
                                                                    'SVI_scaled': 'first',
                                                                    'Factor 1': 'first',
                                                                    'Factor 2': 'first',
                                                                    'PPUNIT': 'first',
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
            if aggregate is True:
                all_arrays[b] = np.add(array, all_arrays[b])

            # calculate lorenz values
            lorenz = array.cumsum() / array.sum()
            lorenz = np.insert(lorenz, 0, 0)
            val = np.arange(lorenz.size)/(lorenz.size-1)

            if method == 'SVI_scaled':
                c_line = 'r'
            elif method == 'Factor 1':
                c_line = 'b'
            elif method == 'Factor 2':
                c_line = 'g'

            # lorenze curve
            lorenz_plot, = axes[row][column].plot(
                val, lorenz, label='tsamp', color=c_line, linewidth=2)

            perc_25 = np.percentile(lorenz, 25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz, 75)

            table_vals[b][0] = f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3] = f'{100-round(perc_75*100)}%'

            # calc gini and save
            ginis[b] = gini(array)

        col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
        # determine appropriate row labels
        row_labels = []
        for method in methods:
            if method == 'SVI_scaled':
                row_labels.append('SVI')
            elif method == 'Factor 1':
                row_labels.append('SS')
            elif method == 'Factor 2':
                row_labels.append('ES')

        # determine appropriate row colors
        row_colors=[]
        for method in methods:
            if method == 'SVI_scaled':
                row_colors.append('r')
            elif method == 'Factor 1':
                row_colors.append('b')
            elif method == 'Factor 2':
                row_colors.append('g')

        # the rectangle is where I want to place the table
        if annotations is True:

            # plot ginis
            textstr = []
            for i, row_label in enumerate(row_labels):
                textstr.append(f'{row_label} Gini: {round(ginis[i],2)}')
            textstr = '\n'.join(textstr)

            # Gini patch properties
            props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
            # text box in upper left in axes coords
            axes[row][column].text(0.05, 0.95, textstr, transform=axes[row][column].transAxes, fontsize=10,
                                   verticalalignment='top', bbox=props)

            # Plot quartile information
            the_table = axes[row][column].table(cellText=table_vals,
                                                colWidths=[0.10]*4,
                                                rowColours=row_colors,
                                                colColours=['lightgrey', 'lightgrey',
                                                            'lightgrey', 'lightgrey'],
                                                rowLabels=row_labels,
                                                colLabels=col_labels,
                                                loc='lower right')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(10)
            the_table.scale(1, 1)
            for i, val in enumerate(methods):
                the_table[i+1,-1].get_text().set_color('white')


        # Equality line
        eq_line, = axes[row][column].plot(
            [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                                     '0%', '25%', '50%', '75%', '100%'])
        axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                                     '0%', '25%', '50%', '75%', '100%'])

        # set x and y limits
        axes[row][column].set_xlim(0, 1)
        axes[row][column].set_ylim(0, 1)

        # set axes title
        axes[row][column].set_title(res_names[idx],style='italic')

        # update row and column names
        column += 1
        if column == columns:
            column = 0
            row += 1

    if aggregate is True:
        # create variable to hold information for percentile table
        table_vals = [[0, 0, 0, 0] for _ in range(len(methods))]
        # create a variable to hold information regarding gini coefficeints
        ginis = [0 for _ in range(len(methods))]
        for b, method in enumerate(methods):
            # calculate lorenz values
            lorenz = all_arrays[b].cumsum() / all_arrays[b].sum()
            lorenz = np.insert(lorenz, 0, 0)
            val = np.arange(lorenz.size)/(lorenz.size-1)

            if method == 'SVI_scaled':
                c_line = 'r'
            elif method == 'Factor 1':
                c_line = 'b'
            elif method == 'Factor 2':
                c_line = 'g'

            # lorenze curve
            lorenz_plot, = axes[row][column].plot(
                val, lorenz, label='tsamp', color=c_line, linewidth=2)

            perc_25 = np.percentile(lorenz, 25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz, 75)

            table_vals[b][0] = f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3] = f'{100-round(perc_75*100)}%'

            # calc gini and save
            ginis[b] = gini(all_arrays[b])

        col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
        # determine appropriate row labels
        row_labels = []
        for method in methods:
            if method == 'SVI_scaled':
                row_labels.append('SVI')
            elif method == 'Factor 1':
                row_labels.append('SS')
            elif method == 'Factor 2':
                row_labels.append('ES')

        # determine appropriate row colors
        row_colors = []
        for method in methods:
            if method == 'SVI_scaled':
                row_colors.append('r')
            elif method == 'Factor 1':
                row_colors.append('b')
            elif method == 'Factor 2':
                row_colors.append('g')

        # the rectangle is where I want to place the table
        if annotations is True:

            # plot ginis
            textstr = []
            for i, row_label in enumerate(row_labels):
                textstr.append(f'{row_label} Gini: {round(ginis[i],2)}')
            textstr = '\n'.join(textstr)

            # Gini patch properties
            props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
            # text box in upper left in axes coords
            axes[row][column].text(0.05, 0.95, textstr, transform=axes[row][column].transAxes, fontsize=10,
                                    verticalalignment='top', bbox=props)

            # Plot quartile information
            the_table = axes[row][column].table(cellText=table_vals,
                                                colWidths=[0.10]*4,
                                                rowColours=row_colors,
                                                colColours=['lightgrey', 'lightgrey',
                                                            'lightgrey', 'lightgrey'],
                                                rowLabels=row_labels,
                                                colLabels=col_labels,
                                                loc='lower right')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(10)
            the_table.scale(1, 1)
            for i, val in enumerate(methods):
                the_table[i+1,-1].get_text().set_color('white')
           

        # Equality line
        eq_line, = axes[row][column].plot(
            [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                                        '0%', '25%', '50%', '75%', '100%'])
        axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                                        '0%', '25%', '50%', '75%', '100%'])

        # set x and y limits
        axes[row][column].set_xlim(0, 1)
        axes[row][column].set_ylim(0, 1)

        # set axes title
        axes[row][column].set_title('Aggregate', style='italic')

# Individual reliability, cost to access resources Lorenz curve, only plot aggregate
def individual_reliability_aggregate(axes,
                                     prefix,
                                     fpath,
                                     times,
                                     columns=1,
                                     methods=['SVI_scaled'],
                                     svi_weight=None):
    """
    Function to plot Lorenz curve of individual reliability for resource accessibility, analyzing the cost to access resources, but only plotting the aggregate results.

    cost_atrs:
        - this is a list of strings associated with the res_parcels geodataframe attributes that are the variables of interest. They should be the attributes that represent the cost to access each resource.
        - each variable is weighted by the PPUNIT Census variable to account for population density differences

    methods:
        -the x-axis variable the Lorenz curve is plotted against
        -acceptable values: 'SVI_scaled', 'Factor 1', and/or 'Factor 2' (refers to attributes of svi geodataframe)
        -Currently, Factor 2 is inversed, because it is the economic factor (inversely related to vulnerability)
        -TODO: should be a seperate function method to allow input on whether variable is inversed or not

    fpath, extension, and times:
        - Want to be able to plot time series data of the aggregates, therefore need to read in multiple residential parcel shapefiles with flow decomposition information. 'fpath' is the folder location of the flow decomposition shapefiles. Times, in the format yyyymmddhh, are the times that are read in. The files that are read in take the format {fpath}/{time in times}{extension}

    outlier masking:
        - Due to the occurance of severe outliers occuring in travel times, (e.g., a residential parcel located at the end of a long roadway that has an extremely high travel time), the cost to access resources are "projected" onto a feasible set. This is done by calculating the interquartile range (IQR). A value is considered an outlier if it is equal to or greater than 3 times the IQR. These values are set to be equal to 3*IQR.


    :param axes: axes of figure being manipulated
    :type axes: a matplotlib.axes element
    :param cost_atrs: the attributes within the res_parcels geodataframe that need to be considered
    :type cost_atrs: list of strings
    :param svi: the social vulnerabilty dataset to aggregate results to
    :type svi: geodataframe 
    :param fpath: the folder location where the flow decompostiion residential parcels are saved
    :type fpath: string
    :param extension: the file ending for each low decompostiion residential parcel that will be loaded
    :type extension: string
    :param times: the times to be loaded, in yyyymmddhh format
    :type times: list of stings
    :param columns: the number of columns in the time series plots
    :type columns: int
    :param methods: the x-axis variable
    :type methods: list of strings
    :param weight: weight to apply to parcels that cannot access the resource. The value will be set to the maximum travel time for that resource times weight
    :type weight: int

    """
    column = 0
    row = 0   

    # for every time to analyze,
    for time in times:
        # create variable to hold information for percentile table
        table_vals = [[0, 0, 0, 0] for _ in range(len(methods))]

        # create a variable to hold information regarding gini coefficeints
        ginis = [0 for _ in range(len(methods))]

        # load aggregate household reliablity data
        res_parcels = gpd.read_file(f'{fpath}/{prefix}_{time}.shp')

        # For each method, 
        for b, method in enumerate(methods):

            # weight by PPUNIT
            res_parcels[f'aggregate_weight_{method}'] = res_parcels['aggregate'] * res_parcels['PPUNIT']

            # Sort by appropriate column 
            if method == 'Factor 2':
                res_parcels.sort_values(method, axis=0, inplace=True, ascending=False)
            else:
                res_parcels.sort_values(method, axis=0, inplace=True)

            # convert to array
            array = res_parcels[f'aggregate_weight_{method}'].values

            # Apply SVI weight
            if svi_weight is None:
                pass
            elif svi_weight == 0:
                pass
            else:
                # create rank array
                ranks = np.arange(0,len(array),1)
                # determine min and max
                metric_min = np.min(ranks)
                metric_max = np.max(ranks)
                # perform calcuation 
                array = array*((1-svi_weight)+((ranks-metric_min)
                               * (2*svi_weight))/(metric_max-metric_min))

            # calculate lorenz values
            lorenz = array.cumsum() / array.sum()
            lorenz = np.insert(lorenz, 0, 0)
            val = np.arange(lorenz.size)/(lorenz.size-1)

            if method == 'SVI_scaled':
                c_line = 'r'
            elif method == 'Factor 1':
                c_line = 'b'
            elif method == 'Factor 2':
                c_line = 'g'

            # lorenze curve
            lorenz_plot, = axes[row][column].plot(val, lorenz, label='tsamp', color=c_line, linewidth=2)

            perc_25 = np.percentile(lorenz, 25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz, 75)

            table_vals[b][0] = f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3] = f'{100-round(perc_75*100)}%'

            # calc gini and save
            ginis[b] = gini(array)

        col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
        # determine appropriate row labels
        row_labels = []
        for method in methods:
            if method == 'SVI_scaled':
                row_labels.append('SVI')
            elif method == 'Factor 1':
                row_labels.append('SS')
            elif method == 'Factor 2':
                row_labels.append('ES')

        # determine appropriate row colors
        row_colors = []
        for method in methods:
            if method == 'SVI_scaled':
                row_colors.append('r')
            elif method == 'Factor 1':
                row_colors.append('b')
            elif method == 'Factor 2':
                row_colors.append('g')

        # plot ginis
        textstr = []
        for i, row_label in enumerate(row_labels):
            textstr.append(f'{row_label} G: {round(ginis[i],2)}')
        textstr = '\n'.join(textstr)

        # Gini patch properties
        props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
        # text box in upper left in axes coords
        axes[row][column].text(0.05, 0.95, textstr, transform=axes[row][column].transAxes, fontsize=8,
                                verticalalignment='top', bbox=props)

        # Plot quartile information
        the_table = axes[row][column].table(cellText=table_vals,
                                            # colWidths=[0.10]*4,
                                            rowColours=row_colors,
                                            colColours=['lightgrey', 'lightgrey',
                                                        'lightgrey', 'lightgrey'],
                                            rowLabels=row_labels,
                                            colLabels=col_labels,
                                            loc='bottom')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)
        the_table.scale(1, 0.8)
        for i, val in enumerate(methods):
            the_table[i+1, -1].get_text().set_color('white')

        # Equality line
        eq_line, = axes[row][column].plot(
            [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

        # set x and y limits
        axes[row][column].set_xlim(0, 1)
        axes[row][column].set_ylim(0, 1)

        # Set x and y tick markers
        axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[' ', ' ', ' ', ' ', ' '])
                                    #labels=[' ', '25%', '50%', '75%', '100%'])
        axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0],
                                    labels=['0%', '25%', '50%', '75%', '100%'])

        # adjust ticks
        axes[row][column].tick_params(axis="x", direction="in")

        # set axes title
        axes[row][column].set_title(
            f'{time[-2:]}:00 {time[6:8]}/{time[4:6]}/{time[:4]}', style='italic')

        # update row and column names
        column += 1
        if column == columns:
            column = 0
            row += 1

    return axes

# Network redundancy
def network_redundancy(axes,
                       fpath,
                       svi,
                       method='SVI_scaled',
                       variable='net_redun',
                       times=[],
                       annotations=True,
                       columns=3,
                       ftypes=['_']):
    """
    Function to plot lorenz curve of network redundancy, showing impacts on compound, fluvial, and/or pluvial flooding.

    This should automate the creation of network redundancy lorenz curves, with a highly customizable plot tool. Some important usage notes:

    fpath and prefix:
        - used to find the network redundancy results that need to be analyzed
        - redundancy files should have the following naming convention: {fpath}/network_redundancy{ftype}{time} matching the 'times' and 'ftypes' inputs. NOTE: Be careful of placement of '_' in names (i.e., in cases where compound inundation results are saved without an ftype parameter, and ftype should only be '_')
    
    method:
        -the x-axis variable the Lorenz curve is plotted against
        -acceptable values: 'SVI_scaled', 'Factor 1', or 'Factor 2' (refers to attributes of svi geodataframe)
        -Currently, Factor 2 is inversed, because it is the economic factor (inversely related to vulnerability)
        -TODO: should be a seperate function method to allow input on whether variable is inversed or not

    variable:
        - The y-axis varaible being plotted. The default is 'net_redun' and is a column of the network results shapefile being analyzed. 'net_redun' is aggregate network redundancy, and individual resource network results can be plotted by using the appropriate column name from the network redundancy result being plotted
        - regardless of the variable, it is also weighted by PPUNIT (people per unit) Census data, to account for population density differences

    times:
        - should be in the format YYYYMMDDHH
        - used to plot subplot titles as well as point to which network needs to be loaded

    columns:
        - the number of columns in the figure, used to iterate through plotting

    ftypes:
        - Which floods should be considered, list of strings
        - Acceptable values are '_' for compound, '_fluvial_' for fluvial, and '_pluvial_' for pluvial


    :param axes: axes of figure being manipulated. If only plotting a single instance, put axes in [[]]
    :type axes: array of matplotlib.axes
    :param fpath: the folder path where the graphs are stored
    :type fpath: string
    :param svi: social vulnerability estimate shapefile
    :type svi: geodataframe
    :param method: the x-axis variable
    :type method: string
    :param variable: y-axis variable
    :type variable: string
    :param times: the times associated with each network redundancy shapefile to be loaded and analyzed 
    :type times: list of strings
    :param annotations: If True, plots annotation box containing quartile burden information and Gini coefficients
    :type annotations: bool
    :param columns: number of columns in figure being manipulated
    :type columns: int
    :param ftypes: which flood types should be plotted 
    :type ftypes: list of strings
    
    """
    column = 0
    row = 0

    # read in appropriate network redundancy shapefile
    for time in times:
        # create variable to hold information for percentile table
        table_vals = [[0, 0, 0, 0] for _ in range(len(ftypes))]
        # create a variable to hold information regarding gini coefficeints
        ginis = [0 for _ in range(len(ftypes))]

        # for each flood type of interest,
        for b, ftype in enumerate(ftypes):
            # read in network redundancy shapefile
            network_redundancy_parcels = gpd.read_file(f'{fpath}/network_redundancy{ftype}{time}.shp')
            
            # Spatial join
            sjoined_data = gpd.sjoin(left_df=svi,
                                    right_df=network_redundancy_parcels,
                                        how='left')

            # # drop all parcels with 0 redundancy 
            # sjoined_data = sjoined_data[sjoined_data[variable]!=0]

            # count the number of residential parcels within each block group THAT HAVE ZERO REDUNDANCY and relate back to geodataframe
            count_dict = sjoined_data[sjoined_data[variable]==0]['GEOID'].value_counts().to_dict()
            sjoined_data["count"] = sjoined_data["GEOID"].apply(lambda x: count_dict.get(x))

            # Aggregate results
            summed_gdf = sjoined_data.dissolve(by='GEOID', aggfunc={variable: 'sum',
                                                                    'SVI_scaled': 'first',
                                                                    'Factor 1': 'first',
                                                                    'Factor 2': 'first',
                                                                    'PPUNIT': 'first',
                                                                    'count': 'first'})

            # TODO: create ppunit and zero redundancy weighted variable
            summed_gdf[variable] = np.where(
                summed_gdf['count'] > 0, summed_gdf[variable]/summed_gdf['count'], summed_gdf[variable])
            summed_gdf[variable] = summed_gdf[variable]*summed_gdf["PPUNIT"]
        
            # SORT BY appropriate column
            # reliability: factor 2 ascending is False
            # redundancy: factor 2 reliability is True
            if method == 'Factor 2':
                summed_gdf.sort_values(method, axis=0, inplace=True, ascending=True)
            else:
                summed_gdf.sort_values(method, axis=0, inplace=True, ascending=False)

            # TODO: HAVE TO MULTIPLY ARRAY by -1 IF FACTOR 2 DUE TO INVERSE WEALTH
            # TODO: I don't actually think this does anything, what actually matters is just the sort order....
            # single resource specific flood conditions
            # if method == 'Factor 2':
            #     array = summed_gdf[variable].values*-1
            # else:
            #     array = summed_gdf[variable].values
            array = summed_gdf[variable].values

            # calculate lorenz values
            lorenz = array.cumsum() / array.sum()
            lorenz = np.insert(lorenz, 0, 0)
            val = np.arange(lorenz.size)/(lorenz.size-1)

            if ftype == '_':
                c_line = '#cc57a4'
            elif ftype == '_pluvial_':
                c_line = '#ff772e'
            elif ftype == '_fluvial_':
                c_line = '#216aa6'

            # Plot lorenz curve
            lorenz_plot, = axes[row][column].plot(val, lorenz, color=c_line, linewidth=2)

            perc_25 = np.percentile(lorenz, 25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz, 75)

            table_vals[b][0] = f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3] = f'{100-round(perc_75*100)}%'

            # calc gini and save
            ginis[b] = gini(array)*1

        col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
        # determine appropriate row labels
        row_labels = []
        for ftype in ftypes:
            if ftype == '_':
                row_labels.append('Cmpd.')
            elif ftype == '_fluvial_':
                row_labels.append('Fluvial')
            elif ftype == '_pluvial_':
                row_labels.append('Pluvial')

        # determine appropriate row colors
        row_colors=[]
        for ftype in ftypes:
            if ftype == '_':
                row_colors.append('#cc57a4')
            elif ftype == '_fluvial_':
                row_colors.append('#216aa6')
            elif ftype == '_pluvial_':
                row_colors.append('#ff772e')

        # the rectangle is where I want to place the table
        if annotations is True:

            # plot ginis
            textstr = []
            for i, row_label in enumerate(row_labels):
                textstr.append(f'{row_label} G: {round(ginis[i],2)}')
            textstr = '\n'.join(textstr)

            # Gini patch properties
            props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
            # text box in upper left in axes coords
            axes[row][column].text(0.05, 0.95, textstr, transform=axes[row][column].transAxes, fontsize=8,
                        verticalalignment='top', bbox=props)

            # Plot quartile information IN GRAPH
            # the_table = axes[row][column].table(cellText=table_vals,
            #                         colWidths=[0.11]*4,
            #                         rowColours=row_colors,
            #                         colColours=['lightgrey', 'lightgrey',
            #                                                 'lightgrey', 'lightgrey'],
            #                         rowLabels=row_labels,
            #                         colLabels=col_labels,
            #                         loc='lower right')
            # Plot quartile information UNDER GRAPH
            the_table = axes[row][column].table(cellText=table_vals,
                                colWidths=[0.25]*4,
                                rowColours=row_colors,
                                colColours=['lightgrey', 'lightgrey',
                                                        'lightgrey', 'lightgrey'],
                                rowLabels=row_labels,
                                colLabels=col_labels,
                                loc='bottom')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(8)
            the_table.scale(1, 1)
            # for i, val in enumerate(methods):
            #     the_table[i+1, -1].get_text().set_color('white')

        # Equality line
        eq_line, = axes[row][column].plot(
            [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

        # Set x and y tick markers 
        # axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
        #     '0%', '25%', '50%', '75%', '100%'])

        axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                                     ' ', ' ', ' ', ' ', ' '])
        axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
            '0%', '25%', '50%', '75%', '100%'])

        # set x and y limits
        axes[row][column].set_xlim(0, 1)
        axes[row][column].set_ylim(0, 1)

        # set axes title
        axes[row][column].set_title(
            f'{time[-2:]}:00 {time[6:8]}/{time[4:6]}/{time[:4]}', style='italic')

        # update row and column names
        column += 1
        if column == columns:
            column = 0
            row += 1



    pass

# Individidaul redundancy
def individual_redundancy(axes,
                          svi,
                          fpath,
                          times,
                          extension,
                          columns,
                          methods=['svi_scaled'],
                          atr='redundancy',
                          annotations=True,
                          svi_weight=None):
    """
    Function to plot Lorenz curves of individual redundancy for resource accessibility, analyzing the cost to access resources. Some important usage notes:

    atr:
        - this is the string associated with the column within the redundancy_parcels shapefile to calculate the lorenz curve. the default, 'redundancy', is the aggregate individual redundancy score
        - individual resources could also be plotted if the appropriate column name is used
        - regardless, the variable weighted by the PPUNIT Census variable to account for population density differences

    methods:
        -the x-axis variable the Lorenz curve is plotted against
        -acceptable values: 'SVI_scaled', 'Factor 1', and/or 'Factor 2' (refers to attributes of svi geodataframe)
        -Currently, Factor 2 is inversed, because it is the economic factor (inversely related to vulnerability)
        -TODO: should be a seperate function method to allow input on whether variable is inversed or not

    :param axes: axes of figure being manipulated
    :type axes: matplotlib.axes
    :param svi: the social vulnerabilty dataset to aggregate results to
    :type svi: geodataframe
    : param fpath: the filepath to the folder where the individiaul redundancy metric parcel shapefiles are saved
    :type fpath: string
    :param times: the times to be loaded and analyzed in yyyymmddhh format
    :type times: list of strings
    :param extension: the file extension used in loading the individual redundancy metric parcel shapefiles 
    :type extension: string
    :param columns: number of columns in the multiplot
    :type columns: int
    :param methods: the x-axis variable, acceptable values include 'SVI_scaled', 'Factor 1', and/or 'Factor 2'
    :type methods: list of strings
    :param atr: name of the column in 'redundancy_parcels' to be used as y-axis variable of lorenz curve.
    :type atr: string
    :param annotations: If True, plots annotation box containing quartile burden information and Gini coefficients
    :type annotations: bool

    """

    column = 0
    row = 0 

    # read in appropriat individual redundancy shapefile and convert to geopandas dataframe
    for time in times:
        redundancy_parcels = gpd.read_file(f'{fpath}/{extension}{time}.shp')

        # create variable to hold information for percentile table
        table_vals = [[0, 0, 0, 0] for _ in range(len(methods))]
        # create a variable to hold information regarding gini coefficeints
        ginis = [0 for _ in range(len(methods))]

        # for each method,
        for b, method in enumerate(methods):

            # Spatial join
            sjoined_data = gpd.sjoin(left_df=svi,
                                    right_df=redundancy_parcels,
                                    how='left')

            # # drop all parcels with 0
            # sjoined_data = sjoined_data[sjoined_data[atr] != 0]

            # count the number of residential parcels within each block group WITH 0 REDUNDANCY and relate back to geodataframe
            count_dict = sjoined_data[sjoined_data[atr] == 0]['GEOID'].value_counts().to_dict()
            sjoined_data["count"] = sjoined_data["GEOID"].apply(lambda x: count_dict.get(x))

            # Aggregate results
            summed_gdf = sjoined_data.dissolve(by='GEOID', aggfunc={atr: 'sum',
                                                                    'SVI_scaled': 'first',
                                                                    'Factor 1': 'first',
                                                                    'Factor 2': 'first',
                                                                    'PPUNIT': 'first',
                                                                    'count': 'first'})

            # TODO: create ppunit and zero redundancy weighted variable
            summed_gdf[atr] = np.where(
                summed_gdf['count'] > 0, summed_gdf[atr]/summed_gdf['count'], summed_gdf[atr])
            summed_gdf[atr] = summed_gdf[atr]*summed_gdf["PPUNIT"]

            # SORT BY appropriate column
            # in reliability, ascending is False for economic factor
            # in redundancy, ascending is True for economic factor
            if method == 'Factor 2':
                summed_gdf.sort_values(method, axis=0, inplace=True, ascending=True)
            else:
                summed_gdf.sort_values(method, axis=0, inplace=True, ascending=False)

            # TODO: HAVE TO MULTIPLY ARRAY by -1 IF FACTOR 2 DUE TO INVERSE WEALTH
            # THIS DOESN'T ACUTALLY DO ANYTHING... what matters is the order above
            # single resource specific flood conditions
            if method == 'Factor 2':
                array = summed_gdf[atr].values*-1
            else:
                array = summed_gdf[atr].values

            # Apply SVI weight
            if svi_weight is None:
                pass
            elif svi_weight == 0:
                pass
            else:
                # create rank array
                ranks = np.arange(0,len(array),1)
                # determine min and max
                metric_min = np.min(ranks)
                metric_max = np.max(ranks)
                # perform calcuation 
                array = array*((1-svi_weight)+((ranks-metric_min)
                               * (2*svi_weight))/(metric_max-metric_min))

            # calculate lorenz values
            lorenz = array.cumsum() / array.sum()
            lorenz = np.insert(lorenz, 0, 0)
            val = np.arange(lorenz.size)/(lorenz.size-1)

            if method == 'SVI_scaled':
                c_line = 'r'
            elif method == 'Factor 1':
                c_line = 'b'
            elif method == 'Factor 2':
                c_line = 'g'

            # lorenze curve
            lorenz_plot, = axes[row][column].plot(val, lorenz, color=c_line, linewidth=2)

            perc_25 = np.percentile(lorenz, 25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz, 75)

            table_vals[b][0] = f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3] = f'{100-round(perc_75*100)}%'

            # calc gini and save
            ginis[b] = gini(array)*1

        col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
        # determine appropriate row labels
        row_labels = []
        for method in methods:
            if method == 'SVI_scaled':
                row_labels.append('SVI')
            elif method == 'Factor 1':
                row_labels.append('SS')
            elif method == 'Factor 2':
                row_labels.append('ES')

        # determine appropriate row colors
        row_colors = []
        for method in methods:
            if method == 'SVI_scaled':
                row_colors.append('r')
            elif method == 'Factor 1':
                row_colors.append('b')
            elif method == 'Factor 2':
                row_colors.append('g')

        # the rectangle is where I want to place the table
        if annotations is True:

            # plot ginis
            textstr = []
            for i, row_label in enumerate(row_labels):
                textstr.append(f'{row_label} G: {round(ginis[i],2)}')
            textstr = '\n'.join(textstr)

            # Gini patch properties
            props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
            # text box in upper left in axes coords
            axes[row][column].text(0.05, 0.95, textstr, transform=axes[row][column].transAxes, fontsize=8,
                    verticalalignment='top', bbox=props)

            # Plot quartile information
            the_table = axes[row][column].table(cellText=table_vals,
                                # colWidths=[0.10]*4,
                                rowColours=row_colors,
                                colColours=['lightgrey', 'lightgrey',
                                            'lightgrey', 'lightgrey'],
                                rowLabels=row_labels,
                                colLabels=col_labels,
                                loc='bottom')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(8)
            the_table.scale(1, 0.8)
            for i, val in enumerate(methods):
                the_table[i+1, -1].get_text().set_color('white')

        # Equality line
        eq_line, = axes[row][column].plot(
            [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[' ', ' ', ' ', ' ', ' '])
                                    #labels=[' ', '25%', '50%', '75%', '100%'])
        axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0],
                                      labels=['0%', '25%', '50%', '75%', '100%'])

        # adjust ticks
        axes[row][column].tick_params(axis="x", direction="in")

        # set x and y limits
        axes[row][column].set_xlim(0, 1)
        axes[row][column].set_ylim(0, 1)

        # set axes title
        axes[row][column].set_title(
            f'{time[-2:]}:00 {time[6:8]}/{time[4:6]}/{time[:4]}', style='italic')

        # update row and column names
        column += 1
        if column == columns:
            column = 0
            row += 1

# Individual recoverability
def individual_recoverability(axes,
                                 svi,
                                 fpath,
                                 init_time,
                                 times,
                                 prefix,
                                 perc=0.5,
                                 columns=1,
                                 methods=['SVI_scaled'],
                                 svi_weight=None,
                                 annotations='in'):
    """
    calculates the amount of time it takes to get within X percentage of initial conditions post peak and creates associated lorenz curve    
    """

    row=0
    column=0

    # read in shapefile of base conditions 
    res_parcels = gpd.read_file(f'{fpath}/{prefix}_{init_time}.shp')

    # CREATE EMPTY ARRAYS
    # create an empty array to hold base data (the initial reliability conditions)
    base = [np.zeros_like(svi.shape[0], dtype='float') for _ in range(len(methods))]
    # create an empty array to hold final time data
    time_series = [np.full(svi.shape[0], fill_value=len(times), dtype='float') for _ in range(len(methods))]
    # create an empty array to keep track of the ppunit weights to apply to the final values
    ppunit_weights = [np.zeros_like(svi.shape[0], dtype='float') for _ in range(len(methods))]
    
    # CALCULATE BASE VARIABLES
    # for each method,
    for b, method in enumerate(methods):

        # sort by appropriate column appropriate column
        if method == 'Factor 2':
            res_parcels.sort_values(method, axis=0, inplace=True, ascending=False)
        else:
            res_parcels.sort_values(method, axis=0, inplace=True)

        # TODO: HAVE TO MULTIPLY ARRAY by -1 IF FACTOR 2 DUE TO INVERSE WEALTH
        # Not sure if this actually does anything???
        # single resource specific flood conditions
        if method == 'Factor 2':
            array = res_parcels['aggregate'].values*-1
            
        else:
            array = res_parcels['aggregate'].values

        # fill array data into base array to hold for later comparison
        base[b] = array

        # fill ppunit_weights data
        ppunit_weights[b] = res_parcels['PPUNIT'].values


    # PERFORM METRIC CALCULATION
    for time_idx, time in enumerate(times):
        # load aggregate household reliability data
        res_parcels = gpd.read_file(f'{fpath}/{prefix}_{time}.shp')

        # For each method, 
        for b, method in enumerate(methods):

            # Sort by appropriate column 
            if method == 'Factor 2':
                res_parcels.sort_values(method, axis=0, inplace=True, ascending=False)
            else:
                res_parcels.sort_values(method, axis=0, inplace=True)

            # TODO: HAVE TO MULTIPLY ARRAY by -1 IF FACTOR 2 DUE TO INVERSE WEALTH
            # I don't think this actually does anything
            # single resource specific flood conditions
            if method == 'Factor 2':
                array = res_parcels['aggregate'].values*-1
            else:
                array = res_parcels['aggregate'].values

            # compare to stored values and replace time_series data when applicable
            time_series[b] = np.where(base[b]/array >= perc, np.minimum(time_series[b], time_idx), time_series[b])
            print(sum(1 for i in time_series[b] if i < time_idx))

    # weight the time_series data by the appropriate PPUNIT values
    for b, method in enumerate(methods):
        time_series[b] *= ppunit_weights[b]

        # Apply SVI weight
        if svi_weight is None:
            pass
        elif svi_weight == 0:
            pass
        else:
            # create rank array
            ranks = np.arange(0, len(time_series[b]), 1)
            # determine min and max
            metric_min = np.min(ranks)
            metric_max = np.max(ranks)
            # perform calcuation 
            time_series[b] = time_series[b]*((1-svi_weight)+((ranks-metric_min) * (2*svi_weight))/(metric_max-metric_min))

    # create variable to hold information for percentile table
    table_vals = [[0, 0, 0, 0] for _ in range(len(methods))]
    # create a variable to hold information regarding gini coefficeints
    ginis = [0 for _ in range(len(methods))]
    for b, method in enumerate(methods):

        # calculate lorenz values
        lorenz = time_series[b].cumsum() / time_series[b].sum()
        lorenz = np.insert(lorenz, 0, 0)
        val = np.arange(lorenz.size)/(lorenz.size-1)

        if method == 'SVI_scaled':
            c_line = 'r'
        elif method == 'Factor 1':
            c_line = 'b'
        elif method == 'Factor 2':
            c_line = 'g'

        # lorenze curve
        lorenz_plot, = axes[row][column].plot(
            val, lorenz, label='tsamp', color=c_line, linewidth=2)

        perc_25 = np.percentile(lorenz, 25)
        perc_50 = np.percentile(lorenz, 50)
        perc_75 = np.percentile(lorenz, 75)

        table_vals[b][0] = f'{round(perc_25*100)}%'
        table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
        table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
        table_vals[b][3] = f'{100-round(perc_75*100)}%'

        # calc gini and save
        ginis[b] = gini(time_series[b])

    col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
    # determine appropriate row labels
    row_labels = []
    for method in methods:
        if method == 'SVI_scaled':
            row_labels.append('SVI')
        elif method == 'Factor 1':
            row_labels.append('SS')
        elif method == 'Factor 2':
            row_labels.append('ES')

    # determine appropriate row colors
    row_colors = []
    for method in methods:
        if method == 'SVI_scaled':
            row_colors.append('r')
        elif method == 'Factor 1':
            row_colors.append('b')
        elif method == 'Factor 2':
            row_colors.append('g')

    # plot ginis
    textstr = []
    for i, row_label in enumerate(row_labels):
        textstr.append(f'{row_label} G: {round(ginis[i],2)}')
    textstr = '\n'.join(textstr)

    # Gini patch properties
    props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
    # text box in upper left in axes coords
    axes[row][column].text(0.05, 0.95, textstr, transform=axes[row][column].transAxes, fontsize=8,
                            verticalalignment='top', bbox=props)


    # for single plots, put annotations in the lower right
    if annotations == 'in':  
        # Plot quartile information
        the_table = axes[row][column].table(cellText=table_vals,
                                            colWidths=[0.11]*4,
                                            rowColours=row_colors,
                                            colColours=['lightgrey', 'lightgrey',
                                                        'lightgrey', 'lightgrey'],
                                            rowLabels=row_labels,
                                            colLabels=col_labels,
                                            loc='lower right')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)
        the_table.scale(1, 1)
        for i, val in enumerate(methods):
            the_table[i+1, -1].get_text().set_color('white')

        # Set x and y tick markers
        axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
            '0%', '25%', '50%', '75%', '100%'])
        axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
            '0%', '25%', '50%', '75%', '100%'])

    # if in a multiplot (e.g., equity plot) put annotations belwo the graph
    elif annotations == 'below':
        # Plot quartile information
        the_table = axes[row][column].table(cellText=table_vals,
                            # colWidths=[0.10]*4,
                            rowColours=row_colors,
                            colColours=['lightgrey', 'lightgrey',
                                        'lightgrey', 'lightgrey'],
                            rowLabels=row_labels,
                            colLabels=col_labels,
                            loc='bottom')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)
        the_table.scale(1, 0.8)
        for i, val in enumerate(methods):
            the_table[i+1, -1].get_text().set_color('white')

        # Equality line
        eq_line, = axes[row][column].plot(
            [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[' ', ' ', ' ', ' ', ' '])
        axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0],
                                      labels=['0%', '25%', '50%', '75%', '100%'])

    # Equality line
    eq_line, = axes[row][column].plot(
        [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

    # set x and y limits
    axes[row][column].set_xlim(0, 1)
    axes[row][column].set_ylim(0, 1)

# individual recoverability cumulative time series
def individual_recoverability_cumulative_ts(svi,
                                            fpath,
                                            init_time,
                                            times,
                                            prefix,
                                            percs=[0.4,0.5,0.6,0.7,0.8,0.9,1]):
    """
    Based on different percent recoverabilities, plot ts of cumulative # of BGs that have reached recovered state
    """
    # read in shapefile of base conditions
    res_parcels = gpd.read_file(f'{fpath}/{prefix}_{init_time}.shp')

    # CREATE EMPTY ARRAYS
    # create an empty array to hold base data (the initial reliability conditions)
    base = [np.zeros_like(svi.shape[0], dtype='float')]
    # create an empty array to hold final time data
    time_series = [np.full(svi.shape[0], fill_value=len(times), dtype='float') for _ in range(len(percs))]

    # create output array
    output = [np.full(len(times), fill_value=0, dtype='float') for _ in range(len(percs)) ]

    # CALCULATE BASE VARIABLES
    base = res_parcels['aggregate'].values


    # PERFORM METRIC CALCULATION
    for time_idx, time in enumerate(times):
        # load aggregate household reliability data
        res_parcels = gpd.read_file(f'{fpath}/{prefix}_{time}.shp')

        # For each method,
        for p, perc in enumerate(percs):

            array = res_parcels['aggregate'].values
            time_series[p] = np.where(base/array >= perc, np.minimum(time_series[p], time_idx), time_series[p])
            output[p][time_idx] = sum(1 for i in time_series[p] if i < time_idx)

    # plot the cumulative results
    # create the plot
    fig, ax = plt.subplots(1, 1, figsize=(180*mm, 90*mm))
    fig.set_facecolor('none')
    # Graph attributes when multiple plots
    # fig.suptitle('Time Series of Cumulative\nHousehold Recovery by Block Group',
    #              fontsize=18, weight='bold', style='italic')
    fig.supxlabel('Time',
                  weight='bold',
                  x=0.5,
                  fontsize=16)
    fig.supylabel(' BGs returned to\nX% functionality',
                  weight='bold',
                  fontsize=16,)

    # plot the appropriate lines
    for idx, data in enumerate(output):
        ax.plot(data,
                linewidth=2,
                markersize=6,
                marker='o',
                label=percs[idx])

    # plot horizontal line for maximum number of block groups
    # plt.axhline(y=len(svi), color='black',linestyle='dashed')

    label = ax.text(2,len(svi)-10, ' Max Number of \n Block Groups ', va='center', ha='center',fontsize=12)
    ax.get_figure().canvas.draw()
    bbox = label.get_window_extent().transformed(ax.transData.inverted())
    # add hlines next to bounding box
    left, right = ax.get_xlim()
    ax.hlines([len(svi)]*2, [left, bbox.x1], [bbox.x0, right],
              color='black', linestyle='dashed')


    # # set xticks
    xticks = []
    for time in times:
        xticks.append(f'{time[-2:]}:00')
    ax.set_xticks(np.arange(0, len(output[0]), 1))
    ax.set_xticklabels(xticks)
    plt.xticks(rotation=45)

    # set legend
    ax.legend(loc='lower right')
    # change legend to multiple columns
    ax.legend(ncol=4)
    plt.grid()
    plt.subplots_adjust(top=0.843,
                        bottom=0.22,
                        left=0.135,
                        right=0.979,
                        hspace=0.2,
                        wspace=0.2)
    ax.set_xlim((-0.5,20))

# Network recoverability
def network_recoverability(axes,
                           fpath,
                           prefix,
                           svi,
                           init_time,
                           method='SVI_scaled',
                           variable='agr_no_cap',
                           times=[],
                           annotations=True,
                           ftypes=['inundation']):
    """
    calculates the amount of time it takes to get within X percentage of initial conditions post peak flood and creates associated lorenz curve
    """
    perc=0.75
    row=0
    column=0
    columns=0

    # CREATE EMPTY ARRAYS TO STORE DATA
    # create an empty array to hold base data (the initial reliability conditions)
    base = [np.zeros_like(svi.shape[0], dtype='float') for _ in range(len(ftypes))]
    # create an empty array to hold final time data
    time_series = [np.full(svi.shape[0], fill_value=len(times), dtype='float') for _ in range(len(ftypes))]
    # create an empty array to keep track of the ppunit weights to apply to the final values
    ppunit_weights = [np.zeros_like(svi.shape[0], dtype='float') for _ in range(len(ftypes))]

    # DETERMINE BASE CONDITIONS
    for idx, ftype in enumerate(ftypes):
        network = mynet.read_graph_from_disk(path=fpath,
                                             name=f'{prefix}_{init_time}_{ftype}')
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
        svi[variable] = svi["GEOID"].apply(lambda x: count_dict.get(x))
        svi[variable] = svi[variable].fillna(0)

        # sort svi dataframe by appropriate column
        if method == 'Factor 2':
            svi.sort_values('Factor 2', axis=0, inplace=True, ascending=False)
        else:
            svi.sort_values(method, axis=0, inplace=True, ascending=True)

        # LORENZ CURVE
        # extract values to use in Lorenz curve
        # DOES THIS ACTUALLY DO ANYTHING?
        if method == 'Factor 2':
            array = svi[variable].values*-1
        else:
            array = svi[variable].values

        # Save the values in the base array
        base[idx] = array

        # update appropriate ppunit weights
        ppunit_weights[idx] = svi['PPUNIT'].array

    # for each time to analyze,
    for time_idx, time in enumerate(times):
    # for each type of flooding,
        for idx, ftype in enumerate(ftypes):
        
            network = mynet.read_graph_from_disk(path=fpath,
                                                    name=f'{prefix}_{time}_{ftype}')
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

            # sort svi dataframe by appropriate column
            if method == 'Factor 2':
                svi.sort_values('Factor 2', axis=0, inplace=True, ascending=False)
            else:
                svi.sort_values(method, axis=0, inplace=True, ascending=True)

            # DOES THIS ACTUALLY DO ANYTHING?
            if method == 'Factor 2':
                array = svi['agr_no_cap'].values*-1
            else:
                array = svi['agr_no_cap'].values

            # compare to stored values and replace time_series data when applicable
            time_series[idx] = np.where(np.divide(base[idx], array, out=np.zeros_like(base[idx]), where=array!=0) >= perc, 
                                        np.minimum(time_series[idx], time_idx),
                                        time_series[idx])
            time_series[idx] = np.where(array == 0, np.minimum(time_series[idx], time_idx), time_series[idx])                      
            
            # time_series[idx] = np.where(base[idx]/array >= perc, np.minimum(time_series[idx], time_idx),time_series[idx])
            
    # weight the time_series data by the appropriate PPUNIT values
    for result in time_series:
        print(np.count_nonzero(result == 19))
    for b, method in enumerate(ftypes):
        time_series[b] *= ppunit_weights[b]

    # create variable to hold information for percentile table
    table_vals = [[0, 0, 0, 0] for _ in range(len(ftypes))]
    # create a variable to hold information regarding gini coefficeints
    ginis = [0 for _ in range(len(ftypes))]
    
    for b, ftype in enumerate(ftypes):

        # calculate lorenz values
        lorenz = time_series[b].cumsum() / time_series[b].sum()
        lorenz = np.insert(lorenz, 0, 0)
        val = np.arange(lorenz.size)/(lorenz.size-1)

        if ftype == 'inundation':
            c_line = '#cc57a4'
        elif ftype == 'inundation_pluvial':
            c_line = '#ff772e'
        elif ftype == 'inundation_fluvial':
            c_line = '#216aa6'

        # lorenze curve
        lorenz_plot, = axes[row][column].plot(val, lorenz, label='tsamp', color=c_line, linewidth=2)

        perc_25 = np.percentile(lorenz, 25)
        perc_50 = np.percentile(lorenz, 50)
        perc_75 = np.percentile(lorenz, 75)

        table_vals[b][0] = f'{round(perc_25*100)}%'
        table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
        table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
        table_vals[b][3] = f'{100-round(perc_75*100)}%'

        # calc gini and save
        ginis[b] = gini(time_series[b])

    col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
    # determine appropriate row labels
    row_labels = []
    for ftype in ftypes:
        if ftype == 'inundation':
            row_labels.append('Cmpd.')
        elif ftype == 'inundation_fluvial':
            row_labels.append('Fluvial')
        elif ftype == 'inundation_pluvial':
            row_labels.append('Pluvial')

    # determine appropriate row colors
    row_colors=[]
    for ftype in ftypes:
        if ftype == 'inundation':
            row_colors.append('#cc57a4')
        elif ftype == 'inundation_fluvial':
            row_colors.append('#216aa6')
        elif ftype == 'inundation_pluvial':
            row_colors.append('#ff772e')

    # the rectangle is where I want to place the table
    if annotations is True:

        # plot ginis
        textstr = []
        for i, row_label in enumerate(row_labels):
            textstr.append(f'{row_label} G: {round(ginis[i],2)}')
        textstr = '\n'.join(textstr)

        # Gini patch properties
        props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
        # text box in upper left in axes coords
        axes[row][column].text(0.05, 0.95, textstr, transform=axes[row][column].transAxes, fontsize=8,
                                verticalalignment='top', bbox=props)

        # Plot quartile information
        the_table = axes[row][column].table(cellText=table_vals,
                                            colWidths=[0.11]*4,
                                            rowColours=row_colors,
                                            colColours=['lightgrey', 'lightgrey',
                                                        'lightgrey', 'lightgrey'],
                                            rowLabels=row_labels,
                                            colLabels=col_labels,
                                            loc='lower right')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)
        the_table.scale(1, 1)

    # Equality line
    eq_line, = axes[row][column].plot(
        [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

    # Set x and y tick markers
    axes[row][column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], 
                                labels=[' ', '25%', '50%', '75%', '100%'])
    axes[row][column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0],
                                    labels=['0%', '25%', '50%', '75%', '100%'])

    # set x and y limits
    axes[row][column].set_xlim(0, 1)
    axes[row][column].set_ylim(0, 1)

    # update row and column names
    column += 1
    if column == columns:
        column = 0
        row += 1

# network recoverability cumulative time series
def network_recovoverability_cumulative_ts(fpath,
                                           prefix,
                                           svi,
                                           init_time,
                                           variable='agr_no_cap',
                                           times=[],
                                           percs=[0,3,0.4,0.5,0.6,0.7,0.8,0.9,1],
                                           ftypes=['inundation']):

    # CREATE EMPTY ARRAYS TO STORE DATA
    # create an empty array to hold base data (the initial reliability conditions)
    base = [np.zeros_like(svi.shape[0], dtype='float')for _ in range(len(ftypes))]
    # create an empty array to hold final time data
    time_series = [[np.full(svi.shape[0], fill_value=len(times), dtype='float') for _ in  range(len(ftypes))] for _ in range(len(percs))]

    # create output array
    output = [[np.full(len(times), fill_value=0, dtype='float')for _ in range(len(ftypes))] for _ in range(len(percs))]

    # DETERMINE BASE CONDITIONS
    for idx, ftype in enumerate(ftypes):
        network = mynet.read_graph_from_disk(path=fpath,
                                             name=f'{prefix}_{init_time}_{ftype}')
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
        svi[variable] = svi["GEOID"].apply(lambda x: count_dict.get(x))
        svi[variable] = svi[variable].fillna(0)

        # Save the values in the base array
        base[idx] = svi[variable].array

    # for each time to analyze,
    for time_idx, time in enumerate(times):
    # for each type of flooding,
        for idx, ftype in enumerate(ftypes):
        
            network = mynet.read_graph_from_disk(path=fpath,
                                                    name=f'{prefix}_{time}_{ftype}')
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
            svi[variable] = svi["GEOID"].apply(lambda x: count_dict.get(x))
            svi[variable] = svi[variable].fillna(0)
            
            array = svi[variable].values
          
            for p, perc in enumerate(percs):
                # compare to stored values and replace time_series data when applicable
                time_series[p][idx] = np.where(np.divide(base[idx], array, out=np.zeros_like(base[idx]), where=array!=0) >= perc, 
                                            np.minimum(time_series[p][idx], time_idx),
                                            time_series[p][idx])
                time_series[p][idx] = np.where(array == 0, np.minimum(time_series[p][idx], time_idx), time_series[p][idx]) 

                output[p][idx][time_idx] = sum(1 for i in time_series[p][idx] if i < time_idx)

    # plot the results
    fig, axes = plt.subplots(3,1,figsize=(180*mm,220*mm))
    fig.set_facecolor('none')
    # fig.suptitle('Time Series of Cumulative\nNetwork Recovery by Block Group',
    #              fontsize=18, weight='bold', style='italic')
    fig.supxlabel('Time',
                  weight='bold',
                  x=0.5,
                  fontsize=16)
    fig.supylabel(' BGs returned to X% functionality',
                  weight='bold',
                  fontsize=16,)

    # create plotable data
    # for each flood type,
    for idx, ftype in enumerate(ftypes):
        # for each percentile,
        for p, perc in enumerate(percs):
            
            axes[idx].plot(output[p][idx],
                            linewidth=2,
                            markersize=6,
                            marker='o',
                            label=percs[p])

            axes[idx].legend(loc='lower right')
            axes[idx].legend(ncol=4)
            # # set xticks
            xticks = []
            for time in times:
                xticks.append(f'{time[-2:]}:00')
            axes[idx].set_xticks(np.arange(0, len(output[p][idx]), 1))
            axes[idx].set_xticklabels(xticks, rotation=45)

            label = axes[idx].text(2,len(svi)-10, ' Max Number of \n Block Groups ', va='center', ha='center',fontsize=10)
            axes[idx].get_figure().canvas.draw()
            bbox = label.get_window_extent().transformed(axes[idx].transData.inverted())
            # add hlines next to bounding box
            left, right = axes[idx].get_xlim()
            axes[idx].hlines([len(svi)]*2, [left, bbox.x1], [bbox.x0, right],
                    color='black', linestyle='dashed')
            axes[idx].set_xlim((-0.5, 20))

    axes[0].set_ylabel('Compound Inundation', style ='italic',fontsize=12)
    axes[1].set_ylabel('Fluvial Inundation', style ='italic',fontsize=12)
    axes[2].set_ylabel('Pluvial Inundation', style ='italic',fontsize=12)
    axes[0].grid()
    axes[1].grid()
    axes[2].grid()

    plt.tight_layout()


# Individual equity
def individual_equity(reliability_fpath,
                      reliability_prefix,
                      redundancy_fpath,
                      redundancy_extension,
                      recoverability_fpath,
                      recoverability_init_time,
                      recoverability_prefix,
                      recoverability_times,
                      svi,
                      time,
                      methods=['SVI_scaled'],
                    #   methods=['SVI_scaled', 'Factor 1', 'Factor 2'],
                      svi_weights=[0, .5 , 1]):
    """
    Create plot of equity for household metrics.

    """
    fig, axes = plt.subplots(3, 3, figsize=(180*mm, 180*mm))
    fig.set_facecolor('none')
    fig.supxlabel('Indicator Percentile',
                weight='bold',
                x=0.5,
                fontsize=16)
    fig.supylabel('Normalized Cumulative Sum of Household Metric',
                weight='bold',
                fontsize=16)

    for idx, svi_weight in enumerate(svi_weights):
        # # first column: redundancy
        individual_redundancy([[axes[0,idx]]],
                              svi=svi,
                              fpath=redundancy_fpath,
                              times=[time],
                              extension=redundancy_extension,
                              columns=1,
                              methods=methods,
                              atr='redundancy',
                              annotations=True,
                              svi_weight=svi_weight)
  
        # second column: reliability
        individual_reliability_aggregate([[axes[1,idx]]],
                                         prefix=reliability_prefix,
                                         fpath=reliability_fpath,
                                         times=[time],
                                         columns=1,
                                         methods=methods,
                                         svi_weight=svi_weight)

        # third column: recoverability
        individual_recoverability(axes=[[axes[2,idx]]],
                                  svi=svi,
                                  fpath=recoverability_fpath,
                                  init_time=recoverability_init_time,
                                  times=recoverability_times,
                                  prefix=recoverability_prefix,
                                  columns=1,
                                  methods=methods,
                                  svi_weight=svi_weight,
                                  annotations='below')
        
    # edit column and row headings 
    for idx, weight in enumerate(svi_weights):
        axes[0][idx].set_title(f'{np.round(weight*100,1)}% Rank Weight', style='italic', fontsize=12)
        axes[1][idx].set_title(' ')
    
    axes[0][0].set_ylabel('Redundancy', style ='italic',fontsize=12)
    axes[1][0].set_ylabel('Reliability', style ='italic',fontsize=12)
    axes[2][0].set_ylabel('Recoverability', style ='italic',fontsize=12)

    return fig

# Flooded roads
def flooded_roads(times,
                  fpath,
                  prefix,
                  suffixes,
                  closed_road_var):
    """
    Plots a time series figure of the number of pluvially, fluivilly, and compound flooded roads.

    
    The graphs should be saved with the following file convention, for each suffix in suffixes and time in times: f'{fpath}/{prefix}_{time}_{suffix}'


    :param times: list of times to include in the plot in YYYYMMDDHH format
    :type times: list of strings
    :param fpath: path where the graph networks are saved
    :type fpath: string
    :param prefix: prefix of the graphs to be loaded
    :type prefix: string
    :param suffixes: the pluvial, fluvial, and compound inundated graphs (IN THAT ORDER) file sufixes
    :type suffixes: list of strings
    :param closed_road_var: the variable within each graph related to the number of closed roads
    :type closed_road_var: string

    """
    # create empty nested list to keep track of pluvial, fluvial, and compound flooding closed roads
    closed_roads=[[],[],[]]
 
    # for each suffix, or each type of flooding,
    for idx, suffix in enumerate(suffixes):
        # for each time to be included,
        for timex in times:
            # read in graph, determine the number of closed roads, and append the list
            network = mynet.read_graph_from_disk(path=fpath, name=f'{prefix}_{timex}_{suffix}')
            # number of edges that are one way
            num_edges = nx.MultiDiGraph([(u,v,d) for u,v,d in network.edges(data=True) if d[closed_road_var] == 0]).number_of_edges()
            closed_roads[idx].append(num_edges)


    # create the plot
    fig, ax = plt.subplots(1, 1, figsize=(180*mm, 90*mm))
    fig.set_facecolor('none')
    # Graph attributes when multiple plots
    # fig.suptitle('Time Series of Closed Roads',
    #             fontsize=18, weight='bold', style='italic')
    fig.supxlabel('Time',
                weight='bold',
                x=0.5,
                fontsize=16)
    fig.supylabel('Number of roads with 0 capacity',
                weight='bold',
                fontsize=16)

    # plot the appropriate lines
    colors = ['#ff772e', '#216aa6', '#cc57a4']
    labels = ['Pluvially inundated roads', 'Fluvially inundated roads', 'Compound Flooded Roads ']

    for idx, data in enumerate(closed_roads):
        print(data)
        ax.plot(data,
                  color=colors[idx],
                  linewidth=2,
                  markersize=6,
                  marker='o',
                  label=labels[idx])
        
    # # set xticks
    xticks=[]
    for time in times:
        xticks.append(f'{time[-2:]}:00')
    ax.set_xticks(np.arange(0, len(closed_roads[0]),1))
    ax.set_xticklabels(xticks)
    plt.xticks(rotation=45)

    # set legend
    ax.legend(loc='upper right')
    plt.grid()
    plt.tight_layout()

# reliability pre-process
def reliability_preprocess(cost_atrs,
                            load_fpath,
                            save_fpath,
                            names,
                            svi,
                            times,
                            extension,
                            methods=['SVI_scaled'],
                            weight=4):
    """
    Taking in a flow decomposition shapefile, calculate the aggregae reliability shapefile.

    This function saves a file of the aggregated household reliability with GEOID information to easily relate back to SVI data for plotting and analysis purpoess
    """

    for time in times:
        res_parcels = gpd.read_file(f'{load_fpath}/{time}{extension}')

        all_arrays = [np.zeros_like(svi.shape[0], dtype='float') for _ in range(len(methods))]
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

        # Spatially join flow decomposition data to SVI 
        sjoined_data = gpd.sjoin(left_df=svi,
                                    right_df=res_parcels,
                                    how='left')
        
        # fill column NaNs with appropriately weighted values
        max_values = sjoined_data[cost_atrs].max().to_dict()
        for key in max_values:
            max_values[key]*=weight
        sjoined_data.fillna(value=max_values, inplace=True)

        # sum across the columns
        sjoined_data['aggregate'] = sjoined_data[cost_atrs].sum(axis=1)

        # disssolve by GEOID
        aggfunc={'SVI_scaled': 'first',
                    'Factor 1': 'first',
                    'Factor 2': 'first',
                    'PPUNIT': 'first',
                    'aggregate': 'sum'}
        for cost_atr in cost_atrs:
            aggfunc[cost_atr]='sum'
        output = sjoined_data.dissolve(by='GEOID', aggfunc=aggfunc)

        # rename columns
        names_dict={}
        for idx, cost_atr in enumerate(cost_atrs):
            names_dict[cost_atr] = names[idx]
            
        output.rename(names_dict, inplace=True, axis='columns')

        # save the aggregate shapefile
        output.to_file(f'{save_fpath}/agg_individual_reliability_{time}.shp')
    

##############################################
# RELIABILITY PRE-PROCESSING
def reliability_preprocess_func():
    cost_atrs = ['cost_of_fl',
                'cost_of__1',
                'cost_of__2',
                'cost_of__3',
                'cost_of__4',
                'cost_of__5',
                'cost_of__6',
                'cost_of__7']
    names=['grocery','er','pharm','police','conven','fire','ems','fuel']
    times = ['2015052603']
    reliability_preprocess(cost_atrs=cost_atrs,
                        load_fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs/Flow_decomp',
                        save_fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Reliability/individual_reliability',
                        svi=svi,
                        times=times,
                        extension='_inundation_res_parcel_flow_decomp.shp',
                        names=names,
                        weight=4)

# reliability_preprocess_func()

##############################################
# PLOT EXAMPLE OF LORENZ CURVE
# example_lorenz()

##############################################
# NUMBER OF FLOODED ROADS
def num_flooded_roads():
    times = ['2015052518', '2015052519',
            '2015052520', '2015052521',
            '2015052522', '2015052523',
            '2015052600', '2015052601',
            '2015052602', '2015052603',
            '2015052604', '2015052605',
            '2015052606', '2015052607',
            '2015052608', '2015052609',
            '2015052610', '2015052611',
            '2015052612', '2015052613',
            '2015052614', '2015052615',
            '2015052616', '2015052617']

    suffixes=['inundation_pluvial' , 'inundation_fluvial', 'inundation']

    flooded_roads(times=times,
                fpath=graphs_fp,
                prefix='AN_Graph',
                suffixes=suffixes,
                closed_road_var='inundation_capacity_agr')

num_flooded_roads()
plt.show()

##############################################
# INDIVIDUAL EQUITY
def ind_equity_func():
    reliability_fpath = '/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Reliability/individual_reliability'
    reliability_prefix='agg_individual_reliability'
    redundancy_fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Redundancy/individual_redundancy/'
    redundancy_extension='individual_redundancy_'
    recoverability_fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Reliability/individual_reliability'
    recoverability_init_time=2015052518
    recoverability_prefix='agg_individual_reliability'
    recoverability_times = ['2015052522', '2015052523',
                            '2015052600', '2015052601',
                            '2015052602', '2015052603',
                            '2015052604', '2015052605',
                            '2015052606', '2015052607',
                            '2015052608', '2015052609',
                            '2015052610', '2015052611',
                            '2015052612', '2015052613',
                            '2015052614', '2015052615',
                            '2015052616', '2015052617']

    individual_equity(reliability_fpath=reliability_fpath,
                    reliability_prefix=reliability_prefix,
                    redundancy_fpath=redundancy_fpath,
                    redundancy_extension=redundancy_extension,
                    recoverability_fpath=recoverability_fpath,
                    recoverability_init_time=recoverability_init_time,
                    recoverability_prefix=recoverability_prefix,
                    recoverability_times=recoverability_times,
                    svi=svi,
                    time='2015052522',
                    methods=['SVI_scaled', 'Factor 1', 'Factor 2'],
                    svi_weights=[.25, .5, 1])

    plt.tight_layout()
    plt.subplots_adjust(top=0.93,
                        bottom=0.12,
                        left=0.159,
                        right=0.954,
                        hspace=0.45,
                        wspace=0.29)

# ind_equity_func()
# plt.show()

##############################################
# NETWORK RECOVERABILITY
def net_recoverability_func():
    times = ['2015052522', '2015052523',
            '2015052600', '2015052601',
            '2015052602',
            '2015052604', '2015052605',
            '2015052606', '2015052607',
            '2015052608', '2015052609',
            '2015052610', '2015052611',
            '2015052612', '2015052613',
            '2015052614', '2015052615',
            '2015052616', '2015052617']
    ftypes = ['inundation' , 'inundation_fluvial', 'inundation_pluvial']

    # Aggregate individual reliability
    fig, axes = plt.subplots(1, 1, figsize=(85*mm, 85*mm))
    fig.set_facecolor('none')
    fig.supxlabel('SVI Percentile',
                weight='bold',
                x=0.5,
                fontsize=12)
    fig.supylabel('Norm. Cumsum. Network Recovery',
                weight='bold',
                fontsize=12)

    network_recoverability(axes=[[axes]],
                            fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs',
                            prefix='AN_Graph',
                            svi=svi,
                            init_time='2015052518',
                            method='SVI_scaled',
                            variable='agr_no_cap',
                            times=times,
                            annotations=True,
                        ftypes=ftypes)
    plt.tight_layout()
    plt.subplots_adjust(top=0.938,
                        bottom=0.162,
                        left=0.214,
                        right=0.902,
                        hspace=0.2,
                        wspace=0.2)

# net_recoverability_func()
# plt.show()

##########################################
# NETWORK RECOVERABILITY CUMULATIVE TIME SERIES
def network_recoverability_ts_func():
    times = ['2015052522', '2015052523',
                '2015052600', '2015052601',
                '2015052602', '2015052603',
                '2015052604', '2015052605',
                '2015052606', '2015052607',
                '2015052608', '2015052609',
                '2015052610', '2015052611',
                '2015052612', '2015052613',
                '2015052614', '2015052615',
                '2015052616', '2015052617']
    ftypes = ['inundation' , 'inundation_fluvial', 'inundation_pluvial']
    percs=[0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    network_recovoverability_cumulative_ts(fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs',
                                           prefix='AN_Graph',
                                           svi=svi,
                                           init_time='2015052518',
                                           variable='agr_no_cap',
                                           times=times,
                                           percs=percs,
                                           ftypes=ftypes)

# network_recoverability_ts_func()
# plt.show()

#############################################
# HOUSEHOLD RECOVERABILITY
def house_recoverability_func():
    times = ['2015052522', '2015052523',
            '2015052600', '2015052601',
            '2015052602', '2015052603',
            '2015052604', '2015052605',
            '2015052606', '2015052607',
            '2015052608', '2015052609',
            '2015052610', '2015052611',
            '2015052612', '2015052613',
            '2015052614', '2015052615',
            '2015052616', '2015052617']
    # Aggregate individual reliability
    fig, axes = plt.subplots(1, 1, figsize=(85*mm, 85*mm))
    fig.set_facecolor('none')
    # Graph attributes when multiple plots

    fig.supxlabel('Indicator Percentile',
                weight='bold',
                x=0.5,
                fontsize=12)
    fig.supylabel('Norm. Cumsum. Residential Recovery',
                weight='bold',
                fontsize=12)

    individual_recoverability(axes=[[axes]],
                            svi=svi,
                            fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Reliability/individual_reliability',
                            init_time=2015052518,
                            times=times,
                            prefix='agg_individual_reliability',
                            columns=1,
                            methods=['SVI_scaled', 'Factor 1', 'Factor 2'],
                            perc=0.95)
    plt.tight_layout()
    plt.subplots_adjust(top = 0.938,
                        bottom = 0.162,
                        left = 0.214,
                        right = 0.902,
                        hspace = 0.2,
                        wspace = 0.2)

# house_recoverability_func()
# plt.show()

###########################################################################################
# HOUSEHOLD RECOVERABILTY CUMULATIVE TIME SERIES
def house_recoverability_ts_func():
    times = ['2015052522', '2015052523',
            '2015052600', '2015052601',
            '2015052602', '2015052603',
            '2015052604', '2015052605',
            '2015052606', '2015052607',
            '2015052608', '2015052609',
            '2015052610', '2015052611',
            '2015052612', '2015052613',
            '2015052614', '2015052615',
            '2015052616', '2015052617']
    individual_recoverability_cumulative_ts(svi=svi,
                                            fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Reliability/individual_reliability',
                                            init_time=2015052518,
                                            times=times,
                                            prefix='agg_individual_reliability',
                                            percs=[0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    
# house_recoverability_ts_func()
# plt.show()

#############################################################################################
# # # INDIVIDUAL REDUNDANCY
def house_redundancy_func():
    methods = ['SVI_scaled', 'Factor 1', 'Factor 2']
    cost_atr = 'redundancy'
    fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Redundancy/individual_redundancy/'

    times = ['2015052518','2015052520','2015052522',
            '2015052600','2015052602','2015052604',
            '2015052606','2015052608','2015052610']

    extension='individual_redundancy_'

    fig, axes = plt.subplots(3,3,figsize=(180*mm,200*mm))
    fig.set_facecolor('none')
    # Graph attributes when multiple plots
    # fig.suptitle('Lorenz Curves of Residential Redundancy',
    #              fontsize=18, weight='bold',style='italic')
    fig.supxlabel('Indicator Percentile',
                weight='bold',
                x=0.5,
                fontsize=16)
    fig.supylabel('Normalized Cumulative Sum of Household Redundancy',
                weight='bold',
                fontsize=16)

    individual_redundancy(axes=axes,
                        svi=svi,
                        fpath=fpath,
                        times=times,
                        extension=extension,
                        columns=3,
                        methods=methods,
                        atr=cost_atr,
                        annotations=True)
    plt.tight_layout()

# house_redundancy_func()    
# plt.show()

#############################################################################################
# NETWORK REDUNDANCY
def net_redundancy_func():
    method = 'SVI_scaled'
    variable='net_redun'
    fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Redundancy/network_redundancy'

    times=['2015052518','2015052522','2015052602']

    fig, axes = plt.subplots(1,3,figsize=(180*mm,85*mm))
    fig.set_facecolor('none')
    # Graph attributes when multiple plots
    # fig.suptitle('Lorenz Curves of Network Redundancy',
    #              fontsize=14, weight='bold')
    fig.supxlabel('SVI Percentile',
                weight='bold',
                x=0.5)
    fig.supylabel('Norm. Cumsum. of Network Red.',
                weight='bold')


    network_redundancy(axes=[axes],
                    fpath=fpath,
                    svi=svi,
                       method=method,
                       variable=variable,
                    times=times,
                    annotations=True,
                    columns=3,
                    ftypes=['_', '_fluvial_', '_pluvial_'])

    plt.tight_layout()
    # layout for single plot
    plt.subplots_adjust(top=0.91,
                        bottom=0.135,
                        left=0.175,
                        right=0.91,
                        hspace=0.2,
                        wspace=0.2)

# net_redundancy_func()
# plt.show()

#############################################################################################
# # INDIVIDUAL RELIABILITY
def house_reliability_func():
    # times
    times = ['2015052518','2015052520','2015052522',
            '2015052600','2015052602','2015052604',
            '2015052606','2015052608','2015052610']

    # Aggregate individual reliability
    fig, axes = plt.subplots(3, 3, figsize=(180*mm, 200*mm))
    fig.set_facecolor('none')
    # Graph attributes when multiple plots
    # fig.suptitle('Lorenz Curves of Residential Accessability',
    #              fontsize=18, weight='bold', style='italic')
    fig.supxlabel('Indicator Percentile',
                weight='bold',
                x=0.5,
                fontsize=16)
    fig.supylabel('Normalized Cumulative Sum of Household Reliability',
                weight='bold',
                fontsize=16)


    individual_reliability_aggregate(axes=axes,
                                    fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Reliability/individual_reliability',
                                    times=times,
                                    prefix='agg_individual_reliability',
                                    columns=3,
                                    methods=['SVI_scaled', 'Factor 1','Factor 2'])


    plt.tight_layout()
    plt.subplots_adjust(top=0.958,
                        bottom=0.11,
                        left=0.119,
                        right=0.954,
                        hspace=0.45,
                        wspace=0.29)

# house_reliability_func()
# plt.show()

#############################################################################################
# # NETWORK RELIABILITY
def net_reliability_func():
    # Impacted roads time series
    fig, axes = plt.subplots(3, 3, figsize=(180*mm, 200*mm))
    fig.set_facecolor('none')
    # fig.suptitle('Lorenz Curves of Impacted Roads',
    #              fontsize=14, weight='bold')
    fig.supxlabel('Social Vulnerability Index Percentile',
                weight='bold',
                x=0.5,
                fontsize=14)
    fig.supylabel('Normalized Cumulative Sum of Network Reliability',
                weight='bold',
                fontsize=14)

    times = ['2015052518','2015052520','2015052522',
            '2015052600','2015052602','2015052604',
            '2015052606','2015052608','2015052610']

    ftypes = ['inundation' , 'inundation_fluvial', 'inundation_pluvial']

    network_reliability(axes=axes,
                        fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs',
                                prefix='AN_Graph',
                                svi=svi,
                                method='SVI_scaled',
                                variable='agr_no_cap',
                                times=times,
                                annotations=True,
                                columns=3,
                                ftypes=ftypes)

    plt.tight_layout()
    plt.subplots_adjust(top=0.958,
                        bottom=0.11,
                        left=0.119,
                        right=0.954,
                        hspace=0.45,
                        wspace=0.29)

# net_reliability_func()
# plt.show()

def net_reliability_func_individual():
    # network reliability single plot
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))
    fig.set_facecolor('none')
    # fig.suptitle('Lorenz Curves of Impacted Roads',
    #              fontsize=14, weight='bold')
    axes.set_xlabel('Social Vulnerability Index Percentile',
                weight='bold',
                x=0.5,
                fontsize=16)
    axes.set_ylabel('Normalized Cumulative Sum of Closed Roads',
                weight='bold',
                fontsize=16)

    times = ['2015052602']
    ftypes = ['inundation', 'inundation_fluvial', 'inundation_pluvial']

    network_reliability(axes=[[axes]],
                        fpath='/home/mdp0023/Desktop/external/Data/Network_Data/Austin_North/AN_Graphs',
                        prefix='AN_Graph',
                        svi=svi,
                        method='SVI_scaled',
                        variable='agr_no_cap',
                        times=times,
                        annotations=True,
                        columns=1,
                        ftypes=ftypes)

    plt.tight_layout()
    plt.subplots_adjust(top=0.95,
                        bottom=0.088,
                        left=0.076,
                        right=0.973,
                        hspace=0.206,
                        wspace=0.255)

# net_reliability_func_individual()
# plt.show()







# functions no longer in use
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
    for a, time in enumerate(times):
        # create variable to hold information for percentile table
        table_vals = [[0, 0, 0, 0] for _ in range(len(ftypes))]
        # create a variable to hold information regarding gini coefficeints
        ginis = [0 * len(ftypes)]
        for b, ftype in enumerate(ftypes):

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
            svi['agr_increased_travel_time'] = svi['agr_increased_travel_time'].fillna(
                0)

            # count the number of roads impacted (increased travel time & 0 capacity)
            svi['agr_impact'] = svi['agr_increased_travel_time'] + svi['agr_no_cap']

            # calc percentage of roads within a BG with no capacity
            svi['perc_no_cap'] = svi['agr_no_cap']/svi['count']*100

            # calc percentage of roads within a BG with an increased travel time
            svi['perc_increased_travel_time'] = svi['agr_increased_travel_time'] / \
                svi['count']*100

            # calc percentage of roads within a BG impacted
            svi['perc_impact'] = svi['agr_impact']/svi['count']*100

            # sort svi dataframe by appropriate column
            if method == 'Factor 2':
                svi.sort_values('Factor 2', axis=0,
                                inplace=True, ascending=False)
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
            lorenz_plot, = axes[row, column].plot(
                val, lorenz, label='tsamp', color=c_line, linewidth=2)

            perc_25 = np.percentile(lorenz, 25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz, 75)

            table_vals[b][0] = f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3] = f'{100-round(perc_75*100)}%'

            # calc gini and save
            ginis[b] = gini(array)

        col_labels = ['Q1', 'Q2', 'Q3', 'Q4']
        row_labels = []
        for ftype in ftypes:
            if ftype == 'inundation':
                row_labels.append('Compound')
            elif ftype == 'inundation_fluvial':
                row_labels.append('Fluvial')
            elif ftype == 'inundation_pluvial':
                row_labels.append('Pluvial')

        # the rectangle is where I want to place the table
        if annotations is True:

            # plot ginis
            textstr = []
            for i, row_label in enumerate(row_labels):
                textstr.append(f'{row_label} Gini: {round(ginis[i],2)}')
            textstr = '\n'.join(textstr)

            # Gini patch properties
            props = dict(boxstyle='square', facecolor='lightgray', alpha=0.5)
            # text box in upper left in axes coords
            axes[row, column].text(0.05, 0.95, textstr, transform=axes[row, column].transAxes, fontsize=14,
                                   verticalalignment='top', bbox=props)

            # Plot quartile information
            the_table = axes[row, column].table(cellText=table_vals,
                                                colWidths=[0.12]*4,
                                                rowColours=[
                                                    '#cc57a4', '#ff772e', '#4d88b8'],
                                                colColours=['lightgrey', 'lightgrey',
                                                            'lightgrey', 'lightgrey'],
                                                rowLabels=row_labels,
                                                colLabels=col_labels,
                                                loc='lower right')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(9)

        # Equality line
        eq_line, = axes[row, column].plot(
            [0, 1], [0, 1], color='k', linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes[row, column].set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                                     '0%', '25%', '50%', '75%', '100%'])
        axes[row, column].set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                                     '0%', '25%', '50%', '75%', '100%'])

        # set x and y limits
        axes[row, column].set_xlim(0, 1)
        axes[row, column].set_ylim(0, 1)

        # set axes title
        axes[row, column].set_title(
            f'{time[-2:]}:00 {time[6:8]}/{time[4:6]}/{time[:4]}', style='italic')

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
    for a, time in enumerate(times):
        ginis = [0, 0, 0]
        table_vals = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
        for b, ftype in enumerate(ftypes):

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
            svi['agr_increased_travel_time'] = svi['agr_increased_travel_time'].fillna(
                0)

            # count the number of roads impacted (increased travel time & 0 capacity)
            svi['agr_impact'] = svi['agr_increased_travel_time'] + svi['agr_no_cap']

            # calc percentage of roads within a BG with no capacity
            svi['perc_no_cap'] = svi['agr_no_cap']/svi['count']*100

            # calc percentage of roads within a BG with an increased travel time
            svi['perc_increased_travel_time'] = svi['agr_increased_travel_time'] / \
                svi['count']*100

            # calc percentage of roads within a BG impacted
            svi['perc_impact'] = svi['agr_impact']/svi['count']*100

            # sort svi dataframe by appropriate column
            if method == 'Factor 2':
                svi.sort_values('Factor 2', axis=0,
                                inplace=True, ascending=False)
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
            lorenz_plot, = axes.plot(
                val, lorenz, label='tsamp', color=c_line, linewidth=2)

            perc_25 = np.percentile(lorenz, 25)
            perc_50 = np.percentile(lorenz, 50)
            perc_75 = np.percentile(lorenz, 75)

            table_vals[b][0] = f'{round(perc_25*100)}%'
            table_vals[b][1] = f'{round(perc_50*100)-round(perc_25*100)}%'
            table_vals[b][2] = f'{round(perc_75*100)-round(perc_50*100)}%'
            table_vals[b][3] = f'{100-round(perc_75*100)}%'

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
                                   rowColours=['#cc57a4',
                                               '#ff772e', '#4d88b8'],
                                   rowLabels=row_labels,
                                   colLabels=col_labels,
                                   colColours=['lightgrey', 'lightgrey',
                                               'lightgrey', 'lightgrey'],
                                   loc='lower right')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(14)
            the_table.scale(1, 1.5)

        # Equality line
        eq_line, = axes.plot([0, 1], [0, 1], color='k',
                             linestyle='--', label='Equality Line')

        # Set x and y tick markers
        axes.set_xticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                        '0%', '25%', '50%', '75%', '100%'], fontsize=14)
        axes.set_yticks(ticks=[0.0, 0.25, 0.5, 0.75, 1.0], labels=[
                        '0%', '25%', '50%', '75%', '100%'], fontsize=14)

        # set x and y limits
        axes.set_xlim(0, 1)
        axes.set_ylim(0, 1)

        # set axes title
        axes.set_title(
            f'{time[-2:]}:00 {time[6:8]}/{time[4:6]}/{time[:4]}', style='italic', fontsize=14)


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
    if aggregate is True:
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
        if aggregate is True:
            all_arrays = np.add(array, all_arrays)

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

    idx += 1
    if aggregate is True:
        # add the final aggregate plot, summing across all of the resources
        res_names.append('Aggregate')
        lorenz_curve_multiplot(X=all_arrays,
                               row=row,
                               column=column,
                               axes=axes,
                               idx=idx,
                               names=res_names,
                               method=method,
                               annotations=annotations)


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
    
