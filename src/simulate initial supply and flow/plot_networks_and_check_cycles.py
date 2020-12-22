# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 15:49:00 2020

@author: yuj5
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

import os
dir_path = os.path.dirname(os.path.abs(__file__))    # go to the directory of the current script

import networkx as nx
from matplotlib import patches 

# set default plot parameters
def set_default_plot_param():
    
    plt.style.use('classic')
    
    plt.rcParams["font.family"] = "Helvetica"
    plt.rcParams['font.weight']= 'normal'
    plt.rcParams['figure.figsize'] = [6, 6*3/4]
   
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    
    plt.rc('axes', titlesize=14, labelsize=12, linewidth=0.75)    # fontsize of the axes title, the x and y labels
    
    plt.rc('lines', linewidth=1.5, markersize=4)
    
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    plt.rcParams['axes.formatter.useoffset'] = False # turn off offset
    # To turn off scientific notation, use: ax.ticklabel_format(style='plain') or
    # plt.ticklabel_format(style='plain')

    
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams["legend.fancybox"] = True
    plt.rcParams["legend.loc"] = "best"
    plt.rcParams["legend.framealpha"] = 0.5
    
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['savefig.dpi'] = 800
    
#    plt.rc('text', usetex=False)
    
set_default_plot_param()


#%%
'''plot networks in which the weights show the restorative importance, flow more.
   ref.: https://qxf2.com/blog/drawing-weighted-graphs-with-networkx/
'''

def min_pos_scale(data_df, scale_max=1, scale_min=-1):
    '''scale data into the desired range
    '''
    for i in np.arange(data_df.shape[1]):
        # for y position, the range should be -0.4, 0.4
        if i==1:
            scale_max, scale_min = 0.4, -0.4            
        data_max = data_df.iloc[:,i].max()
        data_min = data_df.iloc[:,i].min()
        slope = (scale_max-scale_min)/(data_max-data_min)
        intercept = scale_min - slope*data_min
        data_df.iloc[:,i] = intercept + data_df.iloc[:,i]*slope
    
    return data_df

def plot_weighted_graph(node_df, arc_df, node_pos, fig_size = (16, 16*4/11)):
    '''draw a weighted graph according to the defined node positions
    '''
    
    # 0.0 extract nodes label, net_id, etc.
    node_list = node_df['node_id'].tolist()
    net_id = node_df['net_id'].tolist()
    n_node_power = net_id.count(1)
    n_node_total = node_df.shape[0]
    
    # 0.1 extract arcs label
    start_node_id = arc_df['start_node']
    end_node_id = arc_df['end_node']
    n_arc_total = arc_df.shape[0]
 
    # 1.0 draw basic graph and nodes
    G = nx.DiGraph(directed=True)
    for i in np.arange(n_node_total):
        G.add_node(i, pos=(node_pos.iloc[i,0], node_pos.iloc[i,1]))
 
    # 1.1 draw nodes
    pos=nx.get_node_attributes(G,'pos')
    plt.figure(figsize=fig_size) 
    nx.draw_networkx_nodes(G, pos, nodelist=range(0, n_node_power),
                           node_color='royalblue', node_shape='o', node_size=500)
    nx.draw_networkx_nodes(G, pos, nodelist=range(n_node_power, n_node_total),
                           node_color='tab:red', node_shape='o', node_size=500)

    # 1.2 remove edge of nodes      
    ax = plt.gca() # to get the current axis
    for i in np.arange(max(net_id)):
        ax.collections[i].set_edgecolor('none')
    
    # 1.3. add labels to the nodes
    labels = {}
    for j in np.arange(n_node_total):
        labels[j] = node_list[j]
    nx.draw_networkx_labels(G, pos, labels, font_size=9)  #, font_family='serif')
  
    # 2.0 add edges
    for k in np.arange(n_arc_total):
        start_node_num = list(labels.values()).index(start_node_id[k])
        end_node_num = list(labels.values()).index(end_node_id[k])
        G.add_edge(start_node_num, end_node_num, weight=1.5, color='k')
    
    all_weights = []
    #4.0 Iterate through the graph nodes to gather all the weights
    for (node1,node2,data) in G.edges(data=True):
        all_weights.append(data['weight']) #we'll use this when determining edge thickness
 
    #4.1 Get unique weights
    unique_weights = list(set(all_weights))
 
    #4.2 Plot the edges - one by one!
#    testArrow = patches.ArrowStyle.Fancy(head_length=.4, head_width=.4, tail_width=.1)
    for weight in unique_weights:
        #4 d. Form a filtered list with just the weight you want to draw
        weighted_edges = [(n_1,n_2) for (n_1,n_2,edge_attr) in G.edges(data=True) if edge_attr['weight']==weight]
        width = 0.75
        nx.draw_networkx_edges(G, pos, edgelist=weighted_edges, width=width,
                               arrowsize=15)  #, alpha=0.75)
 
    #Plot the graph
    plt.axis('off')
    plt.show()
    
    return G
    

def main():   
    
    # import network data    
    node_df = pd.read_csv('./data/nodes_data.csv')
    arc_df = pd.read_csv('./data/arcs_data.csv')
    
    # import node position data
    node_pos_df = pd.read_csv('./data/node_position.csv')   
    
    # plot
    G = plot_weighted_graph(node_df, arc_df, node_pos=node_pos_df[['pos_x', 'pos_y']]*2)
    
    # check cycles
    cycles = nx.find_cycle(G)
    print(cycles)

    
main()