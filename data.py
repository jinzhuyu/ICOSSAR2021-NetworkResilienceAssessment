# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 17:43:02 2020

@author: 10624
"""

""" This script contains functions for preprocessing the data
"""
import pandas as pd
import numpy as np
import networkx as nx
import os
import os 
dir_path = os.path.dirname(os.path.realpath(__file__)) #Change the directory to the one where the current script is located

nodedata = pd.read_csv('./data/nodes_data.csv')
arcdata = pd.read_csv('./data/arcs_data.csv')

p_nodenum = 24
g_nodenum = 25
p_nodedata = nodedata.iloc[0:p_nodenum, :]
g_nodedata = nodedata.iloc[p_nodenum:(p_nodenum + g_nodenum), :]

p_edgenum = 34
g_edgenum = 24
p2g_edgenum = 4
g2p_edgenum = 1
p_edgedata = arcdata.iloc[0:p_edgenum, :]
g_edgedata = arcdata.iloc[p_edgenum:(p_edgenum + g_edgenum), :]
g2p_edgedata = arcdata.iloc[(p_edgenum + g_edgenum):(p_edgenum + g_edgenum + p2g_edgenum), :]
p2g_edgedata = arcdata.iloc[(p_edgenum + g_edgenum + p2g_edgenum):(p_edgenum + g_edgenum + p2g_edgenum + g2p_edgenum), :]


