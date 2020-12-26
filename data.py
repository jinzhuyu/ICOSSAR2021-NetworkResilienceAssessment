# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 17:43:02 2020

@author: 10624
"""

""" This script contains functions for preprocessing the data
"""
import pandas as pd
#import numpy as np
#import networkx as nx
import os

dir_path = os.path.dirname(os.path.realpath(__file__)) #Change the directory to the one where the current script is located

nodedata = pd.read_csv('./data/nodes_data.csv')
arcdata = pd.read_csv('./data/arcs_data.csv')

p_nodenum = 24
g_nodenum = 25
p_nodedata = nodedata.iloc[0:p_nodenum, :]
g_nodedata = nodedata.iloc[p_nodenum:(p_nodenum + g_nodenum), :]

p_arcnum = 34
g_arcnum = 24
p2g_arcnum = 4
#g2p_arcnum = 1
p_arcdata = arcdata.iloc[0:p_arcnum, :]
g_arcdata = arcdata.iloc[p_arcnum:(p_arcnum + g_arcnum), :]
g2p_arcdata = arcdata.iloc[(p_arcnum + g_arcnum):(p_arcnum + g_arcnum + p2g_arcnum), :]
#p2g_arcdata = arcdata.iloc[(p_arcnum + g_arcnum + p2g_arcnum):(p_arcnum + g_arcnum + p2g_arcnum + g2p_arcnum), :]

# power to gas link is removed for now because the current flow redistribution cannot handle cycles with flow


