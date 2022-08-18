# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import networkx as nx
class network(object):
    """ Define the network class for modeling the power and gas networks
    """
    def __init__(self, node_data, arc_data):
        """
        Input：
            node_data: pandas.dataframe, the data of the nodes
            arc_data: pandas.dataframe, the data of the arcs
        """
        #read the node data
        self.nodenum = len(node_data)
        self.demand = node_data['demand']
        self.supply_cap = node_data['supply_cap']
        self.supply = node_data['supply']
        self.nodeid = node_data['node_id']
        
        
        #read the arc data
        self.arcnum = len(arc_data)
        self.start_node_id = arc_data['start_node']
        self.end_node_id = arc_data['end_node']
        self.flow_cap = arc_data['flow_cap']
        self.conv_rate = arc_data['conv_rate']
        self.flow = arc_data['flow']
        
        self.nodeid2num()
        self.adj_matrix()
        self.adj_list()
        self.networkx_graph()
        self.topo_sort()
        self.flow_matrix()
        # self.flow_check()
        self.centrality()
    def nodeid2num(self):
        """ Mapping the ID of the node to the number of the node in the network
        """
        self.ID2num = {}
        for i in range(self.nodenum):
            self.ID2num[self.nodeid.iloc[i]] = i
    
    def adj_matrix(self):
        """ Create the adjacency matrix of the network, directed graph
        """
        self.adjmatrix = np.zeros((self.nodenum, self.nodenum), dtype = int)
        
        for i in range(self.arcnum):
            self.adjmatrix[self.ID2num[self.start_node_id.iloc[i]], self.ID2num[self.end_node_id.iloc[i]]] = 1
    
    def flow_matrix(self):
        """ Flor matrix of the network, directed graph
        """
        self.flowmatrix = np.zeros((self.nodenum, self.nodenum), dtype = float)
        
        for i in range(self.arcnum):
            self.flowmatrix[self.ID2num[self.start_node_id.iloc[i]], self.ID2num[self.end_node_id.iloc[i]]] = self.flow.iloc[i]
    
    def flow_check(self):
        """ Check whether the initial flow is balanced
        """
        for i in range(self.nodenum):
            print(i + 1, np.sum(self.flowmatrix[:, i]) - np.sum(self.flowmatrix[i, :]) + self.supply.iloc[i] - self.demand.iloc[i])
    
    def adj_list(self):
        """ Create the adjacency list of the network, directed graph
        """
        self.adjlist = {}
        for i in range(self.arcnum):
            if(self.ID2num[self.start_node_id.iloc[i]] not in self.adjlist.keys()):
                self.adjlist[self.ID2num[self.start_node_id.iloc[i]]] = []
            
            self.adjlist[self.ID2num[self.start_node_id.iloc[i]]].append(self.ID2num[self.end_node_id.iloc[i]])
    
    def networkx_graph(self):
        """ Create the networkx object of the network
        """
        self.graph = nx.convert_matrix.from_numpy_matrix(self.adjmatrix, create_using=nx.DiGraph)
    
    def topo_sort(self):
        """ Perform the topological sort of the network, directed graph
        """
        if(nx.algorithms.dag.is_directed_acyclic_graph(self.graph)):
            self.topo_order = list(nx.topological_sort(self.graph))
        else:
            print('The current network is not the DAG')
    
    def centrality(self):
        """ Calculate the centrality of the graph
        Choice: degree, katz, closeness, betweenness
        """
        self.dc = list(nx.algorithms.centrality.degree_centrality(self.graph).values())
        self.kc = list(nx.algorithms.centrality.katz_centrality(self.graph).values())
        self.cc = list(nx.algorithms.centrality.closeness_centrality(self.graph).values())
        self.bc = list(nx.algorithms.centrality.betweenness_centrality(self.graph).values())
