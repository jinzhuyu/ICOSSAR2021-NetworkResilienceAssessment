# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 22:47:00 2020

@author: 10624
"""
import data as dt
from infrasnetwork import network
import failsimulation as fs
import copy
import numpy as np
import networkx as nx
import math
import random

class system(object):
    """ Couple the power and gas networks
    """
    def __init__(self, power, gas, g2p_edgedata): ##Ignore the edge from the power network to gas network
        self.power = power
        self.gas = gas
        self.nodenum = self.power.nodenum + self.gas.nodenum
        self.edgenum = self.power.edgenum + self.gas.edgenum + len(g2p_edgedata)
        
        self.g2p_edgedata = g2p_edgedata
        
        self.adjflow_matrix()
        self.nodeid2num()
        self.networkx_graph()
        self.topo_sort()
        self.centrality()
        self.propertymapping()
    
    def adjflow_matrix(self): #first gas then power
        """ Create the adjacency matrix of the whole system by plugging the adjacency matrix of the power and gas networks 
            with the addition of the interdependent links
        """
        
        self.adjmatrix = np.zeros((self.nodenum, self.nodenum), dtype = int)
        self.flowmatrix = np.zeros((self.nodenum, self.nodenum), dtype = int)
        
        #add the power and gas edges
        self.adjmatrix[0:self.gas.nodenum, 0:self.gas.nodenum] = copy.copy(self.gas.adjmatrix)
        self.flowmatrix[0:self.gas.nodenum, 0:self.gas.nodenum] = copy.copy(self.gas.flowmatrix)
        self.adjmatrix[self.gas.nodenum:(self.gas.nodenum + self.power.nodenum), self.gas.nodenum:(self.gas.nodenum + self.power.nodenum)] = copy.copy(self.power.adjmatrix)
        self.flowmatrix[self.gas.nodenum:(self.gas.nodenum + self.power.nodenum), self.gas.nodenum:(self.gas.nodenum + self.power.nodenum)] = copy.copy(self.power.flowmatrix)
        
        #add the gas2power edges
        for i in range(len(self.g2p_edgedata)):
            self.adjmatrix[self.gas.ID2num[self.g2p_edgedata['start_node'].iloc[i]], self.power.ID2num[self.g2p_edgedata['end_node'].iloc[i]] + self.gas.nodenum] = 1
            self.flowmatrix[self.gas.ID2num[self.g2p_edgedata['start_node'].iloc[i]], self.power.ID2num[self.g2p_edgedata['end_node'].iloc[i]] + self.gas.nodenum] = self.g2p_edgedata['flow'].iloc[i]
    
    def nodeid2num(self):
        """ Mapping the ID of the node to the number of the node in the system
        """
        ##update the power ID
        update_powerID = {x: self.power.ID2num[x]+self.gas.nodenum for x in self.power.ID2num}
        self.ID2num = {**self.gas.ID2num, **update_powerID}
        self.num2ID = dict(map(reversed, self.ID2num.items()))
    
    def networkx_graph(self):
        """ Create the networkx object of the network
        """
        
        self.graph = nx.convert_matrix.from_numpy_matrix(self.adjmatrix, create_using = nx.DiGraph)
    
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
        
    def initial_failure(self, Type, ratio):
        """ Simulate the initial failure sequence
        Input:
            Type - the type of the initial failure sequence, choice:
                   'randomness', 'dc' - degree centrality, 'bc' - betweenness centrality, 'kc' - katz centrality, 'cc': closeness centrality
            ratio - how much percentage of nodes is failed
        Output:
            the initial failure sequence
        """
        fail_num = math.floor(ratio*self.nodenum)
        
        if(Type == 'randomness'):
            self.initial_fail_seq = random.sample(range(self.nodenum), fail_num)
        else:
            exec('self.initial_fail_seq = np.argsort(self.{})[-fail_num:]'.format(Type))
    
    def propertymapping(self):
        """ Mapping all properties in the network to the system
        """
        self.demand = np.concatenate((np.array(self.gas.demand), np.array(self.power.demand)))
        self.supply_cap = np.concatenate((np.array(self.gas.supply_cap), np.array(self.power.supply_cap)))
        self.supply = np.concatenate((np.array(self.gas.supply), np.array(self.power.supply)))
        
        
        self.flowcapmatrix = np.zeros((self.nodenum, self.nodenum), dtype = float)
        self.convratematrix = np.zeros((self.nodenum, self.nodenum), dtype = float)
        for i in range(self.gas.edgenum):
            self.flowcapmatrix[self.ID2num[self.gas.start_node_id.iloc[i]], self.ID2num[self.gas.end_node_id.iloc[i]]] = self.gas.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.gas.start_node_id.iloc[i]], self.ID2num[self.gas.end_node_id.iloc[i]]] = self.gas.conv_rate.iloc[i]
            
        for i in range(self.power.edgenum):
            self.flowcapmatrix[self.ID2num[self.power.start_node_id.iloc[i]], self.ID2num[self.power.end_node_id.iloc[i]]] = self.power.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.power.start_node_id.iloc[i]], self.ID2num[self.power.end_node_id.iloc[i]]] = self.power.conv_rate.iloc[i]
            
        for i in range(len(self.g2p_edgedata)):
            self.flowcapmatrix[self.ID2num[self.g2p_edgedata.start_node.iloc[i]], self.ID2num[self.g2p_edgedata.end_node.iloc[i]]] = self.g2p_edgedata.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.g2p_edgedata.start_node.iloc[i]], self.ID2num[self.g2p_edgedata.end_node.iloc[i]]] = self.g2p_edgedata.conv_rate.iloc[i]
        
    
    
power = network(dt.p_nodedata, dt.p_edgedata)
gas = network(dt.g_nodedata, dt.g_edgedata)
s = system(power, gas, dt.g2p_edgedata)
s.initial_failure('randomness', 0.1) #the type of the initial failure sequence, choice: 'randomness', 'dc' - degree centrality, 'bc' - betweenness centrality, 'kc' - katz centrality, 'cc': closeness centrality




