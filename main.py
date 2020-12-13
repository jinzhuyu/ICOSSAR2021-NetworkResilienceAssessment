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
        
    def propertymapping(self):
        """ Mapping all properties in the network to the system
        """
        self.demand = np.concatenate((np.array(self.gas.demand), np.array(self.power.demand)))
        self.supply_cap = np.concatenate((np.array(self.gas.supply_cap), np.array(self.power.supply_cap)))
        self.supply = np.concatenate((np.array(self.gas.supply), np.array(self.power.supply)))
        self.edgelist = np.empty((self.power.edgenum + self.gas.edgenum + len(self.g2p_edgedata), 2), dtype = int)
        
        
        self.flowcapmatrix = np.zeros((self.nodenum, self.nodenum), dtype = float)
        self.convratematrix = np.zeros((self.nodenum, self.nodenum), dtype = float)
        temp = 0
        for i in range(self.gas.edgenum):
            self.flowcapmatrix[self.ID2num[self.gas.start_node_id.iloc[i]], self.ID2num[self.gas.end_node_id.iloc[i]]] = self.gas.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.gas.start_node_id.iloc[i]], self.ID2num[self.gas.end_node_id.iloc[i]]] = self.gas.conv_rate.iloc[i]
            self.edgelist[temp, 0], self.edgelist[temp, 1] = self.ID2num[self.gas.start_node_id.iloc[i]], self.ID2num[self.gas.end_node_id.iloc[i]]
            temp += 1
            
        for i in range(self.power.edgenum):
            self.flowcapmatrix[self.ID2num[self.power.start_node_id.iloc[i]], self.ID2num[self.power.end_node_id.iloc[i]]] = self.power.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.power.start_node_id.iloc[i]], self.ID2num[self.power.end_node_id.iloc[i]]] = self.power.conv_rate.iloc[i]
            self.edgelist[temp, 0], self.edgelist[temp, 1] = self.ID2num[self.power.start_node_id.iloc[i]], self.ID2num[self.power.end_node_id.iloc[i]]
            temp += 1
            
        for i in range(len(self.g2p_edgedata)):
            self.flowcapmatrix[self.ID2num[self.g2p_edgedata.start_node.iloc[i]], self.ID2num[self.g2p_edgedata.end_node.iloc[i]]] = self.g2p_edgedata.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.g2p_edgedata.start_node.iloc[i]], self.ID2num[self.g2p_edgedata.end_node.iloc[i]]] = self.g2p_edgedata.conv_rate.iloc[i]
            self.edgelist[temp, 0], self.edgelist[temp, 1] = self.ID2num[self.g2p_edgedata.start_node.iloc[i]], self.ID2num[self.g2p_edgedata.end_node.iloc[i]]
            temp += 1
    
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
        self.initial_fail_seq_onehotcode = np.zeros(self.nodenum, dtype = int) #1 - failure, 0 - survive
        
        if(Type == 'randomness'):
            self.initial_fail_seq = random.sample(range(self.nodenum), fail_num)
        else:
            exec('self.initial_fail_seq = np.argsort(self.{})[-fail_num:]'.format(Type))
            
        self.initial_fail_seq_onehotcode[np.array(self.initial_fail_seq)] = 1
        self.initial_fail_link_onehotcode = np.zeros(self.edgenum, dtype = int) #There are no initial failed links
        
    def update_matrix_node_fail(self, adjmatrix, node_fail_seq):
        """ Update the adjacency matrix caused by node failure
        Input:
            adjmatrix - the current adjacency matrix
            node_fail_seq - the failed node number at current time step, numpy1Darray
        """
        adjmatrix[node_fail_seq == 1, :] = 0
        adjmatrix[:, node_fail_seq == 1] = 0
        
        return adjmatrix
    
    def update_matrix_link_fail(self, adjmatrix, link_fail_seq):
        """ Update the adjacency matrix caused by link failure
        Input:
            adjmatrix - the current adjacency matrix
            link_fail_seq - the failed link number at current time step, numpy1Darray
        """
        for i in range(len(link_fail_seq)):
            if(link_fail_seq[i] == 1):
                adjmatrix[self.edgelist[i, 0], self.edgelist[i, 1]] = 0
        
        return adjmatrix
    
    def cascading_failure(self, alpha):
        """ Simulate the cascading failure
        Input:
            alpha: redundancy
        """
        #Update the initial failure scenario: update the adjmatrix and flowmatrix
        self.adjmatrix_init = copy.copy(self.adjmatrix)
        self.flowmatrix_init = copy.copy(self.flowmatrix)
        
        #Update the adjmatrix caused by failed nodes
        self.adjmatrix_init = self.update_matrix_node_fail(copy.copy(self.adjmatrix_init), self.initial_fail_seq_onehotcode)
        
        #Update the adjmatrix caused by failed links
        self.adjmatrix_init = self.update_matrix_link_fail(copy.copy(self.adjmatrix_init), self.initial_fail_link_onehotcode)
        
        #Update the flow matrix
        self.flowmatrix_init = self.flowmatrix*self.adjmatrix_init
        
        self.adjmatrix_evol, self.flowmatrix_evol = [self.adjmatrix, self.adjmatrix_init], [self.flowmatrix, self.flowmatrix_init]
        self.node_fail_evol, self.link_fail_evol = [self.initial_fail_seq_onehotcode], [self.initial_fail_link_onehotcode]
        self.satisfy_node_evol = [self.demand]
        self.performance = [1]
        self.node_fail_evol_track, self.link_fail_evol_track = [self.initial_fail_seq_onehotcode], [self.initial_fail_link_onehotcode]
        
        while(1): #perform the flow redistribution until stable
            satisfynode = copy.copy(self.satisfy_node_evol[-1])
            flowmatrix = copy.copy(self.flowmatrix_evol[-1])
            adjmatrix = copy.copy(self.adjmatrix_evol[-1])
        
            for node in self.topo_order:
                flowin = np.sum(flowmatrix[:, node]*self.convratematrix[:, node]) + self.supply[node]*self.node_fail_evol_track[-1][node] #calculate the total flow going into the node
                print(node, np.sum(flowmatrix[:, node]*self.convratematrix[:, node]))
                #Some flows serve for the node demand value, the remaining part goes to the following distribution process
                if(flowin >= self.demand[node]):
                    
                    satisfynode[node] = self.demand[node]
                    flowin -= self.demand[node]
                else:
                    satisfynode[node] = flowin/2 #if not enough to supply for the demand of the node, supply a half
                    flowin = flowin/2
                
                self.performance.append(np.sum(satisfynode)/np.sum(self.demand)) #Track down the performance
                
                #Redistribute the flow, here we can introduce some randomness to account for the uncertainty
                if(np.sum(adjmatrix[node, :]) == 0 or np.sum(self.flowmatrix_evol[-2][node, :]) == 0):
                    flowout = np.zeros(self.nodenum, dtype = float)
                else:
                    if(np.sum(adjmatrix[node, :]) != np.sum(self.adjmatrix_evol[-2][node, :])): # Some links fail, redistribute evenly with some random noise
                        flowout = 1/np.sum(adjmatrix[node, :])*flowin*adjmatrix[node, :]
                        
                        if(flowin != 0): #The uncertainty only happens when there are multiple out-links and inflow !=0
                            index = np.random.choice(np.argwhere(flowout != 0).reshape(-1)) #flow - beta where beta~U(0, flow/2)
                            if(len(np.argwhere(flowout != 0).reshape(-1)) != 1):
                                noise_flow = np.random.rand()*flowout[index]/2
                                unit_flow = noise_flow/(len(np.argwhere(flowout != 0)) - 1)
    
                                flowout[flowout != 0] = flowout[flowout != 0] + unit_flow
                                flowout[index] = flowout[index] - unit_flow - noise_flow
                        
                        
                    else: #No links fail,redistribute according to the ratio of flow at last time step
                        flowout = self.flowmatrix_evol[-2][node, :]/np.sum(self.flowmatrix_evol[-2][node, :])*flowin
                
                flowmatrix[node, :] = flowout
                
                node_seq_track = np.zeros(self.nodenum, dtype = int) #1 - failure, 0 - survive at current time step
                link_seq_track = np.zeros(self.edgenum, dtype = int) #1 - failure, 0 - survive at current time step
                
                for i in range(self.edgenum):
                    node1, node2 = self.edgelist[i, 0], self.edgelist[i, 1] 
                    if(np.abs(flowmatrix[node1, node2]) > (1 + alpha)*self.flowcapmatrix[node1, node2]):
                        link_seq_track[i] = 1
            
                #node failure caused by flow overload 
                for i in range(self.nodenum):
                    # print(i, np.sum(flowmatrix[:, i]*self.convratematrix[:, i]), np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))
                    if((np.abs(np.sum(flowmatrix[:, i]*self.convratematrix[:, i]))) > (1 + alpha)*np.abs(np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))):
                        node_seq_track[i] = 1
                
                self.node_fail_evol_track.append(node_seq_track)
                self.link_fail_evol_track.append(link_seq_track)
                
            self.satisfy_node_evol.append(satisfynode)
            
            node_seq = np.zeros(self.nodenum, dtype = int) #1 - failure, 0 - survive at current time step
            link_seq = np.zeros(self.edgenum, dtype = int) #1 - failure, 0 - survive at current time step
            #adjacent matrix update
            #link failure caused by flow overload
            for i in range(self.edgenum):
                node1, node2 = self.edgelist[i, 0], self.edgelist[i, 1] 
                if(np.abs(flowmatrix[node1, node2]) > (1 + alpha)*self.flowcapmatrix[node1, node2]):
                    link_seq[i] = 1
            
            #node failure caused by flow overload 
            for i in range(self.nodenum):
                if((np.abs(np.sum(flowmatrix[:, i]*self.convratematrix[:, i]))) > (1 + alpha)*np.abs(np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))):
                    node_seq[i] = 1
            
            self.node_fail_evol.append(node_seq)
            self.link_fail_evol.append(link_seq)
            
            #Update the adjmatrix caused by failed nodes
            adjmatrix = self.update_matrix_node_fail(copy.copy(adjmatrix), node_seq)
            
            #Update the adjmatrix caused by failed links
            adjmatrix = self.update_matrix_link_fail(copy.copy(adjmatrix), link_seq)
            
                    
            #Update the flow matrix
            flowmatrix = adjmatrix*flowmatrix
            
            self.adjmatrix_evol.append(adjmatrix)
            self.flowmatrix_evol.append(flowmatrix)
            
            
            #Check the stability: no newly failed nodes and links
            if(np.sum(link_seq) == 0 and np.sum(node_seq) == 0):
                break

    
    
power = network(dt.p_nodedata, dt.p_edgedata)
gas = network(dt.g_nodedata, dt.g_edgedata)
s = system(power, gas, dt.g2p_edgedata)
s.initial_failure('randomness', 0.1) #the type of the initial failure sequence, choice: 'randomness', 'dc' - degree centrality, 'bc' - betweenness centrality, 'kc' - katz centrality, 'cc': closeness centrality
s.cascading_failure(0.3)



