# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 22:47:00 2020

@author: 10624
"""

import numpy as np
import math
import random

import pandas as pd
import matplotlib.pyplot as plt
from sharefunction import set_default_plot_param
set_default_plot_param()

import copy
import networkx as nx
import pygraphviz


import os
dir_path = os.path.dirname(os.path.realpath(__file__))    # go to the directory of the current script

import data as dt
from infrasnetwork import network
#import failsimulation as fs


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
        self.flowmatrix = np.zeros((self.nodenum, self.nodenum), dtype = float)
        
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
        else:  # case fail_num > self.nodenum needs not be handled because [-fail_num:]= ['all elements']
            # sort in descending order (negate the original) and return the indice of the top fail_num nodes
            exec('self.initial_fail_seq = np.argsort(-np.array(self.{}))[:fail_num]'.format(Type))  # sort in an ascending approach
    
        if(self.initial_fail_seq != []):
            self.initial_fail_seq_onehotcode[np.array(self.initial_fail_seq)] = 1
        self.initial_fail_link_onehotcode = np.zeros(self.edgenum, dtype = int) #There are no initial failed links
        
        return self.initial_fail_seq
        
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
        time = 0
        while(1): #perform the flow redistribution until stable
            satisfynode = copy.copy(self.satisfy_node_evol[-1])
            flowmatrix = copy.copy(self.flowmatrix_evol[-1])
            adjmatrix = copy.copy(self.adjmatrix_evol[-1])
        
            for node in self.topo_order:
                flowin = np.sum(flowmatrix[:, node]*self.convratematrix[:, node]) + self.supply[node]*(1 - self.node_fail_evol_track[-1][node]) #calculate the total flow going into the node
                #Some flows serve for the node demand value, the remaining part goes to the following distribution process
                if(np.round(flowin, 3) >= np.round(self.demand[node], 3)):
                    satisfynode[node] = self.demand[node]
                    flowin -= self.demand[node]
                else:
                    # print(2, node, flowin, self.demand[node])
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
                    # print(node1, node2, flowmatrix[node1, node2], self.flowcapmatrix[node1, node2])
                    link_seq_track[i] = 1
        
            #node failure caused by flow overload 
            for i in range(self.nodenum):
                # print(i, np.sum(flowmatrix[:, i]*self.convratematrix[:, i]), np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))
                if((np.abs(np.sum(flowmatrix[:, i]*self.convratematrix[:, i]))) > (1 + alpha)*np.abs(np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))):
                    # print(time, 'node', i, np.sum(flowmatrix[:, i]*self.convratematrix[:, i]), np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))
                    node_seq_track[i] = 1
            
            self.node_fail_evol_track.append(node_seq_track)
            self.link_fail_evol_track.append(link_seq_track)
            time += 1
                
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
            
            
    # plot graph
    def plot_inter_networks(self, node_list, link_df, is_save=0):
        
        G = pygraphviz.AGraph(strict=False, directed=True)
        
        G.add_nodes_from(node_list)
    
        for i in np.arange(len(link_df)):    
            G.add_edge(link_df[i,0], link_df[i,1], label=link_df[i,3])
    
        G.layout()
        if is_save==1:
            G.draw('inter_networks.pdf')
        else:
            G.draw()

#%%
def compare_attack_types(s=s, attack_types = ['randomness', 'dc', 'bc', 'kc', 'cc'],
                         attack_portions=np.round(np.arange(0,1.001,0.05),2),
                         redun_rate = 0.2, n_repeat_random=50):
    '''
    obtain network performance given different failure types and failure rate
    
    inputs:
        s - an instance of system representing the power-gas systems
        n_repeat_random - int: number of repetition for random attacks to obtain the average performance
    
    returns:
        performance_df - df: performance after all types of attacks over 
        performance_random_attack - df: performance after random attacks for n_repeat_random times
    '''
    performance = np.zeros([len(attack_portions), len(attack_types), n_repeat_random])
    performance_mean = np.zeros([len(attack_portions), len(attack_types)])
    
    for i in np.arange(len(attack_portions)):    
        for j in np.arange(len(attack_types)):           
            # repeat 50 times to obtain the average performance
            for k in np.arange(n_repeat_random):
                # print(k, "-th repetition for random attack")
                s.initial_failure(attack_types[j], attack_portions[i])
                s.cascading_failure(redun_rate)
                performance_final_temp = s.performance[-1]
                performance[i,j,k] = performance_final_temp
            performance_mean[i,j] = np.mean(performance[i,j,:])              
#                if attack_portions[i]>=0.05:
#                    print('Attack type: ', attack_types[j])
#                    print('Attack portion: ', attack_portions[i])
#                    print('Initial failure sequence: ', initial_fail_seq_temp)
#                    print('Number of nodes failed: ', len(initial_fail_seq_temp))                                 

    
    # convert array into pandas df
    performance_mean_df = pd.DataFrame(data=performance_mean, index=attack_portions.tolist(), columns=attack_types)
    print('Cascading failure finished.')
    print('Redundancy rate is: ', redun_rate)
    
    return performance_mean_df
    
   
def plot_performance_different_attack(df, is_save=0):
# plot the performance of different attack types under different attack portions
    # attach both nodes and links? But links do not have a degree.
    plt.figure(figsize=(4, 3))
    
    styles = ['k-','b--','r-.X','c--o','m--s']
    df.plot(style=styles)
    
    plt.xlabel('Percentage of attacked nodes')
    plt.ylabel('Performance of networks')
    
    plt.grid(axis='both')
    
    
    legend_labels = ('Random', 'Degree-based', 'Betweenness-based', 'Katz-based', 'Closeness-based')
    plt.legend(legend_labels)
    
    if is_save==1:
        plt.savefig('compare_attack_types.pdf')
    plt.show()

# change directory
os.chdir('C:/Users/yuj5/Documents/GitHub/ICOSSAR2021')

# create network 
power = network(dt.p_nodedata, dt.p_edgedata) #initialize the power network
gas = network(dt.g_nodedata, dt.g_edgedata) #initialize the gas network
s = system(power, gas, dt.g2p_edgedata)

# simulate failure
#s.initial_failure('randomness', 0.3) #the type of the initial failure sequence, choice: 'randomness', 'dc' - degree centrality, 'bc' - betweenness centrality, 'kc' - katz centrality, 'cc': closeness centrality
#s.cascading_failure(0.5)
attack_portions = np.round(np.arange(0,1.02,0.05), 2)
performance_mean_df = compare_attack_types(attack_portions=attack_portions, redun_rate=0.2, n_repeat_random=100)

# plot
plot_performance_different_attack(df=performance_mean_df, is_save=1)

#%%
#from itertools import product
from mip import Model, minimize, xsum, BINARY, OptimizationStatus

def optimize_restore(y_node_init, y_arc_init, demand, flow_cap, supply_cap):

    '''optimize the restoration of damaged components
    input:
        nodes_damaged - list:
        arcs damaged - dictionary: xxx
        
    output:
        resil_over_time
        components to restore at each time step, i.e. the optimal schedule
    '''
    # x[i,t] - binary: whether or not to restore a node at time t, 0 otherwise.
    # y[k,t] - binary: whether or not to restore a link at time t, 0 otherwise.
    
    # 1. initiate model
    model = Model()
    
    # 2.0 add decision variable
    # 2.1 schedule of repairing components
    time_list = sum(y_arc_init) + sum(y_node_init)
    x_node = [[model.add_var(name="x({},{})".format(i, t), var_type=BINARY) for t in time_list] for i in node_list]
    x_arc = [[model.add_var(name="x({},{})".format(k, t), var_type=BINARY) for t in time_list] for k in node_list]
    
    # 2.2 functional state of nodes and arcs
    y_node = [[model.add_var(name="y({},{})".format(i, t), var_type=BINARY) for t in time_list] for i in node_list]
    y_arc = [[model.add_var(name="y({},{})".format(k, t), var_type=BINARY) for t in time_list] for k in arc_list]
    
    # 2.3 flow, supply, demand, and slack
    flow = [[model.add_var(name="flow({},{})".format(k, t), lb=0) for t in time_list] for k in arc_list]   
    supply = [[model.add_var(name="supply({},{})".format(i, t), lb=0) for t in time_list] for i in node_list] 
    slack = [[model.add_var(name="slack({},{})".format(i, t), lb=0) for t in time_list] for i in node_list] 
     
    # 3.0 obejctive function
    model.objective = minimize(xsum((1- s[i][t]/demand[i]) for i in node_list for t in time_list))
    
    # 4.0 add constraints
    # 4.1 component is restored at one time step
        # These two sets of constraints might be unnecessary.
    for i in node_list:
        model.add_constr(xsum(x_node[i][t] for t in time_list) == 1)
    for k in arc_list:
        model.add_constr(xsum(x_arc[k][t] for t in time_list) == 1)

    # 4.2 restore 2 components at a time period
    for t in time_list:
        model.add_constr(xsum(x_node[i][t] for i in node_list)+ xsum(x_arc[k][t] for k in arc_list) == 2)

    # 4.3 flow conservation: outflow - inflow = supply + slack - demand
    for i in node_list:
        for t in time_list:
            model.add_constr(xsum(f[k][t] for k in arc_list if i==k[0]) - xsum(f[k][t] for k in arc_list if i==k[1]) ==
                             supply[i][t] + slack[i][t] - demand[i])
    
    # 4.4 ub of flow, supply, and slack
    for k in arc_list:
        model.add_constr(flow[k][t] <= y_arc[k][t]*flow_cap[k])
    for i in node_list:
        model.add_constr(slack[i][t] <= y_node[i][t]*demand[i])
        model.add_constr(supply[i][t] <= y_node[i][t]*supply_cap[i])
        
    # 4.5 start node of an arc should be restored before the arc
        # this might not be necessary
    for k in arc_list:
        for i in node_list:
            if i==k[0]:
                model.add_constr(x_arc[k][t] <= x_node[i][t])
    
    # 4.6 components will be functional once repaired
    for i in node_list:
        model.add_constr(y_node[i][t] <= y_node[i][t-1] + x_node[i][t-1])
    for k in arc_list:
        model.add_constr(y_arc[k][t] <= y_arc[k][t-1] + x_arc[i][t-1])
        
    # 4.7 non-deteriorating state of components
    for i in node_list:
        model.add_constr(y_node[i][t-1] <= y_node[i][t])
    for i in arc_list:
        model.add_constr(y_arc[k][t-1] <= y_arc[k][t])
        
    # 5.0 solve the model and check status
    model.max_gap = 0.01
    status = model.optimize(max_seconds=60*10)
    
    if status == OptimizationStatus.OPTIMAL:
        print('optimal solution cost {} found'.format(m.objective_value))
    
    # 6.0 query optimization results
    
    # number of variables and constraints
    print('model has {} vars, {} constraints and {} nzs'.format(model.num_cols, model.num_rows, model.num_nz))
    
#%%    

    
    
    
    
    
    
    