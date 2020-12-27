# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 22:47:00 2020

@author: 10624
"""

import numpy as np
import math
import random

import pandas as pd
from matplotlib import pyplot as plt
from sharefunction import set_default_plot_param
set_default_plot_param(plt)

import copy
import networkx as nx
import pygraphviz


import os
dir_path = os.path.dirname(os.path.realpath(__file__))    # go to the directory of the current script

import data as dt
from infrasnetwork import network
#import failsimulation as fs


class System(object):
    """ Couple the power and gas networks
    """
    def __init__(self, power, gas, g2p_arcdata): ##Ignore the arc from the power network to gas network
        self.power = power
        self.gas = gas
        self.nodenum = self.power.nodenum + self.gas.nodenum
        self.arcnum = self.power.arcnum + self.gas.arcnum + len(g2p_arcdata)
            
        self.g2p_arcdata = g2p_arcdata
        
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
        
        #add the power and gas arcs
        self.adjmatrix[0:self.gas.nodenum, 0:self.gas.nodenum] = copy.copy(self.gas.adjmatrix)
        self.flowmatrix[0:self.gas.nodenum, 0:self.gas.nodenum] = copy.copy(self.gas.flowmatrix)
        self.adjmatrix[self.gas.nodenum:(self.gas.nodenum + self.power.nodenum), self.gas.nodenum:(self.gas.nodenum + self.power.nodenum)] = copy.copy(self.power.adjmatrix)
        self.flowmatrix[self.gas.nodenum:(self.gas.nodenum + self.power.nodenum), self.gas.nodenum:(self.gas.nodenum + self.power.nodenum)] = copy.copy(self.power.flowmatrix)
        
        #add the gas2power arcs
        for i in range(len(self.g2p_arcdata)):
            self.adjmatrix[self.gas.ID2num[self.g2p_arcdata['start_node'].iloc[i]], self.power.ID2num[self.g2p_arcdata['end_node'].iloc[i]] + self.gas.nodenum] = 1
            self.flowmatrix[self.gas.ID2num[self.g2p_arcdata['start_node'].iloc[i]], self.power.ID2num[self.g2p_arcdata['end_node'].iloc[i]] + self.gas.nodenum] = self.g2p_arcdata['flow'].iloc[i]
    
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
        if (nx.algorithms.dag.is_directed_acyclic_graph(self.graph)):
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


#    def get_arc_id_list(self):
#        start_node_id_list = self.gas.start_node_id.to_list() + self.power.start_node_id.to_list()
#        end_node_id_list = self.gas.end_node_id.to_list() + self.power.end_node_id.to_list()
#        self.arc_id_list = [(start_node_id_list[i], end_node_id_list[i]) for i in range(len(start_node_id_list))]
        
    def propertymapping(self):
        """ Mapping all properties in the network to the system
        """
        
        self.demand = np.concatenate((np.array(self.gas.demand), np.array(self.power.demand)))
        self.supply_cap = np.concatenate((np.array(self.gas.supply_cap), np.array(self.power.supply_cap)))
        self.supply = np.concatenate((np.array(self.gas.supply), np.array(self.power.supply)))
        
        # array of list, each element is a list that will contain the start node number and end node number of an arc
        self.arclist = np.empty((self.power.arcnum + self.gas.arcnum + len(self.g2p_arcdata), 2), dtype = int)
        
        self.conv_rate = np.concatenate((np.array(self.gas.conv_rate), np.array(self.power.conv_rate), np.array(self.g2p_arcdata.conv_rate)))
        self.flow_cap = np.concatenate((np.array(self.gas.flow_cap), np.array(self.power.flow_cap), np.array(self.g2p_arcdata.flow_cap)))  
        
        # get start node id and end node id
        self.node_id_list = self.gas.nodeid.to_list() + self.power.nodeid.to_list()
        start_node_id_list = self.gas.start_node_id.to_list() + self.power.start_node_id.to_list() + self.g2p_arcdata.start_node.to_list()
        end_node_id_list = self.gas.end_node_id.to_list() + self.power.end_node_id.to_list() + self.g2p_arcdata.end_node.to_list()
        self.arc_id_list = [(start_node_id_list[i], end_node_id_list[i]) for i in range(len(start_node_id_list))]
#        self.arc_id_list = self.get_arc_id_list
#        
#        print(self.arc_id_list)
        
        self.flowcapmatrix = np.zeros((self.nodenum, self.nodenum), dtype = float)
        self.convratematrix = np.zeros((self.nodenum, self.nodenum), dtype = float)
        temp = 0
        for i in range(self.gas.arcnum):
            self.flowcapmatrix[self.ID2num[self.gas.start_node_id.iloc[i]], self.ID2num[self.gas.end_node_id.iloc[i]]] = self.gas.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.gas.start_node_id.iloc[i]], self.ID2num[self.gas.end_node_id.iloc[i]]] = self.gas.conv_rate.iloc[i]
            self.arclist[temp, 0], self.arclist[temp, 1] = self.ID2num[self.gas.start_node_id.iloc[i]], self.ID2num[self.gas.end_node_id.iloc[i]]
            temp += 1
            
        for i in range(self.power.arcnum):
            self.flowcapmatrix[self.ID2num[self.power.start_node_id.iloc[i]], self.ID2num[self.power.end_node_id.iloc[i]]] = self.power.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.power.start_node_id.iloc[i]], self.ID2num[self.power.end_node_id.iloc[i]]] = self.power.conv_rate.iloc[i]
            self.arclist[temp, 0], self.arclist[temp, 1] = self.ID2num[self.power.start_node_id.iloc[i]], self.ID2num[self.power.end_node_id.iloc[i]]
            temp += 1
            
        for i in range(len(self.g2p_arcdata)):
            self.flowcapmatrix[self.ID2num[self.g2p_arcdata.start_node.iloc[i]], self.ID2num[self.g2p_arcdata.end_node.iloc[i]]] = self.g2p_arcdata.flow_cap.iloc[i]
            self.convratematrix[self.ID2num[self.g2p_arcdata.start_node.iloc[i]], self.ID2num[self.g2p_arcdata.end_node.iloc[i]]] = self.g2p_arcdata.conv_rate.iloc[i]
            self.arclist[temp, 0], self.arclist[temp, 1] = self.ID2num[self.g2p_arcdata.start_node.iloc[i]], self.ID2num[self.g2p_arcdata.end_node.iloc[i]]
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
        
        if (Type == 'randomness'):
            self.initial_fail_seq = random.sample(range(self.nodenum), fail_num)
        else:  # case fail_num > self.nodenum needs not be handled because [-fail_num:]= ['all elements']
            # sort in descending order (negate the original) and return the indice of the top fail_num nodes
            exec('self.initial_fail_seq = np.argsort(-np.array(self.{}))[:fail_num]'.format(Type))  # sort in an ascending approach
    
        if (self.initial_fail_seq != []):
            self.initial_fail_seq_onehotcode[np.array(self.initial_fail_seq)] = 1
        self.initial_fail_link_onehotcode = np.zeros(self.arcnum, dtype = int) #There are no initial failed links
        
        # print(self.initial_fail_seq)
        
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
            if (link_fail_seq[i] == 1):
                adjmatrix[self.arclist[i, 0], self.arclist[i, 1]] = 0
        
        return adjmatrix
    
    def cascading_failure(self, redun_rate=0.2):
        """ Simulate the cascading failure
        Input:
            redun_rate: redundancy
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
                #calculate the total flow going into the node
                flowin = np.sum(flowmatrix[:, node]*self.convratematrix[:, node]) + self.supply[node]*(1 - self.node_fail_evol_track[-1][node])
                #Some flows serve for the node demand value, the remaining part goes to the following distribution process
                if np.round(flowin, 3) >= np.round(self.demand[node], 3):
                    satisfynode[node] = self.demand[node]
                    flowin -= self.demand[node]
                else:
                    satisfynode[node] = flowin/2 #if not enough to supply for the demand of the node, supply a half
                    flowin = flowin/2
#                print('satisfied power demand', satisfynode[-self.power.nodenum:]/self.demand[-self.power.nodenum:])
                    
                if self.initial_fail_seq != []:
                    performance_gas = np.sum(satisfynode[:self.gas.nodenum])/np.sum(self.demand[:self.gas.nodenum])
                    performance_power = np.sum(satisfynode[-self.power.nodenum:])/np.sum(self.demand[-self.power.nodenum:])
                else:
                    performance_gas =1
                    performance_power = 1
                performance_temp = np.mean([performance_gas, performance_power])
                self.performance.append(performance_temp) #Track down the performance
                
                #Redistribute the flow, here we can introduce some randomness to account for the uncertainty
                if (np.sum(adjmatrix[node, :]) == 0 or np.sum(self.flowmatrix_evol[-2][node, :]) == 0):
                    flowout = np.zeros(self.nodenum, dtype = float)
                else:
                    if (np.sum(adjmatrix[node, :]) != np.sum(self.adjmatrix_evol[-2][node, :])): # Some links fail, redistribute evenly with some random noise
                        flowout = 1/np.sum(adjmatrix[node, :])*flowin*adjmatrix[node, :]
                        
                        if (flowin != 0): #The uncertainty only happens when there are multiple out-links and inflow !=0
                            index = np.random.choice(np.argwhere(flowout != 0).reshape(-1)) #flow - beta where beta~U(0, flow/2)
                            if (len(np.argwhere(flowout != 0).reshape(-1)) != 1):
                                noise_flow = np.random.rand()*flowout[index]/2
                                unit_flow = noise_flow/(len(np.argwhere(flowout != 0)) - 1)
    
                                flowout[flowout != 0] = flowout[flowout != 0] + unit_flow
                                flowout[index] = flowout[index] - unit_flow - noise_flow
                                                
                    else: #No links fail,redistribute according to the ratio of flow at last time step
                        flowout = self.flowmatrix_evol[-2][node, :]/np.sum(self.flowmatrix_evol[-2][node, :])*flowin
                
                flowmatrix[node, :] = flowout
                
            node_seq_track = np.zeros(self.nodenum, dtype = int) #1 - failure, 0 - survive at current time step
            link_seq_track = np.zeros(self.arcnum, dtype = int) #1 - failure, 0 - survive at current time step
            
            for i in range(self.arcnum):
                node1, node2 = self.arclist[i, 0], self.arclist[i, 1] 
                if (np.abs(flowmatrix[node1, node2]) > (1 + redun_rate)*self.flowcapmatrix[node1, node2]):
                    # print(node1, node2, flowmatrix[node1, node2], self.flowcapmatrix[node1, node2])
                    link_seq_track[i] = 1
        
            #node failure caused by flow overload 
            for i in range(self.nodenum):
                # print(i, np.sum(flowmatrix[:, i]*self.convratematrix[:, i]), np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))
                if ((np.abs(np.sum(flowmatrix[:, i]*self.convratematrix[:, i]))) > \
                    (1 + redun_rate)*np.abs(np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))):
                    # print(time, 'node', i, np.sum(flowmatrix[:, i]*self.convratematrix[:, i]), np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))
                    node_seq_track[i] = 1
            
            self.node_fail_evol_track.append(node_seq_track)
            self.link_fail_evol_track.append(link_seq_track)
            time += 1
                
            self.satisfy_node_evol.append(satisfynode)
            
            node_seq = np.zeros(self.nodenum, dtype = int) #1 - failure, 0 - survive at current time step
            link_seq = np.zeros(self.arcnum, dtype = int) #1 - failure, 0 - survive at current time step
            #adjacent matrix update
            #link failure caused by flow overload
            for i in range(self.arcnum):
                node1, node2 = self.arclist[i, 0], self.arclist[i, 1] 
                if (np.abs(flowmatrix[node1, node2]) > (1 + redun_rate)*self.flowcapmatrix[node1, node2]):
                    link_seq[i] = 1
            
            #node failure caused by flow overload 
            for i in range(self.nodenum):
                if ((np.abs(np.sum(flowmatrix[:, i]*self.convratematrix[:, i]))) > \
                    (1 + redun_rate)*np.abs(np.sum(self.flowmatrix_evol[0][:, i]*self.convratematrix[:, i]))):
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
            if (np.sum(link_seq) == 0 and np.sum(node_seq) == 0):
                break
    

    # added functions start from here                 
    # to do: plot the networks
    def plot_inter_networks(self, link_df, is_save=True):
        
        G = pygraphviz.AGraph(strict=False, directed=True)
        
        G.add_nodes_from(self.node_id_list)
    
        for i in np.arange(len(link_df)):    
            G.add_arc(link_df[i,0], link_df[i,1], label=link_df[i,3])
    
        G.layout()
        if is_save:
            G.draw('inter_networks.pdf')
        else:
            G.draw()

    #%%
    def compare_attack_types(self, attack_types = ['randomness', 'dc', 'bc', 'kc', 'cc'],
                             attack_portions=np.round(np.arange(0,1.001,0.05),2),
                             redun_rate = 0.5, n_repeat_random=50):
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
                # repeat multiple times to obtain the average performance
                for k in np.arange(n_repeat_random):
                    # print(k, "-th repetition for random attack")
                    self.initial_failure(attack_types[j], attack_portions[i])
                    self.cascading_failure(redun_rate)
                    performance_final_temp = self.performance[-1]
                    performance[i,j,k] = performance_final_temp
                
                performance_mean[i,j] = np.mean(performance[i,j,:])                                            
        
        # convert array into pandas df
        performance_mean_df = pd.DataFrame(data=performance_mean, index=attack_portions.tolist(), columns=attack_types)
        print('Cascading failure ends.')
        
        return performance_mean_df
    
   
    def plot_performance_different_attack(self, attack_types = ['randomness', 'dc', 'bc', 'kc', 'cc'],
                                          attack_portions=np.round(np.arange(0,1.001,0.05),2),
                                          redun_rate = 0.2, n_repeat_random=50, is_save=True):
    # plot the performance of different attack types under different attack portions
        # attach both nodes and links? But links do not have a degree.
        
        performance_df = self.compare_attack_types(attack_types=attack_types, attack_portions=attack_portions,
                                                   redun_rate=redun_rate, n_repeat_random=n_repeat_random)
        
        plt.figure(figsize=(4, 3))
        
        styles = ['k-','b--','r-.X','c--o','m--s'][:performance_df.shape[1]]
        performance_df.plot(style=styles)
        
        plt.xlabel('Percentage of attacked nodes')
        plt.ylabel('Performance of networks')
        
        plt.ylim(bottom=-0.01)
        
        plt.grid(axis='both')
                
        legend_labels = ('Random', 'Degree-based', 'Betweenness-based', 'Closeness-based', 'Katz-based')[:performance_df.shape[1]]
        plt.legend(legend_labels)
        
        if is_save:
            plt.savefig('compare_attack_types.pdf')
        plt.show()


    def get_comp_damage_state_temp(self, type='node'):
        
        # get index of failed components
        if type=='node':
            comp_fail_evol = copy.copy(self.node_fail_evol)
        else:
            comp_fail_evol = copy.copy(self.link_fail_evol)
        fail_comp_idx = []
        for i in np.arange(len(comp_fail_evol)-1):
            fail_comp_idx += np.where(comp_fail_evol[i]==1)[0].tolist()
        
        
        # set the respective y_init to 0 if the component fails at either stage
        y_comp_init_temp = [1]*len(comp_fail_evol[0])
        y_comp_init = [0 if i in fail_comp_idx else item for i,item in enumerate(y_comp_init_temp)]
        
        return y_comp_init

    def get_damage_state(self, attack_types='randomness',
                         attack_portions=0.2,
                         redun_rate=0.2):
        '''
        get final damage state of components from the results of cascading failure at each time step
        feed into the optimization problem as y_init, the initial damage state of components before the restoration is initiated
        '''
        
        # simulate initial attack
        self.initial_failure(attack_types, attack_portions)
        
        # simulate casacading failure
        self.cascading_failure(redun_rate)
        
        # index of components that fail during failure propagation
            # t=0, the components that fail due to direct
        y_node_init, y_link_init = self.get_comp_damage_state_temp(type='node'), self.get_comp_damage_state_temp(type='link')
            
        return y_node_init, y_link_init

    # optimize the repair schedule  
    def optimize_restore(self, attack_types = 'randomness',
                         attack_portions=0.2,
                         redun_rate=0.2,
                         model_type='flow'):
    
        '''optimize the restoration of damaged components
        input:
            y_node_init and y_arc_init - list: =1 if not damaged, and 0 otherwise
            arcs damaged - dictionary: xxx
            
        output:
            resil_over_time
            components to restore at each time step, i.e. the optimal schedule
    
        notes:
            devise a minimal working example/as small as possible
            build up damage scenario to test the program
                scenario 1: no components are damaged
                scenario 2: one node or link is damaged
                scenario 3: all nodes are damaged
                scenario 4: all components are damaged
            make sure the codes mathes the model
            bugs: 
                wrong index or indentation
                differentiate between idx and the real value, e.g. i != node[i]
        refs.: 
            https://www.python-mip.com/
            https://pysal.org/spaghetti/notebooks/transportation-problem.html
        '''
        # x[i,t] - binary: whether or not to restore a node at time t, 0 otherwise.
        # y[i,t] - binary: whether or not to an link functions at time t, 0 otherwise.
        
        # import packages/functions
        from mip import Model, minimize, xsum, BINARY, OptimizationStatus
        
        # 0 get initial damage state
        y_node_init, y_arc_init = self.get_damage_state(attack_types=attack_types, attack_portions=attack_portions, redun_rate=redun_rate)
        
        # 1 decalre and initiate model
        model = Model()
        
        # 2 add decision variable
        # 2.1 schedule of repairing components
        num_node = copy.copy(self.nodenum)
        num_arc = copy.copy(self.arcnum)
        # total restoration time
        num_restore_max = 2
        num_damage_comp = y_arc_init.count(0) + y_node_init.count(0)
        from math import ceil
        time_list = list(range(ceil(num_damage_comp/num_restore_max) + 1))
        x_node = [[model.add_var(name="x({},{})".format(i, t), var_type=BINARY) for t in time_list] for i in np.arange(num_node)]
        x_arc = [[model.add_var(name="x({},{})".format(k, t), var_type=BINARY) for t in time_list] for k in np.arange(num_arc)]
        
        # 2.2 functional state of nodes and arcs
        y_node = [[model.add_var(name="y({},{})".format(i, t), var_type=BINARY) for t in time_list] for i in np.arange(num_node)]
        y_arc = [[model.add_var(name="y({},{})".format(k, t), var_type=BINARY) for t in time_list] for k in np.arange(num_arc)]
        
        # 2.3 flow, supply, demand, and slack
        flow = [[model.add_var(name="flow({},{})".format(k, t), lb=0) for t in time_list] for k in np.arange(num_arc)]   
        supply = [[model.add_var(name="supply({},{})".format(i, t), lb=0) for t in time_list] for i in np.arange(num_node)] 
        slack = [[model.add_var(name="slack({},{})".format(i, t), lb=0) for t in time_list] for i in np.arange(num_node)]
        
        # seudo variable: resilience at each time step
        resil = [model.add_var(name="resilience({})".format(t), lb=0) for t in time_list]

        # 3 obejctive function: min -1* sum of resilience at each time step
        model.objective = minimize(xsum(-1*resil[t] for t in time_list))
        
        
        # 4 add constraints
        # 4.1 component will be restored at one of the time periods
            # These two sets of constraints might not be necessary
        for i in np.arange(num_node):
            if y_node_init[i]==0:
                model.add_constr(xsum(x_node[i][t] for t in time_list) == 1)
        for k in np.arange(num_arc):
            if y_arc_init[k]==0:
                model.add_constr(xsum(x_arc[k][t] for t in time_list) == 1)
    
        # 4.2 number of components restored at a time period is capped
        for t in time_list:
            model.add_constr(xsum(x_node[i][t] for i in np.arange(num_node)) +
                             xsum(x_arc[k][t] for k in np.arange(num_arc)) <= num_restore_max)
    
    
        # 4.3 flow conservation
        # outflow - inflow = supply + slack - demand
        # node in each network: exclude inflow from interdependent links onto that node
        for i in np.arange(num_node):
            for t in time_list:
                model.add_constr(xsum(flow[k][t] for k in np.arange(num_arc) if self.node_id_list[i]==self.arc_id_list[k][0]) -
                                 xsum(flow[k][t] for k in np.arange(num_arc) if \
                                      (self.node_id_list[i]==self.arc_id_list[k][1] and self.conv_rate[k]==1))\
                                 == supply[i][t] + slack[i][t] - self.demand[i])
        
        # end node of interdependent links: supply <= converted supply/inflow to that node
        for i in np.arange(num_node):
            for k in np.arange(num_arc):
                if self.arc_id_list[k][1]==self.node_id_list[i] and self.conv_rate[k]!=1:
                    for t in time_list:
                        model.add_constr(supply[i][t] <= self.conv_rate[k]*flow[k][t])
    
        # 4.4 ub of flow, supply, and slack
        # 4.4.1.1 add auxillary variables to linearize the product of binary variables
        aux_z_arc = [[model.add_var(name="aux_z_arc({},{})".format(k, t), var_type=BINARY) for t in time_list] for k in np.arange(num_arc)] 
        for t in time_list:
            for k in np.arange(num_arc):
                # use auxillary variables to linearize the product of binary variables
                start_node_idx = self.node_id_list.index(self.arc_id_list[k][0])
                end_node_idx = self.node_id_list.index(self.arc_id_list[k][1])
                n_binary_var = 3
                # aux_z <= each binary variable
                model.add_constr(aux_z_arc[k][t] <= y_node[start_node_idx][t])
                model.add_constr(aux_z_arc[k][t] <= y_node[end_node_idx][t])
                model.add_constr(aux_z_arc[k][t] <= y_arc[k][t])
                # aux_z >= sum of binary variables - number of binary variables + 1
                    # active when all binary varialbes == 1
                model.add_constr(aux_z_arc[k][t] >= y_node[start_node_idx][t] + y_node[end_node_idx][t] + 
                                                    y_arc[k][t] - (n_binary_var-1))
                
                # 4.4.1.2 flow cap
                # flow will be zeros unless the start, end node nad the arc itsel are all functional
                if model_type=='flow':
                    model.add_constr(flow[k][t] <= aux_z_arc[k][t]*self.flow_cap[k])
                else:
                    model.add_constr(flow[k][t] <= aux_z_arc[k][t]*1e5)
    
                
            # 4.4.2 slack and supply cap
            for i in np.arange(num_node):
                    model.add_constr(slack[i][t] <= self.demand[i])
                    if model_type=='flow':
                        model.add_constr(supply[i][t] <= y_node[i][t]*self.supply_cap[i])
                    else:
                        model.add_constr(supply[i][t] <= y_node[i][t]*1e5)
                        
                
        # 4.5.0
        for t in time_list:  # start from the second time period
            # 4.5.1 non-deteriorating state of components
            if t==0:
                for i in np.arange(num_node):
                    model.add_constr(y_node_init[i] <= y_node[i][t])
                    model.add_constr(y_node_init[i] >= y_node[i][t])
                for k in np.arange(num_arc):
                    model.add_constr(y_arc_init[k] <= y_arc[k][t])
                    model.add_constr(y_arc_init[k] >= y_arc[k][t]) 
            else:            
                for i in np.arange(num_node):
                    model.add_constr(y_node[i][t-1] <= y_node[i][t])
                for k in np.arange(num_arc):
                    model.add_constr(y_arc[k][t-1] <= y_arc[k][t])
            # 4.5.2 components will be functional once repaired
                for i in np.arange(num_node):
                    model.add_constr(y_node[i][t] <= y_node[i][t-1] + x_node[i][t-1])
                for k in np.arange(num_arc):
                    model.add_constr(y_arc[k][t] <= y_arc[k][t-1] + x_arc[k][t-1])


        # 4.6 calculate flow-based resilience or topology-based resilience
        node_demand_idx = [idx for idx, demand_val in enumerate(self.demand) if demand_val > 0] 
        if model_type == 'flow':
            for t in time_list:
                # proportion of satisfied demand
                demand_satify_rate = np.mean(1-slack[np.where(self.demand!=0), t]/self.demand[np.where(self.demand!=0)])
                model.add_constr(resil[t] >= demand_satify_rate)
                model.add_constr(resil[t] <= demand_satify_rate)
            else:
                # proportion of demand nodes whose slack is lower than demand, i.e. demand node can receive some supply, irrespective of the amount.
                node_demand_receive_supply = [demand_val for idx, demand_val in enumerate(self.demand) if demand_val > slack[idx,t]]
                demand_connect_to_supply_rate = len(node_demand_receive_supply)/len(node_demand_idx)
                model.add_constr(resil[t] >= demand_connect_to_supply_rate)
                model.add_constr(resil[t] <= demand_connect_to_supply_rate)            
            
        # 5.0 solve the model and check status
        model.max_gap = 1e-5
        status = model.optimize(max_seconds=60*5)
        
        # 6.0 query optimization results    
        # 6.1 check solution status
        if status == OptimizationStatus.OPTIMAL:
            # print('optimal solution: {}'.format(model.objective_value))
        
            # 6.2 get objective value and x
            obj_value = model.objective_value
                        
            # # 6.3 print the number of variables and constraints
            # print('model has {} vars, {} constraints and {} nzs'.format(model.num_cols,\
                        # model.num_rows, model.num_nz))
                   
            return obj_value, x_node, x_arc, y_node, y_arc, resil, time_list, y_node_init, y_arc_init
        
        else:
            print('Infeasible or unbounded problem')
    
    # 2 get solution
    # 2.1 convert solution results
    def convert_solu_list_to_arr(self, var):
        '''convert list of solutions to variables in mip entity format to array
        input:
            var - mip solution list: e.g. var=y_var, solution to the functional state of arcs
        '''
        var_arr = np.zeros([len(var), len(var[0])])
        for j in np.arange(len(var)):
            for t in np.arange(len(var[0])):
                var_arr[j,t] = var[j][t].x
        
        return var_arr
    
    # 2.2 solve the model and extract solutions
    def get_solution(self, attack_types = 'randomness',
                     attack_portions=0.2,
                     redun_rate = 0.2,
                     model_type='flow'):
        
        # 1 solve the model
        obj_value, x_node, x_arc, y_node, y_arc, resil, time_list, y_node_init, y_arc_init = \
            self.optimize_restore(attack_types=attack_types, attack_portions=attack_portions, redun_rate=redun_rate, model_type=model_type)
               
        # 2 extract results
        # extract x and y
        x_node_arr = self.convert_solu_list_to_arr(x_node)
        x_arc_arr = self.convert_solu_list_to_arr(x_arc)
        
        #y_arc_arr = convert_solu_list_to_arr(y_arc, time_list)
        #y_node_arr = convert_solu_list_to_arr(y_node, time_list)
        
        # extract resilience  
        resil_arr = self.convert_solu_list_to_arr(resil)  
        #supply_arr = convert_solu_list_to_arr(supply, time_list) 
        
        return x_node_arr, x_arc_arr, resil_arr, time_list, y_node_init, y_arc_init

        
    # 3 visualize results
    # 3.1.1 prepare schedule data
    def get_schedule_df(self, attack_types='randomness',
                        attack_portions=0.2,
                        redun_rate=0.2,
                        model_type='flow'):
        # store scheduling resulst in a df
            # df: index: damaged component; columns: start_time, duration, finish time
        # get damaged component id

        x_node_arr, x_arc_arr, resil_arr, time_list, y_node_init, y_arc_init = \
            self.get_solution(attack_types=attack_types, attack_portions=attack_portions, redun_rate=redun_rate)

        comp_list = self.node_id_list + self.arc_id_list
        # due to randomness in cascading failure, get damage state can only be called once within the function for solving the problem
        # y_node_init, y_arc_init = self.get_damage_state(attack_types=attack_types,attack_portions=attack_portions, redun_rate=redun_rate)
        init_state_list = y_node_init + y_arc_init
        damaged_comp_list =  [comp_list[i] for i, item in enumerate(init_state_list) if item==0]
        
        # get restoration start time point
        x_comp = np.concatenate((x_node_arr, x_arc_arr), axis=0)
        # select restoration start time of damaged components
        x_comp_damage = x_comp[np.amax(x_comp, axis=1)==1]
        # restore start time of each component
        restore_start_time = np.argmax(x_comp_damage, axis=1)
        
        # create df and sort by restore start time
        schedule_df = pd.DataFrame({'restore_start_time':restore_start_time}, index=damaged_comp_list)
        schedule_df = schedule_df.sort_values(by='restore_start_time')
        schedule_df['duration'] = 1
        schedule_df['restore_end_time'] = schedule_df['restore_start_time'] + schedule_df['duration']
        
        return schedule_df
    
    # 3.1.2 plot schedule
    def plot_repair_schedule(self, attack_types = 'randomness',
                             attack_portions=0.2,
                             redun_rate=0.2, is_save=True,
                             model_type='flow'):
        # plot the restoreation schedule
            # refs.: https://towardsdatascience.com/from-the-bridge-to-tasks-planning-build-gannt-chart-in-python-r-and-tableau-7256fb7615f8
                # https://plotly.com/python/gantt/
        
        # get schedule df 
        schedule_df = self.get_schedule_df(attack_types=attack_types, attack_portions=attack_portions, redun_rate=redun_rate)
        
        # plot parameters
        max_time = schedule_df['restore_end_time'].max()
        bar_ht = 0.75
        off_ht = 0.5
        
        fig, ax = plt.subplots(figsize=(1+max_time/1.25, schedule_df.shape[0]/3))         
        for i in np.arange(schedule_df.shape[0]):
            ax.broken_barh([(schedule_df['restore_start_time'].iloc[i], schedule_df['duration'].iloc[i])],
                            yrange = (i+off_ht/4, bar_ht),
                            facecolors = ('tab:blue') if isinstance(schedule_df.index[i], tuple) else ('tab:red'),
                            edgecolor = "none") 
                    
        ax.set_title('Repair schedule after {} attack'.format(attack_types))
        
        ax.set_ylabel('Component')
        ax.set_xlabel('Time period')
        
        ax.set_yticks([i + off_ht for i in np.arange(schedule_df.shape[0])]) 
        y_tick_labels = ["%s->%s" % item if isinstance(item,tuple) else item for item in schedule_df.index] 
        ax.set_yticklabels(y_tick_labels)
        ax.set_ylim(bottom=-off_ht/2, top=schedule_df.shape[0]+off_ht/2)
        
        ax.set_xticks(np.arange(0, max_time+1, 1.0))
        ax.set_xticklabels(np.arange(1, max_time+2, 1))
        ax.set_xlim(left=-off_ht/2, right=max_time+off_ht/2)
      
        ax.grid(True)
        
        # draw legend manually
        colors = {'Node':'tab:red', 'Link':'tab:blue'}         
        labels = list(colors.keys())
        handles = [plt.Rectangle((0,0),0.1,0.1, color=colors[label]) for label in labels]
        plt.legend(handles, labels, loc='lower right')
        
        if is_save:
            plt.savefig('repair_schedule_{}_{}.pdf'.format(attack_types, model_type))
            
        plt.show()
    

    def get_resil_df(self, attack_types = ['randomness'],
                     attack_portions=0.2,
                     redun_rate=0.2,
                     n_repeat_random=50,
                     model_type='flow'):
        '''get df of resilience over time under each attack types
        
        output:
            
        '''
        
        # get solultion to slack at node at time t
        if isinstance(attack_types, list):
            n_attack_types = len(attack_types)
        else:
            n_attack_types = 1
        time_max = 1
        resil_arr_3d = np.ones([self.nodenum+self.arcnum, n_attack_types, n_repeat_random])  # time list is sufficiently large
        for i in np.arange(len(attack_types)):
            for j in np.arange(n_repeat_random):

                x_node_arr, x_arc_arr, resil_arr, time_list, y_node_init, y_arc_init = self.get_solution(attack_types=attack_types[i],
                                                                     attack_portions=attack_portions,
                                                                     redun_rate=redun_rate, model_type=model_type)        
                # resilience at each time point
                resil_temp = resil_arr               
                time_max_temp = time_list[np.where(resil_temp==1)[-1][0]] + 1  # time to full restoration. [0]: get value of array
                resil_arr_3d[:time_max_temp, i, j] = resil_temp[:time_max_temp]                    
                # update number of time periods
                if time_max_temp > time_max:
                    time_max = time_max_temp
        # get mean over n_repeat_random
        resil_arr_mean = np.mean(resil_arr_3d, axis=2)
        # create df
        # add time periods as index
        idx = [item for item in range(1, time_max+1)]
        resil_df = pd.DataFrame(data=resil_arr_mean[:time_max,:], columns=attack_types, index=idx)      
        
        return resil_df 


    def plot_resil(self, attack_types = ['randomness', 'dc'],
                   attack_portions=0.2,
                   redun_rate=0.2,
                   n_repeat_random=50, is_save=True):
        
        resil_df = self.get_resil_df(attack_types=attack_types, attack_portions=attack_portions,
                                     redun_rate=redun_rate, n_repeat_random=n_repeat_random)
        
        plt.figure(figsize=(4, 3))
        
        styles = ['k-','b--','r-.X','c--o','m--s'][:resil_df.shape[1]]
        resil_df.plot(style=styles)
        
        plt.xlabel('Time period')
        plt.ylabel('Resilience')
        
        plt.ylim(top=1.02)
        
        plt.grid(axis='both')
                
        legend_labels = ('Random', 'Degree-based', 'Betweenness-based', 'Closeness-based', 'Katz-based')[:resil_df.shape[1]]
        plt.legend(legend_labels, loc='lower right')
        
        if is_save:
            plt.savefig('resilience.pdf')
        plt.show()

   

# 3.2 restoration rate over time  
        
def main():
    
    # change directory
    os.chdir('C:/Users/yuj5/Documents/GitHub/ICOSSAR2021')
    
    # create network 
    power = network(dt.p_nodedata, dt.p_arcdata) # instantiate the power network
    gas = network(dt.g_nodedata, dt.g_arcdata) # instantiate the gas network
    s = System(power, gas, dt.g2p_arcdata) # instantiate the system class
    
    ## compare performance under different types of attacks
    #attack_portions = np.round(np.arange(0,1.02,0.05), 2)
    attack_types = ['randomness', 'dc', 'bc', 'cc']
    
    # plot
    # performance under different types of attacks
    REDUN_RATE = 0.4
    N_REPEAT = 50
    s.plot_performance_different_attack(attack_types=attack_types,
                                        attack_portions=np.round(np.arange(0,1.01,0.1),2),
                                        redun_rate=REDUN_RATE, n_repeat_random=N_REPEAT, is_save=True)
    
    # restoration schedule after cascading failures under different types of attacks
    attack_types = 'cc'
    s.plot_repair_schedule(attack_types=attack_types, redun_rate=REDUN_RATE, is_save=True)
    
    s.plot_resil(attack_types=['randomness','dc', 'bc', 'cc'], redun_rate=REDUN_RATE,
                 n_repeat_random=N_REPEAT, is_save=True)
    

#main()  

      
    
    
## test optimization model
    # 1 import data
    # 1.1 prepare arc data
#    def get_arc_id_list(self, start_node, end_node):
#        '''get arc list from lists of start node and end node
#        
#        returns:
#            list of tuple elements, e.g. (start node, end node)
#        '''
#        arc_id_list = []
#        for i in np.arange(len(start_node)):
#            arc_id_list.append((start_node[i],end_node[i]))
#            
#        return arc_id_list
    
#    # 1.2 import data    
#    def import_data(self):
#    
#        # 1 import data    
#        node_data = pd.read_csv('./data/case_6_node/node_data.csv')
#        arc_data = pd.read_csv('./data/case_6_node/arc_data.csv')
#        
#        # 2 extract data
#        # 2.1 node
#        demand, supply_cap = node_data.demand.astype(float).tolist(),\
#                             node_data.supply_cap.tolist() # demand should be float format
#        node_id_list = node_data.node_id.tolist()
#        
#        # 2.2 arc
#        # get list of arcs
#        start_node = arc_data.start_node
#        end_node = arc_data.end_node
#        arc_id_list = get_arc_id_list(start_node, end_node)
#        
#        # flow
#        flow_cap = arc_data.flow_cap.astype(float).tolist()
#        conv_rate = arc_data.conv_rate.astype(float).tolist()
#        
#        # initial state of nodes and arcs
#        y_node_init, y_arc_init = node_data.y_node_init.astype(float).tolist(),\
#                                  arc_data.y_arc_init.astype(float).tolist()
#    
#        return node_id_list, arc_id_list, y_node_init, y_arc_init, demand, flow_cap, supply_cap, conv_rate