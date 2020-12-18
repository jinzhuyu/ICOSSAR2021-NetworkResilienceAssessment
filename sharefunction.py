# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 17:43:40 2020

@author: 10624
"""

# set default plot parameters
def set_default_plot_param():
    
    plt.style.use('classic')
    
    plt.rcParams["font.family"] = "Helvetica"
    plt.rcParams['font.weight']= 'normal'
    plt.rcParams['figure.figsize'] = [6, 6*3/4]
   
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    
    plt.rc('axes', titlesize=16, labelsize=15, linewidth=0.9)    # fontsize of the axes title, the x and y labels
    
    plt.rc('lines', linewidth=1.9, markersize=6)
    
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)

    plt.rcParams['axes.formatter.useoffset'] = False # turn off offset
    # To turn off scientific notation, use: ax.ticklabel_format(style='plain') or
    # plt.ticklabel_format(style='plain')

    
    plt.rcParams['legend.fontsize'] = 13
    plt.rcParams["legend.fancybox"] = True
    plt.rcParams["legend.loc"] = "best"
    plt.rcParams["legend.framealpha"] = 0.5
    
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['savefig.dpi'] = 800
    
#    plt.rc('text', usetex=False)
