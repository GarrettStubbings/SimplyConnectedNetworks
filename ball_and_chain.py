""" 

This program created February 2021 to build and run
"Optimal" networks (hand built max longevity Idea)

"""
import os
import re 
import itertools as it
import pylab as pl
from general_functions import *
from network_from_mean_field import *
from visualization import *
import networkx as nx
from default_params import *
import datetime
from evolution_performance import *
date = datetime.date.today().strftime('%e%b%y')
date = date.strip(' ')
import matplotlib.style
import matplotlib as mpl
from scipy.optimize import curve_fit
from evo_plotting import *

mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['lines.markeredgecolor'] = 'k'
mpl.rcParams['lines.markeredgewidth'] = 0.5

if __name__ == '__main__':

    pl.close('all')
    running_stuff = 0
    check_performance = 0
    show_graph = 1
    matrix_folder = 'AdjacencyMatrices/'

    ##### Setting parameters for network
    N = 1000
    Nc = 1
    n_fi_nodes = 10
    avg_k = 4
    ## parameters of running the simulaton
    number = "100000"
    output_folder = 'OptimalNetworkData'

    seed  = '1'
    Lambda = "0.0"
    beta = "100.0"
    power = "1.0"
    numChanges = "0"
    run_hours = "0.0"
    evo_condition = 'none'

    N_values = [100, 200, 500, 1000]
    Nc_values = [1, 2, 5, 10, 20]

    if running_stuff:
        for i, N in enumerate(N_values):
            for j, Nc in enumerate(Nc_values):
                N = int(N)
                avg_k = float(avg_k)
                initial_distribution = "N{0}Nc{1}Optimal".format(N, Nc)
                G, mort_nodes, fi_nodes  = optimal_network(N, Nc, avg_k,
                                                            n_fi_nodes)
                adjacency_matrix = pl.asarray(nx.adjacency_matrix(G).todense())

                N = str(N)
                avg_k = str(avg_k)

                pl.savetxt('{0}{1}.csv'.format(matrix_folder,initial_distribution),
                    adjacency_matrix, fmt = '%i', delimiter = ',')

                ### running stuff
                # compile if needed
                sub.call(["make", "testNetwork"])

                # parameters to run
                param_names = ['number', 'N', 'numChanges', 'avgDeg', 'folder',
                    'singleseed', 'evoCondition', 'InitialDistribution',
                    'lambda', 'beta', 'power']
                params = [number, N, numChanges, avg_k, output_folder, seed,
                    run_hours, evo_condition, initial_distribution, Lambda,
                    beta, power]

                command = ['./testNetwork'] + params 
                sub.call(command)

    if check_performance:
        print('Checkiing performance')
        data_dir = 'OptimalNetworkData/'
        plot_dir = 'Plots/'
        data = optimization(data_dir, plot_dir)
        number = '100000'
        Nc_values = pl.array([1, 5, 10, 20])
        N_values = pl.array([100, 200, 500, 1000])
        data.update_params({'Number': '100000'})
        Nc_mean_death_ages = []
        
        pl.figure(figsize = (8,6))
        for i, Nc in enumerate(Nc_values):
            N_mean_death_ages = []
            N_errors = []
            for j, N in enumerate(N_values):
                N_str = str(N)
                initial_distribution = 'N{0}Nc{1}Optimal'.format(N, Nc)
                data.update_params({'N': N_str,
                                'InitialDistribution': initial_distribution})
                death_ages = data.get_values("DeathAges")
                N_mean_death_ages.append(pl.average(death_ages))
                N_errors.append(pl.std(death_ages)/pl.sqrt(int(number)))
            
            if Nc == 1:
                max_death_age_vs_N = pl.column_stack(
                    [N_values, N_mean_death_ages])
                pl.savetxt('MaxDeathAgeVsN.txt', max_death_age_vs_N)
            pl.errorbar(N_values, N_mean_death_ages, yerr = N_errors, 
                marker = 'o', ls = 'none', label = 'Nc: {}'.format(Nc),
                capsize = 3)
            pl.legend()
            pl.xscale('log')
            pl.xlabel('N')
            pl.ylabel('Mean Death Age')
            Nc_mean_death_ages.append(N_mean_death_ages)

                

    if show_graph:
        Nc = 10

        pl.figure(figsize = (8,6))
        G, mort_nodes, fi_nodes  = optimal_network(N, Nc, avg_k,
                                                    n_fi_nodes)
        pos = nx.spring_layout(G)
        nx.draw(G, pos, edge_color = 'C7', node_size = 3)
        pl.savefig('Plots/MarchMeetingPlots/OptimalNetworkN{0}Nc{1}.pdf'.format(
            N, Nc))
        
    pl.show()
