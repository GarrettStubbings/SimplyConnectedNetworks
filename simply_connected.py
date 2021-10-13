"""
Slash and grab first look at this approach for optimization

Build networks in this python script and push them to an old C++ GNM script to
get performance metrics etc.

"""

import os
import re 
import subprocess as sub
import time
import itertools as it
import pylab as pl
import string
"""
from general_functions import *
from network_from_mean_field import *
from default_params import *
from evolution_performance import *
from evo_plotting import *
from visualization import *
"""
from cluster_analysis import assortativity_extremes
import networkx as nx
import datetime
date = datetime.date.today().strftime('%e%b%y')
date = date.strip(' ')
import matplotlib.style
import matplotlib as mpl
import matplotlib.cm as cmx
from scipy.optimize import curve_fit
import numpy as np

mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['lines.linewidth'] = 1.5
#mpl.rcParams['lines.markeredgecolor'] = 'k'
mpl.rcParams['lines.markeredgewidth'] = 0.5

def gaussian_parameter_change(param, width, limits):
    """
    Change the parameter proportionally using a fraction generated from a 
    gaussian. Disallow values outside of the limits.
    """
    change = 1 - width*pl.normal()

    while change * param < limits[0] or change * param > limits[1]:
        change = 1 - width*pl.normal()
    
    return param * change

def optimal_node_partition(N, Nc, avg_k):
    """ Calculates how many of each node you need for 'optimal' network of size
    N when trying to protect Nc nodes"""

    det = pl.sqrt(9 - 4 * (2*N + 1 - Nc - N*avg_k))

    sign = 1

    if det <= 3:
        sign = -1
    
    Ns = 0.5 * (3 + sign * det)
    Ns = int(Ns)

    N1 = N - Nc - Ns
    
    #N += Nc - (N1 - 1)%Nc

    return N, Nc, N1, Ns

def optimal_network(N, Nc, avg_k, n_fi_nodes):
    """Hypothetical optimal network of size N for protecting the Nc critical
    function nodes."""
    
    N, Nc, N1, Ns = optimal_node_partition(N, Nc, avg_k)

    # initialize Graph
    G = nx.Graph()
    for i in range(N):
        G.add_node(i)
    
    # Critical Node Ids
    critical = pl.arange(Nc) * (N1 + Nc - 1)/Nc + 1
    critical = pl.round_(critical).astype(int)
    critical[0] = 0
    
    # Peripheral node IDs
    peripheral = pl.asarray([i for i in range(N - Ns) if i not in critical])
    
    # Sink Node IDs
    sink = pl.arange(N - Ns, N)

    # Connect Critical Nodes to Peripheral Nodes
    for i, crit in enumerate(critical):
        if i < Nc - 1:
            G.add_edge(crit, critical[i+1])
            connectees = peripheral[(peripheral > crit) & (
                peripheral < critical[i + 1])]
        else:
            connectees = peripheral[(peripheral > crit) & (
                peripheral < pl.inf)]
 
        for connect_id in connectees:
            G.add_edge(crit, connect_id)

    # connect Sink Nodes to Critical Nodes
    if len(sink) > 0:
        G.add_edge(critical[-1], sink[0])

    # Connect sink Nodes together
    for i, ID in enumerate(sink):
        for j, connect_ID in enumerate(sink[sink != ID]):
            G.add_edge(ID, connect_ID)
    
    mort_nodes = []
    # Make up some Mortality Node and FI Node stuff
    if Nc > 1:
        mort_nodes = [critical[0], critical[int(len(critical)/2)]]
    else:
        if len(sink) > 0:
            mort_nodes = [critical[0], sink[0]]
        
    important_nodes = [i for i in list(critical) + list(sink) +
        list(peripheral) if i not in mort_nodes]
    fi_nodes = []
    fi_nodes = important_nodes[:n_fi_nodes]
           
    return G, mort_nodes, fi_nodes

def get_scale_free_degrees(N, temp_folder, alpha, avg_deg, seed, pA = 0.0,
                    return_graph = True):
    """
    Use Spencer's code to get the degree sequence for scale free network of
    given parameters
    """
    # params for generating network
    params = ["0", "0.0", "0.0", "0.0", str(N), "2", str(alpha), str(avg_deg),
                "AND", "Single", temp_folder, "ScaleFree", "0",
                str(seed), str(pA)]
    #sub.call(["make", "backup"])
    command = ['./main'] + params 
    sub.call(command)

    output_files = os.listdir(temp_folder)
    #print(output_files)
    if not return_graph:
        degree_file = output_files[0]
    else:
        degree_file = [f for f in output_files if "Edge" not in f][0]
        edge_list_file = [f for f in output_files if "Edge" in f][0]
    degree_sequence = pl.loadtxt(temp_folder + degree_file)[:,1]
    if return_graph:
        G = nx.read_edgelist(temp_folder + edge_list_file)

    for f in output_files:
        os.remove(temp_folder + f)

    k_min = pl.amin(degree_sequence)
    k_max = pl.amax(degree_sequence)

    degree_bins = pl.arange(k_min, k_max + 2)
    degrees = degree_bins[:-1]
    dk = pl.histogram(degree_sequence, degree_bins)[0]
    degrees = degrees[dk != 0]
    dk = dk[dk != 0]

    
    if return_graph:
        return degrees, dk, degree_sequence, G
    else:
        return degrees, dk, degree_sequence



def power_law_distribution(alpha, k_min, k_max, target_avg_k):
    degrees = pl.arange(k_min, k_max+1)
    pk = np.power(degrees, alpha)
    prefactor = 1 / pl.sum(pk)
    pk = prefactor*pk
    avg_k = pl.dot(degrees, pk)
    cost = pl.absolute(target_avg_k - avg_k)
    return pk, avg_k, cost

def get_scale_free(alpha, k_min, N, target_avg_deg, tolerance = 0.01,
        output_progress = False):
    """
    A function for getting a power law distribution for a network for a given
    average degree and minimum degree.

    Note that given finite size constraints it could be impossible.
    
    Idea is to just use newtons method to get a good one
    """
    # some initial guesses at what the parameters are going to be
    prefactor = 10.0

    # apply the method
    learning_rate = 0.99

    k_max = 0.1 * N
    k_max = int(k_max)
    d_k_max = 1

    last_derivative = 1
    costs = []
    avg_degs = []
    pks = []
    cost = 10
    i = 0
    i_measure = 1
    measure_points = []
    k_maxes = []
    while cost > tolerance:
        # initial calculation
        pk, avg_k, cost = power_law_distribution(alpha, k_min, k_max,
                                target_avg_deg)
        if pl.sum((pk > 1) | (pk < 0)) > 0:
            print("PANIC")
            return "failure"
        degrees = pl.arange(k_min, k_max + 1)


        if k_max <= k_min:
            print("k max <= k min!")
            return "failure"
        if k_max > N:
            print("Need too big of degree")
            return "failure"

        if i == i_measure:
            i_measure *= 2
            if output_progress:
                avg_degs += [avg_k]
                costs += [cost]
                pks += [pk]
                measure_points += [i]
                k_maxes += [k_max]
            if i == 2**15:
                print("Timed Out")
                return "failure"

        # Get the derivatives in the directions
        derivative = (power_law_distribution(alpha, k_min, k_max + d_k_max,
                                        target_avg_deg)[2] - cost)/d_k_max

        k_max +=  - pl.sign(derivative) * d_k_max          
        i += 1

    if output_progress:
        avg_degs += [avg_k]
        costs += [cost]
        measure_points += [i]
        pks += [pk]
        k_maxes += [k_max]

        return (measure_points, pks, degrees, costs, avg_degs, k_maxes)
    else:
        return pk, degrees

def get_colours_from_cmap(values, cmap = 'viridis'):
    """
    Function for getting line / point colours from colormaps
    Better colorschemes to show sequences.
    e.g. when the legend is something like (a = 1, a = 2, a = 3, etc)
    """

    cm = mpl.cm.viridis(pl.linspace(0, 1, 100))
    cm = mpl.colors.ListedColormap(cm)
    c_norm = pl.Normalize(vmin=0, vmax=values[-1])
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap = cm)
    colours = [scalar_map.to_rgba(v) for v in values]

    return colours

def joint_matrix(dk):
    off_diagonal = pl.outer(dk, dk)
    diagonal = pl.array([math.comb(d, 2) for k, d in enumerate(dk)])
    jkk = off_diagonal
    for i in range(len(dk)):
        jkk[i,i] = diagonal[i]
    return jkk

def random_edge_sampling(probabilities, num_edges):
    """
    Assign the number of edges that nodes of a given degree will have with
    nodes of the set of other degrees
    """
    if pl.sum(probabilities != 0) == 1:
        assigned_edges = list(probabilities*num_edges)
        return assigned_edges
    assigned_edges = [1]
    num_sampled_edges = -1
    while num_sampled_edges != num_edges:

        assigned_edges = []
        
        num_sampled_edges = 0
        for i, p in enumerate(probabilities):
            leftovers = 0.0
            naturally_generated = 0
            leftover_generated = 0

            # calculate expected number
            expected_number = num_edges*p
            # some guaranteed number (rounded down integer)
            number = int(expected_number)
            num_sampled_edges += number
            assigned_edges.append(number)
            # probabilistically add another 
            leftovers = expected_number%1.0
            if leftovers >= pl.random():
                assigned_edges[-1] += 1
                num_sampled_edges += 1
        # the diagonal element being an odd number is bad news?
        if assigned_edges[-1]%2 == 1:
            # assign it to one of the other edges
            weights = pl.asarray(assigned_edges[:-1])
            if pl.sum(weights) == 0:
                continue
            probs = weights/pl.sum(weights)
            assigned_edges[-1] -= 1
            index = pl.choice(pl.arange(len(weights)), p = probs)
            assigned_edges[index] += 1

    return assigned_edges

# Now Defunked.
def weighted_edge_assignment(k_index, degrees, touched_degrees, stubs,
                degree_counts, p_assortative, p_dissassortative, p_random):
    """
    A more advanced edge assignment strategy.
    A mix of assortative/dissassortative/random assigments

    The idea is that we only attach in proceeding degree (anticipating that
    degree will be descending, so k >= k' by default)

    """
    assigned_stubs = pl.zeros_like(stubs)
    # case of nothing left to do (likely in dissassortative case)
    if pl.sum(stubs) == 0:
        return stubs, stubs
    # case of only one option left
    if pl.sum(stubs != 0) == 1:
        if pl.where(stubs != 0)[0] == k_index:
            assigned_stubs[k_index] = stubs[k_index]/2
        else:
            assigned_stubs[k_index] = stubs[k_index]
        return assigned_stubs, pl.zeros_like(stubs)
            
    k = degrees[k_index]
    k_stubs = stubs[k_index]
    dk = degree_counts[k_index]
    stubs_distribution = (
            pl.array([p_assortative, p_dissassortative, p_random]) * k_stubs)
    #print("Stubs Distribution:", stubs_distribution)
    if pl.sum(stubs_distribution.astype(int)) != pl.sum(stubs_distribution):
        if p_dissassortative > p_assortative:
            stubs_distribution[1] += 1
        elif p_assortative > 0:
            stubs_distribution[0] += 1
        else:
            stubs_distribution[2] += 1
    assortative_stubs, dissassortative_stubs, random_stubs = (
                                    stubs_distribution.astype(int))
    #print("Stubs Distribution:", stubs_distribution)
    """
    if pl.sum(stubs_distribution) == 0:
        return assigned_stubs, assigned_stubs
    """
    #print(assortative_stubs)
    
    #while assortative_stubs + dissassortative_stubs + random_stubs != k_stubs:
    #    print("Mismatch")
    #    random_stubs += 1

    size = len(stubs) 

    # attach randomly
    # Now, it is likely that some of the degrees are already at their limits.
    # This means that we will be randomly selecting amongst the AVAILABLE
    # nodes, not over all of the nodes (not sure if this is reasonable).
    if  pl.sum(stubs[:k_index+1]) == 0:
        return assigned_stubs, assigned_stubs
    edge_weights = stubs[:k_index+1]/pl.sum(stubs[:k_index+1])
    random_assignments = random_edge_sampling(edge_weights, random_stubs)
    #print("{0} Random Stubs, targets : ".format(random_stubs),
    #    random_assignments)
    #print(random_assignments)
    total_assigned_stubs = 0
    for k_prime_index in range(k_index+1):
        if p_random == 0:
            break
        # the number of stubs we will try to attach (could be that we don't get
        # there)
        target_stubs = random_assignments[k_prime_index]

        k_prime = degrees[k_prime_index]
        k_prime_stubs = stubs[k_prime_index]
        dk_prime = degree_counts[k_prime_index]
        node_constraint = dk * dk_prime 

        # carefull with diagonal elements
        if k_prime_index == k_index:
            # constraint along the diagonal is different
            node_constraint = dk*(dk-1) // 2
            # again check for preexisting stubs
            node_constraint -= assigned_stubs[k_prime_index]

            if node_constraint < 0:
                node_constraint = 0
            # available connections between these nodes
            available_stubs = pl.amin([target_stubs, node_constraint,
                                        k_prime_stubs])
            available_stubs = pl.amax([0, available_stubs])
            if k not in touched_degrees and pl.sum(stubs[k_index:] != 0) > 1:
                # Reduced by 1 to ensure connectedness (a bit wobbly, but
                # negligible in reasonable networks)
                available_stubs -= 1
                random_assignments[k_prime_index + 1] += 1
                # Deal with enforcing even number along diagonal 
                if available_stubs %2 == 1:
                    #print("Doing it")
                    available_stubs += 1

            # We "use" twice as many stubs as we report, since reported
            # stubs will be doubled in the symmetry imposing step
            used_stubs = available_stubs
            available_stubs /= 2
        # off diagonal we only have to consider the node constraint 
        else:
            available_stubs = pl.amin([target_stubs, node_constraint,
                                        k_prime_stubs])
            available_stubs = pl.amax([0, available_stubs])
            used_stubs = available_stubs

        touched_degrees.append(k_prime)
        touched_degrees.append(k)
        # assign the stubs
        assigned_stubs[k_prime_index] += available_stubs
        stubs[k_index] -= available_stubs
        stubs[k_prime_index] -= available_stubs
        total_assigned_stubs += used_stubs
    previous_assigned_stubs = pl.copy(assigned_stubs)
    #print(total_assigned_stubs , " Assigned Stubs after Random:",
    #    assigned_stubs)
    
    #### ASSORTATIVE WIRING PHASE
    # attach going from most similar down
    k_prime_index = k_index
    assigned_assortative_stubs = 0
    #print("Assortative Stubs:", assortative_stubs)
    while assortative_stubs > 0 and k_prime_index >= 0:
        dk_prime = degree_counts[k_prime_index]
        k_prime_stubs = stubs[k_prime_index]
        k_prime = degrees[k_prime_index]
        available_stubs = assortative_stubs

        if k_prime_stubs <= 0:
            k_prime_index -= 1
            continue
        if k_prime > k or available_stubs == 0:
            break


        # the node constraint is from single connection constraint
        node_constraint = dk * dk_prime 

        # carefull with diagonal elements
        if k_prime_index == k_index:
            # constraint along the diagonal is different
            node_constraint = dk*(dk-1) // 2
            # available connections between these nodes
            available_stubs = pl.amin([available_stubs, node_constraint,
                                        k_prime_stubs])
            available_stubs = pl.amax([0, available_stubs])
            # Reduced by 1 to ensure connectedness (a bit wobbly, but
            # negligible in reasonable networks)
            available_stubs -= 1
            # Deal with enforcing even number along diagonal 
            if available_stubs %2 == 1:
                available_stubs -= 1
            # We "use" twice as many stubs as we report, since reported
            # stubs will be doubled in the symmetry imposing step
            used_stubs = available_stubs
            available_stubs /= 2
        # off diagonal we only have to consider the node constraint 
        else:
            available_stubs = pl.amin([available_stubs, node_constraint,
                                        k_prime_stubs])
            available_stubs = pl.amax([0, available_stubs])
            used_stubs = available_stubs
        # assign the stubs
        assigned_stubs[k_prime_index] += available_stubs
        stubs[k_index] -= available_stubs
        stubs[k_prime_index] -= available_stubs
        # reduce the number left to go
        touched_degrees.append(k_prime)
        assortative_stubs -= used_stubs
        total_assigned_stubs += used_stubs
        assigned_assortative_stubs += used_stubs
        #print(used_stubs)
        k_prime_index -= 1
    """
    print(assigned_assortative_stubs,
            "Stubs assigned during assortative_phase:",
            assigned_stubs - previous_assigned_stubs)
    """
    previous_assigned_stubs = pl.copy(assigned_stubs)


    ##### DISSASSORTATIVE PHASE
    # attach going from the least similar up 
    k_prime_index = 0
    assigned_dissassortative_stubs = 0
    while dissassortative_stubs > 0 and k_prime_index < size:
        dk_prime = degree_counts[k_prime_index]
        k_prime_stubs = stubs[k_prime_index]
        k_prime = degrees[k_prime_index]
        available_stubs = dissassortative_stubs
        #print("k: {0}, Stubs left: {1}.".format(k, available_stubs))
        #print("k': {0}, Stubs left: {1}.".format(k_prime, k_prime_stubs))

        if k_prime_stubs == 0:
            k_prime_index += 1
            continue
        if k_prime > k or available_stubs == 0:
            break

        # the node constraint is from single connection constraint
        # there could already be edges between these nodes, have to account for
        # this in the constraints (likely not to be an issue for
        # assortative/dissassortative)
        node_constraint = dk * dk_prime 

        # carefull with diagonal elements
        if k_prime == k:
            #print("Diagonal")
            # constraint along the diagonal is different
            node_constraint = dk*(dk-1) // 2
            # again check for preexisting stubs
            node_constraint -= assigned_stubs[k_prime_index]
            # available connections between these nodes
            available_stubs = pl.amin([available_stubs, node_constraint,
                                        k_prime_stubs])
            available_stubs = pl.amax([0, available_stubs])
            #print(stubs, touched_degrees)
            if k not in touched_degrees and pl.sum(stubs[k_index:] != 0) > 1:
                # Reduced by 1 to ensure connectedness (a bit wobbly, but
                # negligible in reasonable networks)
                available_stubs -= 1
                # Deal with enforcing even number along diagonal 
                if available_stubs %2 == 1:
                    #print("Doing it")
                    available_stubs += 1
            # We "use" twice as many stubs as we report, since reported
            # stubs will be doubled in the symmetry imposing step
            used_stubs = available_stubs
            available_stubs /= 2
        # off diagonal we only have to consider the node constraint 
        else:
            available_stubs = pl.amin([available_stubs, node_constraint,
                                        k_prime_stubs])
            available_stubs = pl.amax([0, available_stubs])
            # additional concern with dissassortative connectinos is that higher
            # degree nodes go all in on degree 1 connections and become
            # disconnected from the rest of the network
            if k_prime not in touched_degrees and (
                        len(touched_degrees) < 1):
                # enforce 2 connections per node to be outside of this degree
                if k_prime_stubs > available_stubs:
                    available_stubs = pl.amax([available_stubs - 2*dk, 0])
            used_stubs = available_stubs

        touched_degrees.append(k_prime)
        #print(touched_degrees)
        # assign the stubs
        assigned_stubs[k_prime_index] += available_stubs
        stubs[k_index] -= available_stubs
        stubs[k_prime_index] -= available_stubs
        #print("Attached {0}.".format(used_stubs))
        # reduce the number left to go
        dissassortative_stubs -= used_stubs
        total_assigned_stubs += used_stubs
        assigned_dissassortative_stubs += used_stubs
        k_prime_index += 1
    """
    print(assigned_dissassortative_stubs,
            "Stubs assigned during dissassortative_phase:",
            assigned_stubs - previous_assigned_stubs)

    print("total Assigned stubs: ", total_assigned_stubs)
    """
    return assigned_stubs, touched_degrees

def sampled_edge_assignment(k_index, degrees, degree_counts, stubs,
                                                                pa, pd, pr):
    """
    Similar to weighted edge assignment, but different approach
    Resolve any bullshit at the end
    """
    do_print = False #(len(degree_counts) - k_index < 5)
    size = len(degrees)
    available_stubs = pl.copy(stubs)
    k_stubs = stubs[k_index]
    k = degrees[k_index]
    dk = degree_counts[k_index]

    #print("Starting with ", k_stubs, "Stubs")
    #print(len(stubs), len(degree_counts), dk)
    node_constraints = dk*degree_counts
    node_constraints[k_index] = dk*(dk-1)//2

    if do_print:
        print("Stubs:",stubs)
        print("Node Constraints:",node_constraints)

    # have to run everything off of the constrained availability of nodes
    # likely to only matter along the diagonal?
    #print(node_constraints, stubs)
    available_stubs = pl.amin(
                        pl.column_stack([stubs, node_constraints]), axis = 1)
    #print(available_stubs)
    #print(available_stubs)

    assigned_stubs = pl.zeros_like(stubs)
    
    assortative_condition = True

    # current target degrees for assortative/dissassortative
    assortative_index = k_index
    dissassortative_index = 0
    k_prime_index = 0
    num_attempts = 0
    while stubs[k_index] > 0:
        num_attempts += 1
        if num_attempts % 100000 == 0:
            print(stubs[k_index], "Left to go")
            print(degrees[available_stubs>0],
                    available_stubs[available_stubs>0],
                    stubs[available_stubs>0])
            break


        if stubs[k_index] > pl.sum(available_stubs):
            print(stubs[k_index], "Left to go")
            print(degrees[available_stubs>0],
                    available_stubs[available_stubs>0],
                    stubs[available_stubs>0])


        min_k_stubs = pl.amin([stubs[k_index], available_stubs[k_index]])

        if pl.sum(available_stubs) == 0:
            print("Ran Out of available stubs!!!")
            break
        # select how to wire everything up
        method = pl.choice(pl.arange(3), p = [pa, pd, pr])

        # assortative choice
        if method == 0:
            k_prime_index = assortative_index
            # if we've run out of stubs there, find the next best target
            while available_stubs[k_prime_index] == 0:
                k_prime_index -= 1
                assortative_index = k_prime_index
            """
            if do_print:
                print("k:", degrees[k_index], "Assortative Index:",
                        assortative_index,
                        "Degree:, ", degrees[assortative_index])
            """

        # dissassortative choice
        elif method == 1:
            k_prime_index = dissassortative_index
            # if the k nodes have only attached to 1 neighbour so far:
            if stubs[k_index] <= 2*dk and dissassortative_index == 0:
                k_prime_index += 1
                dissassortative_index = k_prime_index
            if pl.sum(assigned_stubs > 0) == 1 and stubs[k_index] < 2:
                k_prime_index += 1
                dissassortative_index = k_prime_index
             
            # if we've run out of stubs there, find the next best target
            while available_stubs[k_prime_index] == 0:
                k_prime_index += 1
                dissassortative_index = k_prime_index
                if dissassortative_index == size:
                    k_prime_index = 0
                    dissassortative_index = k_prime_index

        # random choice
        else:
            edge_weights = pl.copy(available_stubs)
            edge_weights[k_index] /= 2.0
            probabilities = edge_weights/pl.sum(edge_weights)

            try:
                k_prime_index = pl.choice(size, p = probabilities)
            except ValueError:
                print("\n Probabilities Contain NaN\n")
                break

        # picking up assortative problems
        if k_prime_index == size:
            break
            

        # if there's only one stub left, you cant assign it on the diagonal
        #print(available_stubs, stubs, assigned_stubs)
        if (min_k_stubs <= 2) and k_index == k_prime_index and (
            pl.sum(stubs != 0) > 1):
            #print("Doin illegal stuff")

            if assortative_condition:
                assortative_index -= 1
                assortative_condition = False
            continue
            
        if available_stubs[k_prime_index] == 1 and (
                pl.sum(available_stubs != 0) == 2) and (min_k_stubs >= 2):
            #print("\n\nOnly One Move Left")
            #print(available_stubs)

            k_prime_index = k_index


        available_stubs[k_prime_index] -= 1
        stubs[k_prime_index] -= 1
        stubs[k_index] -= 1
        if k_index == k_prime_index:
            available_stubs[k_index] -= 1
        assigned_stubs[k_prime_index] += 1
    if k_prime_index == size:
        print(k_index, available_stubs, assigned_stubs)

    if do_print:
        print("Degree:", degrees[k_index], "Assigned Stubs",
            assigned_stubs)

    return assigned_stubs, stubs

def degree_1_edge_assignment(k_index, degrees, degree_counts, stubs,
                                                                pa, pd, pr):
    """
    Degree 1 stuff can get fiddly quickly: immediately doom yourself on the
    connections from 1 to 1 (or fill e.g. all 2s to 1s, etc).

    Really it's the same as sampled edge assignment except you have to count up
    from 2 for assortative, down from size for dissassortative.

    Might also have some wobbly edge cases to deal with due to k = 1
    """
    size = len(degrees)
    available_stubs = pl.copy(stubs)
    k_stubs = stubs[k_index]
    k = degrees[k_index]
    dk = degree_counts[k_index]


    #print(len(stubs), len(degree_counts), dk)
    node_constraints = dk*degree_counts
    node_constraints[k_index] = dk*(dk-1)//2

    
    node_constraints[1] = 2
    """
    pl.amin([
        (degree_counts[1]*(degree_counts[1] - 1)//2 - 1)//2,
        (stubs[1]//2 - 1)//2]) - 1
    print(node_constraints[1])
    """
    # have to run everything off of the constrained availability of nodes
    # likely to only matter along the diagonal?
    #print(node_constraints, stubs)
    available_stubs = pl.amin(
                        pl.column_stack([stubs, node_constraints]), axis = 1)
    #print(available_stubs)

    assigned_stubs = pl.zeros_like(stubs)
    
    assortative_condition = True

    # current target degree for assortative
    assortative_index = k_index + 1
    # target degree for dissassortative
    dissassortative_index = size - 1
    while stubs[k_index] > 0:
        min_k_stubs = pl.amin([stubs[k_index], available_stubs[k_index]])

        if pl.sum(available_stubs) == 0:
            print(assigned_stubs)
            print("Ran Out of available stubs!!!")
            break
        # select how to wire everything up
        method = pl.choice(pl.arange(3), p = [pa, pd, pr])

        # assortative choice
        if method == 0:
            k_prime_index = assortative_index
            #if min_k_stubs == 1:
            #    k_prime_index += 1
            if (assigned_stubs[k_prime_index] >= 0.5*
                degree_counts[k_prime_index]*(degrees[k_prime_index] - 1) - 1):
                k_prime_index += 1
                assortative_index = k_prime_index
                
            # if we've run out of stubs there, find the next best target
            while available_stubs[k_prime_index] == 0:
                k_prime_index += 1
                assortative_index = k_prime_index

        # dissassortative choice
        elif method == 1:
            k_prime_index = dissassortative_index
            # if the k nodes have only attached to 1 neighbour so far:
            if stubs[k_prime_index] <= 2*degree_counts[k_prime_index]:
                k_prime_index -= 1
                dissassortative_index = k_prime_index
            if min_k_stubs == 1 and k_prime_index == size - 1:
                k_prime_index -= 1
            # if we've run out of stubs there, find the next best target
            while available_stubs[k_prime_index] == 0:
                k_prime_index -= 1
                dissassortative_index = k_prime_index

        # random choice
        else:
            edge_weights = pl.copy(available_stubs)
            edge_weights[0] = 0
            probabilities = edge_weights/pl.sum(edge_weights)
            if pl.nan in probabilities or pl.sum(available_stubs) == 0:
                print("\n\n\n", probabilities, "\n\n\n")
                break
            k_prime_index = pl.choice(size, p = probabilities)

        # picking up assortative problems

        # if there's only one stub left, you cant assign it on the diagonal
        #print(available_stubs, stubs, assigned_stubs)
        if (min_k_stubs <= 2) and k_index == k_prime_index and (
            pl.sum(stubs != 0) > 1):
            #print("Doin illegal stuff")
            if assortative_condition:
                assortative_index -= 1
                assortative_condition = False
            continue
            
        if available_stubs[k_prime_index] == 1 and (
                pl.sum(available_stubs != 0) == 2) and (min_k_stubs >= 2):
            #print("\n\nOnly One Move Left")
            #print(available_stubs)
            k_prime_index = k_index


        available_stubs[k_prime_index] -= 1
        stubs[k_prime_index] -= 1
        stubs[k_index] -= 1
        if k_index == k_prime_index:
            available_stubs[k_index] -= 1
        assigned_stubs[k_prime_index] += 1
    #print(available_stubs, assigned_stubs)
    return assigned_stubs, stubs

def weighted_connection_matrix(degrees, degree_counts,
    p_assortative, p_dissassortative, p_random):
    stubs = degrees*degree_counts
    size = len(degree_counts)
    jkk = pl.zeros([size, size])
    reverse_degrees = pl.copy(degrees)[::-1]
    #print("About to start assigning edges", len(degrees), len(degree_counts))
    # if there are degree 1 nodes deal with them first
    
    use_degree_1 = 0 #degrees[0] == 1
    if use_degree_1:
        print("We Usin Degree 1!")
        k_index = 0
        assigned_stubs, stubs = degree_1_edge_assignment(
                            k_index, degrees, degree_counts, stubs,
                            p_assortative, p_dissassortative, p_random)
        jkk[:,k_index] = assigned_stubs
    for i, k in enumerate(reverse_degrees):
        # if we already assigned the degree 1s, continue
        if use_degree_1 and k == 1:
            continue
        k_index = size - 1 - i
        #print("Workin on k = ", degrees[k_index])
        #print("Assigning stubs for k =", degrees[k_index])
        assigned_stubs, stubs = sampled_edge_assignment(
                            k_index, degrees, degree_counts, stubs,
                            p_assortative, p_dissassortative, p_random)
        #print("\n\n", degrees[k_index], stubs, len(degrees), len(stubs),
        #    len(degree_counts))
        #print(assigned_stubs)
        #stubs -= attached_stubs
        jkk[:,k_index] = assigned_stubs
    jkk += jkk.T
    if degrees[0] == 1:
        print("Sum of Jkk:", pl.sum(jkk), "J(1, 1)=", jkk[0,0])
        jkk[0,0] = 0.0
    return jkk

def build_weighted_graph(degrees, degree_counts,
    p_assortative, p_dissassortative, p_random, return_jkk = False):
    n_tries = 10
    N = pl.sum(degree_counts)
    for n in range(n_tries):
        #print(n, "Tries so far")
        jkk = weighted_connection_matrix(degrees, degree_counts,
                            p_assortative, p_dissassortative, p_random)
        G = network_from_jkk(jkk, degrees, N)
        if type(G) != str:
            #print("Built a Graph with {} Nodes?".format(G.number_of_nodes()))
            #print("Odd considering N =", N)
            if G.number_of_nodes() > 0.9 * N:
                if return_jkk:
                    return G, jkk
                else:
                    return G

    print("Couldn't Build the Graph (Valid Jkk)")
    #pl.figure(figsize = (8,6))
    #pl.imshow(jkk, origin = 'lower')
    #pl.colorbar()
    #print(jkk)
    non_1_edges = pl.sum(jkk[:,1:] > 0, axis = 1)
    print(len(non_1_edges), pl.shape(jkk))
    print(non_1_edges - jkk[:,0])
    print(degree_counts * (degrees - 1)/degrees - jkk[:,0])
    #print(pl.sum(jkk, axis = 1)/degrees)
    return "Couldn't Generate Graph", "Hell"

def connection_dict(nkk, degrees):
    size = pl.shape(nkk)[0]
    jkk = {}
    for i, k in enumerate(degrees):
        jkk[k] = {}
        for j, k_prime in enumerate(degrees):
            jkk[k][k_prime] = nkk[i, j].astype(int)
            
    return jkk

def random_hub_node_degrees(degree_distribution, degrees, N):
    degree_counts = []
    degree_sequence = []
    degree_sum = 1
    while len(degree_counts) < 1 or degree_sum %2 == 1:
        degree_sequence = []
        degree_counts = []
        for n, pk in enumerate(degree_distribution):
            k = degrees[n]
            leftovers = 0.0
            naturally_generated = 0
            leftover_generated = 0

            expected_number = N*pk
            number = int(expected_number)
            degree_counts.append(number)
            for i in range(number):
                degree_sequence.append(k)
                naturally_generated += 1
            leftovers = expected_number%1.0
            if leftovers >= pl.random():
                degree_sequence.append(k)
                degree_counts[-1] += 1
                leftover_generated += 1
        degree_sum = pl.sum(degree_sequence)
    return pl.asarray(degree_sequence), pl.asarray(degree_counts)

def check_graphical(degree_sequence):
    N = len(degree_sequence)
    min_constraint = pl.amin(pl.column_stack([degree_sequence, pl.arange(N)]),
                            axis = 1)
    for k, d in enumerate(degree_sequence):
        lhs = pl.sum(degree_sequence[:k + 1])
        rhs = k*(k-1) + pl.sum(min_constraint[k:])
        if lhs > rhs:
            return False
    return True

def change_pk(original_pk, degrees, n1, n2, frac, avg_deg):
    """
    Change a discrete degree distribution by boosting degrees k1 and/or k2
    (depending on the fraction) while maintaining normalization and average
    degree

    Returns a degree distribution and a boolean value which indicates whether
    or not it has been successfully modified
    """
    pk = pl.copy(original_pk)
    change_success = False

    p1 = pk[n1]
    p2 = pk[n2]
    k1 = degrees[n1]
    k2 = degrees[n2]

    # impossible to make this work if k1 == k2
    if k1 == k2:
        #print("k1 = k2")
        return original_pk, change_success
    
    # dont let constants diverge
    if (p1*k1 + p2*k2) == avg_deg or p1 + p2 == 1:
        #print("diverging a or b")
        return original_pk, change_success
    
    # calculate constants to satisfy constraints
    a = 1/(avg_deg - (p1*k1 + p2*k2))
    b = 1/(1-(p1 + p2))

    # dont let the limit on p1 diverge
    if (a*k1 - b) == 0:
        #print("diverging p1")
        return original_pk, change_success
    # calculate the limit
    new_p1_limit = (a * avg_deg - b)/(a*k1 - b)
    if new_p1_limit < 0:
        #print("Negative p1")
        return original_pk, change_success
    
    # select the values that satisfy the relative weighting
    new_p1 = frac*new_p1_limit

    # Calculating p2: dont allow divergence
    if (b - a * k2) == 0:
        #print("Diverging p2")
        return original_pk, change_success
    # calculate p2
    new_p2 = (b - a*avg_deg - new_p1*(b - a * k1))/(b - a * k2)

    # ensure that everything is satisfied (typical issue is p2 < 0)
    if new_p2 > 1 or new_p2 < 0 or new_p2 + new_p1 > 1:
        #print("Illegal")
        return original_pk, change_success
    
    # calculate scale factor to renormalize the degree distribution
    # shouldn't diverge since denominator checked earlier
    scale_factor = (avg_deg - (k1*new_p1 + k2*new_p2))/(
                                                    avg_deg - (k1*p1 + k2*p2))
    # scale factor must be positive
    if scale_factor < 0:
        #print("Negative scale factor")
        return original_pk, change_success

    # All possible failure cases have been successfully passed so we good to go
    change_success = True

    # set values and renormalize pk
    new_pk = scale_factor * pl.copy(original_pk)
    new_pk[n1] = new_p1
    new_pk[n2] = new_p2
    return new_pk, change_success

def jkk_limits(dk, degrees):
    size = len(degrees)
    constraints = pl.full([size, size, 3], 1e6)
    edges = dk*degrees
    # First constraint: Product of node numbers (off diagonal simple constraint)
    node_count_product = pl.outer(dk, dk)
    constraints[:,:,0] = node_count_product
    # First constraint: On the diagonal
    for i in range(size):
        constraints[i,i,0] = dk[i]*(dk[i] - 1)
    # Also don't let the 1s attach for the love of god
    constraints[0,0,0] = 0
    
    # Second constraint: Number of edges of that degree (on one end)
    constraints[:,:,1] = pl.tile(edges, (size, 1))
    # Third constraint: Number of edges of that degree (on other end)
    constraints[:,:,2] = pl.tile(edges, (size, 1)).T
    
    return pl.amin(constraints, axis = 2)

def conditional_from_jkk(jkk, degrees):
    size = len(degrees)
    degrees_tiled = pl.tile(degrees, (size,1))
    dk = pl.sum(jkk, axis = 1)/degrees
    dk_tiled = pl.tile(dk, (size, 1))
    return jkk/dk_tiled/degrees_tiled

def jkk_from_dk(degree_counts, degrees, method = "dissassortative"):
    size = len(degree_counts)
    # reverse these suckers (sort degrees descending)
    degree_counts = degree_counts[::-1]
    degrees = degrees[::-1]
    # indices corresponding to degrees (helpful for non-spanning degrees)
    indices = pl.arange(size)

    edges = degree_counts*degrees
    # this is basically edges but dynamic (will be reduced as wired up)
    stubs = pl.copy(edges)
    # connection matrix (to be filled)
    jkk = pl.zeros([size, size])

    # if doing dissassortative attachment, prime degrees count up (as degrees
    # count down)
    prime_indices = pl.copy(indices)
    prime_degrees = pl.copy(degrees)
    if "dis" in method or "Dis" in method:
        prime_degrees = prime_degrees[::-1]
        prime_indices = prime_indices[::-1]

    # keep track of degrees which are connected to at least 2 types of nodes 
    touched_degrees = []
    for i, k in enumerate(degrees):
        #print("k = {0}, Stubs: {1}".format(k, sk))
        forcing_neighbour = True
        k_index = indices[i]
        dk = degree_counts[k_index]
        for j, k_prime in enumerate(prime_degrees):
            k_prime_index = prime_indices[j]
            dk_prime = degree_counts[k_prime_index]
            # Get stubs of degree k
            sk = stubs[k_index]
            # if no stubs of degree k left go to next k
            if sk == 0:
                break

            # stubs of degree k prime
            sk_prime = stubs[k_prime_index]
            if sk_prime == 0:
                continue
            #print("k = {0}, Dk = {2}, Stubs: {1}".format(k, sk, dk),
            #    "k' = {0}, Dk' = {2}, Stubs: {1}".format(k_prime,
            #                                    sk_prime,
            #                                    dk_prime))
            # Easy constriants: Cant put more edges than abailable stubs
            constraint_list = [sk, sk_prime]


            last_degree = pl.sum(pl.asarray(stubs) != 0) == 1
            #print(pl.asarray(stubs) != 0)
            # this is to make sure that things end up connected (e.g. if there
            # are 3 star formations, there should be at least 3 connections
            # from the stars hubs to the rest of the nodes
            if k not in touched_degrees and k_prime not in touched_degrees:
                #print(touched_degrees)
                # basically don't let the stars go all in (max is all in -
                # number of nodes of that degree (modulo constant)
                if k != k_prime and not last_degree:
                    constraint_list.append(max([sk - 2*dk, 0]))

            #print("k' = {0}, Stubs: {1}".format(l, sl))
            # diagonal constraint
            # NOTE: we also dont let a diagonal element go full death-ball
            # we save at least 1 edge to connect to another degree
            if k_prime == k:
                if k == 1:
                    continue
                constraint_list.append(dk*(dk-1)/2)
                if k in touched_degrees or last_degree:
                    constraint_list.append(sk_prime//2)
                else:
                    constraint_list.append(sk_prime//2 - 1)
                    
            else:
                constraint_list.append(dk*dk_prime)

            #print("Constraint: ", [sl, sk] + constraint_list)
            max_edges = pl.amin(constraint_list)
            #print(constraint_list)
            if max_edges > 0:
                #print("Adding Edges")
                if k not in touched_degrees:
                    touched_degrees.append(k)
                #if k_prime not in touched_degrees:
                #    touched_degrees.append(k_prime)

            #print("Max Edges: ", max_edges)
            stubs[k_index] -= max_edges
            stubs[k_prime_index] -= max_edges
            jkk[k_index,k_prime_index] = max_edges

    jkk += jkk.T
    jkk = jkk[::-1,::-1]
    #print(jkk)
    return jkk

def output_network(G, output_file):
    """
    output the network structure in the form of a CSV Edge list
    """
    nx.write_edgelist(G, output_file, delimiter = ",", data = False)

def build_network(pk, degrees, N, return_pkk = False,
                    method = "dissassortative"):

    # Build a degree sequence that is simply realizable
    max_attempts = 10000
    num_attempts = 0
    while num_attempts < max_attempts:
        num_attempts += 1
        degree_sequence_attempts = 0
        while degree_sequence_attempts < max_attempts:
            degree_sequence_attempts += 1
            degree_sequence, dk = random_hub_node_degrees(pk, degrees, N)
            # Check if it's realizable (if it is we're done here)
            if check_graphical(degree_sequence):
                break
        if degree_sequence_attempts == max_attempts - 1:
            print("Couldn't sample P(k) properly")
            return 'hell', "bubcus"

        # Build up a maximally disassortative connection matrix
        jkk = jkk_from_dk(dk, degrees, method)
        # Convert it to dictionary for nx functions
        jkk_dict = connection_dict(jkk, degrees)
        #print(jkk_dict)
        # check that it works
        #print("Checking Validity")
        if nx.is_valid_joint_degree(jkk_dict):
            # build the graph
            #print("Building Graph")

            # Here we run into some wobbly coding in the networkx implentation
            # of the algorythm
            try:
                G = nx.joint_degree_graph(jkk_dict)
            # The exception seems to be looking for a neighbour to swap to that
            # does not exist (next(whatever) is getting to the end of the
            # iterator)
            except StopIteration:
                continue
            # ensure it's fully connected
            #print("Checking Connectedness")
            if len(list(nx.connected_components(G))) == 1:
                if return_pkk:
                    return G, jkk/pl.sum(jkk)
                else:
                    return G
        
    print("Couldn't sample Network properly")
    return 'hell', "bubcus"

def run_network(G, input_dir, data_dir, params, param_names, 
    measures = ["DeathAge"]):
    """
    Shove the network to the GNM to get whatever measure (e.g. deathages) back
    """

    N = str(G.number_of_nodes())
    # update the number of nodes (it could change by a bit)
    params[1] = N
    initial_distribution = params[-4]
    #print(params)
    network_travel_start = time.time()
    output_network(G, input_dir + "SimplyConnected.csv")
    print("Outputting to " + input_dir + "SimplyConnected.csv")

    network_travel_time = time.time() - network_travel_start
    ### running stuff
    # compile if needed
    sub.call(["make", "testNetwork"])
    command = ['./testNetwork'] + params 
    raw_start_time = time.time()
    sub.call(command)
    raw_simulation_time = time.time()  - raw_start_time
    print("Real Simulation Time: {0:.3f}, Network Travel() time: {1:.3f}".format(
        raw_simulation_time, network_travel_time))

    print("Reading data from " + data_dir)

    output_files = os.listdir(data_dir)
    print(output_files)

    performances = []
    for m in measures:
        performance_file = [data_dir + f for f in output_files if m in f][0]
        performance = pl.loadtxt(performance_file)
        if "Norm" not in m:
            performance /= scale
        performances.append(performance)

    for f in output_files:
        os.remove(data_dir + f)

    if len(performance) == 1:
        return performances[0]
    else:
        return performances

def calculate_entropy(p):
    valid_p = p[p > 0]
    return -pl.sum(valid_p * pl.log(valid_p))

def consolidate_components(connected_components, degree_sequence, G):
    """
    Genius Idea: connect the disconnected components myself
    
    only gonna bother with the k = 1 case.
    General case is hard
    """
    print("We Consolidating")
    giant_component = pl.asarray(list(connected_components[0]))
    #print(giant_component)
    degree_sequence = pl.asarray(degree_sequence)
    giant_degrees = degree_sequence[giant_component]

    #giant_degrees = pl.asarray([k for node_id, k in giant_component.degrees()])
    other_components = connected_components[1:]
    #print(len(other_components), "Other components to sort out")

    # Neva Give up!
    give_up = False
    # fix all the small components
    num_successfull = 0
    for component in other_components:
        # need the degrees
        component = list(component)
        degrees = pl.asarray([degree_sequence[i] for i in component])
        #print("Component Degrees:", degrees)
        # if all the degrees except 1 are equal to 1 it's fixable
        # if it isn't, we still got a chance?
        #if pl.sum(degrees == 1) != len(degrees) - 1:
        #    give_up = True
        #    break
        # gonna swap one of the degree ones with another node of the not 1
        # degree
        if pl.sum(degrees != 1) == 1:
            old_child_index = pl.choice(pl.where(degrees == 1)[0])
        else:
            give_up = True
            break
        old_child = component[old_child_index]
        high_k_index = pl.amax(degrees)
        old_parent_index = pl.where(degrees == high_k_index)[0][0]
        old_parent = component[old_parent_index]
        # yoink a potential swappee from the giant component
        potential_nodes = giant_component[giant_degrees == high_k_index]
        #print("Potential Nodes", potential_nodes)
        # check for suitable parents (basically if they can swap one of their
        # not degree 1 nodes with a degree 1 node)
        successful_adoption = False
        for new_parent in potential_nodes:
            brother_degrees = pl.asarray([degree_sequence[i] for i in
                                G[new_parent]])
            #print("Brother Degrees:", brother_degrees)
            if pl.sum(brother_degrees != 1) < 3:
                continue
            new_child_index = pl.choice(pl.where(brother_degrees != 1)[0])
            new_child = list(G[new_parent])[new_child_index]
            #print(new_parent, new_child)
            #print(old_parent, old_child)
            G.remove_edge(old_parent, old_child)
            G.remove_edge(new_parent, new_child)
            G.add_edge(new_parent, old_child)
            G.add_edge(old_parent, new_child)
            successful_adoption = True
            num_successfull += 1
            #print(num_successfull, "Successful Adoptions")
            break

        if not successful_adoption:
            give_up = True
            break
            

    if give_up:
        return "Too Hard to Fix."
    else:
        connected_components = list(nx.connected_components(G))
        if len(connected_components) == 1:
            print("Smoked Him")
            return G
        else:
            return "Blew it again man"


def network_from_jkk(jkk, degrees, N, return_pkk = False):
    jkk_dict = connection_dict(jkk, degrees.astype(int))
    #print(jkk_dict)
    # check that it works
    #print("Checking Validity")
    tries = 10
    if nx.is_valid_joint_degree(jkk_dict):
        for t in range(tries):
            # build the graph
            #print("Building Graph")

            # Here we run into some wobbly coding in the networkx implentation
            # of the algorythm
            try:
                G = nx.joint_degree_graph(jkk_dict)
                #print("Just Built Graph, N =", G.number_of_nodes())
            # The exception seems to be looking for a neighbour to swap to that
            # does not exist (next(whatever) is getting to the end of the
            # iterator)
            except StopIteration:
                print("Mystery Breakdown: Stop Iteration")
                continue
            #except TypeError:
            #    print("Mystery Breakdwon: Type Error")
            #    continue
            # ensure it's fully connected
            connected_components = list(nx.connected_components(G))
            connected_components.sort(key=len, reverse=True)
            biggest_component = connected_components[0]
 
            num_components = len(connected_components)
            if num_components == 1:
                if return_pkk:
                    return G, jkk/pl.sum(jkk)
                else:
                    return G
            elif len(biggest_component) >= 0.99 * N:
                G = G.subgraph(biggest_component)
                print("Using the Giant Component")
                G = nx.relabel.convert_node_labels_to_integers(G,
                                    first_label=0, ordering='default')
                if return_pkk:
                    return G, jkk/pl.sum(jkk)
                else:
                    return G 

            else:
                degree_sequence = [k for node_id, k in G.degree()]
                G = consolidate_components(connected_components,
                                                        degree_sequence, G)
                if type(G) != str:
                    if return_pkk:
                        return G, jkk/pl.sum(jkk)
                    else:
                        return G

        # ran out of tries
        print("Couldn't get a connected one")
        for sub_G in connected_components[1:]:
            sub_graphs = [degree_sequence[i] for i in sub_G]
            print(sub_graphs)


        if return_pkk:
            return "Invalid G", "Invalid Pkk"
        else:
            return "Invalid G"
 
    else:
        if return_pkk:
            return "Invalid G", "Invalid Pkk"
        else:
            return "Invalid G"

def build_motif_graph(motif, N):
    """
    function for building some motif graphs
    """
    G = nx.empty_graph(N)
    if motif == "star":
        for i in range(1, N):
            G.add_edge(0, i)
    
    elif motif == "ball":
        for i in range(N):
            for j in range(i + 1, N):
                G.add_edge(i, j)
    elif motif == "chain":
        for i in range(N):
            G.add_edge(i, (i+1)%N)

    elif motif == "parents":
        if N == 2:
            G.add_edge(0,1)
        for i in range(2, N):
            G.add_edge(i, 0)
            G.add_edge(i, 1)
    else:
        print("Pick from: star, ball, chain.")
        return "Invalid Motif"

    return G

def show_log_jkk(degrees, dk, degree_sequence, G, title = "Log Jkk"):
    print("Number of Nodes:", G.number_of_nodes(),
        "Number of edges:", G.number_of_edges())
    edge_list = [e for e in G.edges]
    #print("Number of Edges", len(edge_list))
    #print("First 10:", edge_list[:10])
    edge_degrees = [[G.degree(i) for i in e] for e in edge_list]
    #print("Degrees of First 10:", edge_degrees[:10])
    degree_ranks = pl.arange(len(degrees))
    degrees = pl.asarray(sorted(list(set([k for _, k in G.degree()]))))
    #print(degrees)
    #print(degrees)
    #print(pl.arange(len(degrees)))
    jkk = pl.zeros([len(degrees), len(degrees)])
    for n, e in enumerate(edge_list):
        l, m = [int(node_id) for node_id in e]
        k_l, k_m = edge_degrees[n]
        i = pl.where(degrees == k_l)[0][0]
        j = pl.where(degrees == k_m)[0][0]

        jkk[i,j] += 1
        jkk[j,i] += 1
    #print("Sum of Jkk:", pl.sum(jkk))
    print("Row-Wise:", pl.sum(jkk, axis = 1).astype(int)/degrees)
    jkk[jkk == 0] = pl.nan
    #jkk[~pl.isnan(jkk)] = pl.log(jkk[~pl.isnan(jkk)])
    pl.figure(figsize = (8,6))
    pl.imshow(jkk, origin = 'lower')
    degrees = degrees.astype(int)
    pl.xticks(pl.arange(len(degrees)), labels = degrees)
    pl.yticks(pl.arange(len(degrees)), labels = degrees)

    pl.title(title)
    pl.colorbar()

def swap_edges(e0, e1, k0, k1):
    base_r = (pl.diff(k0) + pl.diff(k1))[0]
    rs = [base_r]
    edge_pairs= [[e0, e1]]
    degree_pairs = [[k0, k1]]
    for i in range(2):
        r = k0[0] - k1[i] + k0[1] - k1[(i+1)%2]
        edge_pair = [[e0[0], e1[i]], [e0[1], e1[(i+1)%2]]]
        degree_pair = [[k0[0], k1[i]], [k0[1], k1[(i+1)%2]]]
        edge_pairs.append(edge_pair)
        degree_pairs.append(degree_pair)
        rs.append(r)
    rs = pl.asarray(rs)
    best_index = pl.where(rs == pl.amax(rs))[0][0]
    best_edges = edge_pairs[best_index]
    best_degrees = degree_pairs[best_index]
    best_r = pl.amax(rs)
    difference = r - best_r
    return best_edges, best_degrees, difference



def shuffle_edges(G, swaps_per_edge = 1):
    edge_list = [e for e in G.edges]
    num_edges = len(edge_list)
    total_swaps = int(num_edges * swaps_per_edge)
    edge_degrees = [[G.degree(i) for i in e] for e in edge_list]
    indices = pl.arange(num_edges)
    r_differences = []
    for i in range(total_swaps):
        i0, i1 = pl.choice(indices, replace = False, size = 2)
        e0 = edge_list[i0]
        e1 = edge_list[i1]
        k0 = edge_degrees[i0]
        k1 = edge_degrees[i1]
        new_edges, new_degrees, r_difference = swap_edges(e0, e1, k0, k1)
        new_e0, new_e1 = new_edges
        new_k0, new_k1 = new_degrees
        r_differences.append(r_difference)
        edge_list[i0] = new_e0
        edge_list[i1] = new_e1
        edge_degrees[i0] = new_k0
        edge_degrees[i1] = new_k1
    
    return edge_list, pl.asarray(r_differences)
     

if __name__ == '__main__':

    pl.close('all')
    #pl.seed(2)
    fs = 14
    matrix_folder = 'AdjacencyMatrices/'
    edge_list_folder = "EdgeLists/"
    plots_folder = 'Plots/AssortativityWeightedNetworks/'

    ##### Setting parameters for network
    N = 500
    avg_k = 4
    scale = 0.00183
    ## parameters of running the simulaton
    number = "1000"
    output_folder = 'SimplyConnectedData/'

    seed  = '1'
    Lambda = "0.0"
    beta = "100.0"
    power = "1.0"
    numChanges = "0"
    run_hours = "0.0"
    evo_condition = 'none'
    initial_distribution = "SimplyConnected"

    param_names = ['number', 'N', 'numChanges', 'avgDeg', 'folder',
        'singleseed', 'evoCondition', 'InitialDistribution',
        'lambda', 'beta', 'power']
    params = [number, N, numChanges, avg_k, output_folder, seed,
        run_hours, evo_condition, initial_distribution, Lambda,
        beta, power]

    test_power_law = 0
    if test_power_law:
        plots_dir = "Plots/ScaleFreeChecks/"
        avg_k = 4
        alpha = 2.27
        alphas = [2.001, 2.1, 2.3, 3.0, 5.0]
        N = 10000
        n_samp = 1000
        for alpha in alphas:
            k = pl.arange(0, N)
            reference = k**(-alpha)
            pl.figure(figsize = (8,6))
            pl.plot(k, reference, label = r'Power Law $\alpha = {0}$'.format(alpha))
            pl.xscale('log')
            pl.yscale('log')

            bins = pl.arange(N+1)
            n_max = pl.log2(N).astype(int) + 1

            log_bins = pl.logspace(0, n_max - 1, n_max, base=2)
            log_combined = [] #pl.zeros(n_max-1)
            combined = []#pl.zeros(N)
            for i in range(n_samp):
                degrees, dk, degree_sequence = get_scale_free_degrees(N, alpha,
                                                                        avg_k, i)
                #combined[pl.asarray(degrees).astype(int)] += pl.asarray(dk)
                #log_combined += pl.histogram(degree_sequence, log_bins)[0]
                combined += [pl.histogram(degree_sequence, bins, density=True)[0]]
                log_combined += [pl.histogram(degree_sequence, log_bins)[0]]
                pl.plot(degrees, dk/N, 'C0o', alpha = 0.3)
                if i == 0:
                    pl.plot(degrees, dk/N, 'C0o', alpha = 0.3,
                        label = "Individual Networks")
            means = pl.average(combined, axis = 0)
            errors = pl.std(combined, axis = 0)/pl.sqrt(n_samp - 1)
            pl.errorbar(k, means, fmt='ko', capsize = 3, label = 'Average')
            pl.xlabel('k')
            pl.ylabel('P(k)')
            pl.legend()
            pl.savefig(plots_dir + "CheckScaleFreeAlpha{}.pdf".format(alpha))


            pl.figure(figsize=(8,6))
            pl.xscale('log')
            pl.yscale('log')
            means = pl.average(log_combined, axis = 0)
            errors = pl.std(log_combined, axis = 0)/pl.sqrt(n_samp - 1)
            pl.errorbar(log_bins[:-1], means, yerr = errors, fmt='ko', capsize = 3,
                            label = 'Binned Average')
            ev = []
            for i in range(1, n_max):
                j = int(2**i)
                k = int(2**(i+1))
                print(j, k)
                ev.append(pl.sum(N*reference[j:k]))
            pl.plot(log_bins[1:-2], ev[:-2], 
                        label = r'Power Law $\alpha = {0}$'.format(alpha))
            pl.legend()
            pl.ylabel('Average Number of Sampled Nodes')
            pl.xlabel('k')
            pl.savefig(plots_dir + "BinnedCheckScaleFreeAlpha{}.pdf".format(alpha))



    test_weighted_assignment_assortativity = 0
    if test_weighted_assignment_assortativity:
        N = 500
        original_N = N
        graph_type = "ScaleFree"
        alpha = -2.27

        avg_deg = 4
        k_min = avg_deg//2
        x_scale = 'linear'
        y_scale = 'linear'

        alphas = [2.27] #[2.01, 2.1, 2.5, 3.0, 5.0]
        if graph_type == "Random":
            alphas = ["Bubcus"]
        for alpha in alphas:
            print(alpha)
            folder = plots_folder + graph_type + "N{0}AvgDeg{1}".format(N,
                                                                    avg_deg)
            if "andom" in graph_type:
                G = nx.fast_gnp_random_graph(N, p = avg_deg/N)
                G.remove_nodes_from(list(nx.isolates(G)))
                degree_sequence = [k for node_id, k in G.degree()]
                k_max = pl.amax(degree_sequence)
                k_min = pl.amin(degree_sequence)
                dummy_degrees = pl.arange(k_min, k_max + 2)
                degrees = pl.copy(dummy_degrees)[:-1]
                dk = pl.histogram(degree_sequence, bins = dummy_degrees)[0]
            else:
                x_scale = 'log'
                y_scale = 'log'
                folder += "kMin{0}Alpha{1}".format(k_min, -alpha)
                degrees, dk, degree_sequence = get_scale_free_degrees(
                                            N, alpha, avg_deg)
                """
                details  = get_scale_free(alpha, k_min, N, avg_deg)
                if type(details) == str:
                    continue
                pk, degrees = details
                
                print("Average Degree from P(k):", pl.dot(pk, degrees))
                while 1:
                    degree_sequence, dk = random_hub_node_degrees(pk, degrees, N)
                    while abs(pl.dot(dk, degrees)/N - pl.dot(pk, degrees)) > 0.05:
                        degree_sequence, dk = random_hub_node_degrees(pk, degrees, N)
                    if check_graphical(degree_sequence):
                        break
                """



            degrees = degrees[dk != 0]
            dk = dk[dk!=0]
            print("Average Degree:", pl.dot(dk, degrees)/N)

            # set up the Big figure
            fig, axes = pl.subplots(2,3, figsize = (12,8))
            fig.subplots_adjust(wspace = 0.4)
            fraction = 0.046
            pad = 0.04

            subplot_labels = [a + ')' for a in string.ascii_lowercase]
            plot_index = 0

            ax = axes.ravel()[plot_index]
            im = ax.plot(degrees, dk, 'ko')
            ax.set_xlabel('k')
            ax.set_ylabel('D(k)')
            ax.set_yscale(y_scale)
            ax.set_xscale(x_scale)
            ax.annotate(subplot_labels[plot_index], xy = (-0.1, 1.1),
                        xycoords = 'axes fraction', fontsize = fs)

            """
            N = 714
            dk = N * pl.array([0,2,3,4,3,2,0])
            N = pl.sum(dk)
            print(N)
            degrees = pl.arange(len(dk)) + 1
            """
            edges = degrees*dk
            p_assortative = 0.0
            p_dissassortative = 0.1
            p_random = 1.0 - p_assortative - p_dissassortative

            avg_k = pl.average(degree_sequence) #pl.dot(degrees, pk)
            params[3] = str(avg_k)

            n = 5
            p_values = pl.linspace(0, 1, n)
            p_labels = ["{0:.1f}".format(p) for p in p_values]
            entropies = pl.full([n,n], pl.nan)
            assortativities = pl.full([n,n], pl.nan)
            measures = ["DeathAges", "HealthyAging", "HANorm"]
            measure_labels = ["Mean Death Age", "Mean QALY", "Mean HA"]
            health_data = [pl.copy(pl.full([n,n], pl.nan)) for i in range(3)]

            for i, p_a in enumerate(p_values):
                p_d_limit = 1-p_a
                print("P Assortative =", p_a)
                for j, p_d in enumerate(p_values[p_values <= p_d_limit]):
                    print("P Dissassortative =", p_d)
                    p_random = 1.0 - p_a - p_d
                    # build  graph
                    G, jkk = build_weighted_graph(degrees, dk,
                                                    p_a, p_d, p_random, True)
                    if type(jkk) == str:
                        continue
                    pkk = jkk/pl.sum(jkk)
                    entropies[i,j] = calculate_entropy(pkk)
                    #G = network_from_jkk(jkk, degrees, N)
                    sampled_N = G.number_of_nodes()
                    if sampled_N < 0.9*original_N:
                        print("Not enough nodes?")
                        continue
                    params[1] = str(sampled_N)
                    # get the health data from the network simulation
                    data = run_network(G, edge_list_folder, output_folder,
                                        params, param_names, measures = measures)
                    # record it
                    for k, m in enumerate(data):
                        health_data[k][i,j] = pl.average(pl.copy(m))
                    assortativities[i,j] = nx.degree_assortativity_coefficient(G)

            # Degree Assortativities
            plot_index += 1
            ax = axes.ravel()[plot_index]
            im = ax.imshow(assortativities, origin='lower')
            fig.colorbar(im, ax = ax, fraction=fraction,pad=pad)
            ax.set_xticks(pl.arange(n))
            ax.set_xticklabels(p_labels)
            ax.set_yticks(pl.arange(n))
            ax.set_yticklabels(p_labels)
            ax.set_xlabel('Dissassortative Fraction')
            ax.set_ylabel('Assortative Fraction')
            ax.set_title('Degree Assortativity Coefficient')
            ax.annotate(subplot_labels[plot_index], xy = (-0.1, 1.1),
                        xycoords = 'axes fraction', fontsize = fs)

            # Entropies
            plot_index += 1
            ax = axes.ravel()[plot_index]
            im = ax.imshow(entropies, origin='lower')
            fig.colorbar(im, ax = ax, fraction=fraction,pad=pad)
            ax.set_xticks(pl.arange(n))
            ax.set_xticklabels(p_labels)
            ax.set_yticks(pl.arange(n))
            ax.set_yticklabels(p_labels)

            ax.set_xlabel('Dissassortative Fraction')
            ax.set_ylabel('Assortative Fraction')
            ax.set_title('Pkk Entropy')
            ax.annotate(subplot_labels[plot_index], xy = (-0.1, 1.1),
                        xycoords = 'axes fraction', fontsize = fs)

           
            for i, average_data in enumerate(health_data):
                plot_index += 1
                ax = axes.ravel()[plot_index]
                im = ax.imshow(average_data, origin='lower')
                fig.colorbar(im, ax = ax, fraction=fraction,pad=pad)
                ax.set_xticks(pl.arange(n))
                ax.set_xticklabels(p_labels)
                ax.set_yticks(pl.arange(n))
                ax.set_yticklabels(p_labels)
                ax.set_xlabel('Dissassortative Fraction')
                ax.set_ylabel('Assortative Fraction')
                ax.set_title(measure_labels[i])
                ax.annotate(subplot_labels[plot_index], xy = (-0.1, 1.1),
                            xycoords = 'axes fraction', fontsize = fs)

            asp = 1/2 #pl.diff(axes[0,0].get_xlim())[0] / pl.diff(axes[0,0].get_ylim())[0]
            axes[0,0].set_aspect('auto', 'box')

            fig.savefig(folder + "AssortativitySweep.pdf".format(
                                                original_N, graph_type))

        # Plot in 1 D the slices vs P_a and P_d with other 0
        plot_1d_slices = 0
        if plot_1d_slices:
            d_entropies = entropies[0,:]
            d_assortativities = assortativities[0,:]
            pl.figure(figsize = (8,6))
            pl.plot(p_values, d_entropies, label = 'Dissassortative Entropy')
            pl.plot(p_values, d_assortativities,
                label = 'Dissassortative\nassortativity')
            a_entropies = entropies[:,0]
            a_assortativities = assortativities[:,0]
            pl.plot(p_values, a_entropies, label = 'Assortative Entropy')
            pl.plot(p_values, a_assortativities,
                label = 'Assortative\nAssortativity')
            pl.xlabel('(Diss)Assortative Fraction')

            pl.axhline(1, c = 'k', ls = ':')
            pl.axhline(-1, c = 'k', ls = ':')
            pl.legend()
            pl.savefig(plots_folder + "WeightedNetworks1DSlices.pdf")

    show_weighted_pkk = 0
    if show_weighted_pkk:
        N = 500
        graph_type = "ScaleFree"
        pl.seed(1)
        alpha = 2.27
        avg_deg = 4.0
        if "andom" in graph_type:
            G = nx.fast_gnp_random_graph(N, p = 4/N)
        else:
            G = nx.barabasi_albert_graph(N, 2)
            degrees, dk, degree_sequence = get_scale_free_degrees(
                                        N, alpha, avg_deg)

        #degree_sequence = [k for node_id, k in G.degree()]
        k_max = pl.amax(degree_sequence)
        #degrees = pl.arange(1, k_max + 2)
        #dk = pl.histogram(degree_sequence, bins = degrees)[0]
        #degrees = degrees[:-1]
        degrees = degrees[dk != 0]
        dk = dk[dk!=0]

        edges = degrees*dk
        p_assortative = 0.0
        p_dissassortative = 0.1
        p_random = 1.0 - p_assortative - p_dissassortative

        n_ticks = 10
        if pl.amax(degrees) > 20:
            n_ticks = 6
        show_every = int(pl.round_(len(degrees)/n_ticks)) + 1

        n = 3
        p_values = pl.linspace(0, 1, n)
        
        fig, ax = pl.subplots(n,n, figsize = (8, 8))
        
        ax[0,2].plot(degrees, dk, 'ko')
        ax[0,2].set_ylim(-0.05*pl.amax(dk), pl.amax(dk)*1.05)
        ax[0,2].set_xlim(0, pl.amax(degrees)+2)

        if graph_type != "Random":
            ax[0,2].set_ylim(0.4, pl.amax(dk)*2)
            ax[0,2].set_xlim(0.8, pl.amax(degrees)*1.5)
            ax[0,2].set_yscale('log')
            ax[0,2].set_xscale('log')
        ax[0,2].set_title('Degree Distribution')
        ax[0,2].set_xlabel('k')
        ax[0,2].set_ylabel('D(k)')


        pl.subplots_adjust(hspace=0.05, wspace=0.05)
        
        for i, p_a in enumerate(p_values):
            p_d_limit = 1-p_a
            print("P Assortative =", p_a)
            for j, p_d in enumerate(p_values):#[p_values <= p_d_limit]):
                use_x_ticks = (n - i - 1) == 2
                use_y_ticks = j == 0
                axis = ax[n-i-1,j]
                print("P Dissassortative =", p_d)
                if p_d > p_d_limit:
                    axis.imshow(pl.full([len(degrees), len(degrees)], pl.nan),
                        origin = 'lower')
                    continue
                p_random = 1.0 - p_a - p_d
                # build  graph
                G, jkk = build_weighted_graph(degrees, dk,
                                                p_a, p_d, p_random, True)
                if type(jkk) == str:
                    print("Failure To Build Graph")
                    axis.imshow(pl.full([len(degrees), len(degrees)], pl.nan),
                        origin = 'lower')
                    continue
                pkk = jkk/pl.sum(jkk)
                axis.imshow(pkk, origin='lower')
                axis.annotate(
                    'P-A:{0}\n'
                    'P-D:{1}'.format(p_a, p_d, p_random),
                        xy = (0.5, 0.98), xycoords = 'axes fraction',
                        ha='center', va='top',color = 'w')
                if use_x_ticks:
                    axis.set_xticks(pl.arange(len(degrees))[::show_every])
                    axis.set_xticklabels(degrees[::show_every].astype(int))
                else:
                    axis.set_xticks([])
                if use_y_ticks:
                    axis.set_yticks(pl.arange(len(degrees))[::show_every])
                    axis.set_yticklabels(degrees[::show_every].astype(int))
                else:
                    axis.set_yticks([])


        ax[0,2].set_aspect('auto', 'box')

        ax[0,1].axis('off')
        ax[1,2].axis('off')
        pl.savefig(plots_folder + "{0}N{1}GraphPkkExamples.pdf".format(
                                                    graph_type, N))
                 
    show_network = 0
    if show_network:
        N = 500
        graph_type = "ScaleFree"
        pl.seed(1)
        alpha = 2.27
        avg_deg = 4.0
        if "andom" in graph_type:
            G = nx.fast_gnp_random_graph(N, p = 4/N)
        else:
            G = nx.barabasi_albert_graph(N, 2)
            degrees, dk, degree_sequence = get_scale_free_degrees(
                                        N, alpha, avg_deg)

        #degree_sequence = [k for node_id, k in G.degree()]
        k_max = pl.amax(degree_sequence)
        #degrees = pl.arange(1, k_max + 2)
        #dk = pl.histogram(degree_sequence, bins = degrees)[0]
        #degrees = degrees[:-1]
        degrees = degrees[dk != 0]
        dk = dk[dk!=0]

        edges = degrees*dk
        p_a = 0.0
        p_d = 0.0
        p_random = 1.0 - p_a -p_d 


        G, jkk = build_weighted_graph(degrees, dk,
                                        p_a, p_d, p_random, True)
        pl.figure(figsize = (8,6))
        pos = nx.spring_layout(G)
        nx.draw(G, pos, node_size = 10, edge_color = 'C7')
        pl.savefig(plots_folder + "Pa{0}Pd{1}{2}Visualization.pdf".format(
                            p_a, p_d, graph_type))


    test_build_network = 0
    if test_build_network:
        N = 1000
        alpha = -2.25
        k_min = 2
        avg_deg = 4.0
        p_assortative = 0.9
        p_dissassortative = 0.0
        p_random = 1 - p_assortative - p_dissassortative

        pk, degrees = get_scale_free(alpha, k_min, N, avg_deg)
        print("Average Degree from P(k):", pl.dot(pk, degrees))
        while 1:
            degree_sequence, dk = random_hub_node_degrees(pk, degrees, N)
            while abs(pl.dot(dk, degrees)/N - pl.dot(pk, degrees)) > 0.05:
                degree_sequence, dk = random_hub_node_degrees(pk, degrees, N)
            if check_graphical(degree_sequence):
                break

        degrees = degrees[dk != 0]
        dk = dk[dk!=0]
        print("Sampled Average Degree:", pl.dot(dk, degrees)/N)
        pl.figure(figsize = (8,6))
        pl.plot(degrees, dk, 'ko')
        pl.xscale('log')
        pl.yscale('log')
        pl.ylabel('P(k)')
        pl.xlabel('k')

        # build the graph
        G, jkk = build_weighted_graph(degrees, dk,
                            p_assortative, p_dissassortative, p_random, True)
        if type(G) != str:
            pl.figure(figsize = (8,6))
            pl.imshow(jkk, origin = 'lower')
            pl.colorbar()

    draw_graph = 0
    test_run = 0
    if test_run:
        pk = pl.array([1, 2, 3, 4, 3, 2, 1])
        pk = pk / pl.sum(pk)
        degrees = pl.arange(len(pk)) + 1
        avg_k = pl.dot(degrees, pk)
        params[3] = str(avg_k)
        N = 4000

        G, pkk = build_network(pk, degrees, N, return_pkk = True)
        N = G.number_of_nodes()
        params[1] = str(N)

        if G == "hell":
            print("Damn")
        else:
            if draw_graph:
                pos = nx.spring_layout(G)
                nx.draw(G, pos, node_size = 6)
            print("Assortativity:", nx.degree_assortativity_coefficient(G))
        death_ages = run_network(G, edge_list_folder, output_folder, params,
                            param_names, measure = "DeathAge")
        pl.figure()
        pl.hist(death_ages/scale)
        pl.figure()
        pl.imshow(pkk, origin = 'lower')
        pl.colorbar()
 
    test_optimization = 0
    if test_optimization:
        k_min = 2
        dk = [0,2,3,4,3,2]
        degrees = [i+1 for i in range(len(dk))]
        bonus_degrees = [i+1 for i in range(len(dk), 100)]
        bonus_degrees = [10, 15, 25, 50, 100, 250, 500]#, 1000, 2000]
        k_max = bonus_degrees[-1]
        dk += [0 for i in range(len(bonus_degrees))]
        degrees += bonus_degrees
        degrees = pl.asarray(degrees)
        dk = pl.asarray(dk)
        dk = dk[degrees >= k_min]
        degrees = degrees[degrees >= k_min]
        pk = pl.asarray(dk)/pl.sum(dk)
        #print(pk, degrees)
        """
        pk = pl.array([1, 2, 3, 4, 3, 2, 1, 0, 0, 0, 0])
        pk = pk / pl.sum(pk)
        degrees = pl.arange(len(pk)) + 1
        degrees[-3] = 15
        degrees[-2] = 25
        degrees[-1] = 50
        """
        avg_k = pl.dot(degrees, pk)
        params[3] = str(avg_k)
        N = 1000
        original_N = N

        static_assortativity = False
        if static_assortativity:
            p_assortative = 0.0
            p_dissassortative = 0.0
            p_random = 1.0 - p_assortative - p_dissassortative

        pl.seed(1)
        
        best_pk = pl.copy(pk)
        best_entropy = -0.01
        best_death_age = -0.01

        best_death_ages = []
        best_entropies = []
        best_pkk_entropies = []
        best_merits = []
        best_pks = []
        best_sampled_dks = []
        best_probs = []
        iterations = []
        best_G = 1
        best_pkk = 1

        entropy_weighting = 1.0
        #run_number = 0
        #run_folder = "Run{}/".format(run_number)
        #edge_list_folder += run_folder
        #output_folder += run_folder
        #params[-4] = run_folder + "SimplyConnected"
        #params[4] =  output_folder[:-1]
        if entropy_weighting == 1:
            p_assortative = 0
            p_dissassortative = 0
            p_random = 1
            mean_death_age = best_death_age
        best_merit = entropy_weighting*best_entropy + (
            1-entropy_weighting)*best_death_age


        fraction = 0.1
        max_changes = 100

        start_time = time.time()
        elapsed_time = (time.time() - start_time)/60 # in minutes
        run_time = 60 # In minutes0w
        measure_time = run_time / 10
        progress = 0
        num_changes = max_changes
        i = 1
        fractions = []
        while elapsed_time < run_time: 
            i += 1
            N = original_N
            elapsed_time = (time.time() - start_time)/60 # in minutes
            if elapsed_time > measure_time:
                measure_time += run_time/10
                progress += 10
                print("{0}% Done. Iteration {1}. Entropy {2}".format(progress,
                        i, best_entropies[-1]))
            if i % 2 == 1:
                num_changes -= 1
            if num_changes < 1:
                num_changes = 1#pl.randint(10) + 1
            pk = pl.copy(best_pk)
            successful_changes = 0
            fraction = pl.random()
            while successful_changes < num_changes:
                n1, n2 = pl.choice(pl.arange(len(pk)), size = 2,
                    replace = False)
                pk, success = change_pk(pk, degrees, n1, n2, fraction, avg_k) 
                #print(success, successful_changes)
                if success:
                    successful_changes += 1
            
            if entropy_weighting == 1:
                #mean_death_age = 0.0
                #entropy = calculate_entropy(pk)
                dk = random_hub_node_degrees(pk, degrees, N)[1]
                dk = pl.asarray(dk)
                build_start = time.time()
                #print("building network, iteration", i)
                G, jkk = build_weighted_graph(degrees, dk,
                    p_assortative, p_dissassortative, p_random, True)
                build_time = time.time() - build_start
                if type(G) == str:
                    continue
                pkk = jkk/pl.sum(jkk)
                entropy = calculate_entropy(pkk)

            else:
                dk = random_hub_node_degrees(pk, degrees, N)[1]
                dk = pl.asarray(dk)
                build_start = time.time()
                print("building network, iteration", i)
                if not static_assortativity:
                    p_assortative = pl.random()
                    p_dissassortative = pl.random()*(1-p_assortative)
                    p_random = 1.0 - p_assortative - p_dissassortative

                G, jkk = build_weighted_graph(degrees, dk,
                    p_assortative, p_dissassortative, p_random, True)
                build_time = time.time() - build_start
                if type(G) == str:
                    continue
                pkk = jkk/pl.sum(jkk)


                sampled_N = G.number_of_nodes()
                if sampled_N < 0.95*original_N:
                    continue
                params[1] = str(sampled_N)


                entropy = calculate_entropy(pkk)
                simulation_start = time.time()
                death_ages = run_network(G, edge_list_folder, output_folder,
                                    params, param_names)
                simulation_time = time.time() - simulation_start
                mean_death_age = pl.average(death_ages)

                print("""Network Build Time: {0:.3f},
                        Simulation Time: {1:.3f}""".format(
                                                build_time, simulation_time))
            merit = entropy_weighting*entropy + (
                1-entropy_weighting)*mean_death_age
            r = pl.random()

            #print("Sum of Pk: ", pl.sum(pk), ", Avg k: ", pl.dot(pk, degrees))

            #if merit > best_merit:
            if pl.exp((merit - best_merit)*pl.sqrt(i)*2) > r:
                best_probs.append([p_assortative, p_dissassortative,
                                                    p_random])
                best_pk = pl.copy(pk)
                best_pks.append(pl.copy(pk))
                best_entropies.append(calculate_entropy(best_pk))
                if entropy_weighting != 1:
                    best_G = G
                    degree_sequence = [k for node_id, k in G.degree()]
                    sampled_dk = pl.histogram(degree_sequence,
                            list(degrees) + [pl.inf])[0]
                    best_sampled_dks.append(sampled_dk)
                    best_death_ages.append(mean_death_age)
                best_sampled_dks.append(dk)
                    
                best_pkk = pkk
                best_pkk_entropies.append(calculate_entropy(pkk))
                best_merit = merit
                best_merits.append(merit)
                fractions.append(fraction)
                iterations.append(i)

        if static_assortativity:
            simulation_details = (
                'N{0}KMin{1}KMax{2}Lambda{3}Pa{4}Pd{5}Pr{6}'.format(
                    N, k_min, k_max, entropy_weighting, p_assortative,
                    p_dissassortative, p_random))
        else:
            simulation_details = (
                'N{0}KMin{1}KMax{2}Lambda{3}RandomAssortativities'.format(
                    N, k_min, k_max, entropy_weighting))

    
        if entropy_weighting != 1:
            pl.figure(figsize = (8,6))
            pl.plot(iterations, best_pkk_entropies, 'C0o', label = "Pkk entropy")
            pl.plot(iterations, best_entropies, 'C1o', label = "Pk Entropy")
            pl.ylabel('Entropy (Pkk used in Merit)')
            pl.xlabel('Iteration')
            pl.xscale('log')
            pl.legend()
            pl.savefig(plots_folder + simulation_details + 'Entropies.pdf')

            pl.figure(figsize = (8,6))
            pl.plot(iterations, best_death_ages, 'C0o')
            pl.ylabel('Mean Death Age')
            pl.xlabel('Iteration')
            pl.xscale('log')
            pl.savefig(plots_folder + simulation_details + "MeanDeathAges.pdf")

            pl.figure(figsize = (8,6))
            prob_labels = ['Assortative', 'Dissassortative', 'Random']
            for i in range(3):
                probs = pl.asarray(best_probs)[:,i]
                pl.plot(iterations, probs, marker = 'o', ls = 'none',
                        label = prob_labels[i])
                pl.ylabel('Assortativity Fraction')
                pl.xlabel('Iteration')
                pl.xscale('log')
                pl.legend()
                pl.savefig(plots_folder + simulation_details +
                        "AssortativityFractions.pdf")

            if len(best_merits) > 1:
                best_merits[0] = best_merits[1]
            colours = get_colours_from_cmap(pl.arange(len(best_merits)))
            pl.figure(figsize = (8,6))
            #for i, colour in enumerate(colours):
            #    pl.plot(degrees, best_sampled_dks[i], c = colour, ls = 'none',
            #        marker = 'o')
            pl.plot(degrees, best_sampled_dks[-1], 'C0o')
            pl.yscale('log')
            pl.xscale('log')
            pl.ylabel('Sampled D(k)')
            pl.xlabel('k')
            pl.savefig(plots_folder + simulation_details + "DegreeCounts.pdf")

            pl.figure(figsize = (8,6))
            G, pkk = best_G, best_pkk
            pos = nx.spring_layout(G)
            nx.draw(G, pos, node_size = 8, edge_color = 'C7')
            pl.savefig(plots_folder + simulation_details + "FinalNetwork.pdf")

            pl.figure(figsize = (8,6))
            pl.imshow(pkk, origin = 'lower')
            pl.xticks(pl.arange(len(degrees)), labels = degrees)
            pl.yticks(pl.arange(len(degrees)), labels = degrees)
            pl.colorbar()
            pl.savefig(plots_folder + simulation_details + "FinalPkk.pdf")
        else:
            pl.figure(figsize=(8,6))
            pl.plot(degrees, best_sampled_dks[-1], 'C0o')
            pl.xlabel('k')
            pl.ylabel('Sampled D(k)')
            pl.yscale('log')
            pl.xscale('log')
            pl.savefig(plots_folder + simulation_details + "DegreeCounts.pdf")

            pl.figure(figsize = (8,6))
            pl.plot(iterations, best_pkk_entropies, 'C0o', label = "Pkk entropy")
            pl.plot(iterations, best_entropies, 'C1o', label = "Pk Entropy")
            pl.ylabel('Entropy (Pkk used in Merit)')
            pl.xlabel('Iteration')
            pl.xscale('log')
            pl.legend()
            pl.savefig(plots_folder + simulation_details + 'Entropies.pdf')


            pl.figure(figsize = (8,6))
            pl.imshow(best_pkk, origin = 'lower')
            pl.xticks(pl.arange(len(degrees)), labels = degrees)
            pl.yticks(pl.arange(len(degrees)), labels = degrees)
            pl.colorbar()
            pl.savefig(plots_folder + simulation_details + "FinalPkk.pdf")


    test_motifs = 0
    if test_motifs:
        
        motifs = ["star", "ball", "chain", "parents"]
        N_values = pl.logspace(1, 8, 8, base = 2).astype(int)  
        measures = ["DeathAges", "HealthyAging", "HANorm"]
        measure_labels = ["Mean Death Age", "Mean QALY", "Mean HA"]

        params[3] = "1"
        
        means = []
        errors = []
        for motif in motifs:
            motif_means = []
            motif_errors = []
            for N in N_values:
                run_errors= []
                run_means = []
                G = build_motif_graph(motif, N)
                data = run_network(G, edge_list_folder, output_folder,
                                    params, param_names, measures = measures)
                num = len(data[0])
                for measure in data:
                    run_means.append(pl.average(measure))
                    run_errors.append(pl.std(measure)/pl.sqrt(num-1))
                motif_means.append(run_means)
                motif_errors.append(run_errors)
            means.append(motif_means)
            errors.append(motif_errors)

        colours = ["C0", "C1", "C2", "C3"]
        markers = ["*", "o", "d", "s"]
        lines = ["-", ":", "--", "-."]
        for k, measure_label in enumerate(measure_labels):
            pl.figure(figsize = (8,6))
            for i, motif in enumerate(motifs):
                run_means = []
                run_errors = []
                for j, N in enumerate(N_values):
                    run_mean = means[i][j][k]
                    run_error = errors[i][j][k]
                    if "HA" in measure_label:
                        run_mean *= 100
                        run_error *= 100
                    run_means.append(run_mean)
                    run_errors.append(run_error)
                pl.errorbar(N_values, run_means, yerr = run_errors,
                            fmt = "{0}{1}{2}".format(colours[i], markers[i],
                                                    lines[k]),
                            label = motif, capsize = 3)
            pl.xlabel('N')
            pl.ylabel(measure_label)
            pl.legend()
            pl.xscale('log')
            pl.savefig(plots_folder + "Motifs{}.pdf".format(measure_label))

    visualize_motifs = 0
    if visualize_motifs:
        motifs = ["star", "ball", "chain", "parents"]
        N = 8
        fig, axes = pl.subplots(2,2, figsize = (8,6))
        subplot_labels = [a + ')' for a in string.ascii_lowercase]
        for i, motif in enumerate(motifs):
            p = axes[int(i > 1), i%2]
            G = build_motif_graph(motif, N)
            pos = nx.spring_layout(G)
            nx.draw(G, pos, node_size = 10, edge_color = 'C7', ax = p)
            p.annotate(subplot_labels[i], xy = (0.05, 0.9),
                xycoords = 'axes fraction', fontsize = fs)

            p.set_title(motif + " motif")

        pl.savefig(plots_folder + "MotifExamples.pdf")


    check_death_ball = 0
    if check_death_ball:
        measures = ["DeathAges", "HealthyAging", "HANorm"]
        measure_labels = ["Mean Death Age", "Mean QALY", "Mean HA"]



        #output_folder = "BallAndChainData/"
        #initial_distribution = "SimplyConnected"
        #edge_list_folder = "BallAndChainEdges/"
        number = "10000"
        avg_k = "4"
        param_names = ['number', 'N', 'numChanges', 'avgDeg', 'folder',
            'singleseed', 'evoCondition', 'InitialDistribution',
            'lambda', 'beta', 'power']
        params = [number, N, numChanges, avg_k, output_folder, seed,
            run_hours, evo_condition, initial_distribution, Lambda,
            beta, power]


        #optimal_data_dir = '../CoolingCode/OptimalNetworkData/'
        #files = os.listdir(optimal_data_dir)
        #file_name_format = ("DeathAgesNumber100000N{0}NumChanges0AvgDeg4" +
        #    "Seed1MeritnoneInitialDistributionN{0}Nc{1}optimalLambda0.0" + 
        #    "Beta100.0Power1.0")
        Nc_values = [1, 5, 10, 20]
        N_values = [50, 100, 200, 500, 1000, 5000, 10000]

        errors = []
        means = []
        for Nc in Nc_values:
            N_means = []
            N_errors = []
            for N in N_values:
                
                avg_k = float(avg_k)
                run_means = []
                run_errors = []
                G, mort_nodes, fi_nodes  = optimal_network(N, Nc, avg_k,
                                                            Nc)
                params[1] = str(N)
                if N == 50:
                    pl.figure()
                    pos = nx.spring_layout(G)
                    nx.draw(G, pos, node_size = 10, edge_color = 'C7')
                data = run_network(G, edge_list_folder, output_folder,
                                    params, param_names, measures = measures)
                num = len(data[0])
                
                for measure in data:
                    run_means.append(pl.average(measure))
                    run_errors.append(pl.std(measure)/pl.sqrt(num-1))
                print("N:{0},Nc:{1},<td>={2:.3f}".format(N, Nc, run_means[0]))
                N_means.append(run_means)
                N_errors.append(run_errors)
            means.append(N_means)
            errors.append(N_errors)

        colours = ["C0", "C1", "C2", "C3"]
        markers = ["*", "o", "d", "s"]
        lines = ["-", ":", "--", "-."]
        for k, measure_label in enumerate(measure_labels):
            pl.figure(figsize = (8,6))
            for i, Nc in enumerate(Nc_values):
                run_means = []
                run_errors = []
                for j, N in enumerate(N_values):
                    run_mean = means[i][j][k]
                    run_error = errors[i][j][k]
                    run_means.append(run_mean)
                    run_errors.append(run_error)
                pl.errorbar(N_values, run_means, yerr = run_errors,
                            fmt = "{0}{1}{2}".format(colours[i], markers[i],
                                                    lines[i]),
                            label = "{} Hubs".format(Nc), capsize = 3)
            pl.xlabel('N')
            pl.ylabel(measure_label)
            pl.legend()
            pl.xscale('log')
            pl.savefig(plots_folder + "BallAndChain{}.pdf".format(
                                                            measure_label))

    variational_optimization = 0
    if variational_optimization:
        alpha = 2.27
        best_alpha = alpha
        width = 0.03
        limits = [2.00001, 5.0]
        avg_deg = 4.0
        N = 1000
        seed = 1
        degrees, dk, degree_sequence = get_scale_free_degrees(
                                    N, alpha, avg_deg, seed)
        k_min = 2
        k_max = "UNDEF"

        avg_k = pl.average(degree_sequence)
        params[3] = str(avg_deg)
        original_N = N
        output_folder = "testTemp/"
        edge_list_folder = "testTemp/"
        params[4] = output_folder

        static_assortativity = False
        if static_assortativity:
            p_assortative = 0.0
            p_dissassortative = 0.0
            p_random = 1.0 - p_assortative - p_dissassortative

        pl.seed(1)
        
        best_entropy = -0.01
        best_death_age = -0.01

        best_death_ages = []
        best_entropies = []
        best_pkk_entropies = []
        best_merits = []
        best_pks = []
        best_sampled_dks = []
        best_probs = []
        iterations = []
        best_degrees = []
        best_params = []
        best_G = 1
        best_pkk = 1

        entropy_weighting = 0.0
        if entropy_weighting == 1:
            static_assortativity = True
            p_assortative = 0
            p_dissassortative = 0
            p_random = 1
            mean_death_age = best_death_age
        best_merit = entropy_weighting*best_entropy + (
            1-entropy_weighting)*best_death_age


        fraction = 0.1
        max_changes = 100

        start_time = time.time()
        elapsed_time = (time.time() - start_time)/60 # in minutes
        run_time = 1 # In minutes0w
        measure_time = run_time / 10
        progress = 0
        num_changes = max_changes
        i = 1
        fractions = []
        while elapsed_time < run_time: 
            i += 1
            N = original_N
            elapsed_time = (time.time() - start_time)/60 # in minutes
            if elapsed_time > measure_time:
                measure_time += run_time/10
                progress += 10
                print("{0}% Done. Iteration {1}. Entropy {2}".format(progress,
                        i, best_entropies[-1]))
            alpha = pl.copy(best_alpha)
            alpha = gaussian_parameter_change(alpha, width, limits)
           
            if entropy_weighting == 1:
                degrees, dk, degree_sequence = get_scale_free_degrees(
                                            N, alpha, avg_deg, seed)
                build_start = time.time()
                #print("building network, iteration", i)
                G, jkk = build_weighted_graph(degrees, dk,
                    p_assortative, p_dissassortative, p_random, True)
                build_time = time.time() - build_start
                if type(G) == str:
                    continue
                pkk = jkk/pl.sum(jkk)
                entropy = calculate_entropy(pkk)

            else:
                degrees, dk, degree_sequence = get_scale_free_degrees(
                                            N, alpha, avg_deg, seed)
                build_start = time.time()
                print("building network, iteration", i)
                if not static_assortativity:
                    p_assortative = pl.random()
                    p_dissassortative = pl.random()*(1-p_assortative)
                    p_random = 1.0 - p_assortative - p_dissassortative

                G, jkk = build_weighted_graph(degrees, dk,
                    p_assortative, p_dissassortative, p_random, True)
                build_time = time.time() - build_start
                if type(G) == str:
                    continue
                pkk = jkk/pl.sum(jkk)


                sampled_N = G.number_of_nodes()
                if sampled_N < 0.95*original_N:
                    continue
                params[1] = str(sampled_N)


                entropy = calculate_entropy(pkk)
                simulation_start = time.time()
                death_ages = run_network(G, edge_list_folder, output_folder,
                                    params, param_names)
                simulation_time = time.time() - simulation_start
                mean_death_age = pl.average(death_ages)

                print("""Network Build Time: {0:.3f},
                        Simulation Time: {1:.3f}""".format(
                                                build_time, simulation_time))
            merit = entropy_weighting*entropy + (
                1-entropy_weighting)*mean_death_age
            r = pl.random()

            #print("Sum of Pk: ", pl.sum(pk), ", Avg k: ", pl.dot(pk, degrees))

            #if merit > best_merit:
            if pl.exp((merit - best_merit)*pl.sqrt(i)*2) > r:
                best_probs.append([p_assortative, p_dissassortative,
                                                    p_random])
                best_alpha = alpha
                if entropy_weighting != 1:
                    best_G = G
                    sampled_dk = dk 
                    best_sampled_dks.append(sampled_dk)
                    best_death_ages.append(mean_death_age)
                best_pk = pl.copy(dk/N)
                best_pks.append(best_pk)
                best_entropies.append(calculate_entropy(best_pk))
                best_params.append(alpha)

                best_sampled_dks.append(dk)
                best_degrees.append(degrees)
                    
                best_pkk = pkk
                best_pkk_entropies.append(calculate_entropy(pkk))
                best_merit = merit
                best_merits.append(merit)
                fractions.append(fraction)
                iterations.append(i)

        if static_assortativity:
            simulation_details = (
                'N{0}KMin{1}KMax{2}Lambda{3}Pa{4}Pd{5}Pr{6}'.format(
                    N, k_min, k_max, entropy_weighting, p_assortative,
                    p_dissassortative, p_random))
        else:
            simulation_details = (
                'VariationalN{0}Lambda{3}RandomAssortativities'.format(
                    N, k_min, k_max, entropy_weighting))

        degrees = best_degrees[-1] 
        if entropy_weighting != 1:
            pl.figure(figsize = (8,6))
            pl.plot(iterations, best_pkk_entropies, 'C0o', label = "Pkk entropy")
            pl.plot(iterations, best_entropies, 'C1o', label = "Pk Entropy")
            pl.ylabel('Entropy (Pkk used in Merit)')
            pl.xlabel('Iteration')
            pl.xscale('log')
            pl.legend()
            pl.savefig(plots_folder + simulation_details + 'Entropies.pdf')

            pl.figure(figsize = (8,6))
            pl.plot(iterations, best_params, 'C0o')
            pl.ylabel('Scale Free Exponent')
            pl.xlabel('Iteration')
            pl.xscale('log')
            pl.savefig(plots_folder + simulation_details + "BestAlphas.pdf")


            pl.figure(figsize = (8,6))
            pl.plot(iterations, best_death_ages, 'C0o')
            pl.ylabel('Mean Death Age')
            pl.xlabel('Iteration')
            pl.xscale('log')
            pl.savefig(plots_folder + simulation_details + "MeanDeathAges.pdf")

            pl.figure(figsize = (8,6))
            prob_labels = ['Assortative', 'Dissassortative', 'Random']
            for i in range(3):
                probs = pl.asarray(best_probs)[:,i]
                pl.plot(iterations, probs, marker = 'o', ls = 'none',
                        label = prob_labels[i])
                pl.ylabel('Assortativity Fraction')
                pl.xlabel('Iteration')
                pl.xscale('log')
                pl.legend()
                pl.savefig(plots_folder + simulation_details +
                        "AssortativityFractions.pdf")

            if len(best_merits) > 1:
                best_merits[0] = best_merits[1]
            colours = get_colours_from_cmap(pl.arange(len(best_merits)))
            pl.figure(figsize = (8,6))
            #for i, colour in enumerate(colours):
            #    pl.plot(degrees, best_sampled_dks[i], c = colour, ls = 'none',
            #        marker = 'o')
            pl.plot(best_degrees[-1], best_sampled_dks[-1], 'C0o')
            pl.yscale('log')
            pl.xscale('log')
            pl.ylabel('Sampled D(k)')
            pl.xlabel('k')
            pl.savefig(plots_folder + simulation_details + "DegreeCounts.pdf")

            pl.figure(figsize = (8,6))
            G, pkk = best_G, best_pkk
            pos = nx.spring_layout(G)
            nx.draw(G, pos, node_size = 8, edge_color = 'C7')
            pl.savefig(plots_folder + simulation_details + "FinalNetwork.pdf")

            pl.figure(figsize = (8,6))
            pl.imshow(pkk, origin = 'lower')
            pl.xticks(pl.arange(len(degrees)), labels = degrees)
            pl.yticks(pl.arange(len(degrees)), labels = degrees)
            pl.colorbar()
            pl.savefig(plots_folder + simulation_details + "FinalPkk.pdf")
        else:
            pl.figure(figsize=(8,6))
            pl.plot(best_degrees[-1], best_sampled_dks[-1], 'C0o')
            pl.xlabel('k')
            pl.ylabel('Sampled D(k)')
            pl.yscale('log')
            pl.xscale('log')
            pl.savefig(plots_folder + simulation_details + "DegreeCounts.pdf")

            pl.figure(figsize = (8,6))
            pl.plot(iterations, best_pkk_entropies, 'C0o', label = "Pkk entropy")
            pl.plot(iterations, best_entropies, 'C1o', label = "Pk Entropy")
            pl.ylabel('Entropy (Pkk used in Merit)')
            pl.xlabel('Iteration')
            pl.xscale('log')
            pl.legend()
            pl.savefig(plots_folder + simulation_details + 'Entropies.pdf')


            pl.figure(figsize = (8,6))
            pl.imshow(best_pkk, origin = 'lower')
            pl.xticks(pl.arange(len(degrees)), labels = degrees)
            pl.yticks(pl.arange(len(degrees)), labels = degrees)
            pl.colorbar()
            pl.savefig(plots_folder + simulation_details + "FinalPkk.pdf")

    scale_free_assortativity_check = 1
    if scale_free_assortativity_check:
        pA = 0.0
        N = 1000
        alpha = 2.3
        avg_deg = 4
        seed = 1
        p_d = 0.0
        temp_folder = "EdgeLists/Tests/"

        check_pA_values = 0
        if check_pA_values:
            pAs = pl.linspace(0, 0.999, 11)
            colours = get_colours_from_cmap(pAs)
            degree_assortativities = []
            remake_rs = []
            Ns = pl.zeros_like(pAs)
            remade_Ns = pl.zeros_like(pAs)
            pl.figure(figsize = (8,6))
            for i, pA in enumerate(pAs):
                degrees, dk, degree_sequence, G = get_scale_free_degrees(N,
                                            temp_folder, alpha, avg_deg, seed, pA,
                                            return_graph=True)
                pl.plot(degrees, dk, c = colours[i], ls = "none", marker = "o")
                remade_G = build_weighted_graph(degrees, dk, pA, p_d, 1-pA)
                remade_Ns[i] = dk[0]#remade_G.number_of_nodes()
                r = nx.degree_assortativity_coefficient(G)
                remade_r =nx.degree_assortativity_coefficient(remade_G)
                degree_assortativities.append(r)
                Ns[i] = G.number_of_nodes()
                remake_rs.append(remade_r)
            pl.xscale('log')
            pl.yscale('log')
            pl.figure(figsize=(8,6))
            pl.plot(pAs, degree_assortativities)
            pl.plot(pAs, remake_rs)
            pl.axhline(0, c = "C7")
            pl.plot(pAs, remade_Ns/N, "C0:")
            pl.plot(pAs, Ns/N, "C1:")
        
        check_pkks = 0
        if check_pkks:
            pA = 0.99
            N = 200
            print("Original Network:")
            degrees, dk, degree_sequence, G = get_scale_free_degrees(N,
                                    temp_folder, alpha, avg_deg, seed, pA,
                                    return_graph=True)
            print([v for k, v in G.degree])
            print("Dot Product", pl.dot(dk, degrees))
            print("Number of Self Loops", nx.number_of_selfloops(G))
            show_log_jkk(degrees, dk, degree_sequence, G, "After The Fact")

            print("Remade Network:")
            remade_G, remade_pkk = build_weighted_graph(degrees, dk, pA,
                                        p_d, 1-pA, True)
            pl.figure()
            #remade_pkk[remade_pkk != 0] = pl.log(remade_pkk[remade_pkk != 0])
            remade_pkk[remade_pkk == 0] = pl.nan
            pl.imshow(remade_pkk, origin = 'lower')
            show_log_jkk(degrees, dk, degree_sequence, remade_G,
                "Garrett's Shitty Method")
            print([v for k, v in remade_G.degree])

        check_rewiring = 1
        if check_rewiring:
            pA = 0.0
            N = 200
            swaps_per_edge = 1
            print("Original Network:")
            degrees, dk, degree_sequence, G = get_scale_free_degrees(N,
                                    temp_folder, alpha, avg_deg, seed, pA,
                                    return_graph=True)
            show_log_jkk(degrees, dk, degree_sequence, G, "Original")

            
            old_r = nx.degree_assortativity_coefficient(G)
            print("Old R:", old_r)

            new_edges, r_differences = shuffle_edges(G, swaps_per_edge)
            new_G = nx.from_edgelist(new_edges)
            new_r = nx.degree_assortativity_coefficient(new_G)
            num_components = [len(c) for c in nx.connected_components(new_G)] 
            print("New R:", new_r, num_components,
                "Connected components.")

            show_log_jkk(degrees, dk, degree_sequence, new_G, "Rewired in Post")

            p_a = 1.0
            remade_G = build_weighted_graph(degrees, dk, p_a, p_d, 1-p_a)
            remade_r = nx.degree_assortativity_coefficient(remade_G)
            print("Garrett's R:", remade_r)
            show_log_jkk(degrees, dk, degree_sequence, remade_G, "Garrett's")

            new_Garrett_edges, r_differences = shuffle_edges(remade_G,
                    swaps_per_edge)
            new_Garrett_G = nx.from_edgelist(new_Garrett_edges)
            new_r = nx.degree_assortativity_coefficient(new_Garrett_G)
            print("New R:", new_r)
            show_log_jkk(degrees, dk, degree_sequence, new_Garrett_G,
                "Garrett's Rewired in Post")

        


    pl.show()
