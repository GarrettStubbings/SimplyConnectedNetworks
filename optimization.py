"""
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
import networkx as nx
import datetime
date = datetime.date.today().strftime('%e%b%y')
date = date.strip(' ')
import matplotlib.style
import matplotlib as mpl
import matplotlib.cm as cmx
from scipy.optimize import curve_fit
import numpy as np
import sys

from folder_management import *

mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['lines.markeredgecolor'] = 'k'
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

def get_scale_free_degrees(N, temp_folder, alpha, avg_deg, seed):
    """
    Use Spencer's code to get the degree sequence for scale free network of
    given parameters
    """
    # params for generating network
    params = ["0", "0.0", "0.0", "0.0", str(N), "2", str(alpha), str(avg_deg),
                "AND", "Single", temp_folder, "ScaleFree", "0",
                str(seed), "0.0"]
    #sub.call(["make", "backup"])
    command = ['./main'] + params 
    sub.call(command)

    output_files = os.listdir(temp_folder)
    #print(output_files)
    degree_file = output_files[0]
    degree_sequence = pl.loadtxt(temp_folder + degree_file)[:,1]

    for f in output_files:
        os.remove(temp_folder + f)

    k_min = pl.amin(degree_sequence)
    k_max = pl.amax(degree_sequence)

    degree_bins = pl.arange(k_min, k_max + 2)
    degrees = degree_bins[:-1]
    dk = pl.histogram(degree_sequence, degree_bins)[0]
    degrees = degrees[dk != 0]
    dk = dk[dk != 0]

    return degrees, dk, degree_sequence

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

# Now Defunked?
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
    size = len(degrees)
    available_stubs = pl.copy(stubs)
    k_stubs = stubs[k_index]
    k = degrees[k_index]
    dk = degree_counts[k_index]

    #print("Starting with ", k_stubs, "Stubs")
    #print(len(stubs), len(degree_counts), dk)
    node_constraints = dk*degree_counts
    node_constraints[k_index] = dk*(dk-1)//2

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
            #print(stubs[k_index], "Left to go")
            #print(degrees[available_stubs>0],
            #        available_stubs[available_stubs>0],
            #        stubs[available_stubs>0])
            break


        #if stubs[k_index] > pl.sum(available_stubs):
            #print(stubs[k_index], "Left to go")
            #print(degrees[available_stubs>0],
            #        available_stubs[available_stubs>0],
            #        stubs[available_stubs>0])


        min_k_stubs = pl.amin([stubs[k_index], available_stubs[k_index]])

        if pl.sum(available_stubs) == 0:
            #print("Ran Out of available stubs!!!")
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
                #print("\n Probabilities Contain NaN\n")
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
            #print(assigned_stubs)
            #print("Ran Out of available stubs!!!")
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
        #print("We Usin Degree 1!")
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
        #print("Sum of Jkk:", pl.sum(jkk), "J(1, 1)=", jkk[0,0])
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

    #print("Couldn't Build the Graph (Valid Jkk)")
    #pl.figure(figsize = (8,6))
    #pl.imshow(jkk, origin = 'lower')
    #pl.colorbar()
    #print(jkk)
    non_1_edges = pl.sum(jkk[:,1:] > 0, axis = 1)
    #print(len(non_1_edges), pl.shape(jkk))
    #print(non_1_edges - jkk[:,0])
    #print(degree_counts * (degrees - 1)/degrees - jkk[:,0])
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

def run_network(G, running_folder, params, param_names, 
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
    output_network(G, running_folder + "SimplyConnected.csv")
    #print("Outputting to " + input_dir + "SimplyConnected.csv")

    network_travel_time = time.time() - network_travel_start
    ### running stuff
    # compile if needed
    command = ['./testNetwork'] + params 
    raw_start_time = time.time()
    sub.call(command)
    raw_simulation_time = time.time()  - raw_start_time
    #print("Real Simulation Time: {0:.3f}, Network Travel() time: {1:.3f}".format(
    #    raw_simulation_time, network_travel_time))

    #print("Reading data from " + running_folder)

    output_files = os.listdir(running_folder)
    #print(output_files)

    means = []
    errors = []
    for m in measures:
        performance_file = [running_folder + f for f in output_files if m in
                                                                        f][0]
        performance = pl.loadtxt(performance_file)
        if "eath" in m or "QALY" in m:
            performance /= scale
        means.append(pl.average(performance))
        errors.append(pl.std(performance)/pl.sqrt(len(performance) - 1))

    for f in output_files:
        os.remove(running_folder + f)

    if len(performance) == 1:
        return means[0], errors[0]
    else:
        return means, errors 

def calculate_entropy(p):
    valid_p = p[p > 0]
    return -pl.sum(valid_p * pl.log(valid_p))

def consolidate_components(connected_components, degree_sequence, G):
    """
    Genius Idea: connect the disconnected components myself
    
    only gonna bother with the k = 1 case.
    General case is hard
    """
    #print("We Consolidating")
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
                #print("Using the Giant Component")
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
        #print("Couldn't get a connected one")
        for sub_G in connected_components[1:]:
            sub_graphs = [degree_sequence[i] for i in sub_G]
            #print(sub_graphs)


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



if __name__ == '__main__':

    pl.close('all')
    #pl.seed(2)
    fs = 14
    matrix_folder = 'AdjacencyMatrices/'
    plots_folder = 'Plots/AssortativityWeightedNetworks/'

    ############## parameters of running the simulaton ##############

    avg_k = 4
    scale = 0.00183

    ############### Input parameters ###################
    temp_folder = sys.argv[1]
    output_folder = sys.argv[2]
    details = sys.argv[3]
    method = sys.argv[4]
    run_time = float(sys.argv[5])
    N = int(sys.argv[6])
    number = sys.argv[7]
    health_measure = sys.argv[8]
    seed = int(sys.argv[-2])
    entropy_weighting = float(sys.argv[-1])

    
    print("RUN PARAMS: Method: {4}, Lambda:{0}, Run Time: {1}, N: {2},"\
            " number: {3}, Health Measure: {6}, seed: {5}".format(
            entropy_weighting, run_time, N, number, method, seed,
            health_measure))

    details += "Lambda"+ str(entropy_weighting) + "/" + "Seed" + str(seed) + "/"
    output_folder = output_folder + details
    running_folder = temp_folder + details

    make_directory(output_folder)
    make_directory(running_folder)
    clean_directory(running_folder)

    pl.seed(seed)
    seed  = str(seed)
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
    params = [number, N, numChanges, avg_k, running_folder, seed,
        run_hours, evo_condition, initial_distribution, Lambda,
        beta, power]


    k_min = 2
    dk = [0,2,3,4,3,2]
    degrees = [i+1 for i in range(len(dk))]
    bonus_degrees = [i+1 for i in range(len(dk), 100)]
    bonus_degrees = [10, 15, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500,
                    1000, 2000]
    bonus_degrees = [k for k in bonus_degrees if k < N]
    k_max = bonus_degrees[-1]
    dk += [0 for i in range(len(bonus_degrees))]
    degrees += bonus_degrees
    degrees = pl.asarray(degrees)
    dk = pl.asarray(dk)
    dk = dk[degrees >= k_min]
    degrees = degrees[degrees >= k_min]
    pk = pl.asarray(dk)/pl.sum(dk)

    avg_k = pl.dot(degrees, pk)
    params[3] = str(avg_k)
    original_N = N

    static_assortativity = False
    if static_assortativity:
        p_assortative = 0.0
        p_dissassortative = 0.0
        p_random = 1.0 - p_assortative - p_dissassortative


    best_pk = pl.copy(pk)
    best_entropy = -0.01
    best_health_mean = -0.01
    health_mean = best_health_mean
    best_alpha = 2.27
    alpha = best_alpha

    measures = ["DeathAge", "HealthyAging", "QALY"]
    merit_index = pl.where(pl.asarray(measures) == health_measure)[0][0]
    best_health_means = []
    best_health_errors = []
    health_means = [0,0,0]
    health_errors = [0,0,0]

    best_entropies = []
    best_pkk_entropies = []
    best_merits = []
    best_pks = []
    best_sampled_dks = []
    best_probs = []
    best_alphas = []
    iterations = []
    best_G = 1
    best_pkk = 1

    #print("\nMethod:", method, "\n")

    if "aria" in method:
        #print("USING VARIATIONAL")
        #sub.call(["make","main"])
        avg_deg = 4.0
        width = 0.03
        limits = [2.0, 5.0]
        degrees, dk, degree_sequence = get_scale_free_degrees(
                                    N, running_folder, alpha, avg_deg, seed)
        k_min = 2
        k_max = "UNDEF"

    if entropy_weighting == 1:
        static_assortativity = 1
        p_assortative = 0
        p_dissassortative = 0
        p_random = 1

    best_merit = entropy_weighting*best_entropy + (
        1-entropy_weighting)*best_health_mean

    fraction = 0.1 #polarized changes away from 0.5 (best way to make progress)
    max_changes = 100

    start_time = time.time()
    elapsed_time = (time.time() - start_time)/60 # in minutes
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
            print("{0}% Done. Iteration {1}. Entropy {2}. <health> {3}".format(
                    progress, i, best_entropy, best_health_mean))
            if "aria" in method:
                print("Best Alpha:", best_alpha)

        ################# Getting the degree sequence -----------~~~~~~~~~~~

        # Using variational method (still only scale-free here)
        if "aria" in method:
            alpha = gaussian_parameter_change(alpha, width, limits)
            degrees, dk, degree_sequence = get_scale_free_degrees(
                                    N, running_folder, alpha, avg_deg, seed)
        # Using non-parametric method
        else:
            if i % 2 == 1:
                num_changes -= 1
            if num_changes < 1:
                num_changes = 1#pl.randint(10) + 1
            pk = pl.copy(best_pk)
            successful_changes = 0
            fraction = pl.random()

            #### Changing Degree Distribution ####
            while successful_changes < num_changes:
                n1, n2 = pl.choice(pl.arange(len(pk)), size = 2,
                    replace = False)
                pk, success = change_pk(pk, degrees, n1, n2, fraction, avg_k) 
                #print(success, successful_changes)
                if success:
                    successful_changes += 1
           
            dk = random_hub_node_degrees(pk, degrees, N)[1]
            dk = pl.asarray(dk)

        ###########  Building the Network $$$$%%%%%%%%%%%%
        build_start = time.time()
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
        entropy = calculate_entropy(pkk)

        sampled_N = G.number_of_nodes()
        if sampled_N < 0.95*original_N:
            continue
        params[1] = str(sampled_N)


        ####### Running The Network #######
        if entropy_weighting != 1:
            simulation_start = time.time()
            health_means, health_errors = run_network(G, running_folder,
                                params, param_names, measures)
            simulation_time = time.time() - simulation_start
            health_mean = health_means[merit_index] 

        merit = entropy_weighting*entropy + (
                                1-entropy_weighting)*health_mean
        r = pl.random()

        ########### Merit Function Evaluation and Update #####################
        if pl.exp((merit - best_merit)*pl.sqrt(i)*2) > r:
            best_probs.append([p_assortative, p_dissassortative,
                                                p_random])
            best_pk = pl.copy(pk)
            best_pks.append(pl.copy(pk))
            best_entropies.append(calculate_entropy(best_pk))
            if entropy_weighting != 1:
                best_health_mean = health_mean
                best_health_means.append(health_means)
                best_health_errors.append(health_errors)
            best_G = G
            best_sampled_dks.append(dk)
                
            best_alpha = alpha
            best_alphas.append(best_alpha)
            best_pkk = pkk
            best_entropy = entropy
            best_pkk_entropies.append(entropy)
            best_merit = merit
            best_merits.append(merit)
            fractions.append(fraction)
            iterations.append(i)

    save_output = 1
    if save_output:

        pl.savetxt(output_folder+"BestPkk.csv", best_pkk, delimiter = ",")
        pl.savetxt(output_folder+"BestHealthMeans.txt",best_health_means)
        pl.savetxt(output_folder+"BestHealthErrors.txt",best_health_errors)
        pl.savetxt(output_folder+"BestPkkEntropies.txt", best_pkk_entropies)
        pl.savetxt(output_folder+"IterationNumbers.txt", iterations)
        pl.savetxt(output_folder+"BestAlphas.txt", best_alphas)
        pl.savetxt(output_folder+"BestProportions.txt", best_probs)

        nx.write_edgelist(best_G, output_folder+"BestNetwork.csv",
                            delimiter = ",", data = False)



    plot_output = 0
    if plot_output:
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


        pl.show()

    remove_directory(running_folder)
