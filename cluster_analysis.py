"""

Code for visualizing performance of optimized networks.

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
from optimization import *

import powerlaw
import networkx as nx
import datetime
date = datetime.date.today().strftime('%e%b%y')
date = date.strip(' ')
import matplotlib.style
import matplotlib as mpl
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy.optimize import curve_fit
from folder_management import *
import numpy as np

erms = 14
mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['lines.linewidth'] = 1.5
#mpl.rcParams['lines.markeredgecolor'] = 'k'
#mpl.rcParams['lines.markeredgewidth'] = 0.5
fs = 14
scale = 1/0.00183

def scrape_data(data_dir, N, optimization_measure = "DeathAge",
                plotting_targets = [1,2]):
    """
    Given a base directory, dig through all the lambda, entropy targets, and
    seeds and grab all the peak performance data.
    """
    if optimization_measure == "DeathAge":
        merit_index = 0
    elif optimization_measure == "QALY":
        merit_index = 1
    else:
        merit_index = 2

    measure_files = ["BestHealthMeans.txt"]*3
    measure_names = ["Death Age", "Healthy Aging", "QALY"]
    print(os.listdir(data_dir))
    lambda_folders = sorted([f for f in os.listdir(data_dir) if
                                    "DS_Store" not in f])
    lambda_targets = [float(f.strip("Lambda")) for f in
                                                    lambda_folders]
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]

    all_measures = []
    all_errors = []
    entropies = []
    paths = []
    moments = []
    rs = []
    k_maxes = []

    for i, m in enumerate(measure_files):
        measures = []
        errors = []
        colour_metrics = []

        # Loop over the different Lambdas
        for n, lambda_path in enumerate(lambda_paths):

            entropy_folders = sorted([f for f in os.listdir(lambda_path) if
                                            "DS_Store" not in f])
            entropy_targets = [float(f.strip("EntropyTarget")) for f in
                                                            entropy_folders]
            entropy_paths = [lambda_path + f + "/" for f in entropy_folders]


            for j, path in enumerate(entropy_paths):
                seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                                    in f]
                seeds = [f.strip("Seed") for f in seed_folders]
                seed_paths = [path + f + "/" for f in seed_folders]
                for k, seed_path in enumerate(seed_paths):
                    
                    # taking care of missing data
                    try:
                        # get entropy stuff (only need to do this once)
                        entropy_data = pl.loadtxt(seed_path +
                                                "BestPkkEntropies.txt")
                        # the health stuff is all in the same file, but in different
                        # columns
                        measure_data = pl.loadtxt(seed_path+m)
                        measure_error_data = pl.loadtxt(seed_path +
                            "BestHealthErrors.txt")
                        try:
                            mode = "2d"
                            entropy = entropy_data[-1]
                        except IndexError:
                            mode = "1d"
                            entropy = entropy_data
                        degree_data = pl.loadtxt(seed_path + "BestDegrees.txt")
                        pk_data = pl.loadtxt(seed_path + "BestDks.txt")/float(N)
                        r_data = pl.loadtxt(seed_path +
                            "BestAssortativities.txt")

                        if mode == "2d":
                            max_degree = degree_data[-1, -1]
                            degrees = degree_data[-1,:]
                            best_pk = pk_data[-1,:]
                            r = r_data[-1]
                        else:
                            max_degree = degree_data[-1]
                            degrees = degree_data
                            best_pk = pk_data
                            r = r_data
                        best_pk = best_pk[~pl.isnan(degrees)]
                        degrees = degrees[~pl.isnan(degrees)]
                        k_max = pl.amax(degrees[best_pk != 0]/int(N))
                        second_moment = pl.dot(best_pk, degrees**2)

                        if pl.shape(measure_data) == (0,):
                            break
                        else:
                            if mode == "2d":
                                measure = measure_data[-1,i]
                                measure_error = measure_error_data[-1,i]
                            else:
                                measure = measure_data[i]
                                measure_error = measure_error_data[i]

                            if "Healthy" in measure_names[i]:
                                measure *= 100
                                measure_error *= 100

                    except OSError:
                        if j == 0:
                            print("No Data For Lambda {0}, Seed {1}".format(
                                    entropy_targets[j], seeds[k]))
                            print("Path is", seed_path)
                        continue

                    if i == 0:
                        entropies.append(entropy)
                        paths.append(seed_path)
                        moments.append(second_moment/int(N))
                        rs.append(r)
                        k_maxes.append(k_max)
                    measures.append(measure)
                    errors.append(measure_error)
        # find the appropriate runs to show graphs/pkk
        entropies = pl.asarray(entropies)
        measures = pl.asarray(measures)
        plotting_weight = 0.99
        if i == merit_index:
            plotting_indices = []
            #plotting_targets = [1,2,3,4]
            for pt in plotting_targets:
                merit = 1 - plotting_weight * pl.absolute(entropies  - pt) + (
                                    (1 - plotting_weight) * measures)
                # look for maximum, if there's not enough data return 0
                try:
                    index = pl.where(merit == pl.amax(merit))[0][0]
                except ValueError:
                    index = 0
                plotting_indices.append(index)


        all_measures.append(measures) 
        all_errors.append(errors)
    moments = pl.asarray(moments)
    k_maxes = pl.asarray(k_maxes)
    return (all_measures, all_errors, entropies, moments, rs, k_maxes, paths,
                plotting_indices)

def generate_entropy_plots(data_dir, reference_data_dir, plots_dir, plot_title,
                    save_plots, N,
                    colour_metric = "Moment",
                    plotting_targets = [0,2,3,10],
                    optimization_measure = "DeathAge"):


    measure_names = ["Death Age", "Healthy Aging", "QALY"]
    measures, errors, entropies, moments, rs, k_maxes, paths, plotting_indices = (
                                            scrape_data(data_dir, N,
                                            plotting_targets=plotting_targets))
    scale_free_data = scrape_data(reference_data_dir, N,
                                    plotting_targets=plotting_targets)
    sf_measures, sf_errors, sf_entropies, sf_moments, sf_rs, sf_k_maxes= (
                                                    scale_free_data[:6])
    entropies = list(entropies) + list(sf_entropies) 
    if colour_metric == "Moment":
        colour_metrics = list(moments) + list(sf_moments)
    elif colour_metric == "kMax":
        colour_metrics = list(k_maxes) + list(sf_k_maxes)
    else:
        colour_metrics = rs + sf_rs

    ranks = pl.searchsorted(pl.sort(colour_metrics), colour_metrics)
    colours = get_colours_from_cmap(colour_metrics)
    cmap = mpl.cm.get_cmap("viridis")#(x)#[np.newaxis, :, :3]
    norm = mpl.colors.Normalize(vmin=pl.amin(colour_metrics),
                                    vmax=pl.amax(colour_metrics))
    
    for i, m in enumerate(measures):
        measure = list(m) + list(sf_measures[i])
        error = list(errors[i]) + list(sf_errors[i])
        measure_fig, measure_ax = pl.subplots(1,2,figsize=(8, 6),
            gridspec_kw={"width_ratios": [20, 1]},)
        measure_fig.subplots_adjust(wspace=0.01)

        if i == 0: 
            if int(N) > 500000:
                graph_fig, graph_ax = pl.subplots(2, len(plotting_indices),
                    gridspec_kw={"height_ratios": [1, 8/6],
                                "width_ratios": 
                                    [1 for i in range(len(plotting_indices))]},
                    figsize = (8, 6))
                pl.subplots_adjust(wspace=0.2, hspace = 0.2)

            else:
                graph_fig, graph_ax = pl.subplots(3, len(plotting_indices),
                    gridspec_kw={"height_ratios": [1,1.2, 1],
                                "width_ratios": 
                                    [1 for i in range(len(plotting_indices))]},
                    figsize = (2.7*len(plotting_indices), 6))
                graph_fig.subplots_adjust(hspace=0.4)
            #pl.subplots_adjust(hspace=0)

        for l, c_m in enumerate(colour_metrics):
            if l < len(moments):
                marker = "o"
                alpha = 1.0
            else:
                marker = "*"
                alpha = 0.7
            if l in (plotting_indices):
                path = paths[l]
                
                if i == 0:
                    compare_graphs(l, plotting_indices, path, int(N),
                                        graph_fig,graph_ax,
                                        entropies[l], scale_free_data,
                                        plots_dir)
                marker_size = 15
                mec = "m"
                zorder=10
                mew=1.5
            elif l - len(moments) in scale_free_data[-1]:
                marker_size = 15
                mec = "m"
                zorder=10
                mew=1.5
                alpha = 1.0                
            else:
                marker_size = 9
                mec = "k"
                zorder=1
                mew=0.5

            rank = ranks[l]
            measure_ax[0].errorbar(entropies[l], measure[l], yerr=error[l],
                    ls = "none", markersize = marker_size, mec = mec,
                    marker = marker, color = colours[l], capsize=3, mew=mew,
                    zorder=zorder, alpha = alpha)
            
        cb1 = mpl.colorbar.ColorbarBase(measure_ax[1], cmap=cmap,
                                        norm=norm,
                                        orientation='vertical')
        #cb1.set_label('k max')

        measure_ax[0].set_ylabel(measure_names[i], fontsize = fs)
        measure_ax[0].set_xlabel('Entropy', fontsize = fs)
        if colour_metric == "Moment":
            measure_ax[1].set_ylabel(r"Scaled Second Moment " + 
                r"$\frac{\langle k^{2} \rangle}{N}$", fontsize = fs)
        elif colour_metric == "Assortativity":
            measure_ax[1].set_ylabel('Degree Assortativity Coefficient')
        else:
            measure_ax[1].set_ylabel('kMax / N')
        measure_ax[0].set_title(plot_title)
 
        if save_plots:
            measure_fig.savefig(plots_dir + plot_title.replace(" ","") +
                measure_names[i].replace(" ","") + "VsEntropy" + 
                colour_metric + "Colours.pdf")
            graph_fig.savefig(plots_dir + plot_title.replace(" ","") +
                "GraphComparison.pdf")

    tds = list(measures[0]) + list(sf_measures[0])
    td_error = list(errors[0]) + list(sf_errors[0])
    colours = get_colours_from_cmap(colour_metrics)
    cmap = mpl.cm.get_cmap("viridis")
    norm = mpl.colors.Normalize(vmin=pl.amin(colour_metrics),
                                    vmax=pl.amax(colour_metrics))
    sf_index = len(measures[0]) 
    marker_size=8
    alpha = 1
    for j in range(1,3):
        measure = list(measures[j]) + list(sf_measures[j])
        error = list(errors[j]) + list(sf_errors[j])
        measure_fig, measure_ax = pl.subplots(1,2,figsize=(8, 6),
            gridspec_kw={"width_ratios": [20, 1]},)
        measure_fig.subplots_adjust(wspace=0.01)

        
        for i, e in enumerate(colour_metrics):
            if i >= sf_index:
                marker = "*"
                alpha = 0.7
            else:
                alpha = 1.0
                marker = "o"
            measure_ax[0].errorbar(tds[i], measure[i], yerr=error[i],
                    xerr=td_error[i],
                    ls = "none", markersize = marker_size, mec = mec,
                    marker = marker, color = colours[i], capsize=3, mew=mew,
                    zorder=zorder, alpha = alpha)

        cb1 = mpl.colorbar.ColorbarBase(measure_ax[1], cmap=cmap,
                                        norm=norm,
                                        orientation='vertical')
        measure_ax[0].set_ylabel(measure_names[j], fontsize = fs)
        measure_ax[0].set_xlabel('Death Age', fontsize = fs)
        #measure_ax[1].set_ylabel("Moment", fontsize = fs)
        if colour_metric == "Moment":
            measure_ax[1].set_ylabel(r"Scaled Second Moment " + 
                r"$\frac{\langle k^{2} \rangle}{N}$", fontsize = fs)
        elif colour_metric == "Assortativity":
            measure_ax[1].set_ylabel('Degree Assortativity Coefficient')
        else:
            measure_ax[1].set_ylabel('kMax / N')

        if j == 2:
            line = pl.array([0, 10000])
            measure_ax[0].plot(line, line-20, 'C0')
            measure_ax[0].set_ylim(pl.amin(measure) - 3,pl.amax(measure) + 3)
            measure_ax[0].set_xlim(pl.amin(tds) - 3,pl.amax(tds) + 3)
        if save_plots:
            measure_fig.savefig(plots_dir + plot_title.replace(" ","") +
                measure_names[j].replace(" ","") + "VsDeathAge.pdf")

def convex_hull(x, y):
    """
    Yoinked from WikiBooks with some modification

    Computes the convex hull of a set of 2D points.

    Input: x and y data (converted from the (x,y) points in original)
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    """
    points = [(a,b) for a,b in zip(x,y)]
    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    hull_x = pl.asarray([p[0] for p in upper])
    hull_y = pl.asarray([p[1] for p in upper])

    return hull_x, hull_y

def power_law(k, alpha, c):
    return c*k**(-alpha)

def fit_alpha(degrees, dk):
    #pk = dk/pl.sum(dk)
    #p0 = [5.0, pl.sum(degrees**-5)]
    #popt, pcov = curve_fit(power_law, degrees, pk, p0 = p0, maxfev=10000)
    #alpha, c = popt
    #return alpha
    dk = dk[~pl.isnan(degrees)]
    degrees = degrees[~pl.isnan(degrees)]
    padded_dk = pl.zeros(int(pl.amax(degrees) + pl.amin(degrees)))
    non_zero_indices = degrees.astype(int)
    padded_dk[non_zero_indices] = dk
    results = powerlaw.Fit(padded_dk)
    return results.power_law.alpha

def get_n_by_m_subplots(n, m, size = (12,8)):
    fig, ax = pl.subplots(n, m, figsize = size)
    sequential_axes = []
    for i in range(n):
        for j in range(m):
            p = ax[int(i), j%m]
            sequential_axes.append(p)
    return fig, sequential_axes

def compare_graphs(index, indices, path, N, fig, axes, entropy, scale_free_data,
                    plot_folder):
    """
    Show the graphs (spring layout) if N < 500, and show Pkk
    Will be a set of graphs selected from the vs entropy plot
    """
    print(index, indices, path, N)
    i = pl.where(pl.asarray(indices) == index)[0][0]
    sf_index = scale_free_data[-1][i]
    sf_path = scale_free_data[-2][sf_index]

    use_x_ticks = True
    use_y_ticks = True
    if i != 0:
        use_y_ticks = False
        
    # small graphs: draw the graphs and show the pkk
    G = nx.read_edgelist(path + "BestNetwork.csv",
                            delimiter = ",")
    details_title = "N {0}, Entropy {1:.3f}".format(N, entropy)
    running_folder = "TempData/ClusterReRunning/"
    rerun_network(G, running_folder, params,
                                    param_names, plot_folder,
                                    details_title, 1000)

    degrees = pl.loadtxt(path + "BestDegrees.txt")
    dk = pl.loadtxt(path + "BestDks.txt")
    sf_degrees = pl.loadtxt(sf_path + "BestDegrees.txt")
    sf_dk = pl.loadtxt(sf_path + "BestDks.txt")
    if len(pl.shape(degrees)) != 1:
        degrees = degrees[-1,:]
        dk = dk[-1,:]
    if len(pl.shape(sf_degrees)) != 1:
        sf_degrees = sf_degrees[-1,:]
        sf_dk = sf_dk[-1,:]
    full_degrees = pl.copy(degrees) 
    alpha_data = pl.loadtxt(path + "BestAlphas.txt")
    sf_alpha_data = pl.loadtxt(sf_path + "BestAlphas.txt")
    show_alpha = 0
    if pl.shape(alpha_data) != (0,):
        show_alpha = 1
        if pl.shape(alpha_data) != ():
            alpha = alpha_data[-1]
        else:
            alpha = alpha_data
        if pl.shape(sf_alpha_data) != ():
            sf_alpha = sf_alpha_data[-1]
        else:
            sf_alpha = sf_alpha_data

    dk = dk[~pl.isnan(degrees)]
    degrees = degrees[~pl.isnan(degrees)]
    sf_dk = sf_dk[~pl.isnan(sf_degrees)]
    sf_degrees = sf_degrees[~pl.isnan(sf_degrees)]
    print("Average Degree:", pl.dot(degrees, dk)/N)
 
    show_every = max([int(len(degrees)/12), 1])
    """
    estimated_alpha = fit_alpha(degrees, dk)
    fitted_power_law = degrees**(-estimated_alpha)
    fitted_power_law *= N/fitted_power_law[0]
    axes[0,i].plot(degrees, fitted_power_law, 'C1',
        label = "Fitted\n" + r"$\alpha = {:.3f}$".format(estimated_alpha))
    if show_alpha:
        direct_power_law  = degrees**(-sf_alpha)
        direct_power_law *= N/direct_power_law[0]
        axes[0,i].plot(degrees, direct_power_law, 'C0',
            label = "Variational\n" + r" $\alpha = {:.3f}$".format(sf_alpha))
    """
    axes[0,i].plot(degrees, dk, 'C0o', label = "Lifespan Optimized")
    axes[0,i].plot(sf_degrees, sf_dk, 'C1^', label = "QALY Optimized")
    ref_k = pl.array([10, 100, 1000]).astype(float)
    ref_line = ref_k**(-3)*10e5
    axes[0,1].plot(ref_k*0.5, ref_line, "k:")
    axes[0,1].annotate(r"$\sim \mathrm{k}^{-3}$", xy = (2e1, 1.5e1))
    axes[0,i].set_xscale('log')
    axes[0,i].set_yscale('log')
    axes[0,i].set_xlabel('k', labelpad = -4.0)
    #axes[0,i].set_title('Entropy: {:.2f}'.format(entropy))
    axes[0,i].set_ylim(0.5, N*1.5)
    axes[1,i].set_xlabel("k", labelpad=-8)
    cluster_ref_line = ref_k**(-1)*1
    #axes[1,i].set(adjustable = "box")
    #axes[0,i].set(adjustable = "box")

    if i == 0:
        axes[0,i].set_ylabel('Count')
        axes[1,i].set_ylabel("k'", labelpad=-4)
    
    axes[0,0].legend(bbox_to_anchor = (0, 0.0), loc='lower left', fontsize=8)
    pkk = pl.loadtxt(path + "BestPkk.csv", delimiter = ",")
    axes[0,0].annotate("A)", xy = (0.9, 0.88),
                                xycoords="axes fraction", fontsize=12)
    axes[0,1].annotate("B)", xy = (0.9, 0.88),
                                xycoords="axes fraction", fontsize=12)
    axes[1,0].annotate("C)", xy = (0.85, 0.9),
                                xycoords="axes fraction", fontsize=12)
    axes[1,1].annotate("D)", xy = (0.85, 0.9),
                                xycoords="axes fraction", fontsize=12)

    try:
        #draw_spring_graph(G, axes[2,i])
        axes[2,0].annotate("E)", xy = (0.9, 0.88),
                                    xycoords="axes fraction", fontsize=12)
        axes[2,1].annotate("F)", xy = (0.9, 0.88),
                                    xycoords="axes fraction", fontsize=12)


        fmts = ["C0o", "C1^"]
        for g_id, data_dir in enumerate([path, sf_path]):
            G = nx.read_edgelist(data_dir + "BestNetwork.csv",
                                    delimiter = ",")

            # Doing the clustering as a function of degree in stead
            degree_sequence = pl.asarray([G.degree(i) for i in G.nodes()])
            degrees = pl.asarray(sorted(list(set(degree_sequence))))
            local_ccs = pl.asarray([v for k, v in nx.clustering(G).items()])
            average_local_ccs = []
            local_ccs_error = []
            for k in degrees:
                ccs = local_ccs[degree_sequence == k]
                average_local_ccs.append(pl.average(ccs))
                local_ccs_error.append(pl.std(ccs)/(len(ccs)-1))
            #axes[2,i].errorbar(degrees, average_local_ccs, yerr = local_ccs_error,
            #        fmt = fmts[g_id])
            axes[2,i].plot(degree_sequence, local_ccs, fmts[g_id])
        axes[2,i].set_xscale('log')
        axes[2,1].plot(ref_k, cluster_ref_line, "k:")
        axes[2,1].annotate(r"$\sim \mathrm{k}^{-1}$", xy = (3e1, 2e-3))
        axes[2,i].set_yscale('log')
        axes[2,0].set_ylabel('Clustering Coefficient')
        axes[2,i].set_xlabel('k', labelpad=-4)
        axes[2,i].tick_params(axis="y", which="major", pad = 0)
        axes[2,i].set_ylim(1e-3, 1.8)
        axes[2,i].set_xlim(1.5, int(N)*1.3)
    except IndexError:
        print("Not doing the third row for this plot")






    # catch test-runs that never get any data...
    try:
        #pkk = pkk[dk > 0, dk > 0]
        im = axes[1,i].imshow(pkk, origin = 'lower', cmap = "Greys")
        cbar = fig.colorbar(im,ax = axes[1,i], fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=10)
        #ax_divider = make_axes_locatable(axes[1,0])
        #cbar_ax = ax_divider.append_axes("right", size="5%", pad="2%")
        #cbar = fig.colorbar(im, ax = cbar_ax)
    except TypeError:
        print("No Joint matrix for this data...")
    if 1:
        degrees = pl.copy(full_degrees)
        axes[1, i].set_xticks(pl.arange(len(degrees))[::show_every])
        axes[1, i].set_xticklabels(degrees[::show_every].astype(int),
                                    rotation = -90, fontsize=8)
    else:
        if i == 3:
            print("Killing the ticks")
        axes[1, i].set_xticks([])
    if 1:
        axes[1, i].set_yticks(pl.arange(len(degrees))[::show_every])
        axes[1, i].set_yticklabels(degrees[::show_every].astype(int),
                                    fontsize=8)
    else:
        axes[1, i].set_yticks([])

def check_progress(data_dir, plots_dir, save_plots):
    """
    For a given lambda, plot the data from each run (seeds) for each measure
    (Vs. Iteration type of deal).
    """
    runs = sorted([f for f in os.listdir(data_dir) if "DS_Store" not in f])
    print(runs)
    run_data_dir = [data_dir + f + "/" for f in runs]

    seeds = [float(run.strip("Seed")) for run in runs]
    colours = get_colours_from_cmap(pl.arange(len(seeds)))

    measure_files = ["BestHealthMeans.txt"]*3 + ["BestPkkEntropies.txt",
                        "BestAlphas.txt", "BestProportions.txt"]
    measure_names = ["Death Age", "Healthy Aging", "QALY", "Pkk Entropy",
                        "Alpha", "Proportions"]

    markers = ['o', '*', '^']
    choices = ["Assortative", "Dissassortative", "Random"]

    fig, linear_axes = get_n_by_m_subplots(2,3)
    pl.subplots_adjust(right=0.87)
    # for each measure
    for i, m in enumerate(measure_files):
        ax = linear_axes[i]
        measure_name = measure_names[i]
        
        # for each seed
        for j, r in enumerate(run_data_dir):
            colour = colours[j]
            # taking care of missing data
            try:
                iterations = pl.loadtxt(r + "IterationNumbers.txt")
                # the health stuff is all in the same file, but in different
                # columns
                if i < 3:
                    measure = pl.loadtxt(r+m)
                    if pl.shape(measure) == (0,):
                        pl.close(fig)
                        break
                    else:
                        measure = measure[:,i]
                else:
                    measure = pl.loadtxt(r + m)
            except OSError:
                continue
            if i == 5:
                for k in range(3):
                    ax.plot(iterations, measure[:,k], ls = "none",
                            marker = markers[k], color = colour,
                            label = choices[k])
                    if j == 0:
                        ax.legend(bbox_to_anchor=(0.95,1),
                            loc = "upper left", frameon=False)
            else:
                ax.plot(iterations, measure, ls = "none", marker = "o",
                        color = colour, label = "Seed: {}".format(seeds[j]))
                if i == 2:
                    ax.legend(bbox_to_anchor=(0.95,1),
                            loc = "upper left", frameon=False)

        ax.set_xlabel('Iteration Number')
        ax.set_ylabel(measure_name)
        ax.set_xscale('log')

def versus_entropy(data_dir, plots_dir, plot_title, save_plots,
                    N, plotting_targets = [], colour_metric = "Moment",
                    optimization_measure="DeathAge"):
    if optimization_measure == "DeathAge":
        merit_index = 0
    elif optimization_measure == "QALY":
        merit_index = 1
    else:
        merit_index = 2

    measure_files = ["BestHealthMeans.txt"]*3
    measure_names = ["Death Age", "Healthy Aging", "QALY"]
    lambda_folders = sorted([f for f in os.listdir(data_dir) if
                                    "DS_Store" not in f])
    lambda_targets = [float(f.strip("Lambda")) for f in
                                                    lambda_folders]
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]


    for i, m in enumerate(measure_files):
        measures = []
        errors = []
        colour_metrics = []
        entropies = []
        paths = []


        # Loop over the different Lambdas
        for n, lambda_path in enumerate(lambda_paths):

            entropy_folders = sorted([f for f in os.listdir(lambda_path) if
                                            "DS_Store" not in f])
            entropy_targets = [float(f.strip("EntropyTarget")) for f in
                                                            entropy_folders]
            entropy_paths = [lambda_path + f + "/" for f in entropy_folders]


            for j, path in enumerate(entropy_paths):
                seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                                    in f]
                seeds = [f.strip("Seed") for f in seed_folders]
                seed_paths = [path + f + "/" for f in seed_folders]
                for k, seed_path in enumerate(seed_paths):
                    
                    # taking care of missing data
                    try:
                        # get entropy stuff (only need to do this once)
                        entropy_data = pl.loadtxt(seed_path +
                                                "BestPkkEntropies.txt")
                        # the health stuff is all in the same file, but in different
                        # columns
                        measure_data = pl.loadtxt(seed_path+m)
                        measure_error_data = pl.loadtxt(seed_path +
                            "BestHealthErrors.txt")
                        try:
                            mode = "2d"
                            entropy = entropy_data[-1]
                        except IndexError:
                            mode = "1d"
                            entropy = entropy_data
                        degree_data = pl.loadtxt(seed_path + "BestDegrees.txt")
                        pk_data = pl.loadtxt(seed_path + "BestDks.txt")/float(N)
                        r_data = pl.loadtxt(seed_path +
                            "BestAssortativities.txt")

                        if mode == "2d":
                            max_degree = degree_data[-1, -1]
                            degrees = degree_data[-1,:]
                            best_pk = pk_data[-1,:]
                            r = r_data[-1]
                        else:
                            max_degree = degree_data[-1]
                            degrees = degree_data
                            best_pk = pk_data
                            r = r_data
                        best_pk = best_pk[~pl.isnan(degrees)]
                        degrees = degrees[~pl.isnan(degrees)]
                        second_moment = pl.dot(best_pk, degrees**2)

                        if pl.shape(measure_data) == (0,):
                            break
                        else:
                            if mode == "2d":
                                measure = measure_data[-1,i]
                                measure_error = measure_error_data[-1,i]
                            else:
                                measure = measure_data[i]
                                measure_error = measure_error_data[i]

                            if "Healthy" in measure_names[i]:
                                measure *= 100
                                measure_error *= 100

                    except OSError:
                        if j == 0:
                            print("No Data For Lambda {0}, Seed {1}".format(
                                    entropy_targets[j], seeds[k]))
                            print("Path is", seed_path)
                        continue

                    entropies.append(entropy)
                    measures.append(measure)
                    errors.append(measure_error)
                    paths.append(seed_path)
                    if colour_metric == "Moment":
                        colour_metrics.append(second_moment/int(N))
                    elif colour_metric == "Assortativity":
                        colour_metrics.append(r)

                    # plot the individual run point
                    #ax[0].errorbar(entropy, measure, yerr=measure_error,
                    #        ls = "none", markersize = 8,
                    #        marker = "o", color = colour, capsize=3)

        # find the appropriate runs to show graphs/pkk
        entropies = pl.asarray(entropies)
        measures = pl.asarray(measures)
        plotting_weight = 0.99
        if i == merit_index:
            plotting_indices = []
            for pt in plotting_targets:
                merit = 1 - plotting_weight * pl.absolute(entropies  - pt) + (
                                    (1 - plotting_weight) * measures)
                index = pl.where(merit == pl.amax(merit))[0][0]
                plotting_indices.append(index)

            

        ranks = pl.searchsorted(pl.sort(colour_metrics), colour_metrics)
        colours = get_colours_from_cmap(colour_metrics)

        measure_fig, measure_ax = pl.subplots(1,2,figsize=(8, 6),
            gridspec_kw={"width_ratios": [20, 1]},)
        measure_fig.subplots_adjust(wspace=0.01)

        if i == 0: 
            if int(N) > 500:
                graph_fig, graph_ax = pl.subplots(2, len(plotting_indices),
                                    subplot_kw={"box_aspect": 1},
                                    constrained_layout=True, figsize = (8, 10))
            else:
                graph_fig, graph_ax = pl.subplots(3, len(plotting_indices),
                                    subplot_kw={"box_aspect": 1},
                                    constrained_layout=True, figsize = (12, 8))
            #gridspec_kw={"height_ratios": [1,1, 1]},
            #pl.subplots_adjust(hspace=0)

        for l, c_m in enumerate(colour_metrics):
            if l in plotting_indices:
                path = paths[l]
                
                if i == 0:
                    compare_graphs(l, plotting_indices, path, int(N), 
                                        graph_fig, graph_ax,
                                        entropies[l], plots_dir)
                marker = "*"
                marker_size = 15
                mec = "m"
                zorder=10
                mew=1.5
            else:
                marker = "o"
                marker_size = 8
                mec = "k"
                zorder=1
                mew=0.5

            rank = ranks[l]
            measure_ax[0].errorbar(entropies[l], measures[l], yerr=errors[l],
                    ls = "none", markersize = marker_size, mec = mec,
                    marker = marker, color = colours[l], capsize=3, mew=mew,
                    zorder=zorder)
            
        cmap = mpl.cm.get_cmap("viridis")#(x)#[np.newaxis, :, :3]
        norm = mpl.colors.Normalize(vmin=pl.amin(colour_metrics),
                                        vmax=pl.amax(colour_metrics))
        cb1 = mpl.colorbar.ColorbarBase(measure_ax[1], cmap=cmap,
                                        norm=norm,
                                        orientation='vertical')
        #cb1.set_label('k max')

        measure_ax[0].set_ylabel(measure_names[i], fontsize = fs)
        measure_ax[0].set_xlabel('Entropy', fontsize = fs)
        if colour_metric == "Moment":
            measure_ax[1].set_ylabel(r"Scaled Second Moment " + 
                r"$\frac{\langle k^{2} \rangle}{N}$", fontsize = fs)
        elif colour_metric == "Assortativity":
            measure_ax[1].set_ylabel('Degree Assortativity Coefficient')
        measure_ax[0].set_title(plot_title)
 
        if save_plots:
            measure_fig.savefig(plots_dir + plot_title.replace(" ","") +
                measure_names[i].replace(" ","") + "VsEntropy" + 
                colour_metric + "Colours.pdf")
            graph_fig.savefig(plots_dir + plot_title.replace(" ","") +
                "GraphComparison.pdf")

def network_metrics_vs_entropy(data_dir, plots_dir, plot_title, save_plots):
    """
    Plot various things measured from the network (final network) against
    entropy. Pull stuff from Newtorkx
    """
    lambda_folders = sorted([f for f in os.listdir(data_dir) if "DS_Store" not
                                                                        in f])
    lambdas = [float(f.strip("Lambda")) for f in lambda_folders]
    print(lambdas)
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]


    measure_names = ["Shortest Path Length", "Clustering Coefficient",
                        "Degree Assortativity Coefficient"]
    measure_functions = [nx.algorithms.average_shortest_path_length,
                            nx.algorithms.average_clustering,
                            nx.algorithms.degree_assortativity_coefficient]
    measure_names[0] = "Alpha"
    measure_functions[0] = "None"
    colours = ["C0", "C1", "C2"]
    markers = ['o', 's', '^']


    metric_means = [[] for i in range(len(measure_names))]
    metric_errors = [[] for i in range(len(measure_names))]
    entropy_means = []
    entropy_errors = []
    # loop over lambda values
    for j, path in enumerate(lambda_paths):

        seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                            in f]
        seeds = [f.strip("Seed") for f in seed_folders]
        seed_paths = [path + f + "/" for f in seed_folders]
        entropy_end_points = []
        metric_end_points = [[] for i in range(len(measure_names))]

        # loop over all random seeds for given lambda
        for k, seed_path in enumerate(seed_paths):
            
            # taking care of missing data
            try:
                entropy_data = pl.loadtxt(seed_path +
                                        "BestPkkEntropies.txt")
                entropy_end_points.append(entropy_data[-1])
                print("Loading Network For Lambda {0}, Seed {1}".format(
                        lambdas[j], seeds[k]))
                G = nx.read_edgelist(seed_path + "BestNetwork.csv",
                                        delimiter = ",")

                # measure everything for that network
                for i, m in enumerate(measure_names):
                    print("Calculating", m)
                    if m == "Alpha":
                        metric = pl.loadtxt(seed_path + "BestAlphas.txt")[-1]
                        metric_end_points[i].append(metric)
                    else:
                        metric = measure_functions[i]
                        metric_end_points[i].append(metric(G))

            except OSError:
                if j == 0:
                    print("No Data For Lambda {0}, Seed {1}".format(
                            lambdas[j], seeds[k]))
                    print("Path is", seed_path)
                continue
           
        entropy_means.append(pl.average(entropy_end_points))
        entropy_error = pl.std(entropy_end_points)/pl.sqrt(
                                    len(entropy_end_points) - 1)
        entropy_errors.append(entropy_error)
        for i in range(len(measure_names)):
            metric_data = metric_end_points[i]
            metric_means[i].append(pl.average(metric_data))
            metric_errors[i].append(pl.std(metric_data)/
                                        pl.sqrt(len(metric_data) - 1))

    fig, ax = pl.subplots(1, len(measure_names), figsize = (12, 4))
    for i, m in enumerate(measure_names):
        p = ax[i]
        means = metric_means[i]
        errors = metric_errors[i]
        #print(means, errors)
        #print(entropy_means, entropy_errors)
        p.errorbar(entropy_means, means, yerr=errors, xerr = entropy_errors,
                    capsize = 3, fmt = "ko")
        p.set_title(m)
        p.set_xlabel('Entropy')
 
    fig.suptitle(plot_title)
    if save_plots:
        pl.savefig(plots_dir + plot_title.replace(" ","") +
                    "NetworkMetricsVsEntropy.pdf")

def versus_lambda(data_dir, plots_dir, plot_title, save_plots):

    lambda_folders = sorted([f for f in os.listdir(data_dir) if "DS_Store" not
                                                                        in f])
    lambdas = [float(f.strip("Lambda")) for f in lambda_folders]
    print(lambdas)
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]

    colours = get_colours_from_cmap(pl.arange(len(lambdas)))

    measure_files = ["BestHealthMeans.txt"]*3 + ["BestPkkEntropies.txt",
                        "BestAlphas.txt", "BestProportions.txt"]
    measure_names = ["Death Age", "Healthy Aging", "QALY", "Pkk Entropy",
                        "Alpha", "Proportions"]

    markers = ['o', '*', '^']
    choices = ["Assortative", "Dissassortative", "Random"]

    fig, linear_axes = get_n_by_m_subplots(2,3)
    pl.subplots_adjust(right=0.87)

    for i, m in enumerate(measure_files):
        ax = linear_axes[i]
        means = []
        errors = []

        for j, path in enumerate(lambda_paths):
            colour = colours[j]
            seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                                in f]
            seeds = [f.strip("Seed") for f in seed_folders]
            seed_paths = [path + f + "/" for f in seed_folders]
            end_points = []
            for k, seed_path in enumerate(seed_paths):
                # taking care of missing data
                try:
                    # the health stuff is all in the same file, but in different
                    # columns
                    if i < 3:
                        measure = pl.loadtxt(seed_path+m)
                        if pl.shape(measure) == (0,):
                            print("No Data Here (Lambda = 1)")
                            end_points = [pl.nan, pl.nan]
                            break
                        else:
                            measure = measure[:,i]
                    else:
                        measure = pl.loadtxt(seed_path + m)
                except OSError:
                    if j == 0:
                        print("No Data For Lambda {0}, Seed {1}".format(
                                lambdas[j], seeds[k]))
                        print("Path is", seed_path)
                    continue
                if i == 5:
                    end_points.append([pl.average(measure[-5:,p]) for
                                                            p in range(3)])
                else:
                    end_points.append(measure[-1])
            if i == 5:
                mean = pl.average(end_points, axis = 0)
                error = pl.std(end_points, axis = 0)/pl.sqrt(
                                                        len(end_points) - 1)
                means.append(mean)
                errors.append(error)
                
            else:
                mean = pl.average(end_points)
                error = pl.std(end_points)/pl.sqrt(len(end_points) - 1)
                means.append(mean)
                errors.append(error)
        if i == 5:
            means = pl.asarray(means)
            errors = pl.asarray(errors)
            for p in range(3):
                
                if pl.shape(means)[0] != 3:
                    print("Hllo")
                ax.errorbar(lambdas, means[:,p], yerr=errors[:,p], ls = "none",
                        marker = markers[p], color = 'k', capsize=3,
                        label = choices[p])
            ax.legend()
        else:
            ax.errorbar(lambdas, means, yerr=errors, ls = "none", marker = "o",
                    color = 'k', capsize=3)#, label = "Lambda: {}".format(lambdas[j]))
        #ax.legend()
        ax.set_ylabel(measure_names[i])
        ax.set_xlabel('Lambda')
 
    fig.suptitle(plot_title)
    if save_plots:
        pl.savefig(plots_dir + plot_title.replace(" ","") + "VsLambda.pdf")

def show_degree_distributions(data_dir, plots_dir, plot_title, save_plots):

    lambda_folders = sorted([f for f in os.listdir(data_dir) if "DS_Store" not
                                                                        in f])
    lambdas = [float(f.strip("Lambda")) for f in lambda_folders]
    print(lambdas)
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]

    colours = get_colours_from_cmap(pl.arange(len(lambdas)))

    fig = pl.figure(figsize=(8,6))
    
    bins = pl.logspace(0, 12, 13, base = 2)
    degree_markers = bins[:-1]

    for j, path in enumerate(lambda_paths):
        colour = colours[j]
        seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                            in f]
        seeds = [f.strip("Seed") for f in seed_folders]
        seed_paths = [path + f + "/" for f in seed_folders]
        end_points = []
        for k, seed_path in enumerate(seed_paths):
            # taking care of missing data
            try:
                G = nx.read_edgelist(seed_path + "BestNetwork.csv",
                                        delimiter = ",")
                degrees = [G.degree(n) for n in G.nodes()]
                dk = pl.histogram(degrees, bins)[0]
                end_points.append(dk)

            except OSError:
                if j == 0:
                    print("No Data For Lambda {0}, Seed {1}".format(
                            lambdas[j], seeds[k]))
                continue

        means = pl.average(end_points, axis = 0)
        errors = pl.std(end_points, axis = 0)/pl.sqrt(len(end_points) - 1)
        print(len(degree_markers), len(means), len(errors))
        pl.errorbar(degree_markers, means, yerr=errors,
                ls = "none", marker = "o", color = colour, capsize=3,
                label = "Lambda: {}".format(lambdas[j]))
    pl.xlim(0.8, 1.2*pl.amax(degree_markers[means != 0]))
    pl.legend(loc = "upper left", bbox_to_anchor=(1,1), frameon=False)
    pl.yscale('log')
    pl.xscale('log')
    pl.xlabel('k')
    pl.ylabel('D(k) (Averaged over {} Runs)'.format(len(seeds)))
    pl.subplots_adjust(right=0.8) 
    fig.suptitle(plot_title)
    if save_plots:
        pl.savefig(plots_dir + plot_title.replace(" ","") +
            "DegreeDistributions.pdf")

def compare_degree_distributions(data_dir, plots_dir, plot_title, save_plots):

    lambda_folders = sorted([f for f in os.listdir(data_dir) if "DS_Store" not
                                                                        in f])
    N = 3000
    lambdas = [float(f.strip("Lambda")) for f in lambda_folders]
    print(lambdas)
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]

    colours = get_colours_from_cmap(pl.arange(len(lambdas)))

    fig = pl.figure(figsize=(8,6))
    
    bins = pl.arange(1500)
    degree_markers = bins
    for j, path in enumerate(lambda_paths):
        colour = colours[j]
        seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                            in f]
        seeds = [f.strip("Seed") for f in seed_folders]
        seed_paths = [path + f + "/" for f in seed_folders]
        end_points = []
        for k, seed_path in enumerate(seed_paths):
            # taking care of missing data
            try:
                degrees = pl.loadtxt(seed_path + "BestDegrees.txt")[-1,:]
                dk = pl.loadtxt(seed_path + "BestDks.txt")[-1,:]
                N = pl.sum(dk)
                padded_dk = pl.zeros_like(bins)
                padded_dk[degrees.astype(int)] = dk
                end_points.append(padded_dk)

            except OSError:
                if j == 0:
                    print("No Data For Lambda {0}, Seed {1}".format(
                            lambdas[j], seeds[k]))
                continue

        means = pl.average(end_points, axis = 0)
        errors = pl.std(end_points, axis = 0)/pl.sqrt(len(end_points) - 1)
        valid = (means != 0)
        pl.errorbar(degree_markers[valid], means[valid], yerr=errors[valid],
                ls = "none", marker = "o", color = colour, capsize=3,
                label = r"$\lambda: {}$".format(lambdas[j]))
    #pl.xlim(0.8, 1.2*pl.amax(degree_markers[means != 0]))
    pl.ylim(0.8/len(seeds), N)
    pl.yscale('log')
    pl.xscale('log')
    pl.xlabel('k')
    pl.ylabel('D(k) (Averaged over {} Runs)'.format(len(seeds)))

    pl.plot(degrees,2*N*degrees**(-2), "C3:", label = r"$\alpha = 2$")
    pl.plot(degrees,2*N*degrees**(-2.27), "C3-", label = r"$\alpha = 2.27$")
    pl.plot(degrees,2*N*degrees**(-2.5), "C3--", label = r"$\alpha = 2.5$")
    pl.legend(loc = "upper left", bbox_to_anchor=(1,1), frameon=False)

    pl.subplots_adjust(right=0.8) 
    fig.suptitle(plot_title)
    if save_plots:
        pl.savefig(plots_dir + plot_title.replace(" ","") +
            "DegreeDistributions.pdf")

def assortativity_extremes(G):
    """
    in: a graph (nx undirected graph)
    out: graphs rewired for the assortative, random, and dissassortative
    extremes (via pa, pd, pr)
    """
    degrees = [G.degree(node) for node in G.nodes()]
    bins = sorted(list(set(degrees))) + [1000000]
    dk = pl.histogram(degrees, bins = bins)[0]
    degrees = pl.asarray(bins[:-1])

    graphs = []
    pkks = []
    for i in range(3):
        p_a = float(i == 0)
        p_d = float(i == 1)
        p_r = float(i == 2)
        print("Building: pa:{0}, pd:{1}, pr:{2}".format(p_a, p_d, p_r))
        G, jkk = build_weighted_graph(degrees, dk, p_a, p_d, p_r, True)
        if type(G) != str:
            pkk = jkk/pl.sum(jkk)
        else:
            pkk = "bubcus"
        graphs.append(G)
        pkks.append(pkk)

    return graphs, pkks

def assortativity_extremes_vs_entropy(data_dir, plots_dir, plot_title,
        params, param_names, running_folder, save_plots):
    """
    Goes through all the lambdas and seeds for the given simulation batch
    in: Directories to find data and save data + plot title stuff.
    out: Saves plots (or just displays them), no return.
    """

    lambda_folders = sorted([f for f in os.listdir(data_dir) if "DS_Store" not
                                                                        in f])
    lambdas = [float(f.strip("Lambda")) for f in lambda_folders]
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]

    measures = ["DeathAge", "HealthyAging", "QALY"]
    measure_files = ["BestHealthMeans.txt"]*3 + ["BestPkkEntropies.txt",
                        "None"]
    measure_names = ["Death Age", "Healthy Aging", "QALY", "Pkk Entropy",
                        "Degree Assortativity Coefficient"]

    # stupid renaming function (bad practice??)
    calculate_r = nx.algorithms.degree_assortativity_coefficient
    colours = ["C0", "C1", "C2"]
    markers = ['o', 's', '^']
    assortativity_labels = ["Assortative", "Dissassortative", "Random"]


    metric_means = [[] for i in range(len(measure_names))]
    metric_errors = [[] for i in range(len(measure_names))]
    entropy_means = []
    entropy_errors = []
    rewired_means = []#[[[] for j in range(3)] for i in range(len(measure_names))]
    rewired_errors = []#[[] for j in range(3)] for i in range(len(measure_names))]

    # loop over lambda values
    for j, path in enumerate(lambda_paths):

        seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                            in f]
        seeds = [f.strip("Seed") for f in seed_folders]
        seed_paths = [path + f + "/" for f in seed_folders]
        entropy_end_points = []
        metric_end_points = [[] for i in range(len(measure_names))]
        
        # for each graph type, and each metric, get measurements from each seed
        # (3 graph types, ?? metrics ~5, ?? seeds ~10)
        rewired_metrics = [[[] for i in measure_names] for j in range(3)]

        # loop over all random seeds for given lambda
        for k, seed_path in enumerate(seed_paths):
            
            # taking care of missing data
            try:
                # Get the original data for that network
                entropy_data = pl.loadtxt(seed_path +
                                        "BestPkkEntropies.txt")
                entropy_end_points.append(entropy_data[-1])
                print("Loading Network For Lambda {0}, Seed {1}".format(
                        lambdas[j], seeds[k]))
                G = nx.read_edgelist(seed_path + "BestNetwork.csv",
                                        delimiter = ",")

                # measure everything for that network
                for i, m in enumerate(measure_names):
                    if lambdas[j] == 1:
                        metric_end_points[i].append(pl.nan)
                        continue
                    if i < 4:
                        metric = pl.loadtxt(seed_path + measure_files[i])[-1]
                        if i < 3:
                            if i == 1:
                                metric_end_points[i].append(metric[i]*100)
                            else:
                                metric_end_points[i].append(metric[i])
                        else:
                            metric_end_points[i].append(metric)
                    else:
                        metric_end_points[i].append(calculate_r(G))


                ##### Now: Rewire that sucker and measure stuff
                print("Rewiring graphs")
                rewired_graphs, rewired_pkks = assortativity_extremes(G)
                print("Running Rewired Graphs")
                for i, G in enumerate(rewired_graphs):
                    if type(G) == str:
                        continue
                    # test each graph in GNM 
                    health_means, health_errors = run_network(G,
                                running_folder, params, param_names, measures)
                    # append health means
                    for n, h_m in enumerate(health_means):
                        if n != 1:
                            rewired_metrics[i][n].append(h_m)#*scale)
                        else:
                            rewired_metrics[i][n].append(h_m*100)

                    # also need entropy
                    rewired_metrics[i][3].append(
                                    calculate_entropy(rewired_pkks[i]))
                    # and assortativity coeff.
                    rewired_metrics[i][4].append(calculate_r(G))
            except OSError:
                if j == 0:
                    print("No Data For Lambda {0}, Seed {1}".format(
                            lambdas[j], seeds[k]))
                    print("Path is", seed_path)
                continue
        entropy_means.append(pl.average(entropy_end_points))
        entropy_error = pl.std(entropy_end_points)/pl.sqrt(
                                    len(entropy_end_points) - 1)
        entropy_errors.append(entropy_error)

        for i in range(len(measure_names)):
            # Data from original runs
            metric_data = metric_end_points[i]
            metric_means[i].append(pl.average(metric_data))
            metric_errors[i].append(pl.std(metric_data)/
                                        pl.sqrt(len(metric_data) - 1))
            
        # Rewired network data
        graph_means = []
        graph_errors = []
        for n in range(3):
            measure_means = []
            measure_errors = []
            for i, measure_name in enumerate(measure_names):
                rewired_data = rewired_metrics[n][i]
                measure_means.append(pl.average(rewired_data))
                measure_errors.append(pl.std(rewired_data)/
                                        pl.sqrt(len(rewired_data) - 1))
            graph_means.append(measure_means)
            graph_errors.append(measure_errors)
        rewired_means.append(graph_means)
        rewired_errors.append(graph_errors)
    fig, ax = pl.subplots(1, 3, figsize = (16, 6))
    #print(rewired_means)
    for i, m in enumerate(measure_names):
        if i < 3:
            p = ax[0]
        else:
            p = ax[i-2]
        #print(means, errors)
        #print(entropy_means, entropy_errors)
        p.set_title(m)
        p.set_xlabel('Entropy')

        # add in the rewired stuff
        for n in range(3):
            means = [rewired_means[l][n][i] for l in range(len(lambdas))]
            errors = [rewired_errors[l][n][i] for l in range(len(lambdas))]
            #print(means)
            p.errorbar(entropy_means, means, yerr=errors, xerr = entropy_errors,
                        capsize = 3, fmt = "{0}o".format(colours[n]),
                        label = "{0}".format(assortativity_labels[n]))

        # plop the original stuff on top
        means = metric_means[i]
        errors = metric_errors[i]

        p.errorbar(entropy_means, means, yerr=errors, xerr = entropy_errors,
                    capsize = 3, fmt = "ko", label="Optimal")


    pl.subplots_adjust(right=0.86)
    p.legend(loc="upper left", bbox_to_anchor=(0.98,1))


 
    fig.suptitle(plot_title)
    if save_plots:
        pl.savefig(plots_dir + plot_title.replace(" ","") +
                    "RewiredNetworksVsEntropy.pdf")

def get_max_merit_curve(measure, entropy, lambda_values):
    entropy_contribution = pl.outer(lambda_values, entropy)
    health_contribution = pl.outer((1 - lambda_values), measure)
    max_merits = pl.amax(entropy_contribution + health_contribution, axis = 1)
    return max_merits

def plot_max_merit(data_dir, plots_dir, plot_title, save_plots,
                    optimization_measure="DeathAge"):
    if optimization_measure == "DeathAge":
        merit_index = 0
    elif optimization_measure == "QALY":
        merit_index = 1
    else:
        merit_index = 2

    # Resolution of lambda_points
    x = pl.linspace(0, 1, 101)

    lambda_folders = sorted([f for f in os.listdir(data_dir) if "DS_Store" not
                                                                        in f])
    lambdas = [float(f.strip("Lambda")) for f in lambda_folders]
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]


    measure_files = ["BestHealthMeans.txt"]*3
    measure_names = ["Death Age", "Healthy Aging", "QALY"]
    colours = get_colours_from_cmap(lambdas)
    markers = ['o', 's', '^']


    max_measures = []
    max_entropies = []
    for i, m in enumerate(measure_files):
        best_measure = 0 
        best_entropy = 0
        for j, path in enumerate(lambda_paths):
            seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                                in f]
            seeds = [f.strip("Seed") for f in seed_folders]
            seed_paths = [path + f + "/" for f in seed_folders]
            for k, seed_path in enumerate(seed_paths):
                # taking care of missing data
                try:
                    # get entropy stuff (only need to do this once)
                    entropy_data = pl.loadtxt(seed_path +
                                            "BestPkkEntropies.txt")
                    entropy = pl.amax(entropy_data)
                    if entropy > best_entropy:
                        best_entropy = entropy
                    # the health stuff is all in the same file, but in different
                    # columns
                    measure_data = pl.loadtxt(seed_path+m)

                    # Case of just no successful iterations (unlikely)
                    if pl.shape(measure_data) == (0,):
                        break
                    else:
                        measure = pl.amax(measure_data[:,i])
                        if "Healthy" in measure_names[i]:
                            measure *= 100
                        if measure > best_measure:
                            best_measure = measure

                # Error case: Cant find the data
                except OSError:
                    if j == 0:
                        print("No Data For Lambda {0}, Seed {1}".format(
                                lambdas[j], seeds[k]))
                        print("Path is", seed_path)
                    continue
        max_measures.append(best_measure)
        max_entropies.append(best_entropy)


    entropy_means = []
    entropy_errors = []
    combined_max_merit_curves = []
    for i, m in enumerate(measure_files):
        pl.figure(figsize = (8,6))
        max_measure = max_measures[i]
        max_entropy = max_entropies[i]
        means = []
        errors = []
        maxes = []
        max_errors = []
        best_entropies = []

        for j, path in enumerate(lambda_paths):
            colour = colours[j]
            seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                                in f]
            seeds = [f.strip("Seed") for f in seed_folders]
            seed_paths = [path + f + "/" for f in seed_folders]
            end_points = []
            entropy_end_points = []
            entropy_weight = lambdas[j]
            if entropy_weight == 1.0:
                continue
            best_merit = 0
            best_measure = 0
            best_error = 0
            best_entropy = 0
            best_merit_curves = []
            for k, seed_path in enumerate(seed_paths):
                
                # taking care of missing data
                try:
                    # get entropy stuff (only need to do this once)
                    entropy_data = pl.loadtxt(seed_path +
                                            "BestPkkEntropies.txt")[:]
                    entropy_end_points.append(entropy_data[-1])
                    # the health stuff is all in the same file, but in different
                    # columns
                    measure_data = pl.loadtxt(seed_path+m)[:,i]
                    measure_error_data = pl.loadtxt(seed_path +
                        "BestHealthErrors.txt")[:,i]

                    # Case of just no successful iterations (unlikely)
                    if pl.shape(measure_data) == (0,):
                        break
                    else:
                        measure = measure_data
                        measure_error = measure_error_data
                        if "Healthy" in measure_names[i]:
                            measure *= 100
                            measure_error *= 100
                        measure /= max_measure 
                        entropy_data /= max_entropy

                # Error case: Cant find the data
                except OSError:
                    if j == 0:
                        print("No Data For Lambda {0}, Seed {1}".format(
                                lambdas[j], seeds[k]))
                        print("Path is", seed_path)
                    continue
                # calculate merit curve
                best_merit_curve = get_max_merit_curve(measure, entropy_data,
                                                                            x)
                best_merit_curves.append(best_merit_curve)

                end_points.append(measure[-1])
                # plot the individual run point
                #pl.errorbar(entropy_data[-1], measure, yerr=measure_error,
                #        ls = "none", markersize = 5,
                #        marker = markers[i], color = colours[i], capsize=3,
                #        alpha = 0.5)
            combined_max_curve = pl.amax(best_merit_curves, axis = 0)
            pl.plot(x, combined_max_curve, c = colour)
               
        pl.ylabel(measure_names[i] +"-Based Merit, Scaled by Best Performance")
        pl.xlabel('Lambda')
 
        #pl.title(plot_title)
        if save_plots:
            pl.savefig(plots_dir + plot_title.replace(" ","") +
                "{}MaxMeritCurves.pdf".format(
                    measure_names[i].replace(" ","")))

def draw_networks(data_dir, plots_dir, plot_title, entropy_indices, num_seeds,
            save_plots = 0):

    lambda_folders = sorted([f for f in os.listdir(data_dir) if "DS_Store" not
                                                                        in f])
    lambdas = [float(f.strip("Lambda")) for f in lambda_folders]
    lambda_paths = [data_dir + f + "/" for f in lambda_folders]

    lambdas = [lambdas[i] for i in entropy_indices]
    lambda_paths = [lambda_paths[i] for i in entropy_indices]

    print(lambdas)
    for j, path in enumerate(lambda_paths):
        seed_folders = [f for f in os.listdir(path) if "DS_Store" not
                                                            in f]
        seeds = [f.strip("Seed") for f in seed_folders]
        seed_paths = [path + f + "/" for f in seed_folders]
        for k, seed_path in enumerate(seed_paths[:num_seeds]):
            # taking care of missing data
            try:
                G = nx.read_edgelist(seed_path + "BestNetwork.csv",
                                        delimiter = ",")
            except OSError:
                if j == 0:
                    print("No Data For Lambda {0}, Seed {1}".format(
                            lambdas[j], seeds[k]))
                continue

            #pl.figure(figsize = (8,6))
            #pos = nx.spring_layout(G)
            #nx.draw(G, pos, node_size = 8, edge_color = 'C7')
            draw_spring_graph(G)
            pl.title("Lambda:{0}, Seed {1}".format(lambdas[j], seeds[k]))

    if save_plots:
        pl.savefig(plots_dir + plot_title.replace(" ","") +
            "DegreeDistributions.pdf")

def calc_mortality_rate(tds, bin_years = 5):
    n = len(tds)
    max_age = int(pl.round_(pl.amax(tds)))
    years = pl.arange(0, max_age + bin_years, bin_years)
    mortality_rate = []
    for y in years:
        deaths = pl.sum((tds >= y) & (tds < y+bin_years))
        count = pl.sum(tds >= y)
        if deaths > 0:
            mortality_rate.append(deaths/count)
        else:
            mortality_rate.append(pl.nan)
    return mortality_rate, years

def avg_ttd_fi(fis, qtd = 30):
    included = pl.sum(~pl.isnan(fis), axis = 1) > qtd
    num_included = pl.sum(included)
    fis = fis[included, :]
    ttd_fis = pl.zeros([num_included, qtd])
    for i in range(num_included):
        # grab the last 30 non-nan values from the fi data
        ttd_fis[i,:] = fis[i,~pl.isnan(fis[i,:])][-qtd:]
    
    # convert quarters to years
    years = (pl.arange(qtd) - qtd)/4
    means = pl.average(ttd_fis, axis = 0)
    errors = pl.std(ttd_fis, axis = 0)/pl.sqrt(pl.sum(included) - 1)
    return means, errors, years

def rerun_network(G, running_folder, params, param_names,
                                        plot_dir, plot_title, number = 10000):
    """
    Run the network and show FI vs t, mortality rate, etc
    """
    #G = nx.read_edgelist(graph_file, delimiter = ",")
    params[0] = str(number)
    death_ages, FIs = run_network(G, running_folder, params, param_names,
                                    ["DeathAges", "PopulationFI"],
                                    output = "array")[0]
    mean_fis = pl.nanmean(FIs, axis = 0)
    counts = pl.sum(~pl.isnan(FIs), axis = 0)
    fi_errors = pl.nanstd(FIs, axis = 0)/pl.sqrt(counts-1)
    #mean_fis = mean_fis[counts > 100]
    #fi_errors = fi_errors[counts > 100]
    years = pl.arange(len(mean_fis))

    fig, ax = pl.subplots(1, 3, figsize = (14,5))
    ax[0].plot(years, mean_fis)
    #for i in range(number):
    #    ax[0].plot(years, FIs[i,:], alpha = 0.5, c = 'C3')#mean_fis)
    upper = mean_fis + fi_errors
    lower = mean_fis - fi_errors
    ax[0].fill_between(years, lower, upper, alpha = 0.5)

    mortality_rate, years = calc_mortality_rate(death_ages)
    ax[1].plot(years, mortality_rate, 'C0o', label = "Optimized Network")
    ax[1].set_yscale('log')

    ttd_means, ttd_errors, ttd_years = avg_ttd_fi(FIs, 120)
    upper = ttd_means + ttd_errors 
    lower = ttd_means - ttd_errors 
    ax[2].fill_between(ttd_years, lower, upper, alpha = 0.5)
    ax[2].plot(ttd_years, ttd_means)


    # plop on the reference data
    old_fi_data = pl.loadtxt("ReferenceData/OldMortMeanFI.txt")
    old_mort_data = pl.loadtxt("ReferenceData/OldModelMortalityRates.txt")
    old_ttd_fi_data = pl.loadtxt("ReferenceData/OldModelTTDFIs.txt")

    old_years = old_fi_data[:,0]
    old_mean = old_fi_data[:,1]
    old_error = old_fi_data[:,2]
    upper = old_mean + old_error
    lower = old_mean - old_error
    ax[0].plot(old_years, old_mean, c = "C3")
    ax[0].fill_between(old_years, lower, upper, alpha = 0.5, color="C3")

    old_years = old_mort_data[:,0]
    old_mort_rate = old_mort_data[:,1]
    ax[1].plot(old_years, old_mort_rate, "C3o", label = "Published GNM")

    old_years = old_ttd_fi_data[:,0]
    old_mean = old_ttd_fi_data[:,1]
    old_error = old_ttd_fi_data[:,2]
    upper = old_mean + old_error
    lower = old_mean - old_error
    ax[2].plot(old_years, old_mean, c = "C3")
    ax[2].fill_between(old_years, lower, upper, alpha = 0.5, color="C3")

    ax[0].set_ylabel('FI')
    ax[0].set_xlabel('Age')

    ax[1].set_ylabel('Mortality Rate', labelpad=-2)
    ax[1].set_xlabel('Age')
    ax[1].legend()

    ax[2].set_ylabel('FI')
    ax[2].set_xlabel('Years Until Mortality')
    ax[2].set_xticklabels(pl.arange(0, 36,5)[::-1])
    ax[2].set_ylim((-0.05, 1.02))

    fig.suptitle(plot_title)
    pl.savefig(plot_dir + plot_title.replace(" ", "") + ".pdf")

def simple_entropy_plots(data_dirs, labels, plots_dir, save_plots, N,
                        temp_folder, params, param_names,
                        optimization_measures,
                        add_scale_free=True, tolerance = 0.05,
                        add_alpha = 1, add_metrics=0):
    merit_indices = []
    for optimization_measure in optimization_measures:
        if optimization_measure == "DeathAge":
            merit_indices.append(0)
        elif optimization_measure == "QALY":
            merit_indices.append(2)
        else:
            merit_indices.append(1)
    print(data_dirs, optimization_measures, merit_indices)
    entropy_figs = [pl.subplots(1,1, figsize=(8,6)) for i in range(3)]
    figs = [f[0] for f in entropy_figs]
    axes = [f[1] for f in entropy_figs]
    denom = 1.2
    comp_fig, comp_ax = pl.subplots(1, 1, figsize=(4/denom,3/denom))
    comp_fig.subplots_adjust(right=0.99, top=0.99, left=0.13, bottom=0.15)
    if add_alpha:
        denom = 1.5
        alpha_fig, alpha_ax = pl.subplots(1, 1, figsize=(4/denom,3/denom)) 
        alpha_fig.subplots_adjust(right=0.99, top=0.99, left=0.21, bottom=0.19)
    if add_metrics:
        metric_fig, metric_ax = pl.subplots(1, 4, figsize=(12,3))
        c_bar_ax = metric_fig.add_axes([0.91, 0.15, 0.02, 0.73])
        #metric_fig.set(tight_layout=True)
        metric_fig.subplots_adjust(wspace = 0.33, bottom=0.15, left=0.08)


    measure_names = ["Lifespan", "Healthy Aging", "QALY"]
    colours = ["C{}".format(i) for i in range(len(data_dirs))]
    markers = ["o", "^", "*", "d"]
    print("Number of Data:", len(data_dirs))
    print("Markers:", len(markers), ", Colours:", len(colours))
    combined_colour_metrics = []

    sf_entropy, sf_measures, sf_errors, sf_metrics = measure_spencers_graph(
                                    N, temp_folder, params, param_names)
    er_entropy, er_measures, er_errors, er_metrics = measure_spencers_graph(
                                    N, temp_folder,params, param_names,
                                    graph_type="Random")
    print(data_dirs)
    drew_high = 0
    drew_low = 0
    for i, data_dir in enumerate(data_dirs):
        merit_index = merit_indices[i]
        measures, errors, entropies, moments, rs, k_maxes, paths = (
                                                scrape_data(data_dir, N))[:-1]

        merit_health = measures[merit_index]
        merit_error = pl.asarray(errors[merit_index])
        lambdas = pl.asarray([0]) #pl.linspace(0, 1, 101)
        merits = pl.outer(lambdas, entropies) + pl.outer(1-lambdas, merit_health)
        max_merits = pl.amax(merits, axis = 0)
        entropies = pl.asarray(entropies)
        sorted_entropies = pl.sort(pl.asarray(list(set(entropies))))
        sorted_merits = []
        for e in sorted_entropies:
            max_merit = pl.amax(max_merits[entropies == e])
            sorted_merits.append(max_merit)
        sorted_merits = pl.asarray(sorted_merits)
        
        j = 0
        max_slope = 0
        best_merits = []
        best_entropies = []
        best_slopes = []
        while j < len(sorted_entropies):
            best_merits.append(sorted_merits[j])
            best_entropies.append(sorted_entropies[j])
            best_slopes.append(max_slope)
            max_slope = -1e9
            k = j + 1
            if k >= len(sorted_entropies):
                break
            lower_merit = sorted_merits[j]
            lower_entropy = sorted_entropies[j]
            while k < len(sorted_entropies):
                slope = (sorted_merits[k] - lower_merit)/(
                            sorted_entropies[k] - lower_entropy)
                if slope > max_slope:
                    max_slope = slope
                    best_k = k
                k += 1
            j = best_k
                
        best_merits = pl.asarray(best_merits)
        best_entropies = pl.asarray(best_entropies)
        peak_merit = pl.amax(best_merits)
        entropy_at_peak = best_entropies[
                        pl.where(best_merits == peak_merit)[0][0]]
        best_merits[best_entropies < entropy_at_peak] = peak_merit
        best_slopes = pl.asarray(best_slopes)
        best_slopes[best_entropies < entropy_at_peak] = 0.0


        # DRAWING THE CONVEX HULL PLOT  
        #axes[merit_index].plot(best_entropies, best_merits, "C2", lw = 2,
        #                            label = "Interpolated Max")


        # find distances to these lines
        distances = pl.zeros(len(entropies))
        for j, e in enumerate(entropies):
            merit = max_merits[j]
            best_index = pl.sum(best_entropies < e)
            best_merit = best_merits[best_index] + best_slopes[best_index] * (
                                e - best_entropies[best_index])
            distances[j] = (best_merit - merit)/best_merit

        failures = (distances > tolerance)
        print("Failure Rate:", pl.average(failures))

        # plot alphas if you want
        if add_alpha:
            paths = pl.asarray(paths)[~failures]
            alphas = []
            for p in paths:
                alpha_data = pl.loadtxt(p + "BestAlphas.txt")
                if "NonParametric" in p:
                    alphas.append(pl.nan)
                    continue
                if type(alpha_data) == np.ndarray:
                    try:
                        alphas.append(alpha_data[-1])
                    except IndexError:
                        print("Only One??")
                        alphas.append(alpha_data)
                elif type(alpha_data) == float:
                    alphas.append(alpha_data)
                else:
                    print("Idk what's up", alpha_data)

            """
            alphas = [pl.loadtxt(p + "BestAlphas.txt")[-1] for p in paths
                        if pl.shape(pl.loadtxt(p+"BestAlphas.txt")) != (1,)
                        and  pl.shape(pl.loadtxt(p+"BestAlphas.txt")) != (0,)]
            """
            #alphas = [a - 2 for a in alphas]
            alpha_ax.plot(entropies[~failures], alphas,
                    "{0}{1}".format(colours[i], markers[i]), label = labels[i])

        # plot networks if wanted (will need to calculate clustering coeff)
        if add_metrics:
            # clustering coefficient (cc)
            colour_metric = measures[merit_index]
            colour_metric = pl.asarray(colour_metric)[~failures]
            combined_colour_metrics += list(colour_metric)
            metric_colours = get_colours_from_cmap(colour_metric)
            if not add_alpha:
                paths = pl.asarray(paths)[~failures]
            ccs = []
            for path_num, p in enumerate(paths):
                G = nx.read_edgelist(p + "BestNetwork.csv",
                                        delimiter = ",")
                cc = nx.algorithms.average_clustering(G)
                ccs.append(cc)

            # degree assortativity coefficient
            rs = pl.asarray(rs)[~failures]
            marker_size = 35
            metric_ax[0].scatter(entropies[~failures], rs, c = metric_colours,
                        marker =  markers[i], label = labels[i],
                        ec="k", lw=1,
                        s=marker_size)
            metric_ax[1].scatter(entropies[~failures], k_maxes[~failures],
                        c = metric_colours,
                        marker =  markers[i],# label = labels[i],
                        ec="k", lw=1,
                        s=marker_size)
            metric_ax[2].scatter(entropies[~failures], moments[~failures],
                        c = metric_colours,
                        marker =  markers[i],# label = labels[i],
                        ec="k", lw=1,
                        s=marker_size)
            metric_ax[3].scatter(entropies[~failures], ccs, c = metric_colours,
                        marker =  markers[i],# label = labels[i],
                        ec = "k", lw=1,
                        s=marker_size)





        # plot the Stuff vs entropy for all measures
        for j, measure_name in enumerate(measure_names):
            measure = pl.asarray(measures[j])
            error = pl.asarray(errors[j])
            if j == 0:
                # limits for scale free plots
                #axes[j].set_ylim(81, 114)
                #axes[j].set_xlim(0.5, 5.5)
                # limits for qaly vs td
                axes[j].set_ylim(79, 114)
                axes[j].set_xlim(0.1, 4.6)

            axes[j].errorbar(entropies[~failures], measure[~failures],
                        yerr = error[~failures], marker = markers[i],
                        color = colours[i], ls = "none",
                        capsize = 3, label = labels[i])

        # plot QALY vs Death Age:
        if 1:
            tds = pl.asarray(measures[0])[~failures]
            qalys = pl.asarray(measures[2])[~failures]
            td_errors = pl.asarray(errors[0])[~failures]
            qaly_errors = pl.asarray(errors[2])[~failures]
            comp_ax.errorbar(tds, qalys, yerr = qaly_errors, xerr = td_errors,
                            marker = markers[i],
                            color = colours[i], ls = "none",
                            capsize = 3, label = labels[i])

    # adding scale free point to the QALY vs DeathAge Stuff
    if add_scale_free:
        comp_ax.errorbar(sf_measures[0], sf_measures[2],
                    yerr = sf_errors[2], xerr = sf_errors[0],
                    fmt = 'C6s', label = "Original GNM", capsize=3,
                    ms = erms - 2)
        comp_ax.errorbar(er_measures[0], er_measures[2],
                    yerr = er_errors[2], xerr = er_errors[0],
                    fmt = 'C4*', label = "Random Graph", capsize=3,
                    ms = erms)

        if add_alpha:
            alpha_ax.plot(sf_entropy, 2.35, "C6s", label = "Original GNM",
                            ms = erms-2)

    if add_alpha:
        alpha_ax.set_ylabel('Scale Free Exponent', labelpad=-0.15)
        alpha_ax.set_xlabel('Entropy', labelpad=-0.1)
        #alpha_ax.set_yscale('log')
        if save_plots:
            alpha_fig.savefig(plots_dir + "AlphasTol{0}.pdf".format(tolerance))
                                                
    if add_metrics:
        cmap = mpl.cm.get_cmap("viridis")#(x)#[np.newaxis, :, :3]
        norm = mpl.colors.Normalize(vmin=pl.amin(combined_colour_metrics),
                                        vmax=pl.amax(combined_colour_metrics))
        cb1 = mpl.colorbar.ColorbarBase(c_bar_ax, cmap=cmap,
                                        norm=norm,
                                        orientation='vertical')
        #metric_fig.tight_layout(w_pad=-10)
        c_bar_ax.yaxis.set_ticks_position('right')
        c_bar_ax.yaxis.set_label_position('right')
        c_bar_ax.set_ylabel('Average Lifespan', labelpad=2)
        #cbar = metric_fig.colorbar(
 
        metric_ax[0].set_ylabel('Degree Assortativity Coefficient', labelpad=0)
        metric_ax[0].set_xlabel('Entropy')
        metric_ax[3].set_ylabel('Average Clustering Coefficient', labelpad=0)
        metric_ax[1].set_xlabel('Entropy')
        metric_ax[2].set_xlabel("Entropy")
        metric_ax[1].set_ylabel("Maximum Degree (Units of N)", labelpad=0)
        metric_ax[3].set_xlabel("Entropy")
        metric_ax[2].set_ylabel(r"Scaled Second Moment " + 
                r"$\frac{\langle k^{2} \rangle}{N}$", labelpad=0)
        reordered_ids = [0,2,3,1]
        letters = ['A', 'B', 'C', 'D']
        # Chaos: I've reordered the Ids like a gd fool
        for m_id in reordered_ids: #m_id in range(len(sf_metrics)):
            metric_p = metric_ax[m_id]#int(m_id > 1), int(m_id % 2)]
            metric_p.annotate(letters[m_id] + ")", xy = (0.9, 0.92),
                                xycoords="axes fraction", fontsize=12)
            if m_id == 0:
                metric_p.plot(sf_entropy, sf_metrics[i], 'C6s',
                                label = "Original GNM", ms = erms - 2)
                metric_p.plot(er_entropy, er_metrics[i], 'C4*',
                                ms = erms, label = "Random Graph")
            else:
                metric_p.plot(sf_entropy, sf_metrics[i], 'C6s',
                                ms = erms - 2)
                metric_p.plot(er_entropy, er_metrics[i], 'C4*',
                                ms = erms)

        #    if m_id == 0:
        metric_fig.legend(loc="upper center", ncol=4)
        #                        mode="expand")
        #metric_ax[0,0].set_ylim(-1.08, 0.3)
        #metric_ax[0,0].set_xlim(0.4, 5.57)
        

        #alpha_ax.set_yscale('log')
        if save_plots:
            metric_fig.savefig(plots_dir + "NetworkMetricsTol{0}.pdf".format(
                                tolerance))
 


    comp_ax.set_xlabel('Lifespan', labelpad=0.1)
    comp_ax.set_ylabel("QALY", labelpad=0.1)
    if save_plots:
        comp_fig.savefig(plots_dir + "QALYvsTdTol{0}.pdf".format(tolerance,
                                            measure_name.replace(" ","")))


    for j, measure_name in enumerate(measure_names):
        p = axes[j]
        p.set_xlabel("Entropy")
        p.set_ylabel(measure_name)
        if add_scale_free:
            p.errorbar(sf_entropy, sf_measures[j], yerr = sf_errors[j],
                        fmt = 'C6s', label = "Original GNM", capsize=3,
                        ms = erms-2)
            p.errorbar(er_entropy, er_measures[j], yerr = er_errors[j],
                        fmt = 'C4*', label = "Random Graph", capsize=3,
                        ms = erms)

        
        
        ### LEGEND FOR TD VS ENTROPY
        ### for the td qaly plot
        p.legend(loc='upper right')
        ### for the scale free plot
        #p.legend(loc='lower left')

        if save_plots:
            make_directory(plots_dir)
            figs[j].savefig(plots_dir + "{1}Tol{0}.pdf".format(tolerance,
                                            measure_name.replace(" ","")))

def measure_spencers_graph(N, temp_folder, params, param_names,
                            graph_type = "ScaleFree"):
    """
    Use Spencer's code to get the degree sequence for scale free network of
    given parameters
    """
    # params for generating network
    #sub.call(["make", "backup"])
    graph_params = ["0", "0.0", "0.0", "0.0", str(N), "2", "2.35", "4",
                "AND", "Single", temp_folder, graph_type, "0",
                "1", "0.0"]

    command = ['./main'] + graph_params 
    sub.call(command)

    output_files = os.listdir(temp_folder)
    #print(output_files)
    degree_file = [f for f in output_files if "Initial" in f][0]
    edge_list_file = [f for f in output_files if "Edge" in f][0]
    try:
        degree_sequence = pl.loadtxt(temp_folder + degree_file)[:,1]
    except IndexError:
        print(degree_file)
    G = nx.read_edgelist(temp_folder + edge_list_file)
    r = nx.degree_assortativity_coefficient(G)
    cc = nx.average_clustering(G)
    degrees = pl.asarray(list(set([G.degree(n) for n in G.nodes()])))
    pk = pl.asarray(nx.degree_histogram(G))
    pk = pk[pk != 0]
    scaled_second_moment = pl.dot(pk, degrees**2)/int(N)
    metrics = [r, cc, pl.amax(degrees), scaled_second_moment]
    


    pkk = pkk_from_G(G)
    entropy = calculate_entropy(pkk)
    measures = ["DeathAge", "HealthyAging", "QALY"]
    for f in output_files:
        os.remove(temp_folder + f)

    health_means, health_errors = run_network(G,
                    temp_folder, params, param_names, measures)

    output_files = os.listdir(temp_folder)
    for f in output_files:
        os.remove(temp_folder + f)


    return entropy, health_means, health_errors, metrics
    

if __name__ == '__main__':
    pl.close('all')
    
    base_data_dir = "ClusterData/" #"ClusterData/"
    base_plots_dir = "Plots/SimplePlots/"

    n_bins = "15"
    method = ["NonParametric", "Variational"][0]
    merit = "DeathAge"
    add_metrics=False
    add_alpha = True

    k_min = "2"
    k_max = "1.0"
    N = "128"
    meta_data = "betterkMax2/NBins{}".format(n_bins) 

    k_mins = ["1", "2"][:1]
    k_maxes = ["0.5", "1.0"][1:]
    data_dirs = []
    labels = []
    k_min_k_max = 0
    if k_min_k_max:
        add_metrics = True
        add_alpha = False
        for k_min in k_mins:
            for k_max in k_maxes:
                details = "{0}/{1}/kMin{4}/kMax{5}/N{2}/{3}/".format(meta_data,
                                                    method, N, merit, k_min, k_max)
                data_dir = base_data_dir + details
                data_dirs.append(data_dir)
                labels.append("kMin{0}kMax{1}".format(k_min, k_max))
            brief_description = "kMax0.5/"

    qaly_vs_td = 0
    if qaly_vs_td:
        add_alpha = False
        brief_description = "BigFinalQALYvsDeathAgekMin2/"
        meta_data = "FinalN1024"#"betterkMin2nBins15"
        add_metrics=True
        N = "1024"
        k_min = "2"
        tolerance = 0.05
        k_max = "1.0"
        n_bins = "15"
        optimization_measures = ["DeathAge", "QALY"]
        merit_labels = ["Lifespan", "QALY"]
        for i, merit in enumerate(optimization_measures):
            details = "{0}/NBins{6}/{1}/kMin{4}/kMax{5}/N{2}/{3}/".format(
                    meta_data, method, N, merit, k_min, k_max, n_bins)
            data_dir = base_data_dir + details
            label = merit_labels[i] + " Optimized"
            labels.append(label)
            data_dirs.append(data_dir)
    

    #plots_folder = "N{0}nBins{1}{2}/".format(N, n_bins, brief_description)
    #plots_dir = base_plots_dir + plots_folder
    #make_directory(data_dir)
    #make_directory(plots_dir)

    save_plots = 1
    
    # params for rerunning networks are in there.
    param_fold = 1
    if param_fold:
        number = "1000"
        avg_k = "4"
        running_folder = "TempData/ClusterReRunning/"
        seed = "1"
        Lambda = "0.0"
        alpha = "2.27"
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

    use_scale_free_data = 1
    if use_scale_free_data:
        add_metrics = 1
        Ns = ["128", "512", "1024"][2:]
        merits = ["DeathAge", "QALY", "HealthyAging"][:1]
        optimization_measures = ["DeathAge", "DeathAge"]
        tolerance = 0.05
        meta_data = "ScaleFreeN1024"
        method = "Variational"
        for N in Ns:
            for merit in merits:
                details = "{0}/{1}/N{2}/{3}/".format(meta_data,
                                                    method, N, merit, k_min, k_max)
                data_dir = base_data_dir + details
                data_dirs.append(data_dir)
                label = "Scale-Free" #"N{0}Merit{1}".format(N, merit)
                labels.append(label)
            brief_description = "ScaleFree/"
        #meta_data = ""
        meta_data = "FinalN1024"#"betterkMin2nBins15"
        N = "1024"
        k_min = "2"
        k_max = "1.0"
        n_bins = "15"
        method = "NonParametric"
        for merit in ["DeathAge", "QALY"][:1]:
            details = "{0}/NBins{6}/{1}/kMin{4}/kMax{5}/N{2}/{3}/".format(
                    meta_data, method, N, merit, k_min, k_max, n_bins)
            data_dir = base_data_dir + details
            label = "Non-Parametric"#merit + " Optimized"
            labels.append(label)
            data_dirs.append(data_dir)
         



    use_test_data = 0
    if use_test_data:
        N = 128
        data_dirs = [
            "Data/TopStart/NBins15/NonParametric/kMin2/kMax1.0/N256/DeathAge/",
            "Data/MixedAlphaChange/Variational/N128/DeathAge/"]
        optimization_measures = ["DeathAge", "DeathAge"]
        tolerance = 0.05
        labels = ["NonParametric", "MixedAlphaChange"]
        brief_description = "MixedAlpha/"
        add_metrics = True


    plots_dir = "Plots/SimplePlots/" + brief_description
    do_simple_plots = 1
    if do_simple_plots:
        make_directory(plots_dir)
        simple_entropy_plots(data_dirs, labels, plots_dir, 
                            save_plots, N, running_folder, params, param_names,
                            optimization_measures,
                            add_scale_free=True, add_metrics = add_metrics,
                            add_alpha=add_alpha,
                            tolerance = tolerance)

    do_full_plots = 0
    if do_full_plots:
        bins_dir = data_dirs[1]
        reference_data_dir = data_dirs[0]
        plots_dir = "Plots/PaperPlots/NetworkGraphs/RemadeQALYvsTd/"
        make_directory(plots_dir)
        plot_title = ""
        colour_metric = "Moment"
        plotting_targets = [0.0, 4.0]#, 10]
        generate_entropy_plots(bins_dir, reference_data_dir, plots_dir,
                        plot_title, save_plots, N,
                        colour_metric = colour_metric,
                        plotting_targets = plotting_targets)



    fold_old_stuff = 1
    if fold_old_stuff:
        """ Old plotting stuff down there VVV  





        check_folder = base_data_dir + details + "Lambda{}/".format(entropy_weight)
        test_graph_file = check_folder + "EntropyTarget1.0/Seed1/BestNetwork.csv"
        plot_title = "Low Entropy N128 Example"
        example_plot_folder = "Plots/"
        rerun_network(test_graph_file, running_folder, params,
                                        param_names, example_plot_folder,
                                        plot_title, 1000)


        #plot_title = "{0} {1} Optimization N {2}".format(method, merit, N)
        #check_progress(check_folder, plots_dir, save_plots)
        #versus_lambda(data_dir, plots_dir, plot_title, save_plots)
        #compare_degree_distributions(data_dir, plots_dir, plot_title, save_plots)
               
        bins_dir = base_data_dir + details
        colour_metric = "kMax"
        plot_title = "{2}OptNBin{0}N{1}kMin{3}kMax{4}".format(n_bins, N, merit,
                                                            k_min, k_max)
        reference_data_dir = "ClusterData/{0}/NBins10/{1}/kMin{4}"\
                                "/kMax{5}/N{2}/{3}/".format(
                                "kMinMax", "NonParametric", N, merit, k_min, k_max)


        #generate_entropy_plots(bins_dir, reference_data_dir, plots_dir, plot_title,
        #                save_plots, N,
        #                colour_metric = colour_metric,
        #                plotting_targets = plotting_targets)
        #network_metrics_vs_entropy(data_dir, plots_dir, plot_title, save_plots)
        #versus_entropy(bins_dir, plots_dir, plot_title, save_plots, N,
        #                colour_metric = colour_metric,
        #                plotting_targets = plotting_targets)
     
        #assortativity_extremes_vs_entropy(data_dir, plots_dir, plot_title,
        #    params, param_names, running_folder, save_plots)

        #plot_max_merit(data_dir, plots_dir, plot_title, save_plots)

        entropy_indices = [0, -1]
        num_seeds = 3
        #draw_networks(data_dir, plots_dir, plot_title, entropy_indices,
        #            num_seeds, save_plots)

        """

    pl.show()
