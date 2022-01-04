"""
Building the table.dat file for running META jobs on the cluster
"""

import itertools
import pylab as pl

def output_table(params_list, output_file):
    
    combined_lists = itertools.product(*params_list)
    f = open(output_file, "w")
    for ps in combined_lists:
        line = ""
        for param in ps:
            line += str(param) + " "

        f.write(line[:-1] + "\n")
    f.close()


if __name__ == "__main__":

    output_folder = "./"
    output_file = output_folder + "table.dat"
    
    Ns = [128, 512, 1024][:1]#, 10000]#, 5000, 10000]
    numbers = [1000]
    health_measures = ["DeathAge"]
    entropy_targets = [1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0][:1]
    n_bins = [10]
    k_mins = [1,2]
    k_maxes = [0.25, 0.5, 0.75, 1.0][:1]
    seeds = [1,2,3][:1]
    lambdas = [0.5, 0.85][:1]



    # Need to generate parameters from N onwards for jobs
    params_list = [Ns, numbers, health_measures, entropy_targets, n_bins,
                    k_mins, k_maxes, seeds, lambdas]
    num_tasks = 1
    for p in params_list:
        num_tasks *= len(p)
    print("Total Number of Tasks:", num_tasks)
    output_table(params_list, output_file)
