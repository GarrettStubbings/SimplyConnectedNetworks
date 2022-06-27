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
    
    Ns = [1024]#128, 512, 1024, 2048, 10000]#[:1]#, 10000]#, 5000, 10000]
    numbers = [10000]
    health_measures = ["DeathAge", "QALY"]#, "HealthyAging"]
    entropy_targets = [0.0, 1.33, 1.67, 2.0, 2.33, 2.67, 3.0, 3.25, 3.5, 3.75,
                        4.0, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 100]
    n_bins = [15]
    k_mins = [1]
    k_maxes = [1.0]#[:1]
    seeds = [1,2,3,4,5,6,7,8]#[:1]
    lambdas = [0.9, 0.95, 0.98]#[:1]



    # Need to generate parameters from N onwards for jobs
    params_list = [Ns, numbers, health_measures, entropy_targets, n_bins,
                    k_mins, k_maxes, seeds, lambdas]
    num_tasks = 1
    for p in params_list:
        num_tasks *= len(p)
    print("Total Number of Tasks:", num_tasks)
    output_table(params_list, output_file)
