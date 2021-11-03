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
    
    Ns = [512, 1024, 2048]#, 5000, 10000]
    numbers = [5000]
    health_measures = ["DeathAge"]#, "QALY"]
    entropy_targets = [1.0, 2.0, 3.0, 4.0, 5.0]
    n_bins = [10]
    seeds = [1,2,3,4,5,6,7,8,9]
    lambdas = [0.5, 0.8, 0.9, 0.99]



    # Need to generate parameters from N onwards for jobs
    params_list = [Ns, numbers, health_measures, entropy_targets, n_bins,
                    seeds, lambdas]
    num_tasks = 1
    for p in params_list:
        num_tasks *= len(p)
    print("Total Number of Tasks:", num_tasks)
    output_table(params_list, output_file)
