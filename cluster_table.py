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
    
    Ns = [100, 250, 500, 1000]#, 2000, 5000, 10000]
    numbers = [5000]
    health_measures = ["QALY"]#, "QALY"]
    kmins = [2]
    kmaxs = [10]
    seeds = [1,2,3,4,5,6]
    lambdas = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                0.55, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    # Pretty good death age stuff
    #[0.4, 0.45, 0.5, 0.55, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725,
    #            0.7, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.95, 0.99, 1.0]


    #[0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25,
    #            0.275]
    #[0.3, 0.4, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65,
    #            0.675, 0.7, 0.75, 0.8, 0.9, 1.0]
    # Need to generate parameters from N onwards for jobs
    params_list = [Ns, numbers, health_measures, kmins, kmaxs,
                    seeds, lambdas]
    num_tasks = 1
    for p in params_list:
        num_tasks *= len(p)
    print("Total Number of Tasks:", num_tasks)
    output_table(params_list, output_file)
