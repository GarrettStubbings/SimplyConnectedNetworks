#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --account=def-arutenbe
#SBATCH --mail-user=garrett.stubbings@dal.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09
module load intel/2018.3
module load gsl/2.5
module load gcc/7.3.0
module load boost/1.66.0

make clusterTestNetwork
make clusterMain

date
