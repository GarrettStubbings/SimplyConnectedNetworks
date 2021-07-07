#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --ntasks=2
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

number=9999 # number of runs
N=200 #network size
numChanges=1 # Minimum number of moves (will descend to this quickly)
avgdeg=4 #average degree(must be even for scalefree and smallworld)
Folder="Local"
SingleSeed=2
runHours=0.01
evoCondition="DeathAge"
initialDistribution="MaxEntropyN200"
lambda=0.05
beta=100.0
power=0.5

make evo

srun="srun -n1 -N1 --exclusive"

parallel="parallel -N 1 --delay 1 -j $SLURM_NTASKS --joblog runlogNov24 --resume"

$parallel "$srun ./evo $number $N $numChanges $avgdeg\
    $Folder $SingleSeed \
    $runHours $evoCondition $initialDistribution {1} {2} {3}" \
    :::: smallLambdas.txt :::: betaValues.txt :::: coolingPowers.txt

