#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --nodes=1
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

module load python/3.6
module load scipy-stack
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip==21.0.1+computecanada
pip install --no-index -r requirements.txt


tempFolder="/home/garretts/scratch/SimplyConnected/"
outputFolder="Data/"
runTime="30" # in Minutes RN
N="1000"
number="1000"
method="Variational"
seed="1"
healthMeasure="HealthyAging" # Options: HealthyAging, DeathAge, QALY
details="Sep9Tests/$method/N$N/$healthMeasure/"
entropyWeight=0.5

make clean
make clusterTestNetwork
make clusterMain

my_parallel="parallel --delay=0.2 -j $SLURM_NTASKS"
my_srun="srun --exclusive -N1 -n1 -c1" # --cpus-per-task=1 --cpu-bind=cores"


$my_parallel "$my_srun python optimization.py $tempFolder $outputFolder\
    $details $method $runTime $N $number $healthMeasure {1} $entropyWeight\
    " :::: seeds.txt

date
