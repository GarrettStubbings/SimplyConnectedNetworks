#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=30
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

my_srun="srun --exclusive -N1 -n1 -c1" # --cpus-per-task=1 --cpu-bind=cores"

$my_srun virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
$my_srun pip install --no-index --upgrade pip==21.0.1+computecanada
$my_srun pip install --no-index -r requirements.txt


tempFolder="/home/garretts/scratch/SimplyConnected/"
outputFolder="Data/"
runTime="710" # in Minutes RN
N="100"
number="1000"
method="NonParametric"
seed="1"
healthMeasure="DeathAge" # Options: HealthyAging, DeathAge, QALY
entropyWeight=0.5
kMin=2
kMax=20
if [[ "$method" == "Variational" ]]; then
    details="Sep21Tests/$method/N$N/$healthMeasure/"
else
    details="Sep21Tests/$method/N$N/$healthMeasure/kMin$kMin/kMax$kMax/"
fi


#make clean
make clusterTestNetwork
make clusterMain

my_parallel="parallel --delay=0.2 -j $SLURM_NTASKS"

$my_parallel "$my_srun python optimization.py $tempFolder $outputFolder\
    $details $method $runTime $N $number $healthMeasure $kMin $kMax {1} {2}\
    " :::: seeds1.txt :::: lambdas.txt

date
