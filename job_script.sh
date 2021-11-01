#!/bin/bash
#SBATCH --time=11:59:00
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


# Don't change this line:
task.run
