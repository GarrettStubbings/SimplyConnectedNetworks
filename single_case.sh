#!/bin/bash
# The actual (serial) computation, for a single case No. $2 in the table $1.

TABLE=$1
i=$2

# Total number of cases:
## If the env. variable N_cases is defined, using it, otherwise computing the number of lines in the table:
if test -z $N_cases
  then
  N_cases=`cat "$TABLE" | wc -l`
  fi

# Exiting if $i goes beyond $N_cases (can happen in bundled farms):
if test $i -lt 1 -o $i -gt $N_cases  
  then
  exit
  fi
  
# Extracing the $i-th line from file $TABLE:
LINE=`sed -n ${i}p $TABLE`
# Case id (from the original cases table):
ID=`echo "$LINE" | cut -d" " -f1`
# The rest of the line:
COMM=`echo "$LINE" | cut -d" " -f2-`

# ++++++++++++++++++++++  This part can be customized:  ++++++++++++++++++++++++
#  Here:
#  $ID contains the case id from the original table (can be used to provide a unique seed to the code etc)
#  $COMM is the line corresponding to the case $ID in the original table, without the ID field
#  $SLURM_JOB_ID is the jobid for the current meta-job (convenient for creating per-job files)

echo "Case $ID:"

# Executing the command (a line from table.dat)
# It's allowed to use more than one shell command (separated by semi-columns) on a single line

# Garrett's code starts here:
my_srun="srun --exclusive -N1 -n1 -c1" # --cpus-per-task=1 --cpu-bind=cores"
tempFolder="$SLURM_TMPDIR/" #"/home/garretts/scratch/SimplyConnected/"
outputFolder="Data/"
runTime="5" # in Minutes RN
N="100"
number="1000"
method="NonParametric"
seed="1"
healthMeasure="DeathAge" # Options: HealthyAging, DeathAge, QALY
entropyWeight=0.5
entropyTarget=3.0
nBins=12
details="targetTest"

eval "$my_srun python optimization.py $tempFolder $outputFolder\
    $details $method $runTime $COMM"

# Exit status of the code:
STATUS=$?


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Saving all current metajob statuses to a single file with a unique name:
echo $ID $STATUS >> status.$SLURM_JOB_ID
