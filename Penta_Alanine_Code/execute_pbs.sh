#!/bin/bash

#PBS -N A5
#PBS -e A5.err
#PBS -o A5.out
#PBS -m aeb
#PBS -M sahn1@stanford.edu
#PBS -l nodes=10:ppn=24
#PBS -l walltime=72:00:00
#PBS -q medium
#PBS -V

### --------------------------------------------------------------------------------------------------------
### Bash shell script is specific for PBS and GROMACS. This can be easily modified for other cluster systems
### (e.g. SLURM) and other molecular dynamics simulation programs (e.g. LAMMPS).
### --------------------------------------------------------------------------------------------------------

export MAIN_DIRECTORY=/home/sahn1/Penta_Alanine  # TODO
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sahn1/  # TODO
export GROMACS=/home/sahn1/gromacs/4.6.4/bin  # TODO

echo The master node of this job is `hostname`
echo The working directory is `echo $MAIN_DIRECTORY`
echo This job runs on the following nodes:
echo `cat $PBS_NODEFILE`

cd $MAIN_DIRECTORY

cat $PBS_NODEFILE > nodefilelist.txt
python main.py
exit
