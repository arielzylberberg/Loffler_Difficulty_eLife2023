#!/bin/sh
 
#SBATCH -A zims
 
#SBATCH --time=900
#SBATCH -N 1
#SBATCH --exclusive
 
module load intel-parallel-studio/2017
 
./a.out
 
# End of script