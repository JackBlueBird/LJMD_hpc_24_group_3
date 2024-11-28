#!/bin/bash
#SBATCH --job-name=LJMD_hybrid4    ## Job name
#SBATCH --output=ljmd_hybrid22.datg
#SBATCH -N 1 # Number of nodes
#SBATCH --ntasks-per-node=1                   ## Run on a single CPU
#SBATCH --cpus-per-task=32
#SBATCH --time=00:10:00               ## Time limit hrs:min:sec
#SBATCH -A ICT24_MHPC
#SBATCH -p boost_usr_prod             ## for Leonardo

echo "Running LJMD hybrid version"
module load gcc/
module load openmpi/

echo "num threads = 2 ----------------------2x2------------------------------"
echo "np = 2"
export OMP_NUM_THREADS=2
echo "108 natoms"
mpirun -np 2 ../build/ljmd.x < argon_108.inp
echo "2916 natoms"
mpirun -np 2 ../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
mpirun -np 2 ../build/ljmd.x < argon_78732.inp
