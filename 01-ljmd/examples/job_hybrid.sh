#!/bin/bash
#SBATCH --job-name=LJMD_hybrid    ## Job name
#SBATCH --output=ljmd_hybrid.datg
#SBATCH -N 1 # Number of nodes
#SBATCH --ntasks-per-node=32                   ## Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00               ## Time limit hrs:min:sec
#SBATCH -A ICT24_MHPC
#SBATCH -p boost_usr_prod             ## for Leonardo

echo "Running LJMD hybrid version"
module load gcc/
module load openmpi/

make
cd examples

echo "num threads = 4"
echo "np = 4"
export OMP_NUM_THREADS=4
echo "108 natoms"
mpirun -np 4 ../build/ljmd.x < argon_108.inp
echo "2916 natoms"
mpirun -np 4 ../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
mpirun -np 4 ../build/ljmd.x < argon_78732.int

##echo "num threads = 4"
##echo "np = 2"
##export OMP_NUM_THREADS=4
##echo "108 natoms"
##mpirun -np 2 ../build/ljmd.x < argon_108.inp
##echo "2916 natoms"
##mpirun -np 2 ../build/ljmd.x < argon_2916.inp
##echo "78732 natoms"
##mpirun -np 2 ../build/ljmd.x < argon_78732.int

##echo "num threads = 4"
##echo "np = 4"
##export OMP_NUM_THREADS=4
##echo "108 natoms"
##mpirun -np 2 ../build/ljmd.x < argon_108.inp
##echo "2916 natoms"
##mpirun -np 2 ../build/ljmd.x < argon_2916.inp
##echo "78732 natoms"
##mpirun -np 2 ../build/ljmd.x < argon_78732.int

