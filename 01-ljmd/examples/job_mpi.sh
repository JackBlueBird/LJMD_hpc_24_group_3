#!/bin/bash
#SBATCH --job-name=LJMD_mpi    # Job name
#SBATCH --output=ljmd_mpi.datg
#SBATCH -N 1 # Number of nodes
#SBATCH --ntasks-per-node=32                   # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00               # Time limit hrs:min:sec
#SBATCH -A ICT24_MHPC
#SBATCH -p boost_usr_prod             # For Leonardo

echo "Running LJMD MPI version"
module load gcc/
module load openmpi/

echo "108 natoms"
mpirun -np 1 ../build/ljmd.x < argon_108.inp
echo "2916 natoms"
mpirun -np 1 ../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
mpirun -np 1 ../build/ljmd.x < argon_78732.int

echo "LJMD MPI with np = 2"
echo "108 natoms"
mpirun -np 2 ../build/ljmd.x < argon_108.inp
echo "2916 natoms"
mpirun -np 2 ../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
mpirun -np 2 ../build/ljmd.x < argon_78732.int

echo "LJMD MPI with np = 4"
echo "108 natoms"
mpirun -np 4 ../build/ljmd.x < argon_108.inp
echo "2916 natoms"
mpirun -np 4 ../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
mpirun -np 4 ../build/ljmd.x < argon_78732.int

echo "LJMD MPI with np = 8"
echo "108 natoms"
mpirun -np 8 ../build/ljmd.x < argon_108.inp
echo "2916 natoms"
mpirun -np 8 ../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
mpirun -np 8 ../build/ljmd.x < argon_78732.int