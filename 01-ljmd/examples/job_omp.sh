#!/bin/bash
#SBATCH -A ICT24_MHPC                    #  <account_name>
#SBATCH --job-name=LJMD_omp      # Job name
#SBATCH --output=ljmd_omp.datg
##SBATCH --mail-type=END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=gparedes@ictp.it     # Where to send mail
#SBATCH -N 1                             # Number of nodes
##SBATCH --ntasks-per-node=1             # Run on a single CPU
#SBATCH --cpus-per-task=32               # 2, 4, 8, 16, 32
#SBATCH --time=00:30:00                  # Time limit hrs:min:sec
#SBATCH -p boost_usr_prod               # For Leonardo # for smaller works use --qos=boost_qos_dbg , boost is dedicated to works with GPU


echo "Running LJMD OPENMP version"
make
cd examples

echo "num threads = 1"
export OMP_NUM_THREADS=1

echo "108 natoms"
../build/ljmd.x < argon_108.inp
echo "2916 natoms"
../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
../build/ljmd.x < argon_78732.int

echo "num threads = 2"
export OMP_NUM_THREADS=2

echo "108 natoms"
../build/ljmd.x < argon_108.inp
echo "2916 natoms"
../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
../build/ljmd.x < argon_78732.int

echo "num threads = 4"
export OMP_NUM_THREADS=4
echo "108 natoms"
../build/ljmd.x < argon_108.inp
echo "2916 natoms"
../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
../build/ljmd.x < argon_78732.int

echo "num threads = 8"
export OMP_NUM_THREADS=8
echo "108 natoms"
../build/ljmd.x < argon_108.inp
echo "2916 natoms"
../build/ljmd.x < argon_2916.inp
echo "78732 natoms"
../build/ljmd.x < argon_78732.int