#!/bin/bash
#SBATCH -A ICT24_MHPC                    #  <account_name>
#SBATCH --job-name=LJMD_omp       # Job name
#SBATCH --output=ljmd_omp.datg
##SBATCH --mail-type=END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=gparedes@ictp.it     # Where to send mail
#SBATCH -N 1                             # Number of nodes
##SBATCH --ntasks-per-node=1             # Run on a single CPU
#SBATCH --cpus-per-task=32               # 2, 4, 8, 16, 32
#SBATCH --time=00:15:00                  # Time limit hrs:min:sec
#SBATCH -p boost_usr_prod                # for smaller works use --qos=boost_qos_dbg , boost is dedicated to works with GPU


echo "Running ljmd with omp"
make
cd examples

echo "num threads = 4"
export OMP_NUM_THREADS=4
../ljmd-serial.x < argon_2916_2.inp

echo "num threads = 8"
export OMP_NUM_THREADS=8
../ljmd-serial.x < argon_2916_2.inp

echo "num threads = 12"
export OMP_NUM_THREADS=12
../ljmd-serial.x < argon_2916_2.inp

echo "num threads = 16"
export OMP_NUM_THREADS=16
../ljmd-serial.x < argon_2916_2.inp

echo "num threads = 20"
export OMP_NUM_THREADS=20
../ljmd-serial.x < argon_2916_2.inp

echo "num threads = 24"
export OMP_NUM_THREADS=24
../ljmd-serial.x < argon_2916_2.inp

echo "num threads = 28"
export OMP_NUM_THREADS=28
../ljmd-serial.x < argon_2916_2.inp

echo "num threads = 32"
export OMP_NUM_THREADS=32
../ljmd-serial.x < argon_2916_2.inp

