#!/usr/bin/sh
#SBATCH --job-name=singleCellExtForce
#SBATCH -output /u/jferrari/singleCellExtForce 
#SBATCH --ntasks=56

echo ' running copy of singleCellExtForce'
srun CC3D_RunScript -i singleCellExternalForce.cc3d -o $PWD -f 9999999