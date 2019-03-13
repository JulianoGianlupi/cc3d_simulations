#!/bin/bash 
#PBS -k o 
#PBS -l nodes=1:ppn=1,walltime=10:15:00
#PBS -M jferrari@iu.edu
#PBS -m abe
#PBS -N simpleBoids 
#PBS -j oe

module load compucell3d/3.7.5


paramScan.sh -i /N/u/jferrari/Karst/boids_simple_force/chimeraBoidsV2.cc3d -o /N/u/jferrari/Karst/boids_simple_force/screenShots/ -f 20
