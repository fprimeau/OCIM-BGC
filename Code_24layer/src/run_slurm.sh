#!/bin/bash

#SBATCH --job-name=BGCmodelrun_240418
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --partition=cas3.6,has2.5,brd2.4,ilg2.3,m-c1.9,m-c2.2,nes2.8,sib2.9
#SBATCH --time=720:00:00
#SBATCH --nodes=1
#SBATCH --exclude=c-13-27,c-13-32,c-12-33
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --mail-type=all
#SBATCH --mail-user=hojons1@uci.edu
#
ml purge
ml matlab/R2022b
echo "----------------------------------------------"
echo "Date: $(date)"
echo "JOB ID: "
echo $SLURM_JOB_ID
echo "I ran on:"
echo $SLURM_NODELIST
echo "Partition: "
echo $SLURM_JOB_PARTITION
echo "NOTE: BGC model run for acqurting the parameters of steady-state solution in order to run the time stepping for C isotopes"
echo
#
matlab -nodesktop -nosplash < driver.m
##exit
##EOF
