#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=700GB
#SBATCH --time=12:00:00
#SBATCH --job-name=gpc-box-21cmfast-test-olivia
#SBATCH --account=oz113
#SBATCH --chdir=/fred/oz113/owalters