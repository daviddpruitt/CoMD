#!/bin/bash
#SBATCH -A TG-ASC140011
#SBATCH -J CoMD
#SBATCH -o CoMD%j.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 03:00:00
#SBATCH --mail-user=ddpruitt@miners.utep.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load intel/15.0.2

./stampSub.sh