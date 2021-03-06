#!/bin/bash
#SBATCH -A TG-ASC140011
#SBATCH -J CoMD
#SBATCH -o CoMD%j.o
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 00:20:00
#SBATCH --mail-user=ddpruitt@miners.utep.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load gcc/4.9.1

for i in `seq 1 2` ;
do
    ibrun /usr/bin/time -v ../bin/CoMD-mpi              -n 10 -N 100 -e -x 16 -y 16 -z 16 &>> timing_16.txt
    ibrun /usr/bin/time -v ../bin/CoMD-mpi-neighbor-old -n 10 -N 100 -e -x 16 -y 16 -z 16 &>> timing_16.txt
done

