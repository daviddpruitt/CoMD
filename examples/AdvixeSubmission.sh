#!/bin/bash
#SBATCH -A TG-ASC140011
#SBATCH -J CoMD
#SBATCH -o CoMD%j.o
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 03:00:00
#SBATCH --mail-user=ddpruitt@miners.utep.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load intel/15.0.2

source ../../intel-adv/parallel_studio_xe_2016.0.035/advisor_xe_2016/advixe-vars.sh

echo "*************************************************************************************************************************"
echo "***Running ibrun advixe-cl -collect survey -project-dir ./myAdvisorProj2 ../bin/CoMD-mpi -n 10 -N 100 -e -x 16 -y 16 -z 16"
echo "*************************************************************************************************************************"
ibrun advixe-cl -collect survey -project-dir ./myAdvisorProj3 ../bin/CoMD-mpi -n 10 -N 100 -e -x 16 -y 16 -z 16
echo "*************************************************************************************************************************"
echo "***Running ibrun advixe-cl -collect tripcounts -project-dir ./myAdvisorProj2 ../bin/CoMD-mpi -n 10 -N 100 -e -x 16 -y 16 -z 16"
echo "*************************************************************************************************************************"
ibrun advixe-cl -collect tripcounts -project-dir ./myAdvisorProj3 ../bin/CoMD-mpi -n 10 -N 100 -e -x 16 -y 16 -z 16
echo "*************************************************************************************************************************"
echo "***Running ibrun advixe-cl -collect dependencies -project-dir ./myAdvisorProj2 ../bin/CoMD-mpi -n 10 -N 100 -e -x 16 -y 16 -z 16"
echo "*************************************************************************************************************************"
ibrun advixe-cl -collect dependencies -project-dir ./myAdvisorProj3 ../bin/CoMD-mpi -n 10 -N 100 -e -x 16 -y 16 -z 16
echo "*************************************************************************************************************************"
echo "***Running ibrun advixe-cl -collect map -project-dir ./myAdvisorProj2 ../bin/CoMD-mpi -n 10 -N 100 -e -x 16 -y 16 -z 16"
echo "*************************************************************************************************************************"
ibrun advixe-cl -collect map -project-dir ./myAdvisorProj3 ../bin/CoMD-mpi -n 10 -N 100 -e -x 16 -y 16 -z 16
#ibrun /usr/bin/time -v ../bin/CoMD-mpi              -n 10 -N 100 -e -x 16 -y 16 -z 16 &>> timing_16.txt
