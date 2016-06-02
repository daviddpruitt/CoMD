#!/usr/bin/env bash

module load intel/15.0.2


for i in `seq 1 9` ;
do
    ibrun /usr/bin/time -v ../bin/CoMD-mpi              -n 10 -N 100 -e -x 16 -y 16 -z 16 &>> timing_16.txt
    ibrun /usr/bin/time -v ../bin/CoMD-mpi-neighbor-old -n 10 -N 100 -e -x 16 -y 16 -z 16 &>> timing_16.txt
    ibrun /usr/bin/time -v ../bin/CoMD-mpi-neighbor     -n 10 -N 100 -e -x 16 -y 16 -z 16 &>> timing_16.txt
done

for i in `seq 1 9` ;
do
    ibrun /usr/bin/time -v ../bin/CoMD-mpi              -n 10 -N 100 -e -x 24 -y 24 -z 24 &>> timing_55.txt
    ibrun /usr/bin/time -v ../bin/CoMD-mpi-neighbor-old -n 10 -N 100 -e -x 24 -y 24 -z 24 &>> timing_55.txt
    ibrun /usr/bin/time -v ../bin/CoMD-mpi-neighbor     -n 10 -N 100 -e -x 24 -y 24 -z 24 &>> timing_55.txt
done

for i in `seq 1 9` ;
do
    ibrun /usr/bin/time -v ../bin/CoMD-mpi              -n 10 -N 100 -e -x 40 -y 40 -z 40 &>> timing_256.txt 
    ibrun /usr/bin/time -v ../bin/CoMD-mpi-neighbor-old -n 10 -N 100 -e -x 40 -y 40 -z 40 &>> timing_256.txt
    ibrun /usr/bin/time -v ../bin/CoMD-mpi-neighbor     -n 10 -N 100 -e -x 40 -y 40 -z 40 &>> timing_256.txt
done

