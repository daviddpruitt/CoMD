#!/bin/sh
#PBS -o strong1.out
#PBS -e strong1.out
#PBS -d ~/testCode
#PBS -q micdev.q
#PBS -V
#PBS -l nodes=8
# request 4 hours and 30 minutes of cpu time -l cput=04:30:00 
# mail is sent to you when the job starts and when it terminates or aborts
#PBS -m bea
# specify your email address
#PBS -M ddpruitt@miners.utep.edu

# Simple strong scaling study with eam potential and 256,000  * 4 atoms
mpirun -np 2  ../bin/CoMD-mpi -N 400 -n 40 -e -i 2 -j 1 -k 1 -x 80 -y 80 -z 40

