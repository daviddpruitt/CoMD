#!/bin/sh
#PBS -o strong4.out
#PBS -e strong4.out
#PBS -d ~/testCode
#PBS -q micdev.q
#PBS -V
#PBS -l nodes=8
# request 4 hours and 30 minutes of cpu time -l cput=04:30:00 
# mail is sent to you when the job starts and when it terminates or aborts
#PBS -m bea
# specify your email address
#PBS -M ddpruitt@miners.utep.edu

# Simple strong scaling study with eam potential and 256,000 atoms
mpirun -np 4  /home/ddpruitt/testCode/CoMD/bin/CoMD-mpi -e -i 2 -j 2 -k 1 -x 40 -y 40 -z 40
