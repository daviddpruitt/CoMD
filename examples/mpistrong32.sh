#!/bin/sh
#PBS -o strong32.out
#PBS -e strong32.out
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
mpirun -np 32 /home/ddpruitt/testCode/CoMD/bin/CoMD-mpi -e -i 4 -j 4 -k 2 -x 40 -y 40 -z 40
