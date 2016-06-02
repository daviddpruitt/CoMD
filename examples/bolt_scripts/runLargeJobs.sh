#!/bin/sh
qsub -cwd -V -N "MPI-LARGE-1" -o mpi-1-large.out -e mpi-1-large.err -q micdev.q -pe orte 1 -M ddpruitt@miners.utep.edu -m abe mpilarge1.sh
qsub -cwd -V -N "MPI-LARGE-2" -o mpi-2-large.out -e mpi-2-large.err -q micdev.q -pe orte 2 -M ddpruitt@miners.utep.edu -m abe mpilarge2.sh
qsub -cwd -V -N "MPI-LARGE-4" -o mpi-4-large.out -e mpi-4-large.err -q micdev.q -pe orte 4 -M ddpruitt@miners.utep.edu -m abe mpilarge4.sh
qsub -cwd -V -N "MPI-LARGE-8" -o mpi-8-large.out -e mpi-8-large.err -q micdev.q -pe orte 8 -M ddpruitt@miners.utep.edu -m abe mpilarge8.sh
