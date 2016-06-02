#!/bin/sh
qsub -cwd -V -N "MPI1" -o mpi-1.out -e mpi-1.err -q micdev.q -pe orte 1 -M ddpruitt@miners.utep.edu -m abe mpistrong1.sh
qsub -cwd -V -N "MPI2" -o mpi-2.out -e mpi-2.err -q micdev.q -pe orte 2 -M ddpruitt@miners.utep.edu -m abe mpistrong2.sh
qsub -cwd -V -N "MPI4" -o mpi-4.out -e mpi-4.err -q micdev.q -pe orte 4 -M ddpruitt@miners.utep.edu -m abe mpistrong4.sh
qsub -cwd -V -N "MPI8" -o mpi-8.out -e mpi-8.err -q micdev.q -pe orte 8 -M ddpruitt@miners.utep.edu -m abe mpistrong8.sh
qsub -cwd -V -N "MPI16" -o mpi-16.out -e mpi-16.err -q micdev.q -pe orte 16 -M ddpruitt@miners.utep.edu -m abe mpistrong16.sh
qsub -cwd -V -N "MPI32" -o mpi-32.out -e mpi-32.err -q micdev.q -pe orte 32 -M ddpruitt@miners.utep.edu -m abe mpistrong32.sh
