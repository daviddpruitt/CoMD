#!/bin/bash
set -x

module load intel/15.0.2

make clean; make FRIV="-DNEIGHBOR_LIST -DSEP_LOOPS -g -O3 -march=native -mtune=native"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-neighbor

make clean; make FRIV="-DNEIGHBOR_LIST -g -O3 -march=native -mtune=native"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-neighbor-old

make clean; make FRIV="-g -O3 -march=native -mtune=native"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi

