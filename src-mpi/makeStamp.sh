#!/bin/bash
set -x

module load intel/15.0.2

make clean; make FRIV="-DNEIGHBOR_LIST -DSEP_LOOPS"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-neighbor

make clean; make FRIV="-DNEIGHBOR_LIST"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-neighbor-old

make clean; make
mv ../bin/CoMD-mpi ../bin/CoMD-mpi

