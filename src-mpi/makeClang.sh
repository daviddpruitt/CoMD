#!/bin/bash
set -x

# set to clang first
export OMPI_CXX=clang
export OMPI_CC=clang

make clean; make FRIV="-DNEIGHBOR_LIST -DSEP_LOOPS"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-neighbor-clang

make clean; make FRIV="-DNEIGHBOR_LIST"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-neighbor-old-clang

make clean; make
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-clang

export OMPI_CXX=gcc
export OMPI_CC=cc

make clean; make FRIV="-DNEIGHBOR_LIST -DSEP_LOOPS"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-neighbor

make clean; make FRIV="-DNEIGHBOR_LIST"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-neighbor-old

make clean; make
mv ../bin/CoMD-mpi ../bin/CoMD-mpi

