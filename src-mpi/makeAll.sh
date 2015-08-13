#!/bin/bash
set -x

make clean; make FRIV="-DSKIP_FORCE_COMPUTATIONS -DDO_FRIVOLOUS_PREFETCH"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-noforce-prefetch

make clean; make FRIV="-DDO_FRIVOLOUS_PREFETCH"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-prefetch

make clean; make FRIV=""
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-normal

make clean; make FRIV="-DSKIP_FORCE_COMPUTATIONS -DDO_FRIVOLOUS_PREFETCH -DDO_FRIVOLOUS_TIMING"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-noforce-prefetch-timing

make clean; make FRIV="-DDO_FRIVOLOUS_PREFETCH  -DDO_FRIVOLOUS_TIMING"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-prefetch-timing

make clean; make FRIV=" -DDO_FRIVOLOUS_TIMING"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-normal-timing

make clean; make FRIV=" -DDO_FRIVOLOUS_TIMING -DDO_FRIVOLOUS_PREFETCH  -DDO_TWO_PREFETCH"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-prefetch-timing-double

make clean; make FRIV="-DDO_FRIVOLOUS_PREFETCH  -DDO_TWO_PREFETCH"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-prefetch-double

make clean; make FRIV="-DDO_PREFETCH_POT -DDO_FRIVOLOUS_PREFETCH  -DDO_TWO_PREFETCH"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-prefetch-double-pot

make clean; make FRIV="-DDO_PREFETCH_POT -DDO_FRIVOLOUS_PREFETCH"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-prefetch-pot

make clean; make FRIV="-DDO_PREFETCH_POT"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-pot

make clean; make FRIV="-DSKIP_RCUT"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-nocut

make clean; make FRIVE="-DUNROLL_DIM"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-unroll

make clean; make FRIVE="-RESTRICT_PTR"
mv ../bin/CoMD-mpi ../bin/CoMD-mpi-restrict



