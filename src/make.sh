#!/bin/bash
module --quiet purge
######module load OpenMPI/4.0.3-GCC-9.3.0
module load OpenMPI/4.0.5-GCC-10.2.0
make clean
make -j4
