#!/bin/bash
export OMP_NUM_THREADS=1
make clean
make
./nse2d
