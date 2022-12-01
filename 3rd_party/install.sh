#!/bin/bash
## Install 3rd party libraries
## Only Silo library is needed.
wget https://github.com/LLNL/Silo/releases/download/v4.11/silo-4.11.tar.gz
tar -xf silo-4.11.tar.gz
cd silo-4.11/
./configure CC=gcc FC=gfortran
make install