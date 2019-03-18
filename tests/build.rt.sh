#!/bin/bash

cd bin
sed 's/RT=0/RT=1/g' Makefile |\
    sed 's/NGROUPS = 0/NGROUPS = 3/g' |\
    sed 's/NIONS = 0/NIONS = 3/g' |\
    sed 's/^F90 =.*$/F90 = mpif90 -frecord-marker=4 -ffree-line-length-none -fbacktrace -g -O -fbounds-check -Wuninitialized -Wall/g' |\
    sed 's/^FFLAGS =.*$/FFLAGS = -x f95-cpp-input -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan  $(DEFINES)/g' > Makefile.rt
make -f Makefile.rt clean

make -f Makefile.rt
