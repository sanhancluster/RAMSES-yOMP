#!/bin/bash

cd bin
sed 's/DICE=0/DICE=1/g' Makefile |\
    sed 's/^F90 =.*$/F90 = mpif90 -frecord-marker=4 -ffree-line-length-none -fbacktrace -g -O -fbounds-check -Wuninitialized -Wall/g' |\
    sed 's/^FFLAGS =.*$/FFLAGS = -x f95-cpp-input -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan  $(DEFINES)/g' > Makefile.dice
make -f Makefile.dice clean

make -f Makefile.dice
