FFLAGS=-std=f2008 --warn-all -Ofast -fimplicit-none
FCC=gfortran

all: run_test


%.o: %.f90
	$(FCC) -c $(FFLAGS) $^ -o $@

test_obb: geom_types.o  geom_utils.o test_utils.o geom_distances.o test_obb.o
	$(FCC) $(FFLAGS) $^ -o $@

test_volumes: geom_types.o geom_utils.o test_utils.o geom_distances.o geom_volumes.o test_volume.o
	$(FCC) $(FFLAGS) $^ -o $@

run_%: %
	./$^

run_test: test_obb test_volumes run_test_obb run_test_volumes


clean:
	rm -rf *.o *.mod
