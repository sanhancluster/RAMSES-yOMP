#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@gmail.com
#############################################################################
# Compilation time parameters
NVECTOR = 128
NDIM = 3
NPRE = 8
NENER = 0
NGROUPS = 0
NIONS = 0
METALS = 0
#Extra variables beyond those required by dimensions, or listed above
NVAR_EXTRA = 0
SOLVER = hydro
#############################################################################
# Options
GRACKLE = 0
RT = 0
DICE = 0
REFMASK = 0
LONGINT = 0
QUADHILBERT = 0
# Set to 1 to use old MPI include statements instead of more recent `use mpi`.
USE_OLD_MPI_INCLUDE = 0
EXEC = ramses
# Activate OpenMP
OPENMP = 1
#############################################################################
# Patches to include (without DICE)
#PATCH = ../patch/AGNSFSNpatch/AGNspinRT:../patch/AGNSFSNpatch/sfsn
#PATCH = ../patch/AGNSFSNpatch/sfsn
#PATCH = ../patch/debug
#############################################################################
# Compute the right value for NVAR
NVAR := $(shell echo $$(($(NDIM) + 2 + $(NVAR_EXTRA) + $(NIONS) + $(NENER) + $(REFMASK)+ $(METALS) )) )
#############################################################################
GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)
GITHASH = $(shell git log --pretty=format:'%H' -n 1)
GITREMOTE = $(shell git config --get branch.$(GITBRANCH).remote)
GITREPO = $(shell git config --get remote.$(GITREMOTE).url)
BUILDDATE = $(shell date +"%D-%T")
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNENER=$(NENER) -DNVAR=$(NVAR) \
          -DSOLVER$(SOLVER)
ifeq ($(USE_OLD_MPI_INCLUDE),1)
DEFINES += -DMPI_OLD
endif
#############################################################################
# Get the right compilation options, patches paths, and exec name
# Please keep things in that order
# --- RT --------------------------------------------------------------------
ifeq ($(RT),1)
DEFINES += -DRT -DNIONS=$(NIONS) -DNGROUPS=$(NGROUPS)
CEXEC += _rt$(NGROUPS)
endif
# --- DICE ------------------------------------------------------------------
ifeq ($(DICE),1)
DEFINES += -DDICE
CEXEC += _dice
endif
# --- RefMask ---------------------------------------------------------------
ifeq ($(REFMASK),1)
CEXEC += _refmask
endif
# --- Grackle ---------------------------------------------------------------
ifeq ($(GRACKLE),1)
DEFINES += -Dgrackle
CEXEC += _grackle
endif
# --- OpenMP ----------------------------------------------------------------
ifeq ($(OPENMP),1)
OMPFLAG += -fopenmp
CEXEC += _omp
endif
# --- LONGINT ---------------------------------------------------------------
ifeq ($(LONGINT),1)
DEFINES += -DLONGINT
CEXEC += _lint
endif
# --- QUADHILBERT -----------------------------------------------------------
ifeq ($(QUADHILBERT),1)
DEFINES += -DQUADHILBERT
CEXEC += _qhil
endif
#############################################################################
# Regroup exec suffixes and patches
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
EXECtmp = $(EXEC) $(CEXEC)
EXEC := $(subst $(SPACE),$(EMPTY),$(EXECtmp))
PATCH := $(subst $(SPACE),:,$(PATCH))
#############################################################################
# Fortran compiler options and directives

# --- No MPI, gfortran -------------------------------
#F90 = gfortran -O3 -frecord-marker=4 -fbacktrace -ffree-line-length-none -g
#FFLAGS = -x f95-cpp-input $(DEFINES) -DWITHOUTMPI

# --- No MPI, tau ----------------------------------
#F90 = tau_f90.sh -optKeepFiles -optPreProcess -optCPPOpts=$(DEFINES) -DWITHOUTMPI

# --- No MPI, pgf90 ----------------------------------
#F90 = pgf90
#FFLAGS = -Mpreprocess $(DEFINES) -DWITHOUTMPI

# --- No MPI, xlf ------------------------------------
#F90 = xlf
#FFLAGS = -WF,-DNDIM=$(NDIM),-DNPRE=$(NPRE),-DNVAR=$(NVAR),-DSOLVER$(SOLVER),-DWITHOUTMPI -qfree=f90 -qsuffix=f=f90 -qsuffix=cpp=f90

# --- No MPI, f90 ------------------------------------
#F90 = f90
#FFLAGS = -cpp $(DEFINES) -DWITHOUTMPI

# --- No MPI, ifort ----------------------------------
#F90 = ifort
#FFLAGS = -cpp $(DEFINES) -DWITHOUTMPI

# --- MPI, gfortran syntax ------------------------------
#F90 = mpif90 -frecord-marker=4 -O3 -ffree-line-length-none -g -fbacktrace $(OMPFLAG)
#FFLAGS = -x f95-cpp-input $(DEFINES)

# --- MPI, gfortran DEBUG syntax ------------------------------
#F90 = mpif90 -frecord-marker=4 -ffree-line-length-none -fbacktrace -g -O0 -fbounds-check -Wuninitialized -Wall $(OMPFLAG)
#FFLAGS = -x f95-cpp-input -ffpe-trap=zero,overflow,invalid -finit-real=nan  $(DEFINES)

# --- MPI, pgf90 syntax ------------------------------
#F90 = mpif90 -O3
#FFLAGS = -Mpreprocess $(DEFINES)

# --- MPI, ifort syntax ------------------------------
F90 = mpiifort $(OMPFLAG)
FFLAGS = -cpp -O3 -xAVX -fpp $(DEFINES) -DNOSYSTEM

# --- MPI, ifort syntax, additional checks -----------
#F90 = mpiifort $(OMPFLAG)
#FFLAGS = -warn all -O0 -g -traceback -fpe0 -ftrapuv -check bounds -cpp $(DEFINES) -DNOSYSTEM

# --- MPI, ifort syntax for Xeon phi ------------------------------
#F90 = mpiifort -Wl,-z,noexecstack $(OMPFLAG)
#FFLAGS = -cpp $(DEFINES) -O3 -fPIC -free -xMIC-AVX512 -fma -align array64byte -finline-functions -DNOSYSTEM

# --- MPI, ifort DEBUG syntax for Xeon phi ------------------------------
#F90 = mpiifort -Wl,-z,noexecstack $(OMPFLAG) 
#FFLAGS = -cpp -O0 -g -traceback -fpe0 -ftrapuv -check bounds $(DEFINES) -DNOSYSTEM

# --- MPI, ifort syntax ------------------------------
#F90 = ftn
#FFLAGS = -xAVX -g -traceback -fpp -fast $(DEFINES) -DNOSYSTEM #-DRT

# --- MPI, ifort syntax, additional checks -----------
#F90 = ftn
#FFLAGS = -O3 -g -traceback -fpe0 -ftrapuv -cpp $(DEFINES) -DNOSYSTEM #-DRT

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI =
#LIBMPI = -lfmpi -lmpi -lelan

# --- CUDA libraries, for Titane ---
LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart

ifeq ($(GRACKLE),1)
# Add include and library install path for grackle and hdf5 here
LIBS_GRACKLE = -L$(HOME)/local/lib -lgrackle -lhdf5 -lz -lgfortran -ldl
LIBS_OBJ     = -I$(HOME)/local/include -DCONFIG_BFLOAT_8 -DH5_USE_16_API -fPIC
endif
LIBS = $(LIBMPI) $(LIBS_GRACKLE)
#############################################################################
# Sources directories are searched in this exact order
VPATH := $(PATCH)
ifeq ($(RT),1)
VPATH := $(VPATH):../rt
endif
VPATH := $(VPATH):../$(SOLVER):../aton:../hydro:../pm:../poisson:../amr:../io:../amr/geometry
#############################################################################
# All objects
MODOBJ = mpi_mod.o amr_parameters.o amr_commons.o random.o pm_parameters.o pm_commons.o poisson_parameters.o dump_utils.o
ifeq ($(GRACKLE),1)
MODOBJ += grackle_parameters.o
endif
MODOBJ += poisson_commons.o hydro_parameters.o hydro_commons.o cooling_module.o bisection.o sparse_mat.o \
	 clfind_commons.o gadgetreadfile.o write_makefile.o write_patch.o  write_gitinfo.o
ifeq ($(RT),1)
MODOBJ += rt_parameters.o rt_hydro_commons.o coolrates_module.o \
         rt_spectra.o rt_cooling_module.o rt_flux_module.o
endif
ifeq ($(OPENMP),1)
MODOBJ += omp_utils.o
endif
ifeq ($(DICE),1)
AMROBJ = dice_commons.o
else
AMROBJ =
endif
AMROBJ += read_params.o init_amr.o init_time.o init_refine.o \
	tracer_utils.o adaptive_loop.o utils.o geom_types.o geom_utils.o \
        geom_distances.o geom_volumes.o sink_particle_tracer.o amr_step.o update_time.o \
         output_amr.o flag_utils.o physical_boundaries.o virtual_boundaries.o refine_utils.o nbors_utils.o \
         hilbert.o load_balance.o title.o sort.o cooling_fine.o units.o light_cone.o movie.o
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o           \
        move_fine.o newdt_fine.o particle_tree.o add_list.o           \
        remove_list.o star_formation.o          \
        sink_particle.o feedback.o clump_finder.o clump_merger.o      \
        flag_formation_sites.o init_sink.o output_sink.o              \
        mechanical_commons.o mechanical_fine.o
# Poisson solver objects
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o multigrid_coarse.o multigrid_fine_commons.o \
             multigrid_fine_fine.o multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.o godunov_fine.o \
           uplmde.o umuscl.o interpol_hydro.o godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o \
           boundana.o read_hydro_params.o synchro_hydro_fine.o
# RT objects
ifeq ($(RT),1)
RTOBJ = rt_init_hydro.o rt_init_xion.o rt_init.o rt_init_flow_fine.o rt_output_hydro.o rt_godunov_fine.o \
        rt_interpol_hydro.o rt_godunov_utils.o rt_condinit.o rt_hydro_flag.o rt_hydro_boundary.o rt_boundana.o \
        rt_read_hydro_params.o rt_units.o
endif

# Patch objects
$(foreach PP, $(subst :, ,$(PATCH)), $(wildcard $(shell [ -z $(PP) ] || -sinclude $(PP)/Makefile)))

# All objects
AMRLIB = $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ)
ifeq ($(RT),1)
AMRLIB += $(RTOBJ)
endif

# ATON objects
ATON_MODOBJ = timing.o radiation_commons.o rad_step.o
ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o rad_backup.o ../aton/atonlib/libaton.a
#############################################################################
ramses:	$(MODOBJ) $(AMRLIB) ramses.o
	 $(F90) $^ -o $(EXEC)$(NDIM)d $(LIBS)
	 rm write_makefile.f90
	 rm write_patch.f90
ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o
	 $(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o -o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)
	 rm write_makefile.f90
	 rm write_patch.f90
#############################################################################
write_gitinfo.o: FORCE
	 $(F90) $(FFLAGS) -DPATCH=\'$(PATCH)\' -DGITBRANCH=\'$(GITBRANCH)\' -DGITHASH=\'"$(GITHASH)"\' \
 -DGITREPO=\'$(GITREPO)\' -DBUILDDATE=\'"$(BUILDDATE)"\' -c ../amr/write_gitinfo.f90 -o $@
write_makefile.o: FORCE
	 ../utils/scripts/cr_write_makefile.sh $(MAKEFILE_LIST)
	 $(F90) $(FFLAGS) -c write_makefile.f90 -o $@
write_patch.o: FORCE
	 ../utils/scripts/cr_write_patch.sh $(PATCH)
	 $(F90) $(FFLAGS) -c write_patch.f90 -o $@
%.o:%.F
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
FORCE:
#############################################################################
clean:
	rm -f *.o *.$(MOD)
#############################################################################
print-%  : ; @echo $* = $($*)
