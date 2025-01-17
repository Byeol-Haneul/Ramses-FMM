#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@gmail.com
#############################################################################
# Compilation time parameters
COMPILER = GNU
NHILBERT = 1
NVECTOR = 32
NPRE = 8
NDIM = 3
NVAR = 5
HYDRO = 1
MHD = 0
RT = 0
GRAV = 0
MPI = 0
INIT =
UNITS =
NENER =
NRTGROUPS =
NIONS =
PATCH =
DEBUG = 0
HDF5 = 0
MPIF90 = mpif90
SOLVER = hydro
EXEC = ramses
#############################################################################
GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)
GITHASH = $(shell git log --pretty=format:'%H' -n 1)
GITREPO = $(shell git config --get remote.origin.url)
BUILDDATE = $(shell date +"%D-%T")
BUILDCOMMAND = ${MAKE} ${MAKEFLAGS}
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNVAR=$(NVAR) -DNHILBERT=$(NHILBERT)
ifeq ($(HYDRO),1)
   DEFINES += -DHYDRO
endif
ifeq ($(MHD),1)
   DEFINES += -DMHD
endif
ifeq ($(GRAV),1)
   DEFINES += -DGRAV
endif
ifneq ($(NENER),)
   DEFINES += -DNENER=$(NENER)
endif
ifeq ($(RT),1)
   DEFINES += -DRT
   #DEFINES += -DRT -DNIONS=$(NIONS) -DNGROUPS=$(NGROUPS)
endif
#############################################################################
# Set proper units and initial conditions
ifneq ($(INIT),)
   DEFINES += -DINIT=$(INIT)
endif
ifneq ($(UNITS),)
   DEFINES += -DUNITS=$(UNITS)
endif
#############################################################################
# GNU compiler (gfortran)
ifeq ($(COMPILER),GNU)
   FFLAGS = -cpp $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
      FC = gfortran
   else
      F90 = gfortran
      FC = gfortran
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -ffree-line-length-none -fimplicit-none
   FC += -ffree-line-length-none -fimplicit-none
   ifeq ($(DEBUG),1)
      F90 += -O0 -fcheck=bounds
      FFLAGS += -g -fbacktrace -ffpe-trap=zero,overflow,invalid -finit-real=nan
      F90 += -Wuninitialized -Wall
   else
      F90 += -O3
   endif
endif
#############################################################################
# Intel compiler
ifeq ($(COMPILER),INTEL)
   FFLAGS = -cpp $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
      FC = ifort
   else
      F90 = ifort
      FC = ifort
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -fp-model source
   FC += -fp-model source
   ifeq ($(DEBUG),1)
      F90 += -warn all -O0 -g -traceback -check bounds
      FFLAGS += -fpe0 -ftrapuv -init=zero -init=snan -init=arrays
   else
      F90 += -O3
   endif
endif
#############################################################################
# NVHPC compiler
ifeq ($(COMPILER),NVHPC)
   FFLAGS = -cpp $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
      FC = nvfortran
   else
      F90 = nvfortran
      FC = nvfortran
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -Mfreeform
   FC += -Mfreeform
   ifeq ($(DEBUG),1)
      F90 += -O0 -Mbounds -traceback
      FFLAGS += -O0 -ffpe-trap=zero,overflow,invalid -Minit-real=snan
   else
      F90 += -O4
   endif
endif
#############################################################################
C = ifort
CFLAGS = -g -O3
#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI =
#LIBMPI = -lfmpi -lmpi -lelan
ifeq ($(HDF5),1)
   LIBHDF5 = -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran
endif
LIBS = $(LIBMPI) $(LIBHDF5)
#############################################################################
# Sources directories are searched in this exact order
VPATH = $(PATCH):../$(SOLVER):../feedback:../star:../cooling:../hydro:../pm:../poisson:../amr:../clump:../mdl1:../rt
#############################################################################
# Modules
MODOBJ = amr_parameters.o hydro_parameters.o rt_parameters.o oct_commons.o \
	 hash.o pm_parameters.o pm_commons.o hilbert.o domain_hilbert.o \
	 hydro_commons.o rt_commons.o amr_commons.o poisson_parameters.o rngstream.o \
	 constants.o cooling_module.o cooling_module_ism.o init_cooling.o \
	 coolrates_module.o rt_cooling_module.o init_neq_chem.o \
	 mdl_commons.o mdl.o clfind_commons.o ramses_commons.o \
	 call_back.o cache_commons.o cache.o gadgetreadfile.o \
	 boundaries.o output_hydro.o output_poisson.o output_part.o output_rt.o output_amr.o movie.o \
	 sort.o nbors_utils.o hydro_flag.o multigrid_fine_coarse.o rho_fine.o marshal.o upload.o \
	 rt_flux_module.o rt_godunov_fine.o rt_input_condinit.o rt_init_flow_fine.o rt_units.o  \
	 rt_condinit.o rt_flag.o rt_upload.o interpol_rt.o rt_godunov_utils.o read_rt_params.o \
	 init_hydro.o init_part.o init_amr.o smooth.o flag_utils.o init_time.o read_params.o \
	 init_rt.o input_part_grafic.o input_part_zoom.o input_part_ascii.o \
	 input_part_restart.o input_part_ramses.o input_part_gadget.o input_part.o \
	 interpol_hydro.o update_time.o timer.o init_refine_basegrid.o load_balance.o refine_utils.o \
	 grav_ana.o phi_fine_cg.o interpol_phi.o force_fine.o poisson_flag.o multigrid_fine_commons.o \
	 boundana.o godunov_utils.o write_screen.o courant_fine.o newdt_fine.o godunov_fine.o \
	 umuscl.o synchro_hydro_fine.o source_hydro_fine.o cooling_fine.o \
	 input_hydro_condinit.o input_hydro_grafic.o input_hydro_gadget.o init_flow_fine.o \
	 condinit.o vecpotentialinit.o init_refine_adaptive.o init_refine_restart.o \
	 init_refine_ramses.o title.o units.o move_fine.o rho_ana.o \
	 star_formation.o feedback.o clump_merger.o output_clump.o clump_finder.o \
	 sink_formation.o tree_formation.o \
	 task_manager.o rt_step.o amr_step.o adaptive_loop.o \
	 write_gitinfo.o
#############################################################################
ramses:	$(MODOBJ) ramses.o
	$(F90) $(MODOBJ) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
#############################################################################
write_gitinfo.o: FORCE
	$(FC) -O0 -cpp -DPATCH=\'"$(PATCH)"\' -DGITBRANCH=\'"$(GITBRANCH)"\' -DGITHASH=\'"$(GITHASH)"\' \
                -DGITREPO=\'"$(GITREPO)"\' -DBUILDDATE=\'"$(BUILDDATE)"\' -DBUILDCOMMAND=\'"$(BUILDCOMMAND)"\' -c ../amr/write_gitinfo.f90 -o $@
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
%.o:%.c
	$(C) $(CFLAGS) -c $^ -o $@
FORCE:
#############################################################################
clean :
	rm -f *.o *.$(MOD)

clobber: clean
	-find . -type d -name "NDIM*" -prune -exec rm -rf '{}' \;
#############################################################################
